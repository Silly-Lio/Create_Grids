#!/usr/bin/env python
# coding: utf-8

import geopandas as gpd
import numpy as np
import shapely.geometry as geom
import multiprocessing

from pyproj import Transformer
from osmnx import projection

import os

#获得边界信息
bound_path = 'D:/_DATA/Wuchang.shp'

Boundary = gpd.read_file(bound_path)
Road_type = 'drive'

#对所获取的原始边界数据Boundary进行处理,Boundary 默认为polygon,GeoDataFrame类型
def PreProcessed(Boundary):
    features = dict()
    Boundary = Boundary.to_crs('epsg:4326')
    features['Bd'] = Boundary
    features['Bd_ugeom'] = Boundary.unary_union
    features['Bd_utm_ugeom'], crs_utm = projection.project_geometry(features['Bd_ugeom'])
    features['Bd_utm'] = Boundary.to_crs(crs_utm)
    features['Bd_utm_ugeom'] = features['Bd_utm'].unary_union
    features['Bd_utm_Convex'] = features['Bd_utm_ugeom'].convex_hull
    features['crs_utm'] = crs_utm
    return features


def Polygon_projection(features):
    Boundary = features['Bd']
    crs_utm = features['crs_utm']
    region_polygon = features['Bd_ugeom']
    
    transformer  = Transformer.from_crs(Boundary.crs, crs_utm)
    sw_x, sw_y = transformer.transform(region_polygon.bounds[1], region_polygon.bounds[0])
    ne_x, ne_y = transformer.transform(region_polygon.bounds[3], region_polygon.bounds[2])
    return sw_x, sw_y, ne_x, ne_y#,crs_utm


def Grids(extent,grid_size):
    sw_x,sw_y,ne_x,ne_y = extent
    #这里使用[for]与map或将更快，有时间可重写
    for x in np.arange(sw_x, ne_x, grid_size):
        for y in np.arange(sw_y, ne_y, grid_size):
            grid = geom.Polygon([(x,y), (x+grid_size, y), (x+grid_size, y+grid_size), (x, y+grid_size)])
            yield grid
    
	        
def Hex_Grids(extent, grid_size):
    sw_x, sw_y, ne_x, ne_y = extent
    i = 0
    y_step = grid_size*3**(1/2)
    x_step = grid_size*3
    x = sw_x
    while x < ne_x:
        i += 1
        if i%2 == 1:
            for y in np.arange(sw_y, ne_y, y_step):
                hexgrid = geom.Polygon([(x-x_step/6,y-y_step/2), (x+x_step/6, y-y_step/2), (x+x_step/3, y),
                                     (x+x_step/6,y+y_step/2), (x-x_step/6, y+y_step/2), (x-x_step/3, y)])
                yield hexgrid
        else:
            for y in np.arange(sw_y+y_step/2, ne_y+y_step/2, y_step):
                hexgrid = geom.Polygon([(x-x_step/6,y-y_step/2), (x+x_step/6, y-y_step/2), (x+x_step/3, y),
                                     (x+x_step/6,y+y_step/2), (x-x_step/6, y+y_step/2), (x-x_step/3, y)])
                yield hexgrid
                
        x += x_step/2


#已获得对角坐标
#仅创建格子
def _create_grids(extent, gridsize, gridtype):
    if gridtype == 'grid':
        grids = Grids(extent, gridsize)
    elif gridtype == 'hex':
        grids = Hex_Grids(extent, gridsize)
    return grids


#格子与原多边形的相交判断
def choose_intersection(gridgene, Convex_Hull_utm, Boundary_utm_geom):
    for grid in gridgene:
        if grid.intersection(Convex_Hull_utm) & grid.intersection(Boundary_utm_geom):
            yield grid
        else:
            pass


def Create_Grids(OriginBd, gridtype = 'grid', gridsize = 2000):
    #crs = 'UTM xx'
    features = PreProcessed(OriginBd)
    extent = Polygon_projection(features)
    grids = _create_grids(extent, gridsize, gridtype)
    grids = choose_intersection(grids, features['Bd_utm_Convex'], features['Bd_utm_ugeom'])
    grids_gs = gpd.GeoSeries(grids, crs = features['crs_utm'])
    return grids_gs


def grid_feat(shpname):
    '''
    Boundary = gpd.read_file(f'D:/_DATA/{shpname}.shp')
    grids = Create_Grids(Boundary, 'grid', 3000)
    grids_wgs = grids.to_crs('epsg:4326')
    '''
    grids_wgs = gpd.read_file(f'D:/_DATA/4cities_grid/{shpname}.shp').geometry
    grid_gdf = gpd.GeoDataFrame(geometry = grids_wgs)

    buildings = gpd.read_file(f'D:/_DATA/buildings_{shpname}.shp').set_index('index')
    buildings =  gpd.overlay(buildings, grid_gdf)
    
    #bd4intsec = buildings.copy()
    #bd4intsec.reset_index(inplace = True)

    print(shpname)
    for idx in grids_wgs.index:
        #if not os.path.exists(f'D:/_DATA/segments/{shpname}_{idx}.shp'):
            #agrid = gpd.GeoDataFrame(geometry = [grids_wgs[idx]], crs = grids_wgs.crs)
            #interected = gpd.overlay(agrid, bd4intsec)
            intersected = buildings.within(grids_wgs[idx])
            if intersected.any():
                intersected = buildings[intersected.values]
                intersected.to_file(f'D:/_DATA/segments/{shpname}_{idx}.shp')
                buildings.drop(intersected.index, inplace = True)
                print(f'{shpname}{idx}')
            else:
                print(f'{shpname}{idx} empty')

if __name__ == '__main__':
	pool = multiprocessing.Pool(processes=4)
	pool.map(grid_feat,['BJGrid','SHGrid','SZGrid','WHGrid'])
    #grids_wgs.to_file(f'D:/_DATA/4cities_grid/{city[0:2]}Grid.shp')
	