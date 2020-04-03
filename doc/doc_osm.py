import numpy as np
from shapely.geometry import Point

import mshapely
from osmgmsh import OSM

def doc_example1():
  obj = {
      "name":"example1", 
      "format":"slf", 
      "localFolder":"../data/example1",
      "minDensity":10,
      "maxDensity":10000,
      "shorelineGrowth":1.2,
      "simplification":{"isimplify":10,"buffer":1000,"fsimplify":10},
      "defaultDomain":{"center":[-63.553987,44.627934],"radius":60,"density":[[-63.553987,44.627934,10,1.2]]},
      
      "input":{
        "osm":"../data/water-polygons-split-4326.zip",
        "sosm":"../data/simplified-water-polygons-split-3857.zip",
      },
    }
  
  osm=OSM(obj)
  # osm.densityFineZone.geo
  # osm.densityCoarseZone.geo
  # osm.osmFine.geo
  # osm.osmCoarse.geo
  # osm.osmCoarseZone.geo
  # osm.osmCoarseS.geo
  # osm.osmDomain.geo
  # osm.osmSimplify.geo
  # osm.osmResample.geo
  # osm.osmMesh.geo
  # osm.osmMeshBoundaries.geo
  # osm.osmMeshEdges.geo
  osm.meshmbtiles.geo

def doc_example2():
  obj = {
      "name":"example2", 
      "format":"slf", 
      "localFolder":"../data/example2",
      "minDensity":10,
      "limitFineDensity":1000,
      "maxDensity":100000,
      "shorelineGrowth":1.2,
      
      
      # 0,1000,50=15minutes -> Needs at least 2.6GB
      # 500,1000,50=10minutes
      # 1000,10000,50=30sec
      # 500,5000,50=4msec
      "simplification":{"isimplify":500,"buffer":5000,"fsimplify":100},
      
      "defaultDomain":{"center":[-90,67],"radius":2000,"density":10,"growth":1.2},
      
      "input":{
        "osm":"../data/water-polygons-split-4326.zip",
        "sosm":"../data/simplified-water-polygons-split-3857.zip",
      },
    }
  
  osm=OSM(obj)
  # osm.osmFine.geo #-> min
  # osm.densityCoarseZone.geo
  # osm.osmCoarse.geo #-> 4min
  # osm.osmCoarseZone.geo
  # osm.osmCoarseS.geo #-> 2min
  # osm.osmDomain.geo
  # osm.osmSimplify.geo #->1min
  # osm.osmResample.geo
  
  
if __name__ == "__main__":
  doc_example1()
  # doc_example2()
  