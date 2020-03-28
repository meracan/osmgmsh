import pytest
import os
import numpy as np
from shapely.geometry import Point,LineString,Polygon,MultiPoint,MultiLineString,MultiPolygon,GeometryCollection
import mshapely





def test_togmsh():
  polygon = Point((0,0)).buffer(100)
  hole1 = Point((-50,0)).buffer(20)
  hole2 = Point((50,0)).buffer(20)
  polygon = Polygon(polygon.exterior,[hole1.exterior.coords[::-1],hole2.exterior.coords[::-1]])
  density=polygon.inearest(maxDistance=100,angle=90)
  density[:,2]=density[:,2]*0.1
  density=np.column_stack((density,np.ones(len(density))*1.2))
  
  polygon=polygon.dresample(density,minDensity=1,maxDensity=100)
  density=polygon.inearest(maxDistance=100,angle=90)
  density[:,2]=density[:,2]*0.1
  density=np.column_stack((density,np.ones(len(density))*1.2))
  polygon.msh("test/data/test.msh",density,minDensity=1,maxDensity=100).plot("test/data/test.png")
  
  
  
  # polygon = Polygon([(0, 0), (0, 1),(1,1),(1,0),(0,0)],[LineString([(0.25, 0.25), (0.25, 0.75),(0.75,0.75),(0.75,0.25),(0.25,0.25)])])
  # getGeo(polygon)
  
  # array=np.array([[0,0,11,1.2],[0,0,233,1.2],[0,0,107,1.2],[0,0,1007,1.2],[0,0,12,1.2],[0,0,950,1.2]])
  # getAttractors(array)
  
if __name__ == "__main__":
  test_togmsh()