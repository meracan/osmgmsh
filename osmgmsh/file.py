import sys,os,urllib,subprocess
import inspect
from pyproj import transform,Proj
from mshapely import GIS
from mmesh import MESH,MSLF,MMSH


class File(object):
  def __init__(self,name,ext=".geojson",parent=None,localFolder="",geoPath=None,fproj=None,fgeo=None,pproj="EPSG:3573",pgeo="EPSG:4326"):
    self.parent=parent
    self.name=name
    
    if parent is not None:
      localFolder=getattr(parent, "localFolder")
      pproj=parent.pproj
      pgeo=parent.pgeo
    self.localFolder=localFolder
    self.pproj = pproj
    self.pgeo = pgeo
    
    if geoPath is None:
      geoPath=os.path.join(localFolder,name+ext)
    else:
      ext=os.path.splitext(geoPath)[1]
    
    self.ext=ext
    self.geoPath=geoPath
    
    self.projPath=self._projPath(geoPath)
    
    self._proj = None
    self._geo = None
    
    self._fproj = fproj
    self._fgeo = fgeo
    
    self.dependencies=None
  
    if fproj is not None:
      s = inspect.signature(fproj)
      self.dependencies=s.parameters['dependencies'].default
    if fgeo is not None:
      s = inspect.signature(fgeo)
      self.dependencies=s.parameters['dependencies'].default
  
  @property
  def proj(self):
    if self._proj is None:
      if not os.path.exists(self.projPath):
        if self._fproj is not None:self._fproj(self.projPath)
        else:
          self.geo
          if self._fproj is None and self._geo is None:raise Exception("Data does not exist and does not know how to get it")
          self._toproj(self.geoPath,self.projPath)  
      self._proj = self.read(self.projPath)
    return self._proj
    
  @property
  def geo(self):
    if self._geo is None:
      
      if not os.path.exists(self.geoPath):
        if self._fgeo is not None:self._fgeo(self.geoPath)
        else:
          self.proj
          if self._fproj is None and self._geo is None:raise Exception("Data does not exist and does not know how to get it")
          self._togeo(self.geoPath,self.projPath) 
      self._geo = self.read(self.geoPath)
    return self._geo
  
  def read(self,path):
    if self.ext==".geojson":return GIS.read(path)
    if self.ext==".slf" or self.ext==".msh":return MESH.read(path)
  
  def delete(self):
    """
    Delete file and child dependencies
    
    """
    if os.path.exists(self.geoPath):os.remove(self.geoPath)
    if os.path.exists(self.projPath):os.remove(self.projPath)
    
    # Delete child dependencies
    if self.parent is not None:
      for file_name in self.parent.listFiles:
        file = getattr(self.parent, file_name)
        if file.dependencies is not None and file.dependencies.contains(self.name):
          file.delete()
  
  def _projPath(self,path):
    """
    Get projected path of the output 
    """
    return "{}.proj{}".format(os.path.splitext(path)[0],self.ext) 
  
  def _toproj(self,geographicPath,projectedPath):
    self._transform(geographicPath,projectedPath,self.pgeo,self.pproj)
  
  def _togeo(self,geographicPath,projectedPath):
    self._transform(projectedPath,geographicPath,self.pproj,self.pgeo)
  
  def _transform(self,source,target,sproj,tproj):
    """
    Transform geographic/projection
    """
    if not os.path.exists(source):raise Exception("{} does not exist".format(source))
    if os.path.exists(target):os.remove(target)
    if self.ext==".geojson":File.ogr2ogr(source,target,sproj,tproj)
    elif self.ext==".slf" or self.ext==".msh": 
      mesh=MESH.read(source)
      p1=Proj(init=sproj)
      p2=Proj(init=tproj)
      mesh.setXY(*transform(p1,p2,mesh.x,mesh.y)).write(target)
    else:raise Exception("Method for {} does not exist".format(self.ext))


  @staticmethod
  def ogr2ogr(inputPath,output,s_srs,t_srs,zipLayer=""):
    """
    zipLayer:To determine the layer within a zipfile:
      >>> vim {path}.zip
    """
    if os.path.splitext(inputPath)[1]==".zip":
      basename = os.path.splitext(os.path.basename(inputPath))[0]
      zipname = '{}/{}.shp'.format(basename,zipLayer)
      inputPath = "\"/vsizip/" + inputPath + "/" + zipname + "\""
    
    command = "ogr2ogr -skipfailures -f \"GeoJSON\" {0} -s_srs \"{1}\" -t_srs \"{2}\" {3}".format(output,s_srs,t_srs,inputPath)
    print(command)
    subprocess.call(command, shell=True)