import sys,os,urllib,subprocess
import inspect
from pyproj import transform,Proj
from mshapely import GIS
from mmesh import MESH


class File(object):
  """
  File object handles writing and reading geometric/mesh object to files, and vice-versa.
  To create the file, it either use fproj or fgeo function.
  fproj and fgeo uses projected and geographic data, respectively.
  The files are automatically transformed from proj data to fgeo data and vice-versa after processing the fproj/fgeo
  
  Parameters
  ----------
  name:str
    Must be a unique name
  ext:str
    File extention
    Default is ".geojson
  parent:Object
    Parent object. Parent attributes has priority.
  localFolder:path
  geoPath:path
  fproj:function
    Function to use to create the file in projected
  fgeo:function
    Function to use to create the file in geographic
  pproj:str
    Default is EPSG:3573
  pgeo:str
    Default is EPSG:4326    
  
  Attributes
  ----------
  proj:Object in projected
  geo:Object in geographic
  dependencies:list of attributes
    List of dependencies in the parent object
  """
  def __init__(self,name,ext=".geojson",parent=None,localFolder=None,versioned=False,geoPath=None,fproj=None,fgeo=None,pproj="EPSG:3573",pgeo="EPSG:4326"):
    self.parent=parent
    self.name=name
    self.fixgeoPath=None
    self.versioned=versioned
    
    _localFolder=""
    if parent is not None:
      _localFolder=getattr(parent, "localFolder")
      pproj=parent.pproj
      pgeo=parent.pgeo
    if localFolder is not None:_localFolder=localFolder   
  
    self.localFolder=_localFolder
    self.pproj = pproj
    self.pgeo = pgeo
    
    if geoPath is not None:
      self.fixgeoPath=geoPath
      ext=os.path.splitext(geoPath)[1]
    self.ext=ext
      
  
    
    
    
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
  
  def _getDenpendency(self):
    for name in self.dependencies:
      att=getattr(self.parent,name)
      if isinstance(att,File):
        att.proj
        att.geo
        
  @property
  def proj(self):
    """ Get projected file object
    """
    if self._proj is None:
      if not os.path.exists(self.projPath):
        if self._fproj is not None:
          self._getDenpendency()
          self._fproj(self.projPath)
        else:
          self.geo
          if self._fproj is None and self._geo is None:raise Exception("Data does not exist and does not know how to get it")
          self._toproj(self.geoPath,self.projPath)  
      self._proj = self.read(self.projPath)
    return self._proj
    
  @property
  def geo(self):
    """ Get geographic file object
    """
    
    if self._geo is None:
      
      if not os.path.exists(self.geoPath):
        if self._fgeo is not None:
          self._getDenpendency()
          self._fgeo(self.geoPath)
        else:
          self.proj
          if self._fproj is None and self._geo is None:raise Exception("Data does not exist and does not know how to get it")
          self._togeo(self.geoPath,self.projPath) 
      self._geo = self.read(self.geoPath)
    return self._geo
  
  def read(self,path):
    """ Read file to object
    """
    if self.ext==".geojson":return GIS.read(path)
    if self.ext==".slf" or self.ext==".msh":return MESH.read(path)
  
  def delete(self):
    """ Delete file and child dependencies
    """
    if os.path.exists(self.geoPath):os.remove(self.geoPath)
    if os.path.exists(self.projPath):os.remove(self.projPath)
    
    #Delete child dependencies
    if self.parent is not None:
      for name in self.parent.listFiles:
        file = getattr(self.parent, name)
        if file.dependencies is not None and self.name in file.dependencies:
          file.delete()
  
  @property
  def geoPath(self):
    if self.fixgeoPath is None:
      if self.versioned and self.parent is not None and self.parent.version is not None:geoPath=os.path.join(self.localFolder,"{0}.{1}{2}".format(self.name,self.parent.version,self.ext))
      else:geoPath=os.path.join(self.localFolder,self.name+self.ext)
      return geoPath
    else:
      return self.fixgeoPath
  
  @property
  def projPath(self):
    return self._projPath(self.geoPath)
  
  def _projPath(self,path):
    """ Get projected path of the output 
    """
    return "{}.proj{}".format(os.path.splitext(path)[0],self.ext) 
  
  def _toproj(self,geographicPath,projectedPath):
    """ Transform data to projected
    """
    self._transform(geographicPath,projectedPath,self.pgeo,self.pproj)
  
  def _togeo(self,geographicPath,projectedPath):
    """ Transform data to geographic
    """
    self._transform(projectedPath,geographicPath,self.pproj,self.pgeo)
  
  def _transform(self,source,target,sproj,tproj):
    """ Transform geographic/projection data
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
    """ Transform data usng ogr2ogr
    """
    """
    zipLayer:To determine the layer within a zipfile:
      >>> vim {path}.zip
    """
    if os.path.splitext(inputPath)[1]==".zip":
      basename = os.path.splitext(os.path.basename(inputPath))[0]
      zipname = '{}/{}.shp'.format(basename,zipLayer)
      inputPath = "\"/vsizip/" + inputPath + "/" + zipname + "\""
    
    command = "ogr2ogr -skipfailures -f \"GeoJSON\" {0} -s_srs \"{1}\" -t_srs \"{2}\" {3}".format(output,s_srs,t_srs,inputPath)
    # print(command)
    subprocess.call(command, shell=True)