import sys,os,urllib,subprocess
import numpy as np
import requests
from tqdm import tqdm
from shapely.geometry import Point,LineString,Polygon,MultiPoint,MultiLineString,MultiPolygon,GeometryCollection
from shapely.ops import transform,cascaded_union,unary_union
from functools import partial
import pyproj
from pyproj import Proj

import mshapely
from mshapely import DF
from mshapely.misc import add_method
from .file import File

class OSM(object):
  """
  OSM object prepares the osm coastline using spatial manipulation for gmsh, such as
  downloading,extracting,simplifying and resampling for gmsh.
  It also creates the mesh using gmsh.
  
  Parameters
  ----------
  obj:obj
    name:str,
    format:str,
    localFolder:path
    minDensity:float
    maxDensity:float
    shorelineGrowth:float
    limitFineDensity:float
    simplification:object
      a:
    defaultDomain:object
      a:
      b:
    input:object
      a:
      b:
    output:object
      a:
    pproj:str
    pgeo:str
    
  Note
  ----
  Any spatial manipulation is performed on the projected coordinate system.
  Results are converted back to geographic coordinates.
  
  Attributes
  ----------
  osm:zip, OSM fine resolution coastline
  sosm:zip, Simplified OSM resolution coastline
  domain: Polygon
    Model domain
  density:MultiPoint
    MultiPoint of the density
  
  """
  def __init__(self,obj):
    self.name = obj.get( 'name', "default")
    self.format = obj.get( 'format', "slf")
    self.localFolder = obj.get( 'localFolder', "")
    self.minDensity = obj.get( 'minDensity', 10)
    self.limitFineDensity= obj.get( 'limitFineDensity', 1000)
    self.limitCoarseDensity= obj.get( 'limitCoarseDensity', 2000)
    self.maxDensity = obj.get( 'maxDensity', 10000)
    self.shorelineGrowth = obj.get( 'shorelineGrowth', 1.2)
    self.simplification = obj.get( 'simplification', {})
    self.defaultDomain = obj.get( 'defaultDomain', {})
    self.input=input= obj.get( 'input', {})
    self.output = obj.get( 'output', {})
    self.pproj = pproj = obj.get( 'proj', "EPSG:3573")
    self.pgeo  = pgeo = obj.get( 'proj', "EPSG:4326")
    self.printCommands = obj.get( 'printCommands', False)
    self.version=None
    
    # "+proj=laea +lat_0=90 +lon_0=-100 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m"
    self.togeo = partial(pyproj.transform,Proj(init=pproj),Proj(init=pgeo))
    self.toproj = partial(pyproj.transform,Proj(init=pgeo),Proj(init=pproj))
    
    #
    # Check input or get/create input if required
    #
    self._checkInput(input,'domain')
    self._checkInput(input,'density')
    self._checkInput(input,'osm')
    self._checkInput(input,'sosm')
    
    if not 'domain' in input:input['domain']=self._getDefaultDomain()
    if not 'density' in input:input['density']=self._getDefaultDensity()  
    if not 'osm' in input:input['osm']=OSM.downloadOSM(self.localFolder,"osm")
    if not 'sosm' in input:input['sosm']=OSM.downloadOSM(self.localFolder,"sosm")
    
    #
    # Input and temporay output are saved into a "File"
    #
    # Input
    self.osm=File('osm',parent=self,geoPath=input['osm'])
    self.sosm=File('sosm',parent=self,geoPath=input['sosm'])
    self.domain=File('domain',parent=self,geoPath=input['domain'])
    self.density=File('density',parent=self,geoPath=input['density'])
    
    # Temporary output
    self.densityFineZone=File('densityFineZone',parent=self,fproj=self._getdensityFineZone)
    self.densityCoarseZone=File('densityCoarseZone',parent=self,fproj=self._getdensityCoarseZone)
    self.osmFine=File('osmFine',parent=self,fproj=self._getOSMFine)
    self.osmCoarse=File('osmCoarse',parent=self,fproj=self._getOSMCoarse)
    self.osmCoarseZone=File('osmCoarseZone',parent=self,fproj=self._getOSMCoarseZone)
    self.osmCoarseS=File('osmCoarseS',versioned=True,parent=self,fproj=self._getOSMCoarseS)
    
    self.osmDomain=File('osmDomain',versioned=True,parent=self,fproj=self._getosmDomain)
    self.osmSimplify=File('osmSimplify',versioned=True,parent=self,fproj=self._getosmSimplify)
    
    self.osmResample=File('osmResample',versioned=True,parent=self,fproj=self._getosmResample)  
    self.osmMesh=File('osmMesh',parent=self,versioned=True,ext=".msh",fproj=self._getosmMesh)
    self.osmMeshBoundaries=File('osmMeshBoundaries',versioned=True,parent=self,fgeo=self._getosmMeshBoundaries)
    self.osmMeshEdges=File('osmMeshEdges',versioned=True,parent=self,fgeo=self._getosmMeshEdges)
    self.meshmbtiles=File('mesh',parent=self,versioned=True,ext=".mbtiles",fgeo=self._getmbtiles)
    
    #
    # Get list of File(s) within the OSM class
    # 
    self.listFiles=[p for p in dir(self) if isinstance(getattr(self,p),File)]
    
    # 
    # Documentation
    #
    self.densityFineZone.__doc__=self._getdensityFineZone.__doc__
    self.densityCoarseZone.__doc__=self._getdensityCoarseZone.__doc__
    self.osmFine.__doc__=self._getOSMFine.__doc__
    self.osmCoarse.__doc__=self._getOSMCoarse.__doc__
    self.osmCoarseZone.__doc__=self._getOSMCoarseZone.__doc__
    self.osmCoarseS.__doc__=self._getOSMCoarseS.__doc__
    self.osmDomain.__doc__=self._getosmDomain.__doc__
    self.osmSimplify.__doc__=self._getosmSimplify.__doc__
    self.osmResample.__doc__=self._getosmResample.__doc__
    self.osmMesh.__doc__=self._getosmMesh.__doc__
    self.osmMeshBoundaries.__doc__=self._getosmMeshBoundaries.__doc__
    self.osmMeshEdges.__doc__=self._getosmMeshEdges.__doc__
    self.meshmbtiles.__doc__=self._getmbtiles.__doc__
 
  def compute(self):
    
    
    # print(os.path.exists(self.osmMesh.geoPath))
    for name in self.listFiles:
      file=getattr(self,name)
      if file.versioned==True:
        file._geo=None
        file._proj=None
    for name in self.listFiles:
      file=getattr(self,name)
      file.geo()
        

  
  def setSimplificaton(self,obj,version=None):
    self.simplification=obj
    self.version=version
    self.compute()
    
    # for name in self.listFiles:
    #   file=getattr(self,name)
    #   if file.dependencies is not None and "simplification" in file.dependencies:
    #     file.delete()
    
  def _checkInput(self,input,name):
    """ Check if input exist
    """
    def pathExist(path):
      localFolder = self.localFolder
      return (os.path.exists(path) or os.path.exists(os.path.join(localFolder,path)))
      
    if name in input and not pathExist(input[name]):raise Exception("{} does not exist".format(input[name]))
    
 
  def _getDefaultDomain(self):
    """ Create default domain geojson
    """
    obj = self.defaultDomain
    toproj =self.toproj
    togeo =self.togeo
    output = os.path.join(self.localFolder,"domain.geojson")
    if os.path.exists(output):return output
    
    center = obj.get( 'center', [-63.553987,44.627934])
    radius = obj.get( 'radius', 10)
    
    geo=Point(center).proj(toproj).buffer(radius*1000)
    
    geo.proj(togeo).write(output)
    return output
    

  def _getDefaultDensity(self):
    """ Create default density geojson
    """
    obj = self.defaultDomain
    minDensity = self.minDensity
    maxDensity = self.maxDensity
    shorelineGrowth = self.shorelineGrowth
    
    output=os.path.join(self.localFolder,"density.geojson")
    if os.path.exists(output):return output
    density = np.array(obj.get( 'density', [[-63.553987,44.627934,1.0,1.2]]))
    
    df=DF(density,minDensity=np.min(density[:,2]),maxDensity=maxDensity,minGrowth=np.min(density[:,3]))
    
    df.write(output)
  
    return output

  def _getdensityFineZone(self,projectedPath=None,dependencies = ['density','minDensity','limitFineDensity']):
    """
    Create density zones/area based on density object.
    This zone is used to extract fined osm coastline.
    """
    
    return self.__getdensityZone(projectedPath,self.limitFineDensity)

  def _getdensityCoarseZone(self,projectedPath=None,dependencies = ['density','minDensity','limitCoarseDensity']):
    """
    Create density zones/area based on density object.
    This zone is used to extract fined osm coastline.
    """
    
    return self.__getdensityZone(projectedPath,self.limitCoarseDensity)

  def __getdensityZone(self,projectedPath,_maxDensity):
    """
    Create density zones/area based on density object.
    This zone is used to extract fined osm coastline.
    """
    minDensity=self.minDensity
    density=DF.read(self.density.projPath)
    
    buffers=[]
    for d in density.dp:
      maxDistance = DF.getl_D(d[2],d[3],_maxDensity)
      buffers.append(Point(d[:2]).buffer(maxDistance))

    buffers=cascaded_union(buffers)
    buffers.write(projectedPath)
    return buffers

  def _name(self,path):
    """ Extract basename without extention
    """
    return os.path.splitext(os.path.basename(path))[0]

  def _getOSMFine(self,projectedPath=None,dependencies = ['densityFineZone']):
    """ Extract fine osm coastline within the densityZone.
    """
    densityFineZone=self.densityFineZone
    epsg=self.pproj.split(":")[1]
    
    osmPath = self.input['osm']
    zipname = 'water-polygons-split-4326/water_polygons.shp'
    zipPath = "\"/vsizip/" + osmPath + "/" + zipname + "\""
    
    zoneName = os.path.splitext(os.path.basename(densityFineZone.projPath))[0]
    zone = densityFineZone.projPath
    # pg_sql = "\"With osm AS(SELECT ST_Transform(water_polygons.geometry,{2}) as geo FROM water_polygons,'{0}'.'{1}' zone WHERE ST_Intersects(ST_Transform(water_polygons.geometry,{2}), ST_Envelope(zone.geometry))),osm2 AS(SELECT ST_Simplify(ST_Buffer(ST_Simplify(osm.geo,10),0),10) as geo FROM osm) SELECT osm2.geo FROM osm2,'{0}'.'{1}' zone WHERE osm2.geo is NOT NULL;\"".format(zone,zoneName,epsg)
    pg_sql = "\"With osm AS(SELECT ST_Transform(water_polygons.geometry,{2}) as geo FROM water_polygons),osm2 AS(SELECT ST_Simplify(ST_Buffer(ST_Simplify(osm.geo,10),0),10) as geo FROM osm,'{0}'.'{1}' zone WHERE ST_Intersects(osm.geo, zone.geometry)) SELECT ST_Intersection(osm2.geo,zone.geometry) FROM osm2,'{0}'.'{1}' zone WHERE osm2.geo is NOT NULL;\"".format(zone,zoneName,epsg)
    command = "ogr2ogr -skipfailures -f \"GeoJSON\" {0} -nln \"{3}\" -nlt POLYGON -dialect \"SQLITE\" -sql {2} {1}".format(projectedPath,zipPath,pg_sql,self._name(projectedPath))
    if self.printCommands: print(command)
    t=tqdm(total=1)
    subprocess.call(command, shell=True)
    t.update(1);t.close()

  def _getOSMCoarse(self,projectedPath=None,dependencies = ['domain']):
    """
    Extract osm coastline from zip file.
    This avoids unpacking the zip file.
    """
    domain=self.domain
    osmPath=self.sosm.geoPath
    epsgp=self.pproj.split(":")[1]
    epsgg=self.pgeo.split(":")[1]
    
    domain=domain.geoPath
    zipname = 'simplified-water-polygons-split-3857/simplified_water_polygons.shp'
    zipPath = "\"/vsizip/" + osmPath + "/" + zipname + "\""
    name = os.path.basename(domain)
    name = os.path.splitext(name)[0]

    pg_sql = "\"With one AS(SELECT ST_Buffer(ST_Transform(A.geometry,{2}),0) as geometry FROM simplified_water_polygons A,'{0}'.'{1}' B WHERE ST_Intersects(ST_Transform(A.geometry,{2}), ST_Transform(SetSRID(B.geometry,{3}),{2}))) SELECT ST_Union(one.geometry) from one WHERE one.geometry is not null;\"".format(domain,name,epsgp,epsgg)
    
    
    
    
    command = "ogr2ogr -skipfailures -f \"GeoJSON\" {0} -nln \"{3}\" -dialect \"SQLITE\" -sql {1} {2}".format(projectedPath,pg_sql,zipPath,self._name(projectedPath))
    if self.printCommands: print(command)
    t=tqdm(total=1)
    subprocess.call(command, shell=True)
    t.update(1);t.close()

  def _getOSMCoarseZone(self,projectedPath=None,dependencies = ['osmCoarse','densityCoarseZone']):
    """
    Extract osm coastline from zip file and simplify based on extent.
    This avoids unpacking the zip file.
    """
    osmCoarse=self.osmCoarse
    densityCoarseZone=self.densityCoarseZone
    
    pg_sql = "\"SELECT ST_Intersection(A.geometry,B.geometry) as geometry FROM '{0}'.'{1}' A,'{2}'.'{3}' B;\"".format(osmCoarse.projPath,self._name(osmCoarse.projPath),densityCoarseZone.projPath,self._name(densityCoarseZone.projPath))
    command = "ogr2ogr -skipfailures -f \"GeoJSON\" {0} -nln \"{3}\" -dialect \"SQLITE\" -sql {1} {2}".format(projectedPath,pg_sql,osmCoarse.projPath,self._name(projectedPath))
    if self.printCommands: print(command)
    t=tqdm(total=1)
    subprocess.call(command, shell=True)
    t.update(1);t.close()
    
    
  def _getOSMCoarseS(self,projectedPath=None,dependencies=['osmCoarse','simplification']):
    """
    Extract osm coastline from zip file and simplify based on extent.
    This avoids unpacking the zip file.
    
    Warning: needs RAM
    """
    osmCoarse=self.osmCoarse
    simplification = self.simplification
    
    
    osmPath=osmCoarse.projPath
    isimplify=simplification['isimplify']
    buffering=simplification['buffer']
    fsimplify=simplification['fsimplify']
    
    # ogrinfo 
    
    # Needs at least 2.6GB without simplify
    # buffer1000,s50=15minutes
    
    # Simplify 500,buffer1000,s50=10minutes
    # Simplify 1000,buffer10000,s50=30sec
    # Simplify 500,buffer5000,s50=4msec
    pg_sql = "\"With one AS(SELECT ST_Simplify(ST_Buffer(ST_Buffer(ST_Simplify(A.geometry,{2}),-{3}),{3}),{4}) as geometry FROM '{0}'.'{1}' A) SELECT one.geometry from one WHERE one.geometry is not null;\"".format(osmPath,self._name(osmPath),isimplify,buffering,fsimplify)
    
    
    
    command = "ogr2ogr -skipfailures -f \"GeoJSON\" {0} -nln \"{3}\" -dialect \"SQLITE\" -sql {1} {2}".format(projectedPath,pg_sql,osmPath,self._name(projectedPath))
    if self.printCommands: print(command)
    t=tqdm(total=1)
    subprocess.call(command, shell=True)
    t.update(1);t.close()    
    

  def _getosmDomain(self,projectedPath=None,dependencies = ['domain','osmCoarseS']):
    """
    Extract osm coastline using the domain.
    It will only keep the largest Polygon.
    """
    
    domain = self.domain
    osmCoarseS = self.osmCoarseS
    
    geo = osmCoarseS.proj().geometry
    domain=domain.proj().geometry
    
    t=tqdm(total=1)
    geo=geo.largest().removeHoles(np.pi*np.power(5000,2))
    geo=geo.intersection(domain)
    geo.write(projectedPath).plot().savePlot(os.path.splitext(projectedPath)[0]+".png")
    t.update(1);t.close()
    
    
    return geo


  def _getosmSimplify(self,projectedPath=None,dependencies = ['density','osmDomain','osmFine']):
    """
    Simplify osm shoreline based on density field
    """
    df = DF.read(self.density.projPath)
    osmDomain = self.osmDomain
    osmFine = self.osmFine
    osmCoarseZone = self.osmCoarseZone
    
    geo=osmDomain.proj().geometry
    
    geo=geo.dsimplify(df,limitFineDensity=self.limitFineDensity,limitCoarseDensity=self.limitCoarseDensity,fine=osmFine.proj().geometry,coarse=osmCoarseZone.proj().geometry,progress=True)
    geo=geo.largest()
    geo.write(projectedPath).plot().savePlot(os.path.splitext(projectedPath)[0]+".png")
    return geo
    
  def _getosmResample(self,projectedPath=None,dependencies = ['density','osmSimplify','minDensity','maxDensity','shorelineGrowth']):
    """
    Resample osm shoreline using interior nearest points and density growth field.
    """
    
    # df = DF.read(self.density.projPath)
    osmSimplify = self.osmSimplify
    minDensity = self.minDensity
    maxDensity = self.maxDensity
    shorelineGrowth = self.shorelineGrowth
    
    geo=osmSimplify.proj().geometry
    df=DF(minDensity=minDensity,maxDensity=maxDensity,minGrowth=shorelineGrowth,maxDensitySimplify=10000,progress=True)
    df=df.inearest(geo,progress=True,minDistance=100)
    
    geo=geo.dresample(df,progress=True)
    
    geo.write(projectedPath).plot().savePlot(os.path.splitext(projectedPath)[0]+".png")
    return geo
  
  def _getosmMesh(self,projectedPath=None,dependencies = ['density','osmResample','minDensity','maxDensity','shorelineGrowth']):
    """
    Resample osm shoreline using interior nearest points and density growth field.
    """
    # df = DF.read(self.density.projPath)
    osmResample = self.osmResample
    minDensity = self.minDensity
    maxDensity = self.maxDensity
    shorelineGrowth = self.shorelineGrowth
    
    geo=osmResample.proj().geometry
    df=DF(minDensity=minDensity,maxDensity=maxDensity,minGrowth=shorelineGrowth,maxDensitySimplify=10000,progress=True)
    df=df.inearest(geo,progress=True,minDistance=10,minLength=True)
    
    geo.msh(projectedPath,df).plot().savePlot(os.path.splitext(projectedPath)[0]+".png")
  
  def _getosmMeshBoundaries(self,geographicPath=None,dependencies=['osmMesh']):
    """
    """
    
    mesh=self.osmMesh.geo()
    mesh.boundaries.write(geographicPath)
  
  def _getosmMeshEdges(self,geographicPath=None,dependencies=['osmMesh']):
    """
    """
    mesh=self.osmMesh.geo()
    mesh.geoedges.write(geographicPath)
  
  def _getmbtiles(self,geographicPath,dependencies=['osmMeshBoundaries','osmMeshEdges']):
    edgembitle=os.path.join(os.path.dirname(geographicPath),"edges.mbtiles")
    boundarymbtile=os.path.join(os.path.dirname(geographicPath),"boundaries.mbtiles")
    if os.path.exists(edgembitle):os.remove(edgembitle)
    if os.path.exists(boundarymbtile):os.remove(boundarymbtile)
    command = "tippecanoe -z11 -o {0} -an -l edges {1};tippecanoe  -z11 -o {2} -an -l boundaries {3};tile-join {0} {2} -o {4}".format(edgembitle,self.osmMeshEdges.geoPath,boundarymbtile,self.osmMeshBoundaries.geoPath,geographicPath)
    if self.printCommands: print(command)
    subprocess.call(command, shell=True)
  
  
  @staticmethod
  def transform_geo(project,geo):
    """
    project:proj4
    geo:Shapely object
    """
    return transform(project,geo)
  
  
  @staticmethod
  def ogr2ogrT(inputPath,output,s_srs,t_srs,zipLayer=""):
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
    
    
  @staticmethod
  def downloadOSM(folder,option,overwrite=False):
    """
    Download OSM coastline file (600MB). 
    This geometry is splitted into partions (simarlar to a kdtree).
    option=[osm,sosm,...]
    """
    if option=="osm":
      http = 'https://osmdata.openstreetmap.de/download/water-polygons-split-4326.zip'
    elif option=="sosm":
      http="https://osmdata.openstreetmap.de/download/simplified-water-polygons-split-3857.zip"
    else:
      raise Exception("Not a choice")
    
    name = os.path.basename(http)
    osmPath =os.path.join(folder,name)
    
    if not os.path.exists(osmPath) or overwrite:
      response = requests.get(http, stream=True)
      total_length = int(response.headers.get('content-length', 0))
      t=tqdm(total=total_length, unit='iB', unit_scale=True)
      with open(osmPath, "wb") as f:
        for chunk in response.iter_content(chunk_size=1024):
          if chunk: # filter out keep-alive new chunks
            t.update(len(chunk))
            f.write(chunk)
      t.close()
    return osmPath