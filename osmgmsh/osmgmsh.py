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

# from .io import createGEO,createMSH

# #
# # Create Gmsh
# #
# @add_method(GeometryCollection)
# def msh(self,*args,**kwargs):
#   return self.toShape().msh(*args,**kwargs)


# @add_method([Point,MultiPoint,LineString,MultiLineString])
# def msh(self,*args,**kwargs):
#   return self
  
# @add_method(Polygon)
# def msh(self,path,density,*args,**kwargs):
#   geo = createGEO(self,path,density,*args,**kwargs)
#   return createMSH(geo,path)

# @add_method(MultiPolygon)
# def msh(self,*args,**kwargs):
#   geo=self.largest()
#   geo.msh(*args,**kwargs)
#   return geo


# class File(object):
#   def __init__(self,parent,name,fproj=None,fgeo=None,geoPath=None):
#     self.name=name
#     self.parent=parent
#     if geoPath is None:
#       geoPath=os.path.join(parent.localFolder,name+".geojson")
    
#     self.geoPath=geoPath
#     self.projPath=self._projPath(geoPath)
    
#     self._fproj = fproj
#     self._fgeo = fgeo
#     self._proj = None
#     self._geo = None
    
#     self.dependencies=None
#     if fproj is not None:self.dependencies=fproj(None,None,getDependenciesOnly=True)
#     if fgeo is not None:self.dependencies=fgeo(None,None,getDependenciesOnly=True)
  
#   @property
#   def proj(self):
#     if self._proj is None:
#       if not os.path.exists(self.projPath):
#         if self._fproj is not None:self._fproj(self.geoPath,self.projPath)
#         else:
#           self.geo
#           self._toproj(self.geoPath,self.projPath)  
#       self._proj = mshapely.readGeometry(self.projPath)
#     return self._proj
    
#   @property
#   def geo(self):
#     if self._geo is None:
#       if not os.path.exists(self.geoPath):
#         if self._fgeo is not None:self._fgeo(self.geoPath,self.projPath)
#         else:
#           self.proj
#           self._togeo(self.geoPath,self.projPath) 
#       self._geo = mshapely.readGeometry(self.geoPath)
#     return self._geo
  
  
#   def delete(self):
#     """
#     Delete file and child dependencies
    
#     """
#     if os.path.exists(self.geoPath):os.remove(self.geoPath)
#     if os.path.exists(self.projPath):os.remove(self.projPath)
    
#     # Delete child dependencies
#     for file_name in self.parent.listFiles:
#         file = getattr(self.parent, file_name)
#         if file.dependencies is not None and file.dependencies.contains(self.name):
#           file.delete()
    
  
#   def _projPath(self,path):
#     """
#     Get projected path of the output 
#     """
#     return "{}.proj.geojson".format(os.path.splitext(path)[0]) 
  
#   def _geoPath(self,path):
#     """
#     Get geographic path of the output 
#     """
#     return "{}.geojson".format(os.path.splitext(path)[0])

#   def _toproj(self,geographicPath,projectedPath):
#     """
#     Transform geographic to projection
#     """
#     if not os.path.exists(geographicPath):raise Exception("{} does not exist".format(geographicPath))
#     if os.path.exists(projectedPath):os.remove(projectedPath)
#     OSM.ogr2ogrT(geographicPath,projectedPath,self.parent.geo,self.parent.proj)

#   def _togeo(self,geographicPath,projectedPath):
#     """
#     Transform geographic to projection
#     """
#     if not os.path.exists(projectedPath):raise Exception("{} does not exist".format(projectedPath))
#     if os.path.exists(geographicPath):os.remove(geographicPath)
#     OSM.ogr2ogrT(projectedPath,geographicPath,self.parent.proj,self.parent.geo)




class OSM(object):
  """
  OSM instance prepares the osm coastline using spatial manipulation for gmsh, such as
  downloading,extracting,simplifying and resampling for gmsh.
  
  Parameters
  ----------
  The instance requires one object parameter that contains input and output variables.
  
  obj: object
    path: object
      osm:path, input
      domain:path,input
      density:path,input
      osmDomain:path,output
      osmSimplify:path,output
      osmResample:path,output
    minDensity:float,input
    maxDensity:float,input
    
    
  Note
  ----
  Any spatial manipulation is performed on the LAEA projection (north pole).
  Results are converted back to geographic coordinates.
  """
  def __init__(self,obj):
    self.name = obj.get( 'name', "default")
    self.format = obj.get( 'format', "slf")
    self.localFolder = obj.get( 'localFolder', "")
    self.minDensity = obj.get( 'minDensity', 10)
    self.limitFineDensity= obj.get( 'limitFineDensity', 1000)
    self.maxDensity = obj.get( 'maxDensity', 10000)
    self.shorelineGrowth = obj.get( 'shorelineGrowth', 1.2)
    self.simplification = obj.get( 'simplification', {})
    self.defaultDomain = obj.get( 'defaultDomain', {})
    self.input=input= obj.get( 'input', {})
    self.output = obj.get( 'output', {})
    self.pproj = pproj = obj.get( 'proj', "EPSG:3573")
    self.pgeo  = pgeo = obj.get( 'proj', "EPSG:4326")
    
    
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
    self.osmCoarseS=File('osmCoarseS',parent=self,fproj=self._getOSMCoarseS)
    
    
    self.osmDomain=File('osmDomain',parent=self,fproj=self._getosmDomain)
    self.osmSimplify=File('osmSimplify',parent=self,fproj=self._getosmSimplify)
    self.osmResample=File('osmResample',parent=self,fproj=self._getosmResample)  
    self.osmMesh=File('osmMesh',parent=self,ext=".msh",fproj=self._getosmMesh)
    self.osmMeshBoundaries=File('osmMeshBoundaries',parent=self,fgeo=self._getosmMeshBoundaries)
    self.osmMeshEdges=File('osmMeshEdges',parent=self,fgeo=self._getosmMeshEdges)
    self.meshmbtiles=File('mesh',parent=self,ext=".mbtiles",fgeo=self._getmbtiles)
    
    #
    # Get list of File(s) within the OSM class
    # 
    self.listFiles=[p for p in dir(OSM) if isinstance(getattr(OSM,p),property) and isinstance(getattr(self,p),File)]

 
  def _checkInput(self,input,name):
    """
    Check if input exist
    """
    def pathExist(path):
      localFolder = self.localFolder
      return (os.path.exists(path) or os.path.exists(os.path.join(localFolder,path)))
      
    if name in input and not pathExist(input[name]):raise Exception("{} does not exist".format(input[name]))
    
 
  def _getDefaultDomain(self):
    """
    Create default domain geojson
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
    """
    Create default density geojson
    """
    obj = self.defaultDomain
    minDensity = self.minDensity
    maxDensity = self.maxDensity
    shorelineGrowth = self.shorelineGrowth
    
    output=os.path.join(self.localFolder,"density.geojson")
    if os.path.exists(output):return output
    density = np.array(obj.get( 'density', [[-63.553987,44.627934,1.0,1.2]]))
    
    df=DF(density,minDensity=minDensity,maxDensity=maxDensity,minGrowth=shorelineGrowth)
    
    df.write(output)
  
    return output

  def _getdensityFineZone(self,projectedPath=None,dependencies = ['density','minDensity','limitFineDensity']):
    """
    Create density zones/area based on density object.
    This zone is used to extract fined osm coastline.
    """
    limitFineDensity = self.limitFineDensity
    return self.__getdensityZone(projectedPath,limitFineDensity)

  def _getdensityCoarseZone(self,projectedPath=None,dependencies = ['density','minDensity','maxDensity']):
    """
    Create density zones/area based on density object.
    This zone is used to extract fined osm coastline.
    """
    maxDensity = self.maxDensity
    return self.__getdensityZone(projectedPath,maxDensity)

  def __getdensityZone(self,projectedPath,density):
    """
    Create density zones/area based on density object.
    This zone is used to extract fined osm coastline.
    """
    minDensity=self.minDensity
    self.density.proj
    density=DF.read(self.density.projPath)
    
    buffers=[]
    for d in density.dp:
      maxDistance = DF.getl_D(minDensity,d[3],d[2])
      buffers.append(Point(d[:2]).buffer(maxDistance))

    buffers=cascaded_union(buffers)
    buffers.write(projectedPath)
    return buffers

  def _name(self,path):
      return os.path.splitext(os.path.basename(path))[0]

  def _getOSMFine(self,projectedPath=None,dependencies = ['densityFineZone']):
    """
    Extract fine osm coastline within the densityZone.
    """
    densityFineZone=self.densityFineZone
    densityFineZone.proj
    epsg=self.pproj.split(":")[1]
    
    osmPath = self.input['osm']
    zipname = 'water-polygons-split-4326/water_polygons.shp'
    zipPath = "\"/vsizip/" + osmPath + "/" + zipname + "\""
    
    # 1 point=>1min
    t=tqdm(total=1)
    zoneName = os.path.splitext(os.path.basename(densityFineZone.projPath))[0]
    zone = densityFineZone.projPath
    pg_sql = "\"With osm AS(SELECT ST_Transform(water_polygons.geometry,{2}) as geo FROM water_polygons) SELECT osm.geo FROM osm,'{0}'.'{1}' zone WHERE ST_Intersects(osm.geo, zone.geometry);\"".format(zone,zoneName,epsg)
    command = "ogr2ogr -skipfailures -f \"GeoJSON\" {0} -nln \"{3}\" -nlt POLYGON -dialect \"SQLITE\" -sql {2} {1}".format(projectedPath,zipPath,pg_sql,self._name(projectedPath))
    print(command)
    subprocess.call(command, shell=True)
    t.update(1)
    t.close()


  def _getOSMCoarse(self,projectedPath=None,dependencies = ['domain']):
    """
    Extract osm coastline from zip file and simplify based on extent.
    This avoids unpacking the zip file.
    """
    domain=self.domain
    osmPath=self.sosm.geoPath
    epsgp=self.pproj.split(":")[1]
    epsgg=self.pgeo.split(":")[1]
    
    self.domain.proj
    domain=domain.geoPath
    zipname = 'simplified-water-polygons-split-3857/simplified_water_polygons.shp'
    zipPath = "\"/vsizip/" + osmPath + "/" + zipname + "\""
    name = os.path.basename(domain)
    name = os.path.splitext(name)[0]

    pg_sql = "\"With one AS(SELECT ST_Buffer(ST_Transform(A.geometry,{2}),0) as geometry FROM simplified_water_polygons A,'{0}'.'{1}' B WHERE ST_Intersects(ST_Transform(A.geometry,{2}), ST_Transform(SetSRID(B.geometry,{3}),{2}))) SELECT ST_Union(one.geometry) from one WHERE one.geometry is not null;\"".format(domain,name,epsgp,epsgg)
    
    
    
    t=tqdm(total=1)
    command = "ogr2ogr -skipfailures -f \"GeoJSON\" {0} -nln \"{3}\" -dialect \"SQLITE\" -sql {1} {2}".format(projectedPath,pg_sql,zipPath,self._name(projectedPath))
    print(command)
    subprocess.call(command, shell=True)
    t.update(1)
    t.close()

  def _getOSMCoarseZone(self,projectedPath=None,dependencies = ['domain','osmCoarse']):
    """
    Extract osm coastline from zip file and simplify based on extent.
    This avoids unpacking the zip file.
    """
    osmCoarse=self.osmCoarse
    densityCoarseZone=self.densityCoarseZone
    
    osmCoarse.proj
    densityCoarseZone.proj
    
    
    pg_sql = "\"SELECT ST_Intersection(A.geometry,B.geometry) as geometry FROM '{0}'.'{1}' A,'{2}'.'{3}' B;\"".format(osmCoarse.projPath,self._name(osmCoarse.projPath),densityCoarseZone.projPath,self._name(densityCoarseZone.projPath))
    
    
    t=tqdm(total=1)
    command = "ogr2ogr -skipfailures -f \"GeoJSON\" {0} -nln \"{3}\" -dialect \"SQLITE\" -sql {1} {2}".format(projectedPath,pg_sql,osmCoarse.projPath,self._name(projectedPath))
    print(command)
    subprocess.call(command, shell=True)
    t.update(1)
    t.close()
    
    
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
    
    
    t=tqdm(total=1)
    command = "ogr2ogr -skipfailures -f \"GeoJSON\" {0} -nln \"{3}\" -dialect \"SQLITE\" -sql {1} {2}".format(projectedPath,pg_sql,osmPath,self._name(projectedPath))
    
    print(command)
    subprocess.call(command, shell=True)
    t.update(1)
    t.close()    
    

  def _getosmDomain(self,projectedPath=None,dependencies = ['domain','osmCoarseS']):
    """
    Extract osm coastline using the domain.
    It will only keep the largest Polygon.
    """
    
    domain = self.domain
    osmCoarseS = self.osmCoarseS
    
    geo = osmCoarseS.proj.geometry
    domain=domain.proj.geometry
    
    t=tqdm(total=1)
    geo=geo.largest()#.removeHoles(10000000)
    
    geo=geo.intersection(domain);t.update(1)
    geo.write(projectedPath)
    t.close()
    
    
    return geo


  def _getosmSimplify(self,projectedPath=None,dependencies = ['density','osmDomain','osmFine']):
    """
    Simplify osm shoreline based on density field
    """
    df = DF.read(self.density.projPath)
    osmDomain = self.osmDomain
    osmFine = self.osmFine
    osmCoarseZone = self.osmCoarseZone
    
    geo=osmDomain.proj.geometry
    geo=geo.dsimplify(df,fine=osmFine.proj.geometry,coarse=osmCoarseZone.proj.geometry)
    geo=geo.largest()
    geo.write(projectedPath)
    return geo
    
    
  def _getosmResample(self,projectedPath=None,dependencies = ['density','osmSimplify','minDensity','maxDensity','growth']):
    """
    Resample osm shoreline using interior nearest points and density growth field.
    """
    
    df = DF.read(self.density.projPath)
    osmSimplify = self.osmSimplify
    minDensity = self.minDensity
    maxDensity = self.maxDensity
    shorelineGrowth = self.shorelineGrowth
    
    geo=osmSimplify.proj.geometry
    
    maxDistance = DF.getl_D(minDensity,shorelineGrowth,maxDensity)
    newdensity=geo.inearest(maxDistance=maxDistance,angle=30.0)
    newdensity[:,2]=DF.getD_l(minDensity,shorelineGrowth,newdensity[:,2])
    newdensity=np.column_stack((newdensity,np.ones(len(newdensity))*shorelineGrowth))
    
    df.add(newdensity)
    geo=geo.dresample(df,progress=True)
    geo.write(projectedPath)
    return geo
  
  def _getosmMesh(self,projectedPath=None,dependencies = ['density','osmResample','minDensity','maxDensity','shorelineGrowth']):
    """
    Resample osm shoreline using interior nearest points and density growth field.
    """
    df = DF.read(self.density.projPath)
    osmResample = self.osmResample
    minDensity = self.minDensity
    maxDensity = self.maxDensity
    shorelineGrowth = self.shorelineGrowth
    
    geo=osmResample.proj.geometry
    maxDistance = DF.getl_D(minDensity,shorelineGrowth,maxDensity)
    density=geo.inearest(maxDistance=maxDistance,angle=90)
    density[:,2]=DF.getD_l(minDensity,shorelineGrowth,density[:,2],-1)
    density=np.column_stack((density,np.ones(len(density))*1.2))
    df.add(density)
    geo.msh(projectedPath,df).plot().savePlot(os.path.splitext(projectedPath)[0]+".png")
  
  def _getosmMeshBoundaries(self,geographicPath=None,dependencies=['osmMesh']):
    """
    """
    mesh=self.osmMesh.geo
    mesh.boundaries.write(geographicPath)
  def _getosmMeshEdges(self,geographicPath=None,dependencies=['osmMesh']):
    """
    """
    mesh=self.osmMesh.geo
    mesh.geoedges.write(geographicPath)
  
  def _getmbtiles(self,geographicPath,dependencies=['osmMeshBoundaries','osmMeshEdges']):
    
    edgembitle=os.path.join(os.path.dirname(geographicPath),"edges.mbtiles")
    boundarymbtile=os.path.join(os.path.dirname(geographicPath),"boundaries.mbtiles")
    if os.path.exists(edgembitle):os.remove(edgembitle)
    if os.path.exists(boundarymbtile):os.remove(boundarymbtile)
    command = "tippecanoe -z14 -o {0} -l edges {1};tippecanoe  -z14 -o {2} -l boundaries {3};tile-join {0} {2} -o {4}".format(edgembitle,self.osmMeshEdges.geoPath,boundarymbtile,self.osmMeshBoundaries.geoPath,geographicPath)
    print(command)
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
    print(command)
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
  
  
  # @staticmethod
  # def write(*args,**kwargs):
  #   """
  #   """
  #   writeGeometry(*args, **kwargs)