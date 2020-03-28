# OSM coastline for gmsh
The package downloads, extract, simplify and resample the osm coastline (link)
based on density growth field for gmsh. 

## Installation

## Usage
```python
from mshapely import OSM
```

#### OSM.download(path)
```
  Download osm coastline
  
  Parameters
  ----------
  folder: path
  
  Output
  ------
  {path}/water-polygons-split-4326.zip
```
Example
```python
OSM.downloadOSM("../data")

```
#### OSM.extract(osmPath,geojsonPath,extent)
```
  Extract osm coastline from zip file
  
  Parameters
  ----------
  osmPath: path, osm zip file path
  geojsonPath: path, output geojson path
  extent:float,[xmin,ymin,xmax,ymax]
  
  Output
  ------
  geojson
```
Example

```python
OSM.extractOSM("../data/water-polygons-split-4326.zip","../data/osm.geojson",[-68,43,-62,47])
```

#### OSM(obj)
```
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
    growth:float,input
    
  Note
  ----
  Any spatial manipulation is performed on the LAEA projection (north pole).
  Results are converted back to geographic coordinates.
  
  Ouput
  -----
  osmDomain:geojson
    
  osmSimplify:geojson
    
  osmResample:geojson
```
Example 1

```python
obj = {
    "minDensity":10,
    "maxDensity":10000,
    "growth":1.2,
    "path":{
      "osm":"../data/osm.example1.geojson",
      "domain":"../data/domain.example1.geojson",
      "density":"../data/density.example1.geojson",
      "osmDomain":"../data/osmDomain.example1.geojson",
      "osmSimplify":"../data/osmSimplify.example1.geojson",
      "osmResample":"../data/osmResample.example1.geojson",
    }
  }

extent = [-68,43,-62,47]
OSM.extractOSM("../data/water-polygons-split-4326.zip","../data/osm.example1.geojson",[-68,43,-62,47])

domain=Point((-63.553987,44.627934)).buffer(1) # 1 degree  
domain.write(obj['path']['domain'])

density = np.array([
  [-63.563342,44.634637,10],
  [-63.553987,44.627934,10],
  [-63.495436,44.606528,20],
  ])
xy=density[:,[0,1]]
properties = [{"density":x} for x in density[:,2]]
OSM.write(list(map(Point,zip(xy[:,0],xy[:,1]))),obj['path']['density'],properties=properties)

OSM.osmResample
```


