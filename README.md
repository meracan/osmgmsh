# OSM GMSH
Mesh generation using osm coastline for Marine Energy Resource Assessment Canada.
 

## Installation
This library was only tested using conda.
It requires shapely,fiona,gdal, gmsh and others

```bash
conda create -n osmgmsh python=3.8
conda activate osmgmsh

conda install -c conda-forge numpy scipy fiona pyproj requests shapely gdal geojson tqdm matplotlib gmsh
# git clone https://github.com/meracan/osmgmsh.git
# pip install -e ./mshapely

# On AWS VM
sudo yum install mesa-libGL


```


### Usage
```python
import osmshapely
```
### User Guide and Examples
[Docs](doc/README.md)

###



### Testing

```bash
conda install pytest
mkdir ../data
pytest
```

For developers and debugging:
```bash
mkdir ../data
cd osmgmsh
conda activate osmgmsh
PYTHONPATH="../mshapely/:../osmgmsh/" python3 doc/doc_osm.py

```
###  


[### License](LICENSE)