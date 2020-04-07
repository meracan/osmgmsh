# OSM GMSH
Mesh generation using osm coastline for Marine Energy Resource Assessment Canada.
 
## Installation
This package was developed,tested and built using conda. OSM GMSH uses mshapely,mmesh, shapely, numpy, scipy, fiona, matplotlib,gmsh.
Only tested with python >=3.6

```bash
conda create -n osmgmsh python=3.8
conda activate osmgmsh
conda install -c meracan osmgmsh
```

For developers and debugging:
```bash
conda create -n osmgmsh python=3.8
conda activate osmgmsh
conda install -c conda-forge numpy scipy fiona shapely pyproj requests geojson tqdm matplotlib gmsh
pip install -e ./mshapely
pip install -e ./slfpy
pip install -e ./mmesh
pip install -e ./osmgmsh
```

For develops using VM:
```bash
sudo yum install mesa-libGL
```

### Usage, user guide and examples
[Docs](doc/doc_osmgmsh.ipynb)

### Testing
[Docs](test/README.md)

### License
[License](LICENSE)