import numpy as np
import json
from osmgmsh import OSM
import time
def run_example(path):
  start=time.time()
  with open(path) as f:
    obj = json.load(f)
    osm=OSM(obj)
    # osm.osmMesh.geo
    osm.compute()
  print(time.time()-start)
    
    
if __name__ == "__main__":
  run_example('example/example2.json')