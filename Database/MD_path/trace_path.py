#!/usr/bin/env python

"""
run the metric search
"""
import os
import sys
import math
import array
import numpy as np
lib_path = os.path.abspath('..')
sys.path.append(lib_path)
import database
import util

#####################################################################
## MAIN
#####################################################################
if  __name__ == "__main__":
  db = database.database()
  # parse
  if (len(sys.argv) < 4):
    print "usage: trace_path.py database natoms nsteps"
    sys.exit()
  db_file = sys.argv[1]
  stem = db_file.rstrip(".db")
  natoms = int(sys.argv[2])
  nsteps = int(sys.argv[3])
  ntotal = natoms*nsteps
  npaths = 1
  if (len(sys.argv) > 4):
    npaths = int(sys.argv[4])
  dmax = 1e10
  if (len(sys.argv) > 5):
    dmax = float(sys.argv[5])
   
  paths = []
  db.init(db_file,False)
  for i in range(npaths):
    path = db.compute_path_distances(range(i,ntotal,natoms))
    paths.append(path)

  f = open(stem+".paths","w")
  for i in range(nsteps-1):
    print >>f, i,
    for j in range(npaths):
      d = paths[j][i]
      if d < dmax: print >>f," {0:12g}".format(d),
      else       : print >>f," *",
    print >>f
  
  util.histogram(util.flatten(paths),stem+".phist",[0,0.05])
