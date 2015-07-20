#!/usr/bin/env python

"""
run the metric search
"""
# system
import os
import sys
import numpy as np
import math
import array
# local
lib_path = os.path.abspath('..')
sys.path.append(lib_path)
import database
import util

#####################################################################
## MAIN
#####################################################################
if  __name__ == "__main__":
  if (len(sys.argv) < 3):
    print "usage: force_interpolation.py database sigma dmax"
    sys.exit()
  db_file = sys.argv[1]
  sigma = float(sys.argv[2])
  dmax = float(sys.argv[3])
  Ntrials = 1
  if (len(sys.argv) > 4):
    Ntrials = int(sys.argv[4])
  db = database.database()
  db.init(db_file)
  db.rbf.set_sigma(sigma)
  # run trials
  o = open("radial_distribution.dat","w")
  i = 0
  j = 0
  while i < Ntrials and j < db.N_neighborhoods:
    sys.stdout.flush()
    q = db.Neighborhoods[j]
    f = db.Forces[j]
    neighbors = db.Neighborhood_distances[j]
    ids = []
    error = -1.
    #print >>o,"# trial",i
    for n in range(len(neighbors)):
      [ii,d] = neighbors[n]
      if (d > dmax): break
      ids.append(ii)
      error = 0.
      ff = db.force_interpolation(q,ids)
      error = util.relative_error(ff,f)
      print >>o, d,n,error
    j += 1
    if (error > 0.):
      print ">> trial {0:3d}/{1:3d} min error {2:8g}".format(i+1,Ntrials,error)
      print >>o
      print >>o
      i += 1
      
      
