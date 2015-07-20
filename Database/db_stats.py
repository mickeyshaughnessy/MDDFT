#!/usr/bin/env python
"""
"""
import os
import sys
import database
import util
import numpy as np
import math


#####################################################################
## MAIN
#####################################################################
if  __name__ == "__main__":
  # parse
  if (len(sys.argv) < 2):
    print "usage: db_stats.py database"
    sys.exit()
  dfile =      sys.argv[1]
  nmax =   int(sys.argv[2])
  dmax = float(sys.argv[3])
  drange = (0,dmax)

  print "> distance file:",dfile," samples: ",nmax," distance range: ",drange

  if not os.path.isfile(dfile):
    stem = (dfile.split("."))[0]
    db = database.database()
    db.init(stem)
  if not os.path.isfile(dfile):
    print dfile,"does not exist"
    sys.exit()

  # init
  distances = util.read_distances(dfile,drange,nmax)
  util.save_distances(distances,dfile+".npy")

  # histogram
  print "============ distance histogram ============="
  util.histogram(distances,dfile+".dhist",drange)
  sys.exit()

  # min distances for each neighborhood
  print "============ minimum distances ============="
  min_ds = db.min_distances()
  util.histogram(min_ds,db_file+".min_dhist",drange)

  print "============ mean distances ============="
  mean_ds = db.mean_distances()
  util.histogram(mean_ds,db_file+".mean_dhist",drange)

  print "============ max distances ============="
  max_ds = db.max_distances()
  util.histogram(max_ds,db_file+".max_dhist",drange)

