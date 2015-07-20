#!/usr/bin/env python
"""
"""
import os
import sys
import database
import cluster
import util
import numpy as np
import math


#####################################################################
## MAIN
#####################################################################
if  __name__ == "__main__":
  # parse
  if (len(sys.argv) < 2):
    print "usage: db_cluster.py database"
    sys.exit()
  dfile =      sys.argv[1]
  K     =  int(sys.argv[2])

  if not os.path.isfile(dfile):
    print dfile,"does not exist"
    sys.exit()

  stem = (dfile.split("."))[0]
  db = database.database()
  db.init(stem,"RMSD")
  A = db.adjacency_matrix()

  c = cluster.kcluster(A,K)
  c.cluster()
  c.report()
