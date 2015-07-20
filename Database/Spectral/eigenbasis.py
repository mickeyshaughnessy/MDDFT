#!/usr/bin/env python

"""
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
  if (len(sys.argv) < 2):
    print "usage: eigenbasis.py database dmax"
    sys.exit()
  db_file = sys.argv[1]
  dmax = float(sys.argv[2])
  db = database.database()
  db.init(db_file)
  #A = db.adjacency_matrix(dmax)
  L = db.laplacian_matrix(dmax)
  print "LAPLACIAN\n",L
  [evals,evecs] = np.linalg.eigh(L)
  print "EVALS\n",evals
  print "EVECS\n",evecs
  for i in range(len(evals)):
    print evals[i]
