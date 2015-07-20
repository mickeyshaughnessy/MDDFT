#!/usr/bin/env python

"""
"""
#import fileinput
import sys
import os
import numpy as np
import math
#import array
#import random
#from operator import itemgetter
lib_path = os.path.abspath('..')
sys.path.append(lib_path)
import metric_search
import util


#####################################################################
## MAIN
#####################################################################
if  __name__ == "__main__":
  s = metric_search.metric_search()
  # parse
  if (len(sys.argv) < 2):
    print "usage: database_test.py database [queries]"
    sys.exit()
  db_file = sys.argv[1]
  # init
  N_neighborhoods = s.init(db_file)

  # min distances for each neighborhood
  print "============ minimum distances ============="
  min_ds = s.Db.min_distance(s.Db.Pair_distances,s.Db.N_neighborhoods)
  util.histogram(min_ds,db_file+".min_distances")

  # histogram
  print "============ distance histogram ============="
  util.histogram(s.Db.Pair_distances,db_file+".histogram")

  # graph cluster 
  print "============ cluster ======================="
  #dcut = 1.5*min(min_ds)
  dcut = 1.5*max(min_ds)
  dcut = 3.0
  print "using d_cut=",dcut
  #dcut = 3.
  #dcut = 2.
  #dcut = 2.0
  s.Db.cluster(s.Db.N_neighborhoods,s.Db.Pair_distances,dcut)

  print "============ visualize ======================"
  s.Db.write_dot(dcut)
  
  if (len(sys.argv) < 3) : sys.exit()
  q_file = sys.argv[2]
  queries = s.Db.read_neighborhoods(q_file)
  qforces = s.Db.read_forces(q_file)

  print "============ search ======================"
  q = queries[0]
  N_refs = 1
  refs = s.initial_refs(N_refs,s.init_type)
  tol = 1.0
  [mid,mdist,N_comps,path] = s.search(q,refs,tol)
  print "using tol=",tol
  print "path : "
  for pt in path:
    print pt
## NOTE ingest test_rbf
  s.Db.write_dot_path(3*tol,path)
  
