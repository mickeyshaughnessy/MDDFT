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
  # parse
  if (len(sys.argv) < 2):
    print "usage: database_test.py database [nbins]"
    sys.exit()
  db_file = (sys.argv[1]).rstrip(".db")
  nb  = int(sys.argv[2])
  stem = (db_file.split("/"))[-1]
  for m in ["RMSD","OGTO"]: 
    print "================= METRIC:",m,"==================="
    tag = stem+"_"+m
    s = metric_search.metric_search()
    N_neighborhoods = s.init(db_file,m)
    # min distances for each neighborhood
    #print "============ minimum distances ============="
    min_ds = s.Db.min_distances()
    util.histogram(min_ds,tag+".min_distances",nbins=nb)
    # histogram
    #print "============ distance histogram ============="
    util.histogram(s.Db.Pair_distances,tag+".histogram",nbins=nb)
    print
