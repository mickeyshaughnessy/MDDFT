#!/usr/bin/env python

"""
"""
#import fileinput
import sys
import metric_search
import numpy as np
import math
import util


#####################################################################
## MAIN
#####################################################################
if  __name__ == "__main__":
  s = metric_search.metric_search()
  # parse
  if (len(sys.argv) < 2):
    print "usage: path.py database"
    sys.exit()
  db_file = sys.argv[1]
  dcut = float(sys.argv[2])
  # init
  N_neighborhoods = s.init(db_file)

  s.Db.write_dot(dcut)
  #s.Db.write_dot_path(dcut)
  
