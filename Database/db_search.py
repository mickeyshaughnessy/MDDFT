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
import metric_search
import util

#####################################################################
## MAIN
#####################################################################
if  __name__ == "__main__":
  s = metric_search.metric_search()
  # parse
  if (len(sys.argv) < 2):
    print "usage: run_search.py database <nqueries> <nrefs> <ntrials>"
    sys.exit()
  db_file = sys.argv[1]
  N_refs = 1
  if (len(sys.argv) > 2):
    N_refs =int(sys.argv[2])

  # init
  N_clusters = s.init(db_file)
  ndim = s.Db.ndim
  qids = s.initial_refs(1,"random")
  Queries = []
  for qid in qids:
    Queries.append(s.Db.Neighborhoods[qid])
    
  # find the right answers
  answers = []
  for q in Queries:
    [min_id,dmin,ids] = s.exhaustive_search(q)
    answers.append([min_id,dmin])
    Refs = s.initial_refs(N_refs,s.init_type)
    tol = 1.1*dmin #  set to ensure unique solution
    [mid,mdist,N_comps,path] = s.search(q,tol)
    if (mdist > tol) :
      print "!!! Metric search failed, id ",mid," not equal to ",min_id," and d=",mdist," > tol"
      print >>o, "*"
    else : 
      [exhaustive_dist,R] = s.Db.distance(s.Db.Neighborhoods[min_id],s.Db.Neighborhoods[mid])
