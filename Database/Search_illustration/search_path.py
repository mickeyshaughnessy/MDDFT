#!/usr/bin/env python

"""
run the metric search
"""
import os
import sys
lib_path = os.path.abspath('..')
sys.path.append(lib_path)
import metric_search

#####################################################################
## MAIN
#####################################################################
if  __name__ == "__main__":
  s = metric_search.metric_search()
  # parse
  if (len(sys.argv) < 3):
    print "usage: metric_search.py database queries <nrefs> <ntrials>"
    sys.exit()
  db_file = sys.argv[1]
  q_file  = sys.argv[2]
  dmax = 1.0
  if (len(sys.argv) > 3):
    dmax  = float(sys.argv[3])
  nrefs = 1
  if (len(sys.argv) > 4):
    nrefs = int(sys.argv[4])
  # init
  N_clusters = s.init(db_file)
  Queries = s.Db.read_neighborhoods(q_file)
  q = Queries[0]
  path = []
  # search
  Refs = s.initial_refs(nrefs,s.init_type)
  [min_id,min_dist,N_comps,path] = s.search(q,Refs,1e-8)
  # output
  #nedges = s.Db.prune(dmax)
  #print ">> pruned to",dmax," nedges",nedges
  #s.Db.write_path(dmax,path,"path")
  s.Db.write_neighborhood_path(dmax,path,"path")
  s.Db.histogram()
