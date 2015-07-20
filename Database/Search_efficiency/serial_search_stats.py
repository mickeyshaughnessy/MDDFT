#!/usr/bin/env python

"""
run the metric search
"""
import os
import sys
import math
import array
import numpy as np
#import imp
#metric_search = imp.load_source('metric_search', '../metric_search.py')
lib_path = os.path.abspath('..')
sys.path.append(lib_path)
import metric_search
import util

#####################################################################
## MAIN
#####################################################################
if  __name__ == "__main__":
  sys_stdout = sys.stdout 
  maxPath = 40
  #m_type = "OGTO_1.000"
  m_type = "RMSD"
  ## parse ===========================================================
  if (len(sys.argv) < 2):
    print "usage: run_search.py database <nqueries> <nref1:nref2:...> <ntrials>"
    print "       trials: reset random refs"
    sys.exit()
  db_file = sys.argv[1]
  nrefs = [1]
  if (len(sys.argv) > 2):
    nrefs = map(int, (sys.argv[2]).split(":"))
  nqueries = 1
  N_trials = 1
  if (len(sys.argv) > 2):
    nqueries,N_trials = map(int, (sys.argv[3]).split("x"))
  ## init ===========================================================
  s = metric_search.metric_search()
  N_clusters = s.init(db_file,m_type)
  ndim = s.Db.ndim
  qids = s.initial_refs(nqueries,s.init_type)
  print "*************************************************************"
  Queries = []
  for qid in qids:
    Queries.append(s.Db.Neighborhoods[qid])
  print "query ids:",qids
    
  ## run tests ======================================================
  r = open("size="+str(N_clusters)+"_average.dat","w")
  # find the right answers
  answers = []
  for q in Queries:
    [min_id,dmin,ids] = s.exhaustive_search(q)
    answers.append([min_id,dmin])
  for nr in nrefs:
    N_refs = nr
    o = open("size="+str(N_clusters)+"_refs="+str(N_refs)+"_results.dat","w")
    f = open("size="+str(N_clusters)+"_refs="+str(N_refs)+"_convergence.dat","w")
    SumEvals = 0.
    Sum2Evals = 0.
    Count = 0
    path = []
    for i in range(N_trials):
      print "==== trial:",i,"============================================="
      Refs = s.initial_refs(N_refs,s.init_type)
      n = 1
      for q in Queries: ## NOTE parallelize
        print "---- trial:", i," query:",n, " ----------------------------"
        [min_id,dmin] = answers[n-1]
        print >>o,"# ",i,n
        tol = 1.1*dmin #  set to ensure unique solution
        sys.stdout = f
        [mid,mdist,N_comps,path] = s.search(q,Refs,tol)
        print >>f
        print >>f
        sys.stdout = sys_stdout
        if (N_clusters < maxPath): s.Db.write_path(2.,path,"path"+str(n))
        if (mdist > tol) :
          print "!!! Metric search failed, id ",mid," not equal to ",min_id," and d=",mdist," > tol"
          print >>o, "*"
        else : 
          # triangle check
          [exhaustive_dist,R] = s.Db.distance(s.Db.Neighborhoods[min_id],s.Db.Neighborhoods[mid])
          if (mid != min_id and not s.triangle_check(mid,min_id,q)):
            print "!!! Triangle check failed for id ",mid,"and id",min_id
          # log
          SumEvals += N_comps
          Sum2Evals += N_comps*N_comps
          Count +=1
          print >>o,N_comps,N_clusters
        n += 1
    print "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"
    Ave = SumEvals/Count
    Var = 0.
    if (Count > 2):
      Var = (Sum2Evals-SumEvals*Ave)/(Count-1)
      print "Average number of evaluations {0:8g} +/- {1:8g} ".format(Ave,math.sqrt(Var/Count))
    else:
      print "Average number of evaluations {0:8g}".format(Ave)
    eff = Ave/N_clusters
    print "Percentage explored {0:8g} for size {1:<8d} ".format(eff,N_clusters)
    print >>o,"# ",N_refs,N_clusters,"  average ", Ave,"/",N_clusters
    print >>r,N_clusters,N_refs,Ave,math.sqrt(Var/Count)
