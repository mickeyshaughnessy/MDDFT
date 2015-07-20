#!/usr/bin/env python

"""
run the metric search
"""
import os
import sys
import math
import array
import numpy as np
import multiprocessing as mp
#import imp
#metric_search = imp.load_source('metric_search', '../metric_search.py')
lib_path = os.path.abspath('..')
sys.path.append(lib_path)
import metric_search
import util


#####################################################################
def worker(s,q,r):
  pname = (mp.current_process().name).replace("rocess-","")
  sys_stdout = sys.stdout
  while True :
    qid = q.get()
    if qid is None:
      q.task_done()
      break
    f = open("q="+str(qid)+"_convergence.dat","w")
    sys.stdout = f
    r[qid] = s.search(s.Db.Neighborhoods[qid])
    sys.stdout = sys_stdout
    f.close()
    q.task_done()

#####################################################################
## MAIN
#####################################################################
if  __name__ == "__main__":
  sys_stdout = sys.stdout 
  maxPath = 40
  nSteps = 6
  tol = 1.e-4
  #m_type = "OGTO_1.000"
  m_type = "RMSD"
  ## parse ===========================================================
  if (len(sys.argv) < 2):
    print "usage: search_stats.py db <nref1:nref2:...> <{nqueries}x{ntrials}>"
    print "       trials: reset random refs"
    sys.exit()
  db_file = sys.argv[1]
  nrefs = [1]
  if (len(sys.argv) > 2):
    nrefs = map(int, (sys.argv[2]).split(":"))
  nqueries = 1
  N_trials = 1
  if (len(sys.argv) > 3):
    nqueries,N_trials = map(int, (sys.argv[3]).split("x"))
  print "> database:",db_file
  print "> queries: ",nqueries
  print "> trials:  ",N_trials 
  ## init ===========================================================
  print "*************************************************************"
  s = metric_search.metric_search()
  N_clusters = s.init(db_file,m_type)
  ndim = s.Db.ndim
  nprocs = mp.cpu_count()
  print "*************************************************************"
  print
  ## run tests ======================================================
  rfile = open("size="+str(N_clusters)+"_average.dat","w")
  dfile = open("size="+str(N_clusters)+"_progress.dat","w")
  for N_refs in nrefs:
    tag = "size="+str(N_clusters)+"_refs="+str(N_refs)
    print "> queries:",
    qids = s.initial_refs(nqueries,s.init_type) # refs as queries
    sumEvals = 0.
    sum2Evals = 0.
    count = 0
    Evals = []
    Steps = nSteps*[0]
    for i in range(N_trials):
      q = mp.JoinableQueue()
      for qid in qids:
        q.put(qid)
      for j in range(nprocs):
        q.put(None)
      #print ">",len(qids),"query ids:",qids
      mgr = mp.Manager()
      r = mgr.dict()
      print "> starts: ",
      refs = s.initial_refs(N_refs,s.init_type) # reset refs
      jobs = []
      for j in range(nprocs):
        jobs.append(mp.Process(target=worker,args=(s,q,r,)))
      for job in jobs:
        job.start()
      q.close()
      print "> trial {0:3d}/{1:3d} evals: ".format(i+1,N_trials),
      for job in jobs:
        job.join()
      rs = dict(r)
      count += len(rs)
      j = 0
      for qid in rs:
        result = rs[qid]
        if result[0] != qid:  print "!!! search failed",qid,"!!!"
        evals = result[2]
        Evals.append(evals)
        path = result[3]
        npath = len(path)
        if (npath>0) :
          for k in range(nSteps):
            if (k < npath) : 
              Steps[k] += path[ k][2]
              print >>dfile,path[k][0],path[k][1],N_refs
            else :           
              Steps[k] += path[-1][2]
        if (j < 10): print evals,
        j += 1
        sumEvals += evals
        sum2Evals += evals*evals
      if len(rs) > 9: print "...",
      print
    print >>dfile
    print >>dfile
    util.histogram(Evals,tag+"_histogram.dat")
    sfile = open(tag+"_steps.dat","w")
    sprev = N_refs
    for j in Steps:
      print >>sfile, j, float(j)/count,(j-sprev)
      sprev = j
    sfile.close()
    ave = sumEvals/count
    eff = ave/N_clusters
    var = 0.
    if (count > 2): var = (sum2Evals-sumEvals*ave)/(count-1)
    sd = math.sqrt(var/count)
    sdeff = sd/N_clusters
    print >>rfile,N_clusters,N_refs,ave,sd
    print "> average evaluations {0:8g} +/-{1:<8g} ".format(ave,math.sqrt(var/count))
    print "> percentage explored {0:8.2f} +/-{1:<4.2f} for size {2:<8d} ".format(eff,sdeff,N_clusters)
    print 
  rfile.close()
  dfile.close()
