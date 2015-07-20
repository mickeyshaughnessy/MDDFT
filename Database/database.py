#!/usr/bin/env python

"""
database
TODO: 
* precompute internal representation of neighborhoods
* read/write distances - limited bandwidth
* add to dist file:: save R_AB as p & theta - append to dist file and also compute force diff?
"""
# system
import fileinput
import sys
import os
import math
import array
import re
import random
from operator import itemgetter
#import gzip
#import subprocess
import numpy as np
hasSparse = True
try:
  import scipy.sparse
except ImportError:
  print "!!! need scipy.sparse for adjacency matrix !!!"
  hasSparse = False
import multiprocessing as mp
import Queue
#import pydot
# local 
import util
import rmsd
import ogto
import soap
import interpolation
import cluster

# see: https://www.sqlite.org/limits.html
SQLITE_MAX_COLUMN=2000
SQLITE_MAX_VARIABLE_NUMBER=999
SQLITE_MAX_COLS=min(SQLITE_MAX_COLUMN,SQLITE_MAX_VARIABLE_NUMBER)

#####################################################
class database:
#####################################################
  parallel = False
  cfile = None  # configurations
  dfile = None  # distances
  nfile = None  # nieghborhoods
  metric =  None
  rbf = None
  N_neighborhoods = 0
  N_pairs = 0
  N_distances = 0
  ndim = 3
  N_references = 0
  References = []
  Neighborhoods = [] # id --> neighborhood / cluster
  Forces = []
  maxConnections = 200
  maxConnectionsBuffer = 40
  Pairs = []     # (i,j) pair 
  PairMap = {}   # k --> (i,j)
  Pair_distances = [] # pair --> distance
  Min_distances = [] #  pair --> distance
  Max_min_distance = 0
  dmax = 0.
  dmin = 0.
  Neighborhood_distances = [] # id --> [id,distance] sorted edges of id
  Adjacency_matrix = []
  Laplacian_matrix = []
  Pair_force_differences = [] # pair --> f_a - R f_b
  distance_bandwidth = 0
  Group_ids = []
###################################################################
  def __init__(self,metric_type=None,rbf_type=None): 
    if metric_type != None : self.set_metric(metric_type)
    if rbf_type    != None : self.set_rbf(rbf_type)

#==================================================================
#> BASIC FUNCTIONS
#==================================================================
  def set_metric(self,metric_type):
    metric_type = metric_type.upper()
    if   (metric_type[:4] == 'RMSD'): 
      self.metric = rmsd.rmsd()
    elif (metric_type[:4] == 'OGTO'):
      if len(metric_type) > 4 :
        sigma = float(metric_type[5:])
        self.metric = ogto.ogto(sigma)
      else :
        self.metric = ogto.ogto()
    elif (metric_type[:4] == 'SOAP'):
      self.metric = soap.soap()
    else:
      print "!!! unknown metric type:",metric_type,"!!!"
      sys.exit()
    print "> metric type:",self.metric.tag()
    sys.stdout.flush()

###################################################################
  def set_rbf(self,rbf_type):
    if   (rbf_type == 'shepard'): 
      self.rbf = interpolation.inverse_quadratic()
    elif (rbf_type == 'exponential'): 
      self.rbf = interpolation.inverse_exponential()
    elif (rbf_type == 'gaussian'): 
      self.rbf = interpolation.gaussian()
    else:
      print "ERROR: unknown rbf type:",rbf_type
    print "> rbf    type:",rbf_type
    sys.stdout.flush()

###################################################################
  def init(self,stem,metric_type=None,n_refs=0):
    self.set_metric(metric_type)
    self.stem = stem.rstrip(".db")
    self.cfile = self.stem+".db"
    self.dfile = self.stem+"."+self.metric.tag()
    self.nfile = self.stem+"."+self.metric.tag()+".dsort"
    if not os.path.isfile(self.cfile):
      print "!!! no neighborhood database file",self.cfile,"exists !!!"
      sys.exit()
    else :
      print "> found neighborhoods in",self.cfile
    self.Neighborhoods = self.read_neighborhoods(self.cfile)
    self.Forces = self.read_forces(self.cfile)
    self.N_neighborhoods = len(self.Neighborhoods)
    if (self.N_neighborhoods < 1) :
      print "!!! no neighborhoods !!!"
      sys.exit()
    if self.metric != None:
      if os.path.isfile(self.dfile):
        print "> found distances in",self.dfile
        self.read_distances(self.dfile)
      else :
        self.compute_pair_distances()
    self.distance_bandwidth = self.N_neighborhoods-1
    self.Group_ids = self.N_neighborhoods*[0]
    if n_refs > 0:
      self.N_references = n_refs
      #self.References = self.most_distant(n_refs)
      self.References = self.cluster(n_refs) 
    sys.stdout.flush()
    return self.N_neighborhoods

#####################################################
  def distance(self,a,b):
    return self.metric.distance(a,b)

###################################################################
  def pair(self,k): # flat to pair
    return self.Pairs[k]
  
###################################################################
  def index(self,i,j): # pair to flat
    try :
      if (i<j) : return self.PairMap[(i,j)]
      else     : return self.PairMap[(j,i)]
    except KeyError :
      print "!!! corrupted pair distances",i,",",j," !!!"
      sys.exit()

###################################################################
  def next_distance(self,i,j):
    ds = Neighborhood_distances[i]
    m = 0
    for k,d in ds:
      m += 1
      if (m==len(ds)):
        return d
      elif (k==j): 
        p = ds[m]
        return p[1]
    return -1. # never 

###################################################################
  def write_database(self,filename):
    f = open(filename,'w')
    for i in range(self.N_neighborhoods):
      C = self.Neighborhoods[i]
      F = self.Forces[i]
      print >>f,"{0:2d} {1:10d} {2:18g} {3:18g} {4:18g} ".format(1,11,F[0],F[1],F[2]),
      for x in C:
        print >>f,"{0:18g} {1:18g} {2:18g} ".format(x[0],x[1],x[2]),
      print >>f
    f.close()
    print "> wrote ", self.N_neighborhoods," neighborhoods"

###################################################################
  def read_neighborhoods(self,filename):
    neighborhoods = []
    f = open(filename,'r')
    nb = 0
    ndim = self.ndim
    for line in f:
      cols = line.split()
      nb = (len(cols)-5)/ndim
      c = np.array(map(float,cols[5:]))
      neighborhoods.append(c.reshape(nb,ndim))
    f.close()
    print "> read", len(neighborhoods),str(nb)+"x"+str(ndim),"neighborhoods"
    sys.stdout.flush()
    return neighborhoods

###################################################################
  def read_forces(self,filename):
    forces = []
    f = open(filename,'r')
    for line in f:
      cols = line.split()
      c = np.array(map(float,cols[2:5]))
      forces.append(c)
    f.close()
    print "> read", len(forces),"forces" 
    sys.stdout.flush()
    return forces 

###################################################################
  def read_dump_forces(self,filename):
    f = open(filename,'r')
    natoms_pattern = re.compile("^ITEM: NUMBER ")
    line = f.readline()
    n = 0
    while line:
      if natoms_pattern.search(line):
        line = f.readline()
        n = int((line.split())[0])
        break
      line = f.readline()
    print "> dump file",filename,"has",n,"atoms"
    item_pattern = re.compile("^ITEM: ATOMS ")
    forces = []
    while line:
      if item_pattern.search(line):
        for i in range(n):
          line = f.readline()
          cols = line.split()
          c = np.array(map(float,cols[2:5]))
          forces.append(c)
      line = f.readline()
    f.close()
    print "> read", len(forces),"forces" 
    sys.stdout.flush()
    return forces 

###################################################################
  def read_distances(self,filename):
    npairs = (self.N_neighborhoods)*(self.N_neighborhoods-1)/2
    freq = max(int(npairs/10),1)
    self.Pairs = [] # clear
    self.PairMap = {} # clear
    self.Pair_distances = [] # clear
    if os.path.isfile(filename+".gz"): 
      import gzip
      f = gzip.open(filename)
    else:
      f = open(filename)
    dmin = util._big
    dmax = -1.
    k = 0
    print "> reading ", 
    sys.stdout.flush()
    for line in f:
      cols = line.split()
      i = int(cols[0])
      j = int(cols[1])
      self.Pairs.append([i,j])
      self.PairMap[(i,j)] = k
      d = float(cols[2])
      if (d < dmin): dmin = d
      if (d > dmax): dmax = d
      self.Pair_distances.append(d)
      k += 1
      if (k % freq == 0): 
        sys.stdout.write("*")
        sys.stdout.flush()
    f.close()
    self.N_distances = len(self.Pair_distances) # N_pairs
    self.N_pairs = self.N_neighborhoods*(self.N_neighborhoods-1)/2
    print " ",self.N_distances,"/",self.N_pairs,"distances",
    print "(%.6f:%.6f)"%(dmin,dmax)
    self.create_neighborhood_distances()
    return self.N_distances 

###################################################################
  def write_distances(self,filename):
    if os.path.isfile(filename):
      print "!!! distance file ",filename," exists !!!"
      sys.exit()
    f = open(filename,'w')
    k = 0
    for d in self.Pair_distances:
      p = self.pair(k)
      print >>f, p[0],p[1],d
      k += 1
    f.close()
    print "> wrote ", len(self.Pair_distances)," distances to ",filename

###################################################################
  def set_parallel(self,p=True):
    self.parallel = p

###################################################################
  def compute_pair_distances(self):
    if self.parallel : self.parallel_compute_pair_distances()
    else             :   self.serial_compute_pair_distances()

###################################################################
  def distance_worker(self,q,r): 
    pname = (mp.current_process().name).replace("rocess-","")
    while True : 
      p = q.get()
      if p is None:
        q.task_done()
        break
      [i,j] = p
      idx   = self.index(i,j)
      [d,R] = self.distance(self.Neighborhoods[i],self.Neighborhoods[j])
      if (idx % 1000==0): print "{0:8d}/{1:8d}:{2:5s} [{3:5d},{4:5d}] ==> {5:7.4f}".format(idx,self.N_pairs,pname,p[0],p[1],d)
      q.task_done()
      r[idx] = (d,R)

###################################################################
  def parallel_compute_pair_distances(self):
    k = 0
    self.Pairs = []
    q = mp.JoinableQueue()
    #q = mp.Queue()
    for i in range(0,self.N_neighborhoods):
      for j in range(i+1,self.N_neighborhoods):
        p = [i,j]
        self.Pairs.append(p)
        self.PairMap[(i,j)] = k
        k += 1
        q.put(p)
    self.N_pairs = len(self.Pairs)
    nprocs = mp.cpu_count()
    for i in range(nprocs):
      q.put(None)
    #r = mp.Queue()
    print "> computing", self.N_pairs,"distances with",nprocs,"processes"
    mgr = mp.Manager()
    r = mgr.dict()
    jobs = []
    for i in range(nprocs):
      jobs.append(mp.Process(target=self.distance_worker,args=(q,r,)))
    for job in jobs:
      job.start()
    q.close()
    print "> collecting ",
    for job in jobs:
      job.join()
    self.Pair_distances = self.N_pairs*[0.]
    self.Pair_force_differences = self.N_pairs*[0.]
    self.dmin = util._big
    self.dmax = -1.
    print len(r)," distances ",
    sys.stdout.flush()
    ds = dict(r) # necessary since the mgr version is oddly behaved
    f = open(self.dfile,"w")
    for idx in ds:
      i,j = self.pair(idx)
      d,R = ds[idx]
      self.Pair_distances[idx] = d
      fi = self.Forces[i]
      fj = self.Forces[j]
      df = util.error(fi,util.rotate(R,fj))
      self.Pair_force_differences[self.index(i,j)] = df
      print >>f, i,j,d,df,R[0][0],R[0][1],R[0][2], \
                          R[1][0],R[1][1],R[1][2], \
                          R[2][0],R[2][1],R[2][2]
      if (self.dmin > d) : self.dmin = d
      if (self.dmax < d) : self.dmax = d
    print " %.6f:%.6f"%(self.dmin,self.dmax),"done"
    sys.stdout.flush()
    f.close()
    self.create_neighborhood_distances()
    
###################################################################
  def serial_compute_pair_distances(self):
    self.N_neighborhoods = len(self.Neighborhoods)
    self.N_pairs = self.N_neighborhoods*(self.N_neighborhoods-1)/2
    print "> computing", self.N_pairs,"distances ",
    every = max(4,int(self.N_pairs/10))
    self.Pair_distances = self.N_pairs*[0.]
    self.Pair_force_differences = self.N_pairs*[0.]
    self.dmin = util._big
    self.dmax = -1.
    k = 0
    self.Pairs = []
    f = open(self.dfile,"w")
    for i in range(0,self.N_neighborhoods):
      for j in range(i+1,self.N_neighborhoods):
        self.Pairs.append([i,j])
        self.PairMap[(i,j)] = k
        [d,R] = self.distance(self.Neighborhoods[i],self.Neighborhoods[j])
        self.Pair_distances[self.index(i,j)] = d
        fi = self.Forces[i]
        fj = self.Forces[j]
        df = util.error(fi,util.rotate(R,fj))
        self.Pair_force_differences[self.index(i,j)] = df
        print >>f, i,j,d,df,R[0][0],R[0][1],R[0][2], \
                            R[1][0],R[1][1],R[1][2], \
                            R[2][0],R[2][1],R[2][2]
        if (self.dmin > d) : self.dmin = d
        if (self.dmax < d) : self.dmax = d
        k += 1
        if (k%every==0) : sys.stdout.write("*")
        sys.stdout.flush()
    print " %.6f:%.6f"%(self.dmin,self.dmax),
    print " done"
    self.create_neighborhood_distances()

###################################################################
  def compute_path_distances(self,path):
    print "> computing distances along path "
    dists = []
    p0 = path[0]
    for p in path[1:]:
      [d,R] = self.distance(self.Neighborhoods[p0],self.Neighborhoods[p])
      print "...{0:8d}->{1:<8d}: {2:8g}".format(p0,p,d)
      dists.append(d)
      p0 = p
    print " done"
    return dists

###################################################################
## NOTE make PARALLEL & and write to file? 
###################################################################
  def create_neighborhood_distances(self):
    self.Min_distances = self.N_neighborhoods*[util._big] 
    f = open(self.nfile,"w")
    for i in range(0,self.N_neighborhoods):
      idds = []
      ds = [] 
      for j in range(0,self.N_neighborhoods):
        if (i==j) : continue
        d = self.Pair_distances[self.index(i,j)]
        idds.append([j,d])
        ds.append(d)
        if (self.Min_distances[i] > d) : self.Min_distances[i] = d
        if (self.Min_distances[j] > d) : self.Min_distances[j] = d
      idds.sort(key=itemgetter(1))
      self.Neighborhood_distances.append(idds)
      print >>f, i,
      for p in idds:
        print >>f, p[0],p[1],
      print >>f
    self.Max_min_distance = max(self.Min_distances)
    print "> max (min_i distance_ij) = ", self.Max_min_distance
    f.close()

###################################################################
  def force_interpolation(self,Q,ids):
    n = len(ids)
    D = np.zeros((n,n))
    for i in range(len(ids)):
      for j in range(i+1,len(ids)):
        idx = self.index(ids[i],ids[j])
        D[i,j] = D[j,i] = self.Pair_distances[idx] 
    fs = []
    ds = []
    for i in ids:
      [d,R] = self.distance(self.Neighborhoods[i],Q)
      ds.append(d)
      fi = self.Forces[i]
      Rfi = util.rotate(R,fi)
      fs.append(Rfi)
    w = self.rbf.solve(D,fs)
    return self.rbf.interpolate(ds)

#==================================================================
#> GRAPH FUNCTIONS
#==================================================================

###################################################################
# adjacency matrix for a single cluster
###################################################################
  def cluster_adjacency_matrix(self,C,dmax=1.e10):
    n = len(C)
    a= scipy.sparse.lil_matrix((n,n))
    for i in range(n):
      for j in range(i+1,n):
        d = util.d(C[i],C[j])
        if (d < dmax): a[j,i] = a[i,j] = d
    A= a.tocsr()
    return A.todense()

###################################################################
# eigen-matching
###################################################################
  def trial_permutation(self,C1,C2,tol=1.e-4):
    n = len(C1)
    A1 = self.cluster_adjacency_matrix(C1)
    lambda1,E1 = np.linalg.eig(A1)
    print "++++++++++++++ 1 ++++++++++++"
    print C1
    print A1
    print lambda1
    print E1
    A2 = self.cluster_adjacency_matrix(C2)
    lambda2,E2 = np.linalg.eig(A2)
    print "++++++++++++++ 2 ++++++++++++"
    print C2
    print A2
    print lambda2
    print E2
    i = 0
    match = False
    for j in range(n):
      if (math.fabs(lambda1[j]-lambda2[j]) < tol): 
        match = True
        break
      i += 1
    p = []
    for j in range(n): p.append(j)
    if (match):
      e1 = map(math.fabs,np.array(E1[i])[0].tolist())
      z1 = (zip(p,e1))
      print z1[0]
      print z1[0][1]
      s1=sorted(z1,key=lambda a: a[1]) 
      print s1
      e2 = map(math.fabs,np.array(E2[i])[0].tolist())
      z2 = (zip(p,e2))
      s2=sorted(z2,key=lambda a: a[1]) 
      print s2
      for j in range(n):
        p[s1[j][0]]= s2[j][0]
        p[j]= s2[s1[j][0]][0]
    print "permutation", p
    return p

###################################################################
# adjacency matrix for database
###################################################################
  def adjacency_matrix(self,dmax=1.e20):
    amat = scipy.sparse.lil_matrix((self.N_neighborhoods,self.N_neighborhoods))
    for i in range(0,self.N_neighborhoods):
      for j in range(0,self.N_neighborhoods):
        if (i==j) : continue
        d = self.Pair_distances[self.index(i,j)]
        if (d < dmax): amat[i,j] = d
    self.Adjacency_matrix = amat.tocsr()
    return self.Adjacency_matrix.todense()

###################################################################
#  L = A - D, D_ij = \sum_k d_ik delta_ij
###################################################################
  def laplacian_matrix(self,dmax):
    amat = scipy.sparse.lil_matrix((self.N_neighborhoods,self.N_neighborhoods))
    for i in range(0,self.N_neighborhoods):
      dsum = 0.
      for j in range(0,self.N_neighborhoods):
        if (i!=j) : 
          d = self.Pair_distances[self.index(i,j)]
          if (d < dmax): 
            amat[i,j] = -d
            dsum += d
      amat[i,i] = dsum
    self.Laplacian_matrix = amat.tocsr()
    return self.Laplacian_matrix.todense()

##################################################################
# cluster via k-medoids
##################################################################
  def cluster(self,K):
    A = self.adjacency_matrix()
    c = cluster.kcluster(A,K)
    centers,self.Group_ids = c.cluster()
    max_cluster_size = c.maxsize
    self.distance_bandwidth = len(centers)+max_cluster_size-1
    return centers

##################################################################
# clip connections outside clusters
##################################################################
# def clip(self):

##################################################################
# cluster based connections shorter than dcut
##################################################################
  def clump(self,n,ds,dcut,name=None):
    keys = []
    for i in range (n):
      keys.append(i)
    changed = True
    p = 1
    while (changed):
      print "pass ",p," disconnected regions: ", len(set(keys))
      #print keys
      #if (p > 1) : histogram(keys)
      p += 1
      changed = False
      for i in range(n):
        for j in range(i+1,n):
          d = ds[self.index(i,j)]
          if ( d < dcut ):
            # take lowest equivalent key          
            if   ( keys[i] > keys[j]):
              keys[i] = keys[j]
              changed = True
            elif ( keys[j] > keys[i]):
              keys[j] = keys[i]
              changed = True
    print "database size",n,"has",len(set(keys)),"disjoint regions, using cutoff", dcut
    self.compact(keys)

##################################################################
#
##################################################################
  def compact(self,keys):
    n = len(keys)
    # compact keys
    rset = set(keys)
    nset = len(rset)
    counts = nset*[0]
    k = 0
    for i in range(n):
      if (i in rset):
        #print "reassigning ",i," to ",k
        for j in range(n):
          if (keys[j] == i):
            counts[k] += 1
            keys[j] = k
        k += 1
    cmax = max(counts)
    for i in range(nset):
      s = int(10*counts[i]/cmax)
      print '%6d:' % (counts[i]), s*"*"

###################################################################
# pick n neighborhoods with largest neighborhood-neighborhood distances
###################################################################
  def most_distant(self,n):
    sorted_indices= list(np.argsort(self.Pair_distances))
    sorted_indices.reverse()
    c = set([])
    for i in sorted_indices:
      p = self.Pairs[i]
      c.add(p[0])
      if (len(c) == n): break
      c.add(p[1])
      if (len(c) == n): break
    return list(c)

###################################################################
# return a set of n indices
###################################################################
  def winnow(self,n):
    c = set(range(N_neighborhoods))
    sorted_indices= list(np.argsort(Pair_distances))
    for i in sorted_indices:
      p = pair(i)
      j = p[0]
      k = p[1]
      if   (j not in c) : c.discard(k)
      elif (k not in c) : c.discard(j)
      else :
        if (next_distance(j,k) < next_distance(k,j)) : c.discard(j)
        else                                         : c.discard(k)
      if (len(c) == n): break
    return list(c)

##################################################################
# prune all connections longer than dmax
###################################################################
  def prune(self,dmax):
    pairs = []
    pairMap = {}
    pair_distances = [] # list, pair --> distance
    ndist = 0
    for k in range(self.N_distances):
      d = self.Pair_distances[k]
      [i,j] = self.Pairs[k]
      if (d < dmax):
        pairs.append([i,j])
        pairMap[(i,j)] = ndist
        pair_distances.append(d)
        ndist += 1
    # copy temp to class
    self.N_distances = ndist
    self.Pairs = pairs
    self.PairMap = pairMap
    self.Pair_distances = pair_distances
    return ndist

#==================================================================
#> CHARACTERIZATION FUNCTIONS
#==================================================================

###################################################################
  def histogram(self):
    print "> histogram of distances"
    util.histogram(self.Pair_distances)
    print "> histogram of min distances"
    util.histogram(self.Min_distances)


###################################################################
  def write_graph(self,dmax,file="db"):
    f = open(file+".gv","w")
    print >>f,"digraph {"
    print >>f,"ratio=1.0"
    print >>f,"  node [shape=circle,color=skyblue,style=filled]"
    for i in range(0,self.N_neighborhoods):
      print >>f, "  C"+str(i+1),";";
    print >>f,"subgraph Dist { edge [dir=none,len=1]"
## NOTE sort on distances
    for i in range(0,self.N_neighborhoods):
      for j in range(i+1,self.N_neighborhoods):
        d = self.Pair_distances[self.index(i,j)]
        if (d < dmax): 
          l = "%.2f" % (d/dmax)
          print >>f, "C"+str(i+1)," -> C"+str(j+1) + " [label = " +str(l)+" ]"
    print >>f,"  }"
    print >>f,"}"
    f.close()
    self.convert_gv(file)

##################################################################
  def write_path(self,dmax,path,file="path"):
    N = self.N_neighborhoods
    M = self.N_distances
    f = open(file+".gv","w")
    print >>f,"digraph {"
    print >>f," ratio=1.0"
    print >>f," node [shape=circle,color=red,style=filled]"
    for pt in path:
      id = pt[0]
      print >>f, "  C"+str(id),";";
    print >>f," node [shape=circle,color=skyblue,style=filled]"
    for i in range(0,N):
      print >>f, "  C"+str(i+1),";";
    print >>f," subgraph Dist { edge [dir=none]"
    for k in range(M):
      d = self.Pair_distances[k]
      if (d < dmax): 
        [i,j] = self.Pairs[k]
        l = "%.2f" % (d)
        #l = "%.2f" % (d/dmax)
        print >>f, "  C"+str(i+1)," -> C"+str(j+1) + " [label = " +str(l)+" ]"
    print >>f,"  }"
    print >>f,' subgraph Path { edge [style="bold",color="red"] '
    pt = path[0]
    lastid = pt[0]
    lastd  = pt[1]
    for pt in path[1:] :
      id = pt[0]
      d  = pt[1]
      inc = lastd-d
      l = "%.2f" % (inc/dmax)
      print >>f, "  C"+str(lastid)," -> C"+str(id) + " [label = " +str(l)+" ]"
      lastid = id
      lastd = d
    # add last to Q
    print >>f," }"
    print >>f,"}"
    f.close()
    self.convert_gv(file)

##################################################################
  def write_neighborhood_path(self,dmax,path,file="path"):
    print "> writing path"
    N = self.N_neighborhoods
    M = self.N_distances
    f = open(file+".gv","w")
    print >>f,"digraph {"
    print >>f," ratio=1.0"
    print >>f," node [shape=circle,color=red,style=filled]"
    for pt in path:
      i = pt[0]
      print >>f, "  C"+str(i+1),";";
    print >>f," node [shape=circle,color=skyblue,style=filled]"
    points = set()
    edges = {}
    dlast = dmax
    for pt in path:
      i  = pt[0]
      points.add(i)
      for c in self.Neighborhood_distances[i]:
        [j,d] = c
        #if d < dmax: 
        if d < 2*dlast: 
          points.add(j)
          edges[(i,j)] = d
          print >>f, "  C"+str(j+1),";";
      print "...",i,"nodes:",len(points),"edges:",len(edges)
      dlast = pt[1]
    print >>f," subgraph Dist { edge [dir=none]"
    for e in edges:
      (i,j) = e
      d = edges[e]
      l = "%.2f" % (d)
      print >>f, "  C"+str(i+1)," -> C"+str(j+1) + " [label = " +str(l)+" ]"
    print >>f,"  }"
    print >>f,' subgraph Path { edge [style="bold",color="red"] '
    pt = path[0]
    lastid = pt[0]
    lastd  = pt[1]
    for pt in path[1:] :
      id = pt[0]
      d  = pt[1]
      inc = lastd-d
      l = "%.2f" % (inc)
      print >>f, "  C"+str(lastid+1)," -> C"+str(id+1) + " [label = " +str(l)+" ]"
      lastid = id
      lastd = d
    # add last to Q
    print >>f," }"
    print >>f,"}"
    f.close()
    self.convert_gv(file)

##################################################################
  def convert_gv(self,file):
    # convert
    c = " -Tpdf "+file+".gv -o "+file+".pdf"
    if (self.N_neighborhoods > 0):
      os.system("dot "+c)
    else:
      os.system("sfdp -Goverlap=prism "+c)
      #os.system("sfdp -Goverlap=false "+c)

##################################################################
  def min_distances(self):
    # assuming ds in a flattened trinagular matrix i.e. complete
    nrows = self.N_neighborhoods
    mins = nrows*[util._big]
    k = 0
    for i in range(nrows):
      for j in range(i+1,nrows):
        d = self.Pair_distances[k]
        if (d > util._small)  :
          if (mins[i] > d) : mins[i] = d
          if (mins[j] > d) : mins[j] = d
        k += 1
    print "min distance: ", min(mins), ":", max(mins)
    return mins

##################################################################
  def max_distances(self):
    # assuming ds in a flattened trinagular matrix i.e. complete
    nrows = self.N_neighborhoods
    maxs = nrows*[0.]
    k = 0
    for i in range(nrows):
      for j in range(i+1,nrows):
        d = self.Pair_distances[k]
        if (maxs[i] < d) : maxs[i] = d
        if (maxs[j] < d) : maxs[j] = d
        k += 1
    print "max distance: ", min(maxs), ":", max(maxs)
    return maxs

##################################################################
  def mean_distances(self):
    # assuming ds in a flattened trinagular matrix
    nrows = self.N_neighborhoods
    means = nrows*[0.]
    k = 0
    for i in range(nrows):
      for j in range(i+1,nrows):
        d = self.Pair_distances[k]
        means[i] += d
        means[j] += d
        k += 1
    s = 1./nrows
    for i in range(nrows):
      means[i] *= s
    print "mean distance: ", min(means), ":", max(means)
    return means


##################################################################
  def is_complete(self):
    return N_distances == N_pairs

#=================================================================
#> OUTPUT
#=================================================================

##################################################################
  def write_sql(self,maxwidth=(SQLITE_MAX_COLS/2-1)):
    try:
      import sqlite3
    except:
      print "!!! no sqlite3 module !!!"
      return
    sqlname = self.stem+".sql"
    if os.path.exists(sqlname): os.unlink(sqlname)
    nb,ndm = self.Neighborhoods[0].shape
    print "> creating",str(nb)+"-cluster SQL table",
    sys.stdout.flush()
    conn = sqlite3.connect(sqlname)
    c = conn.cursor()
    ## cluster table
    cmd = "CREATE TABLE clusters (cid INTEGER, ctype INTEGER, cluster_type INTEGER, fx TEXT, fy TEXT, fz TEXT "
    for i in range(nb):
      ii = str(i+1)
      cmd += ", c"+ii+"x TEXT, c"+ii+"y TEXT, c"+ii+"z TEXT, c"+ii+"r TEXT"
    cmd += ")"
    c.execute(cmd)
    chunk_size = 1000
    nrec = 3+ndm+nb*(ndm+1)
    cmd = 'INSERT INTO clusters VALUES (?'+(nrec-1)*",?"+')'
    entries = []
    i = 0
    ctype = "1"
    cluster_type = nb*ctype
    for C in self.Neighborhoods:
      fx,fy,fz = map(str,self.Forces[i])
      i += 1
      cid = str(i)
      entry = [cid,ctype,cluster_type,fx,fy,fz]
      for row in C: 
        for v in row:
          entry.append(str(v))
        r = np.linalg.norm(row)
        entry.append(str(r))
      entries.append(entry)
      if len(entries) == chunk_size:
        c.executemany(cmd, entries)
        sys.stdout.write("*")
        sys.stdout.flush()
        entries = []
    if len(entries) > 0:
      c.executemany(cmd, entries)
    print "done"
    sys.stdout.flush()
    ## distance table
    #self.maxConnections = min(self.maxConnections,self.N_neighborhoods-1)
    self.maxConnections = min(maxwidth,self.distance_bandwidth+self.maxConnectionsBuffer) # NOTE buffer
    if self.maxConnections < self.distance_bandwidth:
      print "!!! insufficient bandwidth to write a complete SQL distance table"
    print "> creating",str(self.maxConnections)+"-distance SQL table",
    sys.stdout.flush()
    cmd = "CREATE TABLE cluster_distances (cid INTEGER";
    for i in range(self.maxConnections):
      ii = str(i+1)
      cmd += ", c"+ii+" INTEGER, d"+ii+" TEXT"
    cmd += ")"
    c.execute(cmd)
    nrec = 2*self.maxConnections
    cmd = 'INSERT INTO cluster_distances VALUES (?'+(nrec)*",?"+')'
    gids = self.Group_ids
    for i in self.References:
      gids[i] = -(gids[i]+1)
    nrefs = len(self.References)
    maxn = self.maxConnections-nrefs
    entries = []
    for i in range(len(self.Neighborhood_distances)):
      cid = str(i+1)
      entry = [cid]
      gid = gids[i]
      n = 0
      m = 0
      if gid < 0: 
        gid = -(gid+1)
        n = -1
        m = 1
      for j in range(len(self.Neighborhood_distances[i])):
        #print j,m,n
        [k,d] = self.Neighborhood_distances[i][j];
        if gids[k] < 0: # always add reference
          entry.append(str(k+1)) # 1-based
          entry.append(str(d))
          m += 1
        elif n < maxn: # add all regular upto quota
        #elif n < maxn and gids[k] == gid : # add all regular upto quota
          entry.append(str(k+1)) # 1-based
          entry.append(str(d))
          n += 1
        if (n+m == self.maxConnections): break
      #print i,gid,len(entry)
      entries.append(entry)
      if len(entries) == chunk_size:
        sys.stdout.flush()
        c.executemany(cmd, entries)
        sys.stdout.write("*")
        sys.stdout.flush()
        entries = []
    if len(entries) > 0:
      c.executemany(cmd, entries)
    print "done"
    sys.stdout.flush()
    ## reference table
    print "> creating",str(self.N_references)+"-reference SQL table",
    cmd = "CREATE TABLE reference_clusters (cid INTEGER); ";
    c.execute(cmd)
    cmd = 'INSERT INTO reference_clusters VALUES (?)'
    entries = []
    for i in self.References:
      entry = [i+1] # 1-based
      entries.append(entry)
      if len(entries) == chunk_size:
        c.executemany(cmd, entries)
        sys.stdout.write("*")
    if len(entries) > 0:
      c.executemany(cmd, entries)
    print "done"
    sys.stdout.flush()
    ## finalize
    conn.commit()
    conn.close()

#####################################################################
## MAIN
#####################################################################
if  __name__ == "__main__":
  print "!!! TEST !!!" 
  db = database()
  c1 = [[1.1, 1, 1],
        [1, -1, -1],
        [-1, 1, -1],
        [-1, -1, 1]]
  C1 = np.array(c1)
  c2 = [[1.1, 1, 1],
        [-1, 1, -1],
        [1, -1, -1],
        [-1, -1, 1]]
  C2 = np.array(c2)
  ## NOTE move to rmsd?
  db.trial_permutation(C1,C2)


  sys.exit()
  db_file = sys.argv[1]
  mtype = "SOAP"
  db = database(mtype)
  db.init(db_file)
  db.histogram()
  d_file = ((db_file.split("."))[0])+"."+mtype+"COPY"
  db.write_distances(d_file)
  db.read_distances(d_file)

  if (hasSparse) :
    amat = db.adjacency_matrix(3.0)
    #print amat
    #graph = pydot.graph_from_adjacency_matrix(amat.tolist())
    #print "> writing adjacency matrix to graph.pdf"
    #graph.write("./graph.pdf",prog="dot",format="pdf")

  
