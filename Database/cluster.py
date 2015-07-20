#!/usr/bin/env python
import random

#######################################################################
class kcluster:
  # per point
  cluster_ids = [] 
  # per cluster
  clusters = [] # cluster of point ids
  centers  = []
  def __init__(self,A,K):
    self.distance_matrix = A
    self.K = K
    self.N = A.shape[0]
    self.cluster_ids = self.N*[0]
    self.centers     = self.K*[-1]

  ####################################################################
  def cluster(self):
    print "> making",self.K,"clusters of",self.N,"points"
    self.centers = random.sample(range(self.N), self.K)
    previous_centers = self.K*[-2]
    while not self.converged(self.centers, previous_centers):
      previous_centers = [c for c in self.centers]
      self.cluster_points() # given centers
      self.update_centers() # given clusters
      self.report()
    return (self.centers,self.cluster_ids)

  ####################################################################
  def closest(self,a,bs):
    dlist = [(i, self.distance_matrix[a,bs[i]]) for i in range(len(bs))]
    return min(dlist, key=lambda t:t[1])[0]

  ####################################################################
  def min_total_distance(self,cs):
    n = len(cs)
    ds = n*[0.]
    for i in range(n):
      for j in range(i+1,n):
        d = self.distance_matrix[cs[i],cs[j]]
        ds[i] += d
        ds[j] += d
    dlist = [(i, ds[i]) for i in range(n)]
    idx = min(dlist, key=lambda t:t[1])[0]
    return cs[idx]

  ####################################################################
  def distance_range(self,cs):
    n = len(cs)
    dmin = 1.e20
    dmax = 0.
    for i in range(n):
      for j in range(i+1,n):
        d = self.distance_matrix[cs[i],cs[j]]
        if d > dmax: dmax = d
        if d < dmin: dmin = d
    return (dmin,dmax)

  ####################################################################
  def cluster_points(self):
    self.clusters = []
    for i in range(self.K):
      self.clusters.append([])
    for i in range(self.N):
      c = self.closest(i,self.centers)
      self.cluster_ids[i] = c
      self.clusters[c].append(i)
    return self.clusters

  ####################################################################
  def update_centers(self):
    self.centers = []
    for c in self.clusters:
      self.centers.append(self.min_total_distance(c))
    return self.centers

  ####################################################################
  def converged(self, list1, list2): 
    return (set(list1) == set(list2))

  ####################################################################
  def report(self):
    sizes = map(len,self.clusters)
    self.minsize = min(sizes)
    self.maxsize = max(sizes)
    print "> cluster size center d_min d_max:",
    print " (size ",self.minsize,"-",self.maxsize,")"
    for i in range(self.K):
      (dmin,dmax) = self.distance_range(self.clusters[i])
      print " {0:4d}: {1:4d} {2:4d} {3:10.4f} {4:10.4f}".format(i,sizes[i],self.centers[i],dmin,dmax)

  ####################################################################
  def write(self,filename="clusters.dat"):
    print "> writing to "+filename
    f = open(filename,"w")
    print >>f, "# index cluster"
    for i in range(self.N):
      print >>f, "{0:5d} {1:3d}".format(i,self.cluster_ids[i]),
    print >>f, "\n"
    print >>f, "# centers"
    for i in range(self.K):
      c = self.centers[i]
      print >>f, "{0:5d}".format(c),
      print >>f

###############################################################################
if __name__ == "__main__":
  import sys
  import numpy as np
  K = 4 
  N = 40
  A = np.zeros((N,N))
  for i in range(N):
    for j in range(i+1,N):
       A[i,j] = A[j,i] = abs(i-j)
  print A
  c = kcluster(A,K)
  c.cluster()
  c.write()
  c.report()
