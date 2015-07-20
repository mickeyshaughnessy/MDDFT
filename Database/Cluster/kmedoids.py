#!/usr/bin/env python
import random
import sys

#######################################################################
class kmedoids:
  # per point
  points = [] # positions
  cluster_ids = [] 
  # per cluster
  clusters = [] # cluster of point ids
  centers  = []
  def __init__(self,X,K,metric):
    self.points = X
    self.K = K
    self.metric = metric
    self.N = len(X)
    self.cluster_ids = self.N*[0]
    self.centers     = self.K*[-1]

  ####################################################################
  def cluster(self):
    self.centers = random.sample(range(self.N), self.K)
    previous_centers = self.K*[-2]
    while not self.converged(self.centers, previous_centers):
      previous_centers = [c for c in self.centers]
      self.cluster_points() # given centers
      self.update_centers() # given clusters
      self.report()
    return self.centers

  ####################################################################
  def distance(self,X,Y):
    return self.metric.d(X,Y)

  ####################################################################
  def closest(self,x,Xs):
    dlist = [(i, self.distance(x,Xs[i])) for i in range(len(Xs))]
    return min(dlist, key=lambda t:t[1])[0]

  ####################################################################
  def min_total_distance(self,cs):
    n = len(cs)
    ds = n*[0.]
    for i in range(n):
      for j in range(i+1,n):
        d = self.distance(self.points[cs[i]],self.points[cs[j]])
        ds[i] += d
        ds[j] += d
    dlist = [(i, ds[i]) for i in range(n)]
    idx = min(dlist, key=lambda t:t[1])[0]
    return cs[idx]

  ####################################################################
  def cluster_points(self):
    self.clusters = []
    centers = []
    for i in range(self.K):
      self.clusters.append([])
      centers.append(self.points[self.centers[i]])
    for i in range(self.N):
      x = self.points[i]
      c = self.closest(x,centers)
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
    print "> cluster size, center:"
    for i in range(self.K):
      print " {0:4d}: {1:4d} {2:4d}".format(i,len(self.clusters[i]),self.centers[i])

  ####################################################################
  def write(self,filename="clusters.dat"):
    print "> writing to "+filename
    f = open(filename,"w")
    print >>f, "# points"
    i = 0
    for X in self.points:
      for x in X:
        print >>f, "{0:14.9f}".format(x),
      print >>f, self.cluster_ids[i]
      i += 1
    print >>f, "\n"
    print >>f, "# centers"
    for i in range(self.K):
      c = self.centers[i]
      print >>f, "{0:5d}".format(c),
      for x in self.points[c]:
        print >>f, "{0:14.9f}".format(x),
      print >>f

###############################################################################
import numpy as np
class euclidean:
  def distance(self, A, B):
    return np.linalg.norm(A-B)
  def d(self, A, B):
    return np.linalg.norm(A-B)
#######################################################################
def uniform_points(N):
    X = np.array([(random.uniform(-1, 1), random.uniform(-1, 1)) for i in range(N)])
    return X

###############################################################################
if __name__ == "__main__":
  K = 4 
  X = uniform_points(100)
  m = euclidean()
  c = kmedoids(X,K,m)
  c.cluster()
  c.write()
  c.report()
