#!/usr/bin/env python
import numpy as np
import random

# TODO: change to arb distance metric
# TODO: find centers closest to means

#######################################################################
### Lloyd's algorithm
#######################################################################
class kmeans:
  # per point
  points = [] # positions
  cluster_ids = [] 
  # per cluster
  clusters = {} # id-->points
  means    = [] # position
  indices  = [] # 
  centers  = []
  def __init__(self,X,K):
    self.points = X
    self.N = len(X)
    self.K = K
    self.cluster_ids = self.N*[0]
    self.centers     = self.K*[-1]

  ####################################################################
  def cluster(self):
    self.means = random.sample(self.points, self.K)
    previous_centers = self.K*[-2]
    while not self.converged(self.centers, previous_centers):
      previous_centers = [c for c in self.centers]
      self.cluster_points()
      self.update_means()
      print self.means
      self.find_centers()
      print self.centers, previous_centers
    return self.means

  ####################################################################
  def distance(self,X,Y):
    return np.linalg.norm(X-Y)

  ####################################################################
  def average(self,Xs):
    return np.mean(Xs, axis = 0)

  ####################################################################
  def closest(self,x,Xs):
    dlist = [(i, self.distance(x,Xs[i])) for i in range(len(Xs))]
    return min(dlist, key=lambda t:t[1])[0]

  ####################################################################
  def cluster_points(self):
    self.clusters = []
    self.indices  = []
    for i in range(self.K):
      self.clusters.append([])
      self.indices.append([])
    for i in range(self.N):
      x = self.points[i]
      c = self.closest(x,self.means)
      self.cluster_ids[i] = c
      self.clusters[c].append(x)
      self.indices[c].append(i)
    return self.clusters

  ####################################################################
  def update_means(self):
    self.means = []
    for c in self.clusters:
      self.means.append(self.average(c))
    return self.means

  ####################################################################
  def converged(self, list1, list2): 
    return (set(list1) == set(list2))

  ####################################################################
  def find_centers(self):
    for i in range(self.K):
      self.centers[i] = self.closest(self.means[i],self.points)
    return self.centers

  ####################################################################
  def write(self,filename="clusters.dat"):
    f = open(filename,"w")
    print >>f, "# points"
    i = 0
    for X in self.points:
      for x in X:
        print >>f, "{0:14.9f}".format(x),
      print >>f, self.cluster_ids[i]
      i += 1
    print >>f, "\n"
    print >>f, "# clusters"
    for c in self.indices:
      for i in c:
        print >>f, "{0:5d}".format(i),
      print >>f
    print >>f, "\n"
    print >>f, "# means"
    for i in range(self.K):
      for x in self.means[i]:
        print >>f, "{0:14.9f}".format(x),
      c = self.centers[i]
      print >>f, "{0:5d}".format(c),
      for x in self.points[c]:
        print >>f, "{0:14.9f}".format(x),
      print >>f

#######################################################################
def init_board_uniform(N):
    X = np.array([(random.uniform(-1, 1), random.uniform(-1, 1)) for i in range(N)])
    return X

#######################################################################
def init_board_gauss(N, k):
    n = float(N)/k
    X = []
    for i in range(k):
        c = (random.uniform(-1, 1), random.uniform(-1, 1))
        s = random.uniform(0.05,0.5)
        x = []
        while len(x) < n:
            a, b = np.array([np.random.normal(c[0], s), np.random.normal(c[1], s)])
            # Continue drawing points from the distribution in the range [-1,1]
            if abs(a) < 1 and abs(b) < 1:
                x.append([a,b])
        X.extend(x)
    X = np.array(X)[:N]
    return X

###############################################################################
if __name__ == "__main__":
  K = 4
  #X = init_board_uniform(100)
  X = init_board_gauss(100,K)
  km = kmeans(X,K)
  km.cluster()
  km.write()
