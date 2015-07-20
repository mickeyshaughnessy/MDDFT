#!/usr/bin/env python

"""
metric database
TODO: precompute internal representation of neighborhoods
  read/write distances - limited bandwidth
"""
# system
import sys
import os
import math
import array
import re
import random
import numpy as np
import util
hasSparse = True
try:
  import scipy.sparse
except ImportError:
  print "!!! need scipy.sparse for adjacency matrix !!!"
  hasSparse = False


#####################################################
class graph:
#####################################################
  A = [] # adjacency matrix
  L = [] # laplacian matrix
  def __init__(self):
    self.exponent = 1
    self.init()

###################################################################
  def init(self):
    self.nnodes = 0
    self.edges = []
    self.A = []
    self.L = []
    
###################################################################
  def from_points(self,c1,c2,dmax=1.e10):
    self.init()
    ### NOTE assume len(c1) == len(c2)
    self.nnodes = n = len(c1)
    for i in range(n):
      x = c1[i]
      for j in range(n):
        y = c2[i]
        d = util.d(x,y)
        if (d < dmax): 
          self.edges.append([[i,j],d])

###################################################################
  def adjacency_matrix(self):
    if (len(self.A) > 0): return self.A.todense()
    n = self.nnodes
    A= scipy.sparse.lil_matrix((n,n))
    for e in self.edges:
      ii = e[0][0] ### NOTE might need to compact
      jj = e[0][1]
      d = e[1]
      A[jj,ii] = A[ii,jj] = d
    self.A = A.tocsr()
    return self.A.todense()

###################################################################
  def laplacian_matrix(self,dmax):
    if (len(self.L) > 0): return self.L.todense()
    n = self.nnodes
    L= scipy.sparse.lil_matrix((n,n))
    dsum = n*[0.]
    for e in self.edges:
      ii = e[0][0] ### NOTE might need to compact
      jj = e[0][1]
      d = e[1]
      L[jj,ii] = L[ii,jj] = -d
      dsum[ii] += d
    self.L = L.tocsr()
    return self.L.todense()

###################################################################
  def wg(self,d):
    if (math.fabs(d) < util._small) : return 0.0
    return d**self.exponent

###################################################################
  def norm(self):
    if (self.exponent ==1): return np.linalg.norm(self.adjacency_matrix())
    A = self.adjacency_matrix()
    print A
    L2 = 0.0
    n = self.nnodes
    for i in range(n):
      for j in range(i+1,n):
        d = float(A[i,j])
        print i,j,d
        L2 += self.wg(d)
    return math.sqrt(2*L2)

###################################################################
  def eigensystem(self):
    return np.linalg.eig(self.adjacency_matrix())

###################################################################
  def write(self,tag="G"):
    f = open(tag+".gv","w")
    print >>f,"digraph {"
    print >>f,"ratio=1.0"
    print >>f,"  node [shape=circle,color=skyblue,style=filled]"
    for i in range(0,self.nnodes):
      print >>f, "  C"+str(i+1),";";
    print >>f,"subgraph Dist { edge [dir=none,len=1]"
    for e in self.edges:
      i = e[0][0] ### NOTE might need to compact
      j = e[0][1]
      d = e[1]
      l = "%.2f" % d
      print >>f, "C"+str(i+1)," -> C"+str(j+1) + " [label = " +str(l)+" ]"
    print >>f,"  }"
    print >>f,"}"
    f.close()
    os.system("dot -Tpdf "+tag+".gv -o "+tag+".pdf")

#####################################################################
## MAIN
#####################################################################
if  __name__ == "__main__":
  print "!!! TEST !!!" 
  c1 = [[1.1, 1, 1],
        [1, -1, -1],
        [-1, 1, -1],
        [-1, -1, 1]]
  C1 = np.array(c1)
  c2 = [[1.1, 1, 1],
        [-1, 1, -1],
        [1, -1, -1],
        [-1, -1, 1]]
  #c2 = c1
  C2 = np.array(c2)
  g = graph()
  
  g.from_points(c1,c2)
  g.write("c12")
  print "# theta d"
  g.exponent = -2
  for i in range(19):
    theta = math.pi*i/6
    R = util.rotation([1,0,0],theta)
    #print R
    Rc2 = util.rotate_all(R,c2)
    #print Rc2
    g.from_points(c1,Rc2)
    print "{0:10.6f} {1:10.6f}".format(theta,g.norm())
