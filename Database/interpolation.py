#!/usr/bin/env python

"""
interpolation
"""
import sys
import math
import numpy as np

###################################################################
class radial_basis_function:
###################################################################
  sigma = 1.0
  n = 0
  values = []
  weights = []
  def __init__(self,s=1.):
    self.sigma = 1./s
  ###################################################################
  def set_sigma(self,s):
    self.sigma = 1./s
  ###################################################################
  def solve(self,dij,vals):
    self.values = vals
    A = []
    self.n = len(vals)
    #print self.n
    for i in range(self.n):
      row = []
      for j in range(self.n):
        row.append(self.value(dij[i][j]))
      A.append(row)
    #print A
    self.weights = np.linalg.solve(A, vals)
    return self.weights
  ###################################################################
  def interpolate(self,ds):
    v = 0.
    for i in range(self.n):
      v += self.value(ds[i])*self.weights[i]
    return v

###################################################################
class inverse_quadratic(radial_basis_function):
###################################################################
  def __init__(self,s=1.):
    radial_basis_function.__init__(self,s) 
  #################################################################
  def value(self,r):
    y = radial_basis_function.sigma*r
    return 1./(1.+y*y)

###################################################################
class inverse_exponential(radial_basis_function):
###################################################################
  def __init__(self,s=1.):
    radial_basis_function.__init__(self,s) 
  #################################################################
  def value(self,r):
    return math.exp(-r*radial_basis_function.sigma)

###################################################################
class gaussian(radial_basis_function):
###################################################################
  def __init__(self,s=1.):
    radial_basis_function.__init__(self,s) 
  #################################################################
  def value(self,r):
    y = radial_basis_function.sigma*r
    return math.exp(-y*y)

###################################################################
class interpolation:
###################################################################
  def __init__(self,type="shepard",s=1.):
    if   (type == "shepard"):
      rbf = inverse_quadratic(s)
    elif (type == "exponential"):
      rbf = inverse_exponential(s)
    elif (type == "gaussian"):
      rbf = gaussian(s)

######################################################################
  def rbf(self,r): # radial basis function
    return self.rbf.value(r)

###################################################################
if  __name__ == "__main__":
 print "!!! TEST !!!"
 c1 = np.array([[1.1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]])
 f1 = np.array([1,0,0])
 c2 = np.array([[1.1, 1.1, 1.1], [1.1, -1.1, -1.1], [-1.1, 1.1, -1.1], [-1.1, -1.1, 1.1]])
 f2 = np.array([0,1,0])
 c3 = np.array([[1, 1, 1], [1, -1, -1.01], [-1, 1, -1], [-1, -1, 1]])
 f3 = np.array([0,0,1])
 c4 = np.array([[1, 2, 1], [1, -1, -1.01], [-1, 1, -1], [-1, -1, 1]])
 f4 = np.array([1,0,1])
 print "c1\n",c1
 print "f1 ",f1
 print "c2\n",c2
 print "f2 ",f2
 print "c3\n",c3
 print "f3 ",f3
 print "c4\n",c4
 print "f4 ",f4

 import rmsd
 m= rmsd.rmsd()
 [d12,R12] = m.distance(c1,c2)
 [d13,R13] = m.distance(c1,c3)
 [d14,R14] = m.distance(c1,c4)
 [d23,R23] = m.distance(c2,c3)
 [d24,R24] = m.distance(c2,c4)
 [d34,R34] = m.distance(c3,c4)

 Ds = np.array([[0,d12,d13,d14],[d12,0,d23,d24],[d13,d23,0,d34],[d14,d24,d34,0]])

 print "using an inverse quadratic"
 rbf = inverse_quadratic()
 print "relative distances/adjacency matrix\n",Ds
 fs = np.array([f1,f2,f3,f4])
 print "(rotated) force values\n",fs
 w = rbf.solve(Ds,fs)
 print "weights:\n", w

 f = open("interpolation.dat","w")
 n = 10
 dalpha = 1./n
 for j in range(n+1):
   alpha = j*dalpha
   cQ = alpha*c1+(1.-alpha)*c2
   [d1,R1] = m.distance(c1,cQ)
   [d2,R2] = m.distance(c2,cQ)
   [d3,R3] = m.distance(c3,cQ)
   [d4,R4] = m.distance(c4,cQ)
   ds = [d1,d2,d3,d4]
   print >>f, alpha,
   for d in ds:
     print >>f,d,
   fs = rbf.interpolate(ds)
   for fi in fs:
     print >>f,fi,
   print >>f
