#!/usr/bin/env python

"""
computes RMS-D distance
"""
import numpy as np
import math


###################################################################
class rmsd:
  def __init__(self) :
    pass
###################################################################
  def tag(self):
    return "rmsd"
###################################################################
  def dist(self,A, B):
    n = np.shape(A)[0]
    rmsd_sq = (sum(sum(A*A)) + sum(sum(B*B)) - 2.*sum(sum(A*B)))/float(n)
    rmsd_sq = max([rmsd_sq, 0.0])
    R = np.eye(3)
    return [np.sqrt(rmsd_sq),R]

###################################################################
  def distance(self,A, B):
    n = np.shape(A)[0]
    correlation_matrix = np.dot(np.transpose(A), B)
    #print "C=",correlation_matrix
    U,S,V = np.linalg.svd(correlation_matrix)
    is_reflection = (np.linalg.det(U) * np.linalg.det(V)) < 0.0
    if is_reflection: S[-1] = - S[-1]
    E0 = sum(sum(A * A)) + sum(sum(B * B))
    #print "E0 ", E0, sum(sum(A * A)) , sum(sum(B * B))
    rmsd_sq = (E0 - 2.0*sum(S)) / float(n)
    rmsd_sq = max([rmsd_sq, 0.0])
    R = np.dot(U,V)
    #print "D ",s
    #print "U ",np.transpose(U)
    #print "V ",V
    #print "R ",R
    return [np.sqrt(rmsd_sq),R]

###################################################################
  def d(self,A, B):
    dAB,R = self.distance(A,B)
    return dAB

###################################################################
class wrmsd:
  sigma = 1.
  def __init__(self,s=1) :
    self.sigma = 1./s   # need a function with a cutoff
  
###################################################################
  def r2(self,x):
    return x[0]*x[0]+x[1]*x[1]+x[2]*x[2] 

###################################################################
  def r(self,x):
    return math.sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2])

###################################################################
  def w(self,r):
    return math.exp(-r*self.sigma)

###################################################################
  def weight(self,a): # note a separable weight could be included in the db
    n = len(a)
    m = len(a[0])
    rs = map(self.r,a)
    #print rs
    ws = map(self.w,rs)
    wsum = sum(ws)
    #print ws
    wa =  []
    for i in range(n):
      w = ws[i]/wsum
      wa.append(w*a[i])
    return np.array(wa)

###################################################################
  def distance(self,A, B):
    n = np.shape(A)[0]
    wA = self.weight(A)
    wB = self.weight(B)
    #print "WA=", wA
    #print "WB=", wB
    correlation_matrix = np.dot(np.transpose(wA), wB)
    #print "C=",correlation_matrix
    U,S,V = np.linalg.svd(correlation_matrix)
    is_reflection = (np.linalg.det(U) * np.linalg.det(V)) < 0.0
    if is_reflection: S[-1] = - S[-1]
    E0 = sum(sum(wA * wA)) + sum(sum(wB * wB))
    #print "E0 ", E0, sum(sum(wA * wA)) , sum(sum(wB * wB))
    rmsd_sq = (E0 - 2.0*sum(S)) / float(n)
    rmsd_sq = max([rmsd_sq, 0.0])
    R = np.dot(U,V)
    #print "D ",s
    #print "U ",np.transpose(U)
    #print "V ",V
    #print "R ",R
    return [np.sqrt(rmsd_sq),R]

#####################################################################
## MAIN
#####################################################################
if  __name__ == "__main__":
  r = rmsd()
  print "!!! TEST !!!"
  c1 = [[1.1, 1, 1],
        [1, -1, -1],
        [-1, 1, -1],
        [-1, -1, 1]]
  c1 = [[ 0.26289221,  0.23899291,  0.23899291],
        [ 0.25366903, -0.25366903, -0.25366903],
        [-0.25366903,  0.25366903, -0.25366903],
        [-0.25366903, -0.25366903,  0.25366903]]
  C1 = np.array(c1)
  theta = math.pi/3.
  c = math.cos(theta)
  s = math.sin(theta)
  R = [[c,-s,0],
       [s, c,0],
       [0, 0,1]]
  print "R=\n",np.array(R)
  RC1 = np.transpose(np.dot(R,np.transpose(C1)))
  print "  C1 = \n",C1
  print "R.C1 = \n",RC1
  print "\n----- Self    -----"
  [d11,R11] = r.distance(C1,C1)
  print "self d=",d11,"\n R=",R11
  print "\n----- Rotated -----"
  [d1r,R1r] = r.distance(C1,RC1)
  print "rot  d=",d1r,"\n R=",R1r

  print "\n====== Weighted ======"
  c1 = [[1.1, 1, 1],
        [1, -1, -1],
        [-1, 1, -1],
        [-1, -1, 1]]
  C1 = np.array(c1)
  RC1 = np.transpose(np.dot(R,np.transpose(C1)))
  r = wrmsd(1.)
  print "\n----- Self    -----"
  [d11,R11] = r.distance(C1,C1)
  print "self d=",d11," R=\n",R11
  print "\n----- Rotated -----"
  [d1r,R1r] = r.distance(C1,RC1)
  print "rot  d=",d1r," R=\n",R1r
