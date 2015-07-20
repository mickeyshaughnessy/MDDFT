#!/usr/bin/env python

"""
computes gaussian overlap distance
"""
import numpy as np
import math
try: 
  from scipy.optimize import minimize
except:
  print "!!! need scipy.optimize !!!"
  raise

###################################################################
class ogto:
  sigma = 1.0 # 0.5
  A = []
  B = []
  def __init__(self,s=sigma) :
    self.init(s)

###################################################################
  def init(self,s):
    self.sigma = s
    #self.a = 0.5*math.sqrt(math.pi)*s 
    self.a = (2.*math.sqrt(math.pi)*s)**3 # normalized
    #e4 = math.exp(0.25)
    #self.a = 8.*(e4-1.)*e4*s*s
    self.b = -0.25/(s*s)

###################################################################
  def tag(self):
    return "ogto_{0:5.3f}".format(self.sigma)

###################################################################
  def s(self,x, y):
    dx = x[0]-y[0]
    dy = x[1]-y[1]
    dz = x[2]-y[2]
    r2 = dx*dx+dy*dy+dz*dz
    return self.a*math.exp(self.b*r2)

###################################################################
  def norm2(self,A):
    nA = np.shape(A)[0]
    AA = 0.
    for i in range(nA):
      ab = self.s(A[i],A[i])
      AA += ab
      for j in range(i+1,nA):
        ab = self.s(A[i],A[j])
        AA += 2.*ab
    return AA
###################################################################
  def inner_product(self,A, B):
    nA = np.shape(A)[0]
    nB = np.shape(B)[0]
    AB = 0.
    for i in range(nA):
      for j in range(nB):
        ab = self.s(A[i],B[j])
        AB += ab
    return AB

###################################################################
  def rotation(self,angles):
    a = angles[0]
    c = math.cos(a)
    s = math.sin(a)
    R1 = [[c,-s,0],[s,c,0],[0,0,1.]]
    a = angles[1]
    c = math.cos(a)
    s = math.sin(a)
    R2 = [[c,0,s],[0,1.,0],[-s,0,c]]
    a = angles[2]
    c = math.cos(a)
    s = math.sin(a)
    R3 = [[c,-s,0],[s,c,0],[0,0,1.]]
    R = np.dot(R1,np.dot(R2,R3))
    return R

###################################################################
  def rot(self,R,A):
    RA  = np.transpose(np.dot(R,np.transpose(A)))
    return RA

###################################################################
  def r(self,angles):
    R = self.rotation(angles)
    RB = self.rot(R,self.B)
    ARB = -self.inner_product(self.A,RB)
    return ARB

###################################################################
  def optimal_rotation(self,A,B):
    ### minimize
    angles = [0,0,0] # many minima --> need to start close
    self.A = A
    self.B = B
    ntry = 0
    while True:
      #result = minimize(self.r, angles)
      #result = minimize(self.r, angles, method='nelder-mead',options={'xtol': 1e-8})
      #result = minimize(self.r, angles, method='BFGS', jac=Dr, options={'disp': True})
      #result = minimize(self.r, angles, method='CG',options={'xtol': 1e-8})
      result = minimize(self.r, angles, method='CG')
      ntry += 1
      if ntry > 0 : break
      if result.success : break
      else : angles = [angles[0]+math.pi,angles[1],angles[2]]
    if not result.success : print "!!! finding optimal rotation failed !!!"
    angles = result.x
    AB     = -result.fun
    #print "angles: ",angles 
    #print " AB ", AB
    R = self.rotation(angles)
    #RB = np.dot(B,R)
    #print B
    #print RB
    #AB = self.inner_product(A,RB)
    #print " AB ", AB
    return [R,AB]

###################################################################
  def distance(self,A, B,s=0.): # optimal
    if (s > 0. and s != self.sigma): self.init(s)
    #AA = self.inner_product(A,A)
    #BB = self.inner_product(B,B)
    AA = self.norm2(A)
    BB = self.norm2(B)
    [R,AB] = self.optimal_rotation(A,B)
    d2 = AA + BB - 2.*AB
    #print AA,BB,AB, d2
    if d2 < 0. :
      if d2 > -1.e-8 : d2 =0
      else          : print "!!! distance less than zero", d2,"!!!"
    d = math.sqrt(d2)
    return [d,R]

###################################################################
  def dist(self,A, B,s=0.):
    if (s > 0. and s != self.sigma): self.init(s)
    AA = self.inner_product(A,A)
    BB = self.inner_product(B,B)
    AB = self.inner_product(A,B)
    d = math.sqrt(AA + BB - 2*AB)
    R = np.eye(3)
    return [d,R]

#####################################################################
## MAIN
#####################################################################
if  __name__ == "__main__":
  import sys
  np.set_printoptions(suppress=True) # "chop"
  r = ogto()
  sValue = 0.5
  print "!!! TEST !!!"
  c1 = [[0, 1, 0]] 
  c2 = [[0, 2, 0]] 
  C1 = np.array(c1)
  C2 = np.array(c2)
  print "========== 2 point ====================="
  print "C1:",C1
  print "C2:",C2
  [d12,R12] = r.dist    (C1,C2,sValue)
  print "self d=",d12," unrotated"
  [d12,R12] = r.distance(C1,C2,sValue) # underdetermined
  print "self d=",d12,"\n R=\n",R12
  print
  print "======== permutation ==================="
  c1 = [[ 0.26289221,  0.23899291,  0.23899291],
        [ 0.25366903, -0.25366903, -0.25366903],
        [-0.25366903,  0.25366903, -0.25366903],
        [-0.25366903, -0.25366903,  0.25366903]]
  c1p= [[ 0.26289221,  0.23899291,  0.23899291],
        [-0.25366903,  0.25366903, -0.25366903],
        [ 0.25366903, -0.25366903, -0.25366903],
        [-0.25366903, -0.25366903,  0.25366903]]
  C1 = np.array(c1)
  C2 = np.array(c1p)
  [d12,R12] = r.distance(C1,C2)
  print "self d=",d12,"\n R=\n",R12
  print  
  print "========= rotation ====================="
  c1 = [[1.1, 1, 1],
        [1, -1, -1],
        [-1, 1, -1],
        [-1, -1, 1]]
  C1 = np.array(c1)
  print "C1:\n",C1
  theta = math.pi/3.
  print "theta=",theta
  R = np.array(r.rotation([theta,0,0]))
  print "R=\n",R
  RC1 = np.transpose(np.dot(R,np.transpose(C1)))
  print "R.C1 = \n",RC1
  #print "\n----- Self    -----"
  #[d1, R1]  = r.dist(C1,C1)
  #print "self d=",d1, " unrotated"
  #[d11,R11] = r.distance(C1,C1)
  #print "self d=",d11,"\n R=\n",R11
  #print "\n----- Rotated -----"
  [d1, R1]  = r.dist(C1,RC1)
  print "self d=",d1, " unrotated"
  [d1r,R1r] = r.distance(C1,RC1)
  print "rot  d=",d1r,"\n R=\n",R1r

  ##############################################################
  c1 = [[ 0,  0,  0]]
  C1 = np.array(c1)
  n = 40
  f = open("2pt.dat","w")
  print >>f,"# sigma=",r.sigma,r.a,r.b
  dx = 0.1
  sigmas = [0.1,0.5,1.0,2.0]
  for i in range(n+1):
    dr = dx*i
    c2 = [[ c1[0][0]+dr,  0,  0]]
    C2 = np.array(c2)
    ds = []
    for s in sigmas:
      r.init(s)
      [d12, R12]  = r.dist(C1,C2)
      #d = math.sqrt(r.s(c1[0],c1[0])+r.s(c2[0],c2[0])-2*r.s(c1[0],c2[0]))
      #print d-d12,
      ds.append(d12)
    #print
    print >>f,"{0:7.4f}".format(dr),
    for d in ds:
      print >>f,"{0:7.4f}".format(d),
    print >>f
  f.close()
  sys.exit()
  ##############################################################
  c1 = [[ 1,  1,  1],
        [ 1, -1, -1],
        [-1,  1, -1],
        [-1, -1,  1.1]]
  C1 = np.array(c1)
  n = 40
  f = open("ogto.dat","w")
  for i in range(n+1):
    theta = 2.*math.pi*i/n
    c = math.cos(theta)
    s = math.sin(theta)
    R = [[c,-s, 0],
         [s, c, 0],
         [0, 0, 1]]
    RC1 = np.transpose(np.dot(R,np.transpose(C1)))
    print i,theta
    #print RC1
    [d1, R1]  = r.dist(C1,RC1)
    [d1r,R1r] = r.distance(C1,RC1)
    C1r = np.transpose(np.dot(R1r,np.transpose(RC1)))
    #print C1r
    print >>f,"{0:7.4f} {1:7.4f} {2:7.4f}".format(theta,d1,d1r)
  f.close()
