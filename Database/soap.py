#!/usr/bin/env python

"""
metric/s from Bartok's SOAP similarity measure
"""
import sys
import math
import array
try:
  import numpy as np
except ImportError:
  print "!!! need numpy !!!"
  sys.exit(1)

try:
  import scipy.special
except ImportError:
  print "!!! need scipy.special !!!"
  sys.exit(1)
#from scipy.special import sph_harm
#from scipy.special import iv

###################################################################
class soap:
###################################################################
  pi = math.pi
  sigma = 1.0
  pmax = 4
  def __init__(self):
    self.pi = math.pi
    s = 1./math.sqrt(2.)
    self.A = np.array([[-s,    0., s   ],
                       [-s*1j, 0.,-s*1j],
                       [  0.,  1., 0.  ]])
    self.At= self.A.conjugate().transpose()
###################################################################
  def tag(self):
    return "soap_{0:d}_{1:5.3f}".format(self.pmax,self.sigma)

###################################################################
  def cart_to_sph(self,a):
    x = a[0]
    y = a[1]
    z = a[2]
    r = math.sqrt(x*x+y*y+z*z)
    theta = math.atan(y/x)
    phi   = math.acos(z/r)
    return [r,theta,phi]

###################################################################
  def sph_to_cart(self,a):
    r     = a[0]
    theta = a[1]
    phi   = a[2]
    x = r*math.sin(phi)*math.cos(theta)
    y = r*math.sin(phi)*math.sin(theta)
    z = r*math.cos(phi)
    return [x,y,z]

###################################################################
  def gaussian(self,a,sigma=1.0):
    x = a[0]
    y = a[1]
    z = a[2]
    r2 = x*x+y*y+z*z
    exp = math.exp(-r2/(4*sigma*sigma))
    return exp

###################################################################
  def G(self,r,sigma=1.0):
    exp = math.exp(-r*r/(4.*sigma*sigma))
    return exp

###################################################################
  def Y(self,l,m,theta,phi):
    # theta in [0,2pi]  azimuthal/longitudinal
    # phi   in [0,pi]   polar/colatitude
    return scipy.special.sph_harm(m,l,theta,phi) 

###################################################################
  def W(self,l,x,sigma=1.0):
    y = x/(2.*sigma*sigma)
    iota = scipy.special.iv(l,y)  # besselI
    return 4.*self.pi*iota

###################################################################
  def upsilon(self,l,m,x,sigma=1.):
    r     = x[0]
    theta = x[1]
    phi   = x[2]
    exp = self.G(r,sigma)
    y   = self.Y(l,m,theta,phi)
    return exp*y

###################################################################
### NOTE inner_product
  def compute_I(self,c1, c2, p=0, sigma=1):
    n = 2*p+1
    I = np.zeros((n,n),complex)
    for m1 in range(-p,p+1):
      i = m1+p
      for m2 in range(-p,p+1):
        j = m2+p
        Iab = 0
        for a1 in range(len(c1)): # atoms in cluster 1
          x1 = self.cart_to_sph(c1[a1])
          r1     = x1[0]
          u1 = self.upsilon(p,m1,x1,sigma)
          for a2 in range(len(c2)): # atoms in cluster 2
            #x2 = c2[a2]
            x2 = self.cart_to_sph(c2[a2])
            r2     = x2[0]
            u2 = self.upsilon(p,m2,x2,sigma)
            w = self.W(p,r1*r2,sigma)
            Iab += u1.conjugate()*w*u2
        print "I_"+str(i+1)+","+str(j+1),"=",Iab
        I[i,j] = (Iab.conjugate()*Iab)
    return I

###################################################################
### NOTE inner_product
  def compute_S(self,c1, c2, p=0, sigma=1):
    S = 0
    for l in range(p+1):
      for m1 in range(-l,l+1):
        for m2 in range(-l,l+1):
          I= 0.
          for a1 in range(len(c1)): # atoms in cluster 1
            #x1 = c1[a1]
            x1 = self.cart_to_sph(c1[a1])
            r1     = x1[0]
            u1 = self.upsilon(l,m1,x1,sigma)
            for a2 in range(len(c2)): # atoms in cluster 2
              #x2 = c2[a2]
              x2 = self.cart_to_sph(c2[a2])
              r2     = x2[0]
              u2 = self.upsilon(l,m2,x2,sigma)
              w = self.W(l,r1*r2,sigma)
              I += u1.conjugate()*w*u2
          S += (I.conjugate()*I)
    return S

###################################################################
  def distance(self,c1, c2, p=0, sigma=1.): # optimal
    S11 = self.compute_S(c1,c1,p,sigma)
    S12 = self.compute_S(c1,c2,p,sigma)
    S22 = self.compute_S(c2,c2,p,sigma)
    #print S11,S12,S22
    #d = math.sqrt(S11.real*S22.real)/S12.real - 1. ### NOTE
    #d = math.sqrt(S11.real+S22.real-S12.real)
    d = math.sqrt((S11+S22-2.*S12).real) ## #NOTE
    R = np.eye(3) ### HACK
    return [d,R]

###################################################################
  def dist(self,c1, c2, p=0, sigma=1.):
    S11 = self.compute_S(c1,c1,p,sigma)
    S12 = self.compute_S(c1,c2,p,sigma)
    S22 = self.compute_S(c2,c2,p,sigma)
    d = math.sqrt((S11+S22-2.*S12).real) ## #NOTE
    R = np.eye(3) 
    return [d,R]
  
###################################################################
if  __name__ == "__main__":
 import util
 d = soap()
 print "BENCHMARK against mathematica"
 print "!!! TEST !!!"
 c1 = [[1.1, 1, 1], 
       [1, -1, -1], 
       [-1, 1, -1], 
       [-1, -1, 1]]
 c2 = [[1.1, 1.1, 1.1], 
       [1.1, -1.1, -1.1], 
       [-1.1, 1.1, -1.1], 
       [-1.1, -1.1, 1.1]]
 eps = 1.e-4
 c1 = [[eps,   eps, 0]]
 c2 = [[eps, 2*eps, 0]]
 s1 = map(d.cart_to_sph,c1)
 s2 = map(d.cart_to_sph,c2)
 print "c1 "
 for i in range(len(c1)):
   print c1[i]," ==> ",s1[i]
 print "c2 "
 for i in range(len(c2)):
   print c2[i]," ==> ",s2[i]
 print 

 I0 = d.compute_I(c1,c2,0)
 S0 = d.compute_S(c1,c2,0)
 print "I_0=", I0
 print "S_0=", S0
 I1 = d.compute_I(c1,c2,1)
 S1 = d.compute_S(c1,c2,1)-S0
 print "I_1=\n", I1
 print "S_1=", S1
 print

 p = 4
 print "# permutation test order=",p
 c2p = [c2[1],c2[2],c2[3],c2[0]]
 d12,R12 = d.distance(c1,c2,p) 
 d12p,R12p = d.distance(c1,c2,p) 
 print "{0:10f}".format(d12)
 print R12
 print "-- permutation --"
 print "{0:10f}".format(d12p)
 print R12p
 print 

 sigmas = [0.1,0.5,1.0]
 print "# triangle test: alpha, d1 d2, ... for sigmas=",sigmas
 n = 10
 for i in range(n+1):
   a = 1.0*i/n
   c = util.interpolate(c1,c2,a)
   print "{0:4g}".format(a),
   for sig in sigmas:
     d1,R1 = d.distance(c1,c,p,sig) 
     d2,R2 = d.distance(c2,c,p,sig)
     print "{0:10f}{1:10f}".format(d1,d2),
   print 
 print 
   

 print "# d vs p: p d12 .. for sigmas=",sigmas
 for p in range(10):
   print "{0:2d}".format(p),
   for sig in sigmas:
     #d11 = compute_distance(s1,s1,p)
     d12,R12 = d.distance(c1,c2,p,sig)
     #d22 = compute_distance(s2,s2,p)
     print "{0:10f}".format(d12),
   print
 print
