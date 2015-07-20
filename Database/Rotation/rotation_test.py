#!/usr/bin/env python

"""
metric/s from Bartok's SOAP similarity measure
"""
import sys
import math
import os
import numpy as np
lib_path = os.path.abspath('..')
sys.path.append(lib_path)
import util
import soap
import rmsd
import ogto

###################################################################
r = rmsd.rmsd()
o = ogto.ogto()
c1 = [[2.2, 2, 2], 
      [2, -2, -2], 
      [-2, 2, -2], 
      [-2, -2, 2]]
c2 = [[1.1, 1.1, 1.1], 
      [1.1, -1.1, -1.1], 
      [-1.1, 1.1, -1.1], 
      [-1.1, -1.1, 1.1]]
c2 = [[2.2, 2, 2], 
      [2, -2, -2], 
      [-2, 2, -2], 
      [-2, -2, 2]]
C1 = np.array(c1)
C2 = np.array(c2)

sigmas = [0.1,0.5,1.0,2.0]

######################################################################
def rotation_test():
   print ">> running rotation test ",
   sys.stdout.flush()
   f = open("rotation.dat","w")
   p = 4
   print >>f,"# rotation test "
   n = 40
   alphas = []
   cs = []
   p = [math.sqrt(2.),math.sqrt(2.),0]
   p = [0,1,0]
   for i in range(n+1):
     a = 2.0*math.pi*i/n
     R = util.rotation(p,a)
     c = np.array(util.rotate_all(R,c1))
     alphas.append(a)
     cs.append(c)
   g0,R0 = r.distance(C1,C2)
   print >>f,"# RMSD ref =",g0
   print >>f,"# OGTO ref =",
   o0s = []
   for sig in sigmas: 
     d0,R0 = o.distance(c1,c2,sig)
     o0s.append(d0)
     print >>f,d0,
   print >>f
   print >>f,"# alpha d_rmsd, d_ogto @ sig=",
   for sig in sigmas: print >>f, sig,
   print >>f
   for i in range(len(alphas)):
     sys.stdout.write("*")
     sys.stdout.flush()
     a = alphas[i]
     # optimal rot
     d1,R1r = r.distance(C1,cs[i]) 
     print >>f,"{0:8.5f} ".format(a),
     print >>f,"{0:10.7f} ".format(d1),
     for j in range(len(sigmas)):
       sig = sigmas[j]
       d0 = o0s[j]
       d1,R1 = o.distance(c1,cs[i],sig) 
       print >>f,"{0:10.7f} ".format(d1),
#      if (d1 < 0.01):
#        print 
#        print R1r
#        print R1
     # no rot
     d1,R1 = r.dist(C1,cs[i]) 
     print >>f,"{0:10.7f} ".format(d1),
     for j in range(len(sigmas)):
       sig = sigmas[j]
       d0 = o0s[j]
       d1,R1 = o.dist(c1,cs[i],sig) 
       print >>f,"{0:10.7f} ".format(d1),
     print >>f
   print "done"
####################################################
if  __name__ == "__main__":
 rotation_test()
