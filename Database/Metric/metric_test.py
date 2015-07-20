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
d = soap.soap()
g = rmsd.rmsd()
o = ogto.ogto()
c1 = [[1.1, 1, 1], 
      [1, -1, -1], 
      [-1, 1, -1], 
      [-1, -1, 1]]
c1 = [[2.2, 2, 2], 
      [2, -2, -2], 
      [-2, 2, -2], 
      [-2, -2, 2]]
c2 = [[1.1, 1.1, 1.1], 
      [1.1, -1.1, -1.1], 
      [-1.1, 1.1, -1.1], 
      [-1.1, -1.1, 1.1]]
#-------------------------------
c1 = [[ 1.23,  1.43,  1.43], 
      [ 1.43, -1.43, -1.43], 
      [-1.43,  1.41, -1.43], 
      [-1.63, -1.43,  1.43]]
c2 = [[ 1.43,  1.43,  1.45], 
      [ 1.45, -1.43, -1.43], 
      [-1.43,  1.43, -1.43], 
      [-1.43, -1.53,  1.43]]
#-------------------------------
c1 = [[1.1, 1, 1], 
      [1, -1, -1], 
      [-1, 1, -1], 
      [-2, -2, 2]]
c2 = [[2.2, 2, 2], 
      [2, -2, -2], 
      [-2, 2, -2], 
      [-4, -4, 4]]
#-------------------------------
eps = 1.e-4
#c1 = [[eps,   eps, 0]] 
#c2 = [[eps, 2*eps, 0]] 
## NOTE this generated the figure
#c1 = [[1, 0, 0]] 
#c2 = [[2, 0, 0]] 
#-------------------------------
C1 = np.array(c1)
C2 = np.array(c2)
print ">> writing xyz files for c1 & c2"
util.write_xyz2(c1,c2,"cluster_overlap")
util.write_mathematica(c1,"c1")
util.write_mathematica(c2,"c2")

sigmas = [0.1,0.5,1.0,2.0]
#sigmas = [0.5]
#sigmas = [0.0001,1000.0]
maxorder = 8

######################################################################
def order_convergence():
   print ">> running order convergence ",
   sys.stdout.flush()
   f = open("order.dat","w")
   print >>f,"# order convergence p, soap:sigma=",sigmas
   ds = []
   for sig in sigmas:
     d12,R12 = o.distance(c1,c2,sig)
     ds.append(d12)
   for p in range(maxorder):
     sys.stdout.write("*")
     sys.stdout.flush()
     print >>f,"{0:4d}".format(p),
     sds = []
     i = 0
     for sig in sigmas:
       d12,R12 = d.distance(c1,c2,p,sig)
       print >>f,"{0:10f}".format(d12),
       sds.append(d12/ds[i])
       i += 1
     for sd in sds:
       print >>f,"{0:10f}".format(sd),
     print >>f
   print >>f
   print >>f
   print >>f,"# ogto:\n inf",
   for od in ds:
     print >>f,"{0:10f}".format(od),
   print >>f
   f.close()
   print "done"

######################################################################
def inner_product_test():
   print ">> running inner product test ",
   n = 10
   alphas = []
   cs = []
   for i in range(n+1):
     a = 1.0*i/n
     c = np.array(util.interpolate(c1,c2,a))
     alphas.append(a)
     cs.append(c)
   for i in range(len(alphas)):
     a = alphas[i]
     p1 = o.inner_product(C1,cs[i]) 
     p2 = o.inner_product(C2,cs[i]) 
     print "{0:8f} {1:10f} {2:10f}".format(a,p1,p2)

######################################################################
def triangle_test():
   print ">> running triangle test ",
   sys.stdout.flush()
   f = open("triangle.dat","w")
   p = 4
   print >>f,"# triangle test "
   #print >>f,"# sph harmonic order ",p
   n = 10 # 40
   alphas = []
   cs = []
   for i in range(n+1):
     a = 1.0*i/n
     c = np.array(util.interpolate(c1,c2,a))
     alphas.append(a)
     cs.append(c)
   g0,R0 = g.distance(C1,C2)
   print >>f,"# RMSD ref =",g0
   print >>f,"# OGTO ref =",
   o0s = []
   for sig in sigmas: 
     d0,R0 = o.distance(c1,c2,sig)
     #d0,R0 = o.dist(c1,c2,sig)
     o0s.append(d0)
     print >>f,d0,
   print >>f
#  print >>f,"# SOAP ref =",
#  d0s = []
#  for sig in sigmas: 
#    d0,R0 = d.distance(c1,c2,p,sig)
#    d0s.append(d0)
#    print >>f,d0,
#  print >>f
   #print >>f,"# alpha (d1,d2)_rmsd, _ogto @, _soap @ sig=",
   print >>f,"# alpha (d1,d2)_rmsd, _ogto @ sig=",
   for sig in sigmas: print >>f, sig,
   print >>f
   for i in range(len(alphas)):
     sys.stdout.write("*")
     sys.stdout.flush()
     a = alphas[i]
     d1,R1 = g.distance(C1,cs[i]) 
     d2,R2 = g.distance(C2,cs[i]) 
     print "RMSD\n",R1,"\n",R2
     print >>f,"{0:8f} {1:10f} {2:10f}".format(a,d1/g0,d2/g0),
     for j in range(len(sigmas)):
       sig = sigmas[j]
       d0 = o0s[j]
       d1,R1 = o.distance(c1,cs[i],sig) 
       d2,R2 = o.distance(c2,cs[i],sig)
       #d1,R1 = o.dist(c1,cs[i],sig) 
       #d2,R2 = o.dist(c2,cs[i],sig)
       print "OGTO\n",R1,"\n",R2
       print >>f,"{0:10f} {1:10f}".format(d1/d0,d2/d0),
#    for j in range(len(sigmas)):
#      sig = sigmas[j]
#      d0 = d0s[j]
#      d1,R1 = d.distance(c1,cs[i],p,sig) 
#      d2,R2 = d.distance(c2,cs[i],p,sig)
#      print >>f,"{0:10f} {1:10f}".format(d1/d0,d2/d0),
     print >>f
   print "done"
   
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
   for i in range(n+1):
     a = 2.0*math.pi*i/n
     R = util.rotation([1,0,0],a)
     c = np.array(util.rotate_all(R,c1))
     alphas.append(a)
     cs.append(c)
   g0,R0 = g.distance(C1,C2)
   print >>f,"# RMSD ref =",g0
   print >>f,"# OGTO ref =",
   o0s = []
   for sig in sigmas: 
     d0,R0 = o.distance(c1,c2,sig)
     #d0,R0 = o.dist(c1,c2,sig)
     o0s.append(d0)
     print >>f,d0,
   print >>f
   print >>f,"# alpha (d1,d2)_rmsd, _ogto @ sig=",
   for sig in sigmas: print >>f, sig,
   print >>f
   for i in range(len(alphas)):
     sys.stdout.write("*")
     sys.stdout.flush()
     a = alphas[i]
     d1,R1 = g.distance(C1,cs[i]) 
     d2,R2 = g.distance(C2,cs[i]) 
     print "RMSD\n",R1,"\n",R2
     print >>f,"{0:8.5f} {1:10.7f} {2:10.7f}".format(a,d1/g0,d2/g0),
     for j in range(len(sigmas)):
       sig = sigmas[j]
       d0 = o0s[j]
       d1,R1 = o.distance(c1,cs[i],sig) 
       d2,R2 = o.distance(c2,cs[i],sig)
       #d1,R1 = o.dist(c1,cs[i],sig) 
       #d2,R2 = o.dist(c2,cs[i],sig)
       print "OGTO\n",R1,"\n",R2
       print >>f,"{0:10.7f} {1:10.7f}".format(d1/d0,d2/d0),
     print >>f
   print "done"
####################################################
def permutation_test():
   print "\n# permutation test "
   s2p = [s2[1],s2[2],s2[3],s2[0]]
   d12,R12   = d.distance(s1,s2,p)
   dp,Rp = d.distance(s1,s2,p)
   print "d_12=",d12
   print "d_12=",dp
   print "R_12=\n",R12
   print "R_12=\n",Rp

   print "# d vs p "
   for sig in sigmas:
     print "# sigma=",sig
     for p in range(10):
       #d11,R = d.distance(s1,s1,p)
       d12,R = d.distance(s1,s2,p,sig)
       print p,sig,d12
       #d22 = compute_distance(s2,s2,p)
       #print p,":d_12= ",dR[0] #, " d_11 ", d11, " d_22 ",d22
     print
   print "done"
####################################################
if  __name__ == "__main__":
 inner_product_test()
 triangle_test()
 rotation_test()
 #order_convergence()
