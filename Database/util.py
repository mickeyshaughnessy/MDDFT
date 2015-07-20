#!/usr/bin/env python

"""
utilities
"""
import sys
import numpy as np
import math
import array
import random
import gzip
from operator import itemgetter

###################################################################
_big = 1.e20
_small = 1.e-20

###################################################################
def fopen(name,mode="r"):
  if name[-3:] == ".gz":
    return gzip.open(name,mode)
  else:
    return open(name,mode)

###################################################################
def write_mathematica(c,name="c"):
  f = open(name+".mma","w")
  print >>f,to_mathematica(c,name)
  f.close()

###################################################################
def to_mathematica(c,name="c"):
  s = name+"= {",
  for x in c[:-1]:
   s+="{"+str(x[0])+", "+str(x[1])+", "+str(x[2])+"},",
  x = c[-1]
  s+="{"+str(x[0])+", "+str(x[1])+", "+str(x[2])+"}};",
  return s
  
###################################################################
def write_xyz(c,name="c"):
  f = open(name+".xyz","w")
  print >>f,len(c)+1
  print >>f,name
  print >>f,"C {0:10f} {1:10f} {2:10f}".format(0,0,0)
  for x in c:
    print >>f,"C {0:10f} {1:10f} {2:10f}".format(x[0],x[1],x[2])
  f.close()

###################################################################
def write_xyz2(c1,c2,name="c"):
  f = open(name+".xyz","w")
  print >>f,len(c1)+len(c2)+1
  print >>f,name
  print >>f,"O {0:10f} {1:10f} {2:10f}".format(0,0,0)
  for x in c1:
    print >>f,"C {0:10f} {1:10f} {2:10f}".format(x[0],x[1],x[2])
  for x in c2:
    print >>f,"P {0:10f} {1:10f} {2:10f}".format(x[0],x[1],x[2])
  f.close()

###################################################################
def interpolate(c1, c2, a):
  c = []
  for i in range(len(c1)):
    x = a*c1[i][0]+(1.-a)*c2[i][0]
    y = a*c1[i][1]+(1.-a)*c2[i][1]
    z = a*c1[i][2]+(1.-a)*c2[i][2]
    c.append([x,y,z])
  return c

###################################################################
def kernel_density_estimate(ys,name=None,nx=0,xrange=None):
  if nx == 0: nx = 100
  if xrange == None:
    xrange = (min(ys),1.1*max(ys))
  n = len(ys)
  w = n**(-1./(1+4))
  h = 1./w
  s = 1./(math.sqrt(2.*math.pi)*n*w)
  xs = []
  zs = []
  print "kde range:",xrange, " width:",w," nsamples: ",n
  for x in np.arange(xrange[0],xrange[1],(xrange[1]-xrange[0])/nx):
    z = 0.
    for y in ys:
       dx = (x-y)*h
       z += math.exp(-0.5*dx*dx)
    z *= s
    zs.append(z)
    xs.append(x)
  if (name!=None):
    file = open(name,"w")
    for i in range(len(xs)):
      print >>file, xs[i],zs[i]
  return [xs,zs]

###################################################################
def flatten(lists):
  return [item for sublist in lists for item in sublist]

###################################################################
def histogram(ds,name=None,drange=None,nbins=0): # drange = [min,max]
  #counts,partitions = np.histogram(ds,new=True)
  if (nbins==0): 
    if (len(ds) > 1000) : 
      nbins = min(int(len(ds)/100),100)
    else:
      nbins = 10
  if (drange == None):
    counts,partitions = np.histogram(ds,bins=nbins,normed=True)
  else :
    counts,partitions = np.histogram(ds,bins=nbins,normed=True,range=drange)
  max_count = max(counts)
  if (name==None):
    for i in range(len(counts)):
      print '%3f to %3f :' % (partitions[i],partitions[i+1]),str(counts[i]).ljust(5),
      s = int(10*counts[i]/max_count)
      print s*"*",
      print
  else:
    file = open(name,"w")
    sum = 0
    for i in range(len(counts)):
      sum += counts[i]
    s = 1./sum
    for i in range(len(counts)):
      print >>file,"%10f %10f %10d %10f %10f"% ( 0.5*(partitions[i]+partitions[i+1]),counts[i]*s,counts[i],partitions[i],partitions[i+1] )

###################################################################
def read_distances(filename,drange,nmax):
  dmin = drange[0]
  dmax = drange[1]
  freq = int(nmax/10)
  distances = nmax*[0] # clear
  f = open(filename)
  k = 0
  print "> reading ",
  sys.stdout.flush()
  for line in f:
    cols = line.split()
    d = float(cols[2])
    if (d >= dmin and d <= dmax): 
      distances[k] =  float(cols[2])
      k += 1
      if (k % freq == 0):
        print "*",
        sys.stdout.flush()
      if k >= nmax:
        break
  print k," distances"
  f.close()
  ds = np.array(distances[0:k])
  return ds

###################################################################
def load_distances(ds,filename):
  return np.save(filename)

###################################################################
def save_distances(ds,filename):
  np.save(filename,ds)

#####################################################################
def rotation(p,theta):
  # rodroguez formula: R = I + sin K + (1-cos) K^2
  R = np.eye(3)
  axis = normalize(p)
  k1 = axis[0]
  k2 = axis[1]
  k3 = axis[2]
  K = np.array([[0.,-k3,k2],[k3,0.,-k1],[-k2,k1,0.]])
  s = math.sin(theta)
  R += s*K
  K2 = np.dot(K,K)
  c = 1.-math.cos(theta)
  R += c*K2
  return R

#####################################################################
def rotate(R,f):
  Rf = np.dot(R,f)
  return Rf

#####################################################################
def rotate_all(R,C):
  RC = []
  for c in C:
    RC.append(np.dot(R,c))
  return RC

#####################################################################
def add(a,b):
  c = a
  for i in range(len(a)):
    c[i] += b[i]
  return c

#####################################################################
def normalize(v):
  s = 1./norm(v)
  return scale(s,v)

#####################################################################
def scale(s,v):
  w = v
  for i in range(len(v)):
    w[i] *= s
  return w

#####################################################################
def error_sq(a,b):
  diff =  len(a)*[0.] 
  for i in range(len(a)):
    diff[i] = a[i]-b[i]
  error2 = np.dot(diff,diff)
  return error2

#####################################################################
def error(a,b):
  return math.sqrt(error_sq(a,b))

#####################################################################
def norm(a):
  sq = 0.
  for v in a:
    sq += v*v
  return math.sqrt(sq)

#####################################################################
def d(a,b):
  sq = 0.
  for i in range(len(a)):
    v = a[i]-b[i]
    sq += v*v
  return math.sqrt(sq)

#####################################################################
def relative_error(a,b):
  return error(a,b)/norm(b)

#####################################################################
def moments(C,ms): ## WIP
  n = len(C)
  moments = len(ms)*[0.]
  for x in C:
    for m in ms:
      if m == "0":
        moments[0] += 1
  return moments

#####################################################################
def central_moments(C,ms): ## WIP
  n = len(C)
  moments = len(ms)*[0.]
  for x in C:
    for m in ms:
      if m == "0":
        moments[0] += 1
  return moments

#####################################################################
#def principal_axes(C,ms): ## WIP

#####################################################################
## MAIN
#####################################################################
if  __name__ == "__main__":
  print "!!! TEST !!!" 
