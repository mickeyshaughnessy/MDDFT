#!/usr/bin/env python

"""
metric search
"""
# system
import fileinput
import sys
import numpy as np
import math
import array
import random
from operator import itemgetter
# local
import database
import util
import interpolation

#####################################################
class metric_search:
#####################################################
  init_type = "cluster" # "random" 
  N_refs = 0
  Done = [] # flag for already compared
  Refs = [] # index --> id
  Path = []

  def __init__(self): 
    self.Db = database.database()

###################################################################
  def init(self,filename,metric_type="RMSD"):
    self.Db.init(filename,metric_type)
    self.Done = self.Db.N_neighborhoods*[False]
    return self.Db.N_neighborhoods

###################################################################
# initial search : a version of closer & new_refs
# - find a new min and then use triangle inequality
###################################################################
  def closest(self,ref_ids,query_clust):
    d_min = util._big
    id_min = -1
    n = 0
    #for ref in ref_ids:
    while len(ref_ids) > 0:
      ref = ref_ids.pop()
      if not self.Done[ref] :
        [d,R] = self.Db.distance(self.Db.Neighborhoods[ref],query_clust)
        self.Done[ref] = True
        n += 1
        #print ref,":",d
        if (d < d_min):
          id_min = ref
          d_min = d
          ref_ids = self.winnow(id_min,ref_ids,d_min)
    #print id_min,"==>",d_min
    return [id_min,d_min,n]

###################################################################
  def closer(self,d_min,ref_ids,query_clust):
    id_min = -1
    n = 0
    for ref in ref_ids:
      if not self.Done[ref] :
        [d,R] = self.Db.distance(self.Db.Neighborhoods[ref],query_clust)
        self.Done[ref] = True
        n += 1
        #print " ",ref,":",d
        if (d < d_min):
          id_min = ref
          d_min = d
          return [id_min,d_min,n]
    return [id_min,d_min,n]

###################################################################
  def winnow(self, min_id, ref_ids, d_min):
    r = []
    for i in ref_ids:
      idx = self.Db.index(min_id,i)
      d = self.Db.Pair_distances[idx]
      if   (d < 2*d_min):  r.append(i)
    unique_refs = self.uniquify(r)
    #print "winnow:",len(ref_ids),"=>",len(unique_refs)
    return unique_refs

###################################################################
  def new_refs(self, min_id, d_min):
    gid = self.Db.Group_ids[min_id]
    r = []
    for i in range(len(self.Db.Pair_distances)):
      p = self.Db.pair(i)
      d = self.Db.Pair_distances[i]
      #print p,"::",d
      p0,p1 = p
      if   ((p1 == min_id) and d < 2*d_min):  
        if (self.Db.Group_ids[p0] == gid):
          r.append(p0)
      elif ((p0 == min_id) and d < 2*d_min):  
        if (self.Db.Group_ids[p1] == gid):
          r.append(p1)
    unique_refs = self.uniquify(r)
    unique_refs.append(min_id) # match from previous gets to compete 
    return unique_refs

###################################################################
  def uniquify(self,seq): # not order preserving / alternate set(list)
    keys = {}
    for e in seq:
      keys[e] = 1
    return keys.keys()

###################################################################
  def initial_refs(self,N_refs,init_type="random"):
    self.N_refs = N_refs
    refs = []
    if (init_type == "random"):
      r = set() 
      while len(r) < N_refs:
        r.add(random.randint(0,self.Db.N_neighborhoods-1))
      while len(r) > 0:
        refs.append(r.pop())
    elif (init_type == "most_distant"):
      refs = self.Db.most_distant(N_refs)
    elif (init_type == "uniformly_distant"):
      refs = self.Db.winnow(N_refs)
    elif (init_type == "cluster"):
      refs = self.Db.cluster(N_refs)
    # output
    print ">",N_refs,"reference ids:",
    n = min(N_refs,8)
    for i in range(n):
      print refs[i],
    if n < len(refs):
      print "..."
    else:
      print 
    self.Refs = refs
    return refs

######################################################################
  def exhaustive_search(self,Q,tol=0): 
    print "> exhaustive search for ball with radius ",tol
    min_id = -1
    every = 1
    dmin = util._big
    ids= []
    for i in range(self.Db.N_neighborhoods):
      [d,R] = self.Db.distance(self.Db.Neighborhoods[i],Q)
      if (d < dmin):
        ids.append([i,d])
        min_id = i
        dmin = d
        print "{0:4d} id:{1:4d} distance:{2:8g}".format(i,min_id,dmin)
      sys.stdout.flush()
    print '* Exhaustive search: id', min_id, ' distance', dmin,
    print 'with', self.Db.N_neighborhoods, 'operations'
    return [min_id,dmin,ids]

######################################################################
  def search(self,Q,tol=0):  # default
    print "> metric search with tolerance",tol
    refs = self.Refs
    N_comps = 0 
    Path = []
    dmin = util._big
    min_id = -1
    last_id = -2
    self.Done = self.Db.N_neighborhoods*[False]
    # select best ref via an exhaustive search
    [min_id,dmin,evals] = self.closest(refs,Q) 
    N_comps+=evals
    refs = self.new_refs(min_id,dmin)
    print "{0:6d} id:{1:4d}:{4:<4d} distance:{2:8.6f} refs: {3:6d}".format(N_comps,min_id,dmin,len(refs),self.Db.Group_ids[min_id])
    # metric search
    while (dmin > tol and min_id != last_id and len(refs) > 1):
      last_id = min_id
      [min_id,dmin,evals] = self.closer(dmin,refs,Q) 
      N_comps+=evals
      if min_id < 0 : 
        min_id = last_id
        break
      Path.append([min_id,dmin,N_comps])
      refs = self.new_refs(min_id,dmin)
      print "{0:6d} id:{1:4d}:{4:<4d} distance:{2:8.6f} refs: {3:6d}".format(N_comps,min_id,dmin,len(refs),self.Db.Group_ids[min_id])
      sys.stdout.flush()
    print "{0:6d}/{1:<6d} id:{2:4d} distance:{3:8.6f}".format(N_comps,self.Db.N_neighborhoods,min_id,dmin)
    return [min_id,dmin,N_comps,Path]

######################################################################
  def triangle_check(self,id1,id2,q): 
    c1 = self.Db.Neighborhoods[id1]
    c2 = self.Db.Neighborhoods[id2]
    [d12,R] = self.Db.distance(c1,c2)
    [d1q,R] = self.Db.distance(c1,q)
    [d2q,R] = self.Db.distance(c1,q)
    if (id1 != id2):
      if (d12 >= (d1q+d2q)) : return False
    return True

######################################################################
  def interpolate_force(self,q,ids):
    force = np.array([0.,0.,0.])
    sum = 0.
    for idd in ids:
      id = idd[0]
      [d,R] = self.Db.distance(self.Db.Neighborhoods[id],q) # recompute
      w = interpolate.rbf(d)
      f = np.array(self.Db.Forces[id]) 
      Rf = util.rotate(R,f) 
      force += w*np.array(Rf)
      #print f, Rf, w, force
      print w, Rf
      sum += w
    force /= sum
    return force 

#####################################################################
## MAIN
#####################################################################
if  __name__ == "__main__":
  print "!!! TEST !!!" 
  search = metric_search()
  c1 = [[1.1, 1, 1],
        [1, -1, -1],
        [-1, 1, -1],
        [-1, -1, 1]]
  C1 = np.array(c1)
  f = [1,2,3]
  F = np.array(f)
  theta = math.pi/3.
  c = math.cos(theta)
  s = math.sin(theta)
  R = [[c,-s,0],
       [s, c,0],
       [0, 0,1]]
  print "R=\n",np.array(R)
  RC1 = np.transpose(np.dot(R,np.transpose(C1)))
  print "A=  C1 = \n",C1
  print "B=R.C1 = \n",RC1
  RF  = np.dot(R,F)
  print "F @A=  F  = \n",F
  print "F @B=R.F  = \n",RF
  [d,OptRot] = search.Db.distance(C1,RC1)
  print " distance = ",d
  print " Optimal rotation =\n",OptRot
  OptRotTF = util.rotate(np.transpose(OptRot),F)
  print " Force at A rotated to B=\n",OptRotTF,"\n error =",util.error(RF,OptRotTF)
  OptRotRF = util.rotate(OptRot,RF)
  print " Force at B rotated to A=\n",OptRotRF,"\n error =",util.error(F,OptRotRF)
