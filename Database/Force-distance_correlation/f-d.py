#!/usr/bin/env python

"""
metric/s from Bartok's SOAP similarity measure
"""
import sys
import math
import os
import glob
import numpy as np
lib_path = os.path.abspath('..')
sys.path.append(lib_path)
import database

####################################################
if  __name__ == "__main__":
  root = "../Databases/200-clusters/"
  #root = "../Databases/2000-clusters"
  parallel = True
  for d in glob.glob(root+"*.db"):
    stem = ((d.split("/"))[-1]).rstrip(".db")
    out = stem+"_rmsd.fd_corr"
    if not os.path.isfile(out):
      db = database.database()
      if (parallel) : db.set_parallel()
      db.init(d,"RMSD")
      cfile = db.write_force_distance_correlation()
      dest = (cfile.split("/"))[-1]
      os.rename(cfile,dest)
    else:
      print out,"exists"
    out = stem+"_ogto_.fd_corr"
    if not os.path.isfile(out):
      db = database.database()
      if (parallel) : db.set_parallel()
      db.init(d,"OGTO") 
      cfile = db.write_force_distance_correlation()
      dest = (cfile.split("/"))[-1]
      os.rename(cfile,dest)
    else:
      print out,"exists"
