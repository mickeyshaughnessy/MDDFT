#!/usr/bin/env python 
"""
  function: 
  output: 
  usage: 
  to do: 
    handle steps in unique writes by tagging
    step selection
"""
import sys
import os
import math
import re
import numpy as np

def ADF(self):
    """ prints ADF function to stdout"""
    angs = {}  #histogram of angles
    for i in xrange(180):
        angs[i] = 0
    read_dump(sys.argv[1])
    for atom in self.positions:
        print atom
    
    print 'Hello'
    # 1. for each atom, compute its adf as
    #    sum_j cos^-1(x1 dot x2 / |x1||x2|) 

###################################################################
def read_dump(filename):
    """ reads a dump"""
    natoms = 0
    positions = []
    forces = []
    type_counts = []
    f = open(filename)
    lines = f.readlines()
    i = 0
    pattern = re.compile("ITEM: NUMBER OF ATOMS")
    while not pattern.search(lines[i]) : i+=1
    i += 1
    natoms = int((lines[i].split())[0])
    pattern = re.compile("ITEM: BOX BOUNDS")
    while not pattern.search(lines[i]) : i+=1
    i += 1
    cols = (lines[i].split())
    xlo = float(cols[0])
    xhi = float(cols[1])
    i += 1
    cols = (lines[i].split())
    ylo = float(cols[0])
    yhi = float(cols[1])
    i += 1
    cols = (lines[i].split())
    zlo = float(cols[0])
    zhi = float(cols[1])
    Lx = xhi-xlo
    Ly = yhi-ylo
    Lz = zhi-zlo
    b1=[Lx,0.,0.]
    b2=[0.,Ly,0.]
    b3=[0.,0.,Lz]
    f.close()
    # read configurations
    nsteps = 0
    pattern = re.compile("ITEM: ATOMS")
    f = open(filename)
    atoms = False
    n = 0
    freq = 1
    for line in f:

    if atoms :
    cols = line.split()
    type = cols[1]
    self.types.append(type)
    x = float(cols[2])-xlo
    y = float(cols[3])-ylo
    z = float(cols[4])-zlo
    self.positions.append([x,y,z])
    fx = float(cols[2])-xlo
    fy = float(cols[3])-ylo
    fz = float(cols[4])-zlo
    self.forces.append([fx,fy,fz])
    n += 1
    if n == self.natoms :
      self.nsteps += 1 
      self.write()
      # reset
      atoms = False
      n = 0
      self.types = []
      self.positions = []
      self.forces = []
      if (self.nsteps % freq == 0): 
        print "* read step ",self.nsteps
  elif pattern.search(line):
    atoms = True
f.close()


### main ########################################################
if __name__ == '__main__':
    
    infile = sys.argv[1]
    outfile = sys.argv[2]


