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
def read_dump(filename, outfile):
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
            x = float(cols[2])-xlo
            y = float(cols[3])-ylo
            z = float(cols[4])-zlo
            positions.append([x,y,z])
            fx = float(cols[2])-xlo
            fy = float(cols[3])-ylo
            fz = float(cols[4])-zlo
            forces.append([fx,fy,fz])
            n += 1
            if n == natoms :
                compute_ADF(positions, outfile)
                nsteps += 1 
                # reset
                atoms = False
                n = 0
                types = []
                positions = []
                forces = []
        elif pattern.search(line):
            atoms = True
    f.close()

def compute_ADF(positions, outfile):
    print 'Computing ADF'
    """ prints ADF function to file""" 
    angs = {}  #histogram of angles
    for i in xrange(180):
        angs[i] = 0
    for i, atom in enumerate(positions):
        print 'atom #:' % i
        for i in positions:
            v1 = atom[0] - i[0], atom[1] - i[1], atom[2] - i[2]
            for j in positions:
                v2 = atom[0] - j[0], atom[1] - j[1], atom[2] - j[2]
                if v1 !=(0,0,0) and v2 !=(0,0,0):
                    vects = (dot(v1,v2) / (mag(v1) * mag(v2)))
                    angle = np.arccos(trim(vects))
                    angle = int(angle * 57.2957795) # convert to degrees
                    angs[angle] += 1
    with open(outfile, 'w') as f:
        for ang in angs:
            f.write('%s %s \n' % (ang, angs[ang])) 
        f.write('\n')

def mag(x):
    return np.sqrt(sum([i*i for i in x]))
def dot(x,y):
    return sum([X*Y for X,Y in zip(x,y)])
def trim(x):
    if x < -1: return -1
    elif x > 1.0: return 1.0
    else: return x
### main ########################################################
if __name__ == '__main__':
    
    infile = sys.argv[1]
    outfile = sys.argv[2]
    read_dump(infile, outfile) 

