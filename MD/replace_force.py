#!/usr/bin/env python
"""
replaces force in cluster db with force from dmp file 
needed for 3-body potentials where "newton on" is required in LAMMPS

"""
import sys
import re

dbfile   = open(sys.argv[1])
dmpfile  = open(sys.argv[2])

step = 0
cluster = 1
for line in dbfile:
  line2 = dmpfile.readline()
  if (re.search("ITEM",line2)):
    #print "step ",step, " cluster ", cluster
    step += 1
    for i in range(8): dmpfile.readline()
    line2 = dmpfile.readline()
  cols = line2.split()
  fx = cols[2]
  fy = cols[3]
  fz = cols[4]
  cols = line.split()
  cols[2] = fx
  cols[3] = fy
  cols[4] = fz
  #print "   cluster ", cluster
  cluster += 1
  for c in cols:
    print c,
  print
