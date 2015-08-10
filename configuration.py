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
#import util

#################################################################
class configuration:
#################################################################
  atomic_masses = {
    "H": 1.0079 ,
    "He": 4.0026 ,
    "Li": 6.941 ,
    "Be": 9.0122 ,
    "B": 10.811 ,
    "C": 12.0107 ,
    "N": 14.0067 ,
    "O": 15.9994 ,
    "F": 18.9984 ,
    "Ne": 20.1797 ,
    "Na": 22.9897 ,
    "Mg": 24.305 ,
    "Al": 26.9815 ,
    "Si": 28.0855 ,
    "P": 30.9738 ,
    "S": 32.065 ,
    "Cl": 35.453 ,
    "Ar": 39.948 ,
    "K": 39.0983 ,
    "Ca": 40.078 ,
    "Sc": 44.9559 ,
    "Ti": 47.867 ,
    "V": 50.9415 ,
    "Cr": 51.9961 ,
    "Mn": 54.938 ,
    "Fe": 55.845 ,
    "Co": 58.9332 ,
    "Ni": 58.6934 ,
    "Cu": 63.546 ,
    "Zn": 65.39 ,
    "Ga": 69.723 ,
    "Ge": 72.64 ,
    "As": 74.9216 ,
    "Se": 78.96 ,
    "Br": 79.904 ,
    "Kr": 83.8 ,
    "Rb": 85.4678 ,
    "Sr": 87.62 ,
    "Y": 88.9059 ,
    "Zr": 91.224 ,
    "Nb": 92.9064 ,
    "Mo": 95.94 ,
    "Tc": 98. ,
    "Ru": 101.07 ,
    "Rh": 102.9055 ,
    "Pd": 106.42 ,
    "Ag": 107.8682 ,
    "Cd": 112.411 ,
    "In": 114.818 ,
    "Sn": 118.71 ,
    "Sb": 121.76 ,
    "Te": 127.6 ,
    "I": 126.9045 ,
    "Xe": 131.293 ,
    "Cs": 132.9055 ,
    "Ba": 137.327 ,
    "La": 138.9055 ,
    "Ce": 140.116 ,
    "Pr": 140.9077 ,
    "Nd": 144.24 ,
    "Pm": 145 ,
    "Sm": 150.36 ,
    "Eu": 151.964 ,
    "Gd": 157.25 ,
    "Tb": 158.9253 ,
    "Dy": 162.5 ,
    "Ho": 164.9303 ,
    "Er": 167.259 ,
    "Tm": 168.9342 ,
    "Yb": 173.04 ,
    "Lu": 174.967 ,
    "Hf": 178.49 ,
    "Ta": 180.9479 ,
    "W": 183.84 ,
    "Re": 186.207 ,
    "Os": 190.23 ,
    "Ir": 192.217 ,
    "Pt": 195.078 ,
    "Au": 196.9665 ,
    "Hg": 200.59 ,
    "Tl": 204.3833 ,
    "Pb": 207.2 ,
    "Bi": 208.9804 ,
    "Po": 209. ,
    "At": 210. ,
    "Rn": 222. ,
    "Fr": 223. ,
    "Ra": 226. ,
    "Ac": 227. ,
    "Th": 232.0381 ,
    "Pa": 231.0359 ,
    "U": 238.0289 
    }

  # box data
  alat = 1.
  b1 = [0.,0.,0.]
  b2 = [0.,0.,0.]
  b3 = [0.,0.,0.]
  # per atom data
  natoms = 0
  positions = []
  types = []
  masses = []
  forces = []
  # per type data
  ntypes = 0
  type_counts = []
  element_map = []
  type_map = []
  mass_map = []
  # per type data
  nsteps = 0
  source = ""
  outfile = ""
  neighborhoods = []
  rs = []

  def __init__(self):
    pass

###################################################################
  def init(self):
    # construct missing data
    unique_types = set()
    for t in self.types:
      unique_types.add(t)
    self.ntypes = len(unique_types)
    type_counts = self.ntypes*[0]
    #print self.types
    #print self.ntypes
    for t in self.types:
      type_counts[int(t)-1] += 1
    ## NOTE use a type to amu map
    if (len(self.masses) == 0):
      n = self.natoms
      self.masses = n*[1.]

###################################################################
  def read_xyz(self,filename):
    f = open(filename)
    self.natoms  = int(((f.readline()).split())[0])
    print "> reading ",self.natoms,"atoms",
    f.readline()
    self.unique_types.clear()
    self.positions = []
    self.types = []
    for i in range(n):
      cols = (f.readline()).split()
      type = cols[0]
      self.types.append(type)
      x = float(cols[1])
      y = float(cols[2])
      z = float(cols[3])
      self.positions.append([x,y,z])
    self.natoms = len(self.positions)
    self.write()
  
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
  def write_xyz(self,filename):
    f = open(filename,"w")
    print >>f,self.natoms
    print >>f,filename,self.alat
    k = 0
    for xyz in self.positions:
      type = self.types[k]
      x = xyz[0]
      y = xyz[1]
      z = xyz[2]
      k += 1
      print >>f, type,x,y,z
    f.close()

#################################################################
  def read_poscar(self, filename="POSCAR" ) :
    f = open(filename)
    lines = f.readlines()
    # lattice oonstant
    self.alat = float(lines[1])
    # periodic cell
    self.b1 = [0,0,0]
    items = (lines[2]).split()
    self.b1[0] = self.alat*float(items[0])
    self.b1[1] = self.alat*float(items[1])
    self.b1[2] = self.alat*float(items[2])
    items = (lines[3]).split()
    self.b2[0] = self.alat*float(items[0])
    self.b2[1] = self.alat*float(items[1])
    self.b2[2] = self.alat*float(items[2])
    items = (lines[4]).split()
    self.b3[0] = self.alat*float(items[0])
    self.b3[1] = self.alat*float(items[1])
    self.b3[2] = self.alat*float(items[2])
    # types
    items = (lines[5]).split()
    self.ntypes = len(items)
    self.type_counts = self.ntypes*[0]
    self.natoms = 0
    t = 1
    for i in range(self.ntypes):
      n = int(items[i])
      self.type_counts[i] = n
      for j in range(self.natoms,self.natoms+n):
        self.type.append(t)
      t += 1
      self.natoms += n
    # atoms
    self.positions = []
    style = lines[6].rstrip("\n")
    alines = lines[7:]
    for i in range(self.natoms):
      items = (alines[i]).split()
      s1 = float(items[0])
      s2 = float(items[1])
      s3 = float(items[2])
      if (style == "Cartesian"):
        x = s1
        y = s2
        z = s3
      else :
        x = s1*self.b1[0]+s2*self.b2[0]+s3*self.b3[0]
        y = s1*self.b1[1]+s2*self.b2[1]+s3*self.b3[1]
        z = s1*self.b1[2]+s2*self.b2[2]+s3*self.b3[2]
      self.positions.append([x,y,z])
    f.close()
    self.write()

#################################################################
  def write_poscar(self, filename="POSCAR" ) :
    f = open(filename,"w")
    print >>f, self.source
    # lattice oonstant
    print >>f, self.alat 
    # periodic cell
    print >>f, self.b1[0]/self.alat, \
               self.b1[1]/self.alat, \
               self.b1[2]/self.alat
    print >>f, self.b2[0]/self.alat, \
               self.b2[1]/self.alat, \
               self.b2[2]/self.alat
    print >>f, self.b3[0]/self.alat, \
               self.b3[1]/self.alat, \
               self.b3[2]/self.alat
    # types
    for c in self.type_counts :
       print >>f, c,
    print >>f
    #print >>f,"Direct"
    print >>f,"Cartesian"
    # atoms
    for xyz in self.positions:
      print >>f, xyz[0],xyz[1],xyz[2]
    f.close()

#################################################################
  def read_data(self, filename) :
    """ reads LAMMPS data files"""
    f = open(filename)
    lines = f.readlines()
    i = 0
    pattern = re.compile("atoms")
    while not pattern.search(lines[i]) : i+=1
    self.natoms = int((lines[i].split())[0])
    pattern = re.compile("atom types")
    while not pattern.search(lines[i]) : i+=1
    self.ntypes = int((lines[i].split())[0])
    pattern = re.compile("xlo")
    while not pattern.search(lines[i]) : i+=1
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
    self.b1=[Lx,0.,0.]
    self.b2=[0.,Ly,0.]
    self.b3=[0.,0.,Lz]
    pattern = re.compile("Masses")
    while not pattern.search(lines[i]) : i+=1
    i += 2
    for l in lines[i:i+self.ntypes]: # NOTE FRAGILE use re
      cols = l.split()
      self.masses.append(float(cols[1]))
      i += 1
    pattern = re.compile("Atoms")
    while not pattern.search(lines[i]) : i+=1
    i += 2
    for l in lines[i:i+self.natoms]:
      cols = l.split()
      type = cols[1] # NOTE str or int
      x = float(cols[2])-xlo
      y = float(cols[3])-ylo
      z = float(cols[4])-zlo
      self.types.append(type)
      self.positions.append([x,y,z])
    self.type_counts.append(self.natoms) # HACK
    f.close()
    self.write()

#################################################################
  def write_data(self, filename) :
    """ write LAMMPS data file format to file <filename>"""
    f = open(filename,"w")
    print >>f,filename
    print >>f
    print >>f, self.natoms," atoms"
    print >>f 
    print >>f, self.ntypes," atom types"
    print >>f 
    print >>f, 0,self.b1[0],"xlo xhi"
    print >>f, 0,self.b2[1],"ylo yhi"
    print >>f, 0,self.b3[2],"zlo zhi"
    print >>f 
    print >>f, "Masses"
    print >>f 
    for i in range(self.ntypes):
      print >>f, (i+1), 1.
    print >>f 
    print >>f, "Atoms"
    print >>f 
    i = 0
    for j in range(self.ntypes):
      for k in range(self.type_counts[j]):
        print >>f, (i+1),(j+1),self.positions[i][0], \
                               self.positions[i][1], \
                               self.positions[i][2]
        i+=1
    f.close()

#################################################################
  def write_gnuplot(self, filename) :
    status = "w"
    if self.nsteps > 0 : status = "a" # append
    f = open(filename,status)
    hasForces = False
    if (len(self.forces) > 0): hasForces = True
    if self.nsteps == 0:
      print >>f,"#",self.alat
      print >>f,"#",self.b1
      print >>f,"#",self.b2
      print >>f,"#",self.b3
    print >>f,"# step:",self.nsteps
    for i in range(self.natoms):
      type = self.types[i]
      print >>f,"{0:3s} ".format(str(type)),
      xyz = self.positions[i]
      x = xyz[0]
      y = xyz[1]
      z = xyz[2]
      print >>f,"{0:18g} {1:18g} {2:18g} ".format(x,y,z),
      if hasForces:
        force = self.forces[i]
        fx = force[0]
        fy = force[1]
        fz = force[2]
        print >>f,"{0:18g} {1:18g} {2:18g} ".format(fx,fy,fz),
      print >>f
    print >>f
    f.close()

###################################################################
  def read_dump(self,filename,step=-1):
    """ reads a dump"""
    self.natoms = 0
    self.positions = []
    self.forces = []
    self.type_counts = []
    f = open(filename)
    lines = f.readlines()
    i = 0
    pattern = re.compile("ITEM: NUMBER OF ATOMS")
    while not pattern.search(lines[i]) : i+=1
    i += 1
    self.natoms = int((lines[i].split())[0])
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
    self.b1=[Lx,0.,0.]
    self.b2=[0.,Ly,0.]
    self.b3=[0.,0.,Lz]
    f.close()
    # read configurations
    self.nsteps = 0
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


###################################################################
  def read_outcar(self,filename="OUTCAR",step=-1):
    f = open(filename)
    self.natoms = 0
    self.positions = []
    self.forces = []
    self.type_counts = []
    type_pattern = re.compile("ions per type =")
    force_pattern = re.compile("TOTAL-FORCE")
    box1_pattern = re.compile("A1 = ")
    box2_pattern = re.compile("A2 = ")
    box3_pattern = re.compile("A3 = ")
    read_next = False
    read_forces = False
    n = 0
    self.nsteps = 0
    freq = 1
    for line in f:
      if read_forces :
        read_next = False
        cols = line.split()
        x = float(cols[0])
        y = float(cols[1])
        z = float(cols[2])
        fx = float(cols[3])
        fy = float(cols[4])
        fz = float(cols[5])
        self.positions.append([x,y,z])
        self.forces.append([fx,fy,fz])
        n += 1
        if (n == self.natoms): 
          self.nsteps +=1
          self.write()
          # reset
          read_forces = False
          n = 0
          self.positions = []
          self.forces = []
          if (self.nsteps % freq == 0): 
            print "* read step ",self.nsteps
      if (read_next): 
        read_forces = True
        config = []
      if type_pattern.search(line):
        cols = line.split()
        self.type_counts = map(int,cols[4:])
        type = 1
        for c in self.type_counts:
          self.natoms += c
          for i in range(c):
            self.types.append(type)
          type += 1
      elif box1_pattern.search(line):
        cols = line.split()
        Lx = float(cols[3].replace(",",""))
        self.b1 = [Lx,0.,0.]
      elif box2_pattern.search(line):
        cols = line.split()
        Ly = float(cols[4].replace(",",""))
        self.b2 = [0.,Ly,0.]
      elif box3_pattern.search(line):
        cols = line.split()
        Lz = float(cols[5].replace(")",""))
        self.b3 = [0.,0.,Lz]
      elif force_pattern.search(line):
        read_next = True
 
###################################################################
  def neighbor(self, nb,R=0):
    L = [self.b1[0],self.b2[1],self.b3[2]]
    l = [0.5*L[0],0.5*L[1],0.5*L[1]]
    if not R > 0: 
      Lmin = 1.e20
      for i in range(3): 
        if L[i] < Lmin: Lmin = L[i]
      R = Lmin
    #print "> using R=",R," in neighboring"
    self.neighborhoods = []
    R2 = R*R 
    print "> neighboring",
    sys.stdout.flush()
    k = 0
    freq = len(self.positions)/10
    for x in self.positions:
      if (k % freq == 0):
        print "*",
        sys.stdout.flush()
      k += 1
      c = []
      for y in self.positions:
        dx = [0.,0.,0.]
        r2 = 0.
        for i in range(3):
          z = y[i]-x[i]
          if (math.fabs(z) > l[i]): 
            if (z > 0) : z -= L[i] 
            else       : z += L[i]
          dx[i] = z 
          r2 += z*z
        if r2 < R2:
          c.append([r2,dx])
        self.rs.append(math.sqrt(r2))
      c.sort(key=lambda z:z[0])
      clist = []
      for i in range(1,nb+1): # omit self
        #print i,c[i][0],c[i][1]
        clist.append(c[i][1])
      cmat = np.array(clist)
      self.neighborhoods.append(cmat.reshape(nb,3))
    print "done"
    #util.kernel_density_estimate(rs,"rdf.dat")
    #util.histogram(rs,"rdf.dat")
    return self.rs


###################################################################
  def write_database(self, filename):
    (tag,nn,natoms,nsteps) = self.parse_db_name(filename)
    self.neighbor(nn)
    if (len(self.forces) == 0): 
      self.forces = self.natoms*[[0.,0.,0.]]
    status = "w"
    if self.nsteps > 0 : status = "a" # append
    f = open(filename,status) 
    k = 0
    for n in self.neighborhoods:
      F = self.forces[k]
      k += 1
      print >>f, "1 11",F[0],F[1],F[2], # HACK type string
      for v in n.flatten():
        print >>f,v,
      print >>f
    f.close()

###################################################################
  def write_vtk(self, filename):
    xlo = str(0.)
    ylo = str(0.)
    zlo = str(0.)
    xhi = str(self.b1[0])
    yhi = str(self.b2[1])
    zhi = str(self.b3[2])
    npts = len(self.positions)
    # box
    o = open("box_"+filename,"w")
    o.write("# vtk DataFile Version 2.0\n")
    o.write(filename+" box\n")
    o.write("ASCII\n")
    o.write("DATASET POLYDATA\n")
    o.write("POINTS 8 float\n")
    o.write(xlo+" "+ylo+" "+zlo+"\n")
    o.write(xhi+" "+ylo+" "+zlo+"\n")
    o.write(xhi+" "+yhi+" "+zlo+"\n")
    o.write(xlo+" "+yhi+" "+zlo+"\n")
    o.write(xlo+" "+ylo+" "+zhi+"\n")
    o.write(xhi+" "+ylo+" "+zhi+"\n")
    o.write(xhi+" "+yhi+" "+zhi+"\n")
    o.write(xlo+" "+yhi+" "+zhi+"\n")
    o.write("DATASET POLYDATA\n")
    o.write("POLYGONS 6 30\n")
    o.write("4 0 1 2 3\n")
    o.write("4 4 5 6 7\n")
    o.write("4 0 1 5 4\n")
    o.write("4 2 3 7 6\n")
    o.write("4 0 4 7 3\n")
    o.write("4 1 2 6 5\n")
    o.close()
    # points
    o = open(filename,"w")
    o.write("# vtk DataFile Version 2.0\n")
    o.write(filename+" atoms\n")
    o.write("ASCII\n")
    o.write("DATASET POLYDATA\n")
    o.write("POINTS "+str(npts)+" float\n")
    for xyz in self.positions:
      x = xyz[0]
      y = xyz[1]
      z = xyz[2]
      o.write(str(x)+" "+str(y)+" "+str(z)+"\n")
    o.write("VERTICES "+str(npts)+" "+str(2*npts)+"\n")
    for i in range(npts):
      o.write("1 "+str(i)+"\n")
    o.write("POINT_DATA "+str(npts)+"\n")
    o.write("SCALARS type int 1\n")
    o.write("LOOKUP_TABLE default\n")
    for t in self.types:
      o.write(str(t)+"\n")

###################################################################
  def parse_db_name(self, filename):
    (tag,remainder) = filename.split("-")
    (nn,remainder)  = remainder.split("_")
    (natoms,remainder)  = remainder.split("x")
    (nsteps,remainder)  = remainder.split(".")
    return (tag,int(nn),int(natoms),int(nsteps))

#################################################################
  def has_steps(self, filename) :
    if self.style(filename) == "OUTCAR": return True
    if self.style(filename) == "dmp":    return True
    return False

#################################################################
  def can_have_steps(self, filename) :
    if self.style(filename) == "db":    return True
    return False

#################################################################
  def append_step(self, filename) :
    (prefix,suffix) = (self.outfile).split(".")
    name = prefix+"_"+str(self.nsteps)+"."+suffix
    return name

#################################################################
  def style(self, filename) :
    if   filename == "POSCAR": 
      return "POSCAR"
    elif filename == "OUTCAR": 
      return "OUTCAR"
    else :
      suffix = (filename.split("."))[-1]
      if suffix == "data" :
        return "data"
      elif suffix == "xyz" :
        return "xyz"
      elif suffix == "vtk" :
        return "vtk"
      elif suffix == "dat" :
        return "dat"
      elif suffix == "dmp" :
        return "dmp"
      elif suffix == "db" :
        return "db"
      else  :
        return "UNKNOWN"

#################################################################
  def read(self, filename) :
     self.source = filename
     if   (self.style(filename) == "POSCAR"): 
       self.read_poscar(filename)
     elif (self.style(filename) == "OUTCAR"): # forces & steps
       self.read_outcar(filename)
     elif (self.style(filename) == "dmp"):    # forces & steps
       self.read_dump(filename)
     elif (self.style(filename) == "data"):
       self.read_data(filename)
     elif (self.style(filename) == "xyz"):
       self.read_xyz(filename)
     else :
       print "ERROR unknown style:",filename
       #sys.exit()

#################################################################
  def write(self, filename=None) :
     if filename == None: filename = self.outfile
     if (not self.can_have_steps(filename)):
       filename = self.append_step(filename) 
     self.init()
     if   (self.style(filename) == "POSCAR"):
       self.write_poscar(filename)
     elif (self.style(filename) == "data"):
       self.write_data(filename)
     elif (self.style(filename) == "xyz"):
       self.write_xyz(filename)
     elif (self.style(filename) == "vtk"):
       self.write_vtk(filename)
     elif (self.style(filename) == "dat"):
       self.write_gnuplot(filename)
     elif (self.style(filename) == "db"):
       self.write_database(filename)
     else :
       print "ERROR unknown style:",filename
       #sys.exit()

#################################################################
  def convert(self, infile, outfile):
    self.outfile = outfile
    self.read(infile)
 
 
### main ########################################################
if __name__ == '__main__':
  c = configuration()
  if (len(sys.argv) < 3):
    print "usage: configuration.py in.suffix out.suffix"
    print "  input  types: data, xyz, POSCAR, dmp"
    print "  output types: data, xyz, POSCAR, dat, db, vtk"
    #sys.exit()
    
  infile = sys.argv[1]
  outfile = sys.argv[2]

  if (infile[-3:] == ".gz"): 
    os.system("gunzip "+infile)
    infile = infile[:-3]
  print "> converting",infile,"to",outfile 
  c.convert(infile,outfile)

  if c.style(outfile) == "db":
    #print "> radial distribution"
    #util.histogram(c.rs,"rdf.dat")
    nbins = max(min(int(len(c.rs)/100),100),10)
    counts,partitions = np.histogram(c.rs,bins=nbins)
    f = open("rdf.dat","w")
    for i in range(len(counts)):
      print >>f,0.5*(partitions[i]+partitions[i+1]),counts[i]
    f.close()
    print "> wrote rdf.dat"

  #if len(sys.argv) > 3:
  #  """ run adf """
  c.ADF() 
