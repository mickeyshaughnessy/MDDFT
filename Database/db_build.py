#!/usr/bin/env python
"""
"""
import os
import sys
import database
import util
import numpy as np
import math


#####################################################################
## MAIN
#####################################################################
if  __name__ == "__main__":
  # parse
  if (len(sys.argv) < 2):
    print "usage: db_build.py run.db <run.dmp>"
    sys.exit()
  dbfile = sys.argv[1]
  nrefs = int(sys.argv[2])
  db = database.database()
  db.init(dbfile,"RMSD",nrefs)
  if (len(sys.argv) > 3):
    ffile = sys.argv[3]
    db.Forces = db.read_dump_forces(ffile)
  db.write_sql()
  #sqdb = dbfile[:-2]+"sql"
  #sqdmp = sqdb+"_dmp"
  #os.system("sqlite3 "+sqdb+" .dump > "+sqdmp);
