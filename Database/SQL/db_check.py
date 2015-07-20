#!/usr/bin/env python
import sys
import os
try:
  import sqlite3
except:
  print "!!! no sqlite3 module"
  sys.exit()


conn = sqlite3.connect(sys.argv[1])
c = conn.cursor()

for s in  c.execute("SELECT name FROM sqlite_master WHERE type='table'"):
  print s

tables = c.execute("SELECT name FROM sqlite_master WHERE type='table'").fetchall()
print len(tables),"tables"
for table in tables:
  name = table[0] 
  size = c.execute("SELECT count(*) FROM "+name).fetchone()
  print name,size[0]
  
