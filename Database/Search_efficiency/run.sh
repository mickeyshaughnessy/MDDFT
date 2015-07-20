#!/bin/bash
#refs="4:8:16:32:64"
refs="4:8:16"
nqueries=8
ntrials=8
db0="test.db"
dbs=""
for i in 100 200 400 800; do
  db="$i.db"
  head -n$i $db0 > $i.db
  dbs="$dbs $db"
done
#echo $dbs
for db in $dbs; do
  #q=${db//\.db/.q}
  ./run_search.py $db $nqueries $refs $ntrials 
done

cat size=*_average.dat | sort -n -k 2 -k 1 > evaluations.dat
