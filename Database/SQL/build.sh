#!/bin/bash
#g++ -lsqlite3 -I/usr/local/ -I../src database_builder.cpp -DVERBOSE -o db_builder
if [[ database_builder.cpp -nt db_builder ]]; then
  echo "building db_builder"
  g++ -O -lsqlite3 -I/usr/local/ -I../../src database_builder.cpp -o db_builder
else
  echo "db_builder is up to date"
fi 
exit

if [[ generate_FvsD.cpp -nt fvsd ]]; then
  echo "building fvsd"
  g++ -O           -I/usr/local/ -I../../src generate_FvsD.cpp    -o fvsd
else
  echo "fvsd is up to date"
fi 

