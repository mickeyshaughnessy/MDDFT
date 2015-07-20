#!/bin/bash
#g++ -lsqlite3 -I /usr/local/include/eigen3/ /usr/local/lib/libgsl.a /usr/local/lib/libgslcblas.a database_builder.cpp -o db_builder
g++ -lsqlite3 -I /usr/local/ database_builder.cpp -o db_builder
