#!/bin/bash
cd ..
make mac-ext -j 4
cd EXT/Database
./build.sh
