#!/bin/bash
ln -sf ../../DFT/AIMD/Silicon/300K/Si-4_216x200.db 300K-Si-4_216x200.db
ln -sf ../../DFT/AIMD/Silicon/300K/Si-16_216x200.db 300K-Si-16_216x200.db
ln -sf ../../DFT/AIMD/Silicon/2500K/Si-4_216x200.db 2500K-Si-4_216x200.db
ln -sf ../../DFT/AIMD/Silicon/2500K/Si-16_216x200.db 2500K-Si-16_216x200.db

for db in *.db; do 
  ./trace_path.py $db  216 200 10 0.04
done
