#!/bin/bash
#tail -n216 ../../DFT/AIMD/Silicon/300K/Si-4_216x200.db   >  300K-Si-4.db
#tail -n216 ../../DFT/AIMD/Silicon/2500K/Si-4_216x200.db  > 2500K-Si-4.db
#tail -n216 ../../DFT/AIMD/Silicon/300K/Si-16_216x200.db   >  300K-Si-16.db
#tail -n216 ../../DFT/AIMD/Silicon/2500K/Si-16_216x200.db  > 2500K-Si-16.db
#ln -sf ../../DFT/AIMD/Silicon/300K/Si-4_216x200.db     300K-Si-4_216x200.db
#ln -sf ../../DFT/AIMD/Silicon/300K/Si-16_216x200.db   300K-Si-16_216x200.db
#ln -sf ../../DFT/AIMD/Silicon/2500K/Si-4_216x200.db   2500K-Si-4_216x200.db
#ln -sf ../../DFT/AIMD/Silicon/2500K/Si-16_216x200.db 2500K-Si-16_216x200.db

#ln -sf ../../DFT/AIMD/Silicon/300K/Si-4_216x200.rmsd     300K-Si-4_216x200.rmsd
#ln -sf ../../DFT/AIMD/Silicon/300K/Si-16_216x200.rmsd   300K-Si-16_216x200.rmsd
#ln -sf ../../DFT/AIMD/Silicon/2500K/Si-4_216x200.rmsd   2500K-Si-4_216x200.rmsd
#ln -sf ../../DFT/AIMD/Silicon/2500K/Si-16_216x200.rmsd 2500K-Si-16_216x200.rmsd

n=1
../path_db.py ../../DFT/AIMD/Silicon/300K/Si-4_216x200.db   $n >  300K-Si-4.db
../path_db.py ../../DFT/AIMD/Silicon/300K/Si-16_216x200.db  $n >  300K-Si-16.db
../path_db.py ../../DFT/AIMD/Silicon/2500K/Si-4_216x200.db  $n >  2500K-Si-4.db
../path_db.py ../../DFT/AIMD/Silicon/2500K/Si-16_216x200.db $n >  2500K-Si-16.db

dmax=0.5
dat="force-distance.dat"
rm $dat; touch $dat
for db in 300K-Si-4.db 300K-Si-16.db 2500K-Si-4.db 2500K-Si-16.db ; do 
  ../fd_corr.py $db $dmax
  echo
  echo
  out=${db%.db}
  cor="${out}.fdcorr"
  echo "# $out" >> $dat
  cat $cor >> $dat
  echo >> $dat
  echo >> $dat
done
