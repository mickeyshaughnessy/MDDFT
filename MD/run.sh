#!/bin/bash
#echo "!!! cleaning db, sql & data !!!";rm -f *.db *.sql *.data *.rmsd*
## system
os=`uname`
if   [ $os == "Linux" ]; then
  lammps="${HOME}/repo/lammps-atc/src/lmp_linux-ext"
elif [ $os == "Darwin" ]; then
  lammps="${HOME}/repo/lammps-atc/src/lmp_mac-ext"
fi
if [ ! -e $lmps ]; then
  echo "!!! no $lmps exists"
fi
ln -sf $lammps lmps

## usage
if [ $# -eq  0 ]; then
  #echo "usage: run.sh number_of_neighbors init:cell_width init:number_of_steps run:cell_width run:number_of_steps"
  echo "usage: run.sh type nneighbors cell_width nsteps nrefs"
  echo "example: run.sh 200"
  echo "example: run.sh Ar 12 2 200 80"
  echo "example: run.sh Ar 12 3 10 20"
  echo "example: run.sh Si  4 3 10 20"
  echo "FCC: 1st <110> 3x4 = 12"
  echo "   : 2nd <100> 3x2 =  6"
  echo "DC : 1st <111> 8/2 =  4"
  echo "   : 2nd <110> 3x4 = 12"
  exit
fi

## parameters
TYPE="Ar"
NN=12
W=2  #3
NSTEPS=10
NREF=10
if [[ $# -eq 1 ]]; then
  NSTEPS=$1
else
  TYPE=$1
  NN=$2
  W=$3
  NSTEPS=$4
  NREF=$5
fi
NATOMS=1
if [ $TYPE == "Ar" ]; then
  NATOMS=4
fi
if [ $TYPE == "Si" ]; then
  NATOMS=8
fi
NATOMS=$((NATOMS*W*W*W))
if [[ $# -eq 1 ]]; then
  NREF=$(echo "sqrt($NSTEPS*$NATOMS)" | bc)
fi
echo "* type: $TYPE"
echo "* number of neighbors: $NN"
echo "* cell size: $W [uc], $NATOMS atoms"
echo "* number of steps: $NSTEPS"
echo "* number of references: $NREF"

## files
STEM=${TYPE}-${NN}_${W}x${NSTEPS}
DB=${STEM}.db
SQL=${STEM}.sql
DATA="${TYPE}_$W.data"
TOL="2.0e-3"
MAX_LOCAL_CLUSTERS=-12

## reference run
if [ ! -e $SQL ]; then 
  echo -n "> creating database ${DB} ..."
  DMP=${TYPE}-force.dmp
  if [ -e $DMP ]; then rm -f $DMP; fi
  if [ -e init.log ]; then rm -f init.log ; fi
  if [ -e init.dat ]; then rm -f init.dat ; fi
  sed "s/NN/$NN/g;s/W/$W/g;s/NSTEPS/$NSTEPS/g;s/DB/${DB}/g" in.${TYPE}_init.tmpl > in.init
  $lammps < in.init > init.stdout
  if [ ! -s $DB ]; then 
    echo "!!! database ${DB} is zero size"
    exit
  fi
  if [ -e $DMP ]; then
    ../Database/db_build.py $DB ${NREF} $DMP >& db_build.stdout
  else 
    ../Database/db_build.py $DB ${NREF}      >& db_build.stdout
  fi
  echo "done"
else
  echo "> using existing ${DB}"
fi
if [ ! -s $SQL ]; then 
  echo "!!! SQL database ${SQL} is zero size"
  cat db_builder.stdout
  echo "> removing ${SQL}"; rm $SQL;
  echo "> removing ${DB}";  rm $DB;
  echo "> removing distances";  rm -f ${STEM}.rmsd*
  exit
fi 
if [ ! -e init.dat ]; then 
  echo "!!! no output from init run"
  exit
fi

## database run
if [ -e $DATA ]; then
  if [ -e db.log ]; then rm -f db.log ; fi
  if [ -e db.dat ]; then rm -f db.dat ; fi
  echo -n "> running dynamics ..."
  sed "s/NN/$NN/g;s/W/$W/g;s/NSTEPS/$NSTEPS/g;s/DB/${SQL}/g;s/TOL/${TOL}/g;s/MAX_LOCAL_CLUSTERS/${MAX_LOCAL_CLUSTERS}/g" in.${TYPE}_db.tmpl > in.run
  $lammps < in.run > run.stdout
  echo "done"
  if [ ! -e db.dat ]; then 
    echo "!!! no output from db run"
    exit
  fi
  echo "# system temperatures"
  paste init.dat db.dat | tail -n +2 > comp.dat
  cat comp.dat
else 
  echo "!!! no initial data file $DATA"
fi
