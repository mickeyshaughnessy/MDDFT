# Install/unInstall EXT classes in LAMMPS
# called via : sh -f Install.sh

CP="ln -sf"
#SRC=`echo src/*.cpp src/*.h modsrc/*.h modsrc/*.cpp`
#SRC="src/pair_database.h src/pair_database.cpp modsrc/pair_lj_cut.cpp"
SRC=`echo src/*.cpp src/*.h`
HOME=EXT
#echo $SRC

if (test $1 = 1) then
  cd ..
  # link source files
  for src in $SRC ; do
    $CP $HOME/$src .
    #echo $src
  done
  cd MAKE; $CP ../EXT/Makefile.* .;cd .
elif (test $1 = 0) then
  # remove links
  cd ..
  for src in $SRC ; do
    src=${src##*/}
    echo "removing $src"
    rm -f $src
  done
  rm -f MAKE/Makefile.*-ext
fi

