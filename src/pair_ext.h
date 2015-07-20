/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(ext,PairExt)

#else

#ifndef LMP_PAIR_EXT_H
#define LMP_PAIR_EXT_H

#include "pair.h"
#include <fstream>
#include <utility>
#include <vector>
#include <string>

using std::fstream;
using std::vector;
using std::pair;
using std::string;

#include "database_utility.h"
using namespace SQL_interface;


namespace LAMMPS_NS {

class PairExt : public Pair {
 public:
  PairExt(class LAMMPS *);
  virtual ~PairExt();
  virtual void compute(int, int);
  void init_list(int, class NeighList *);
  void init_style();
  double init_one(int, int);
  void settings(int, char**);
  void coeff(int, char**);

  // database
  string database_;
  int n_neighbs_;

 protected:
  void allocate();
  void write_cluster(INDICES,DATA,int, string, double *);
  int ** cut;
  fstream cluster_file;
  int nclusters;
};

}

#endif
#endif

