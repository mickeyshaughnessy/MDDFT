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

PairStyle(database,PairDatabase)

#else

#ifndef LMP_PAIR_DATABASE_H
#define LMP_PAIR_DATABASE_H

#include "pair.h"
#include "database_utility.h"
using namespace SQL_interface;

namespace LAMMPS_NS {

class PairDatabase : public Pair {
 public:
  PairDatabase(class LAMMPS *);
  virtual ~PairDatabase();
  virtual void compute(int, int);
  void init_list(int, class NeighList *);
  void init_style();
  double init_one(int, int);
  void settings(int, char**);
  void coeff(int, char**);

 protected:
  void allocate();
  // database 
  string database_;
  class ClusterDatabase * db_;
  // parameters
  double rmsdTol_;
  int maxLocalClusters_;
  int n_neighbs_;
  // workspace
  int nDistComps_;
  int clusterCount_;
  //int *  prev_ids_;
  // other
  int ** cut;
};
}
#endif
#endif
