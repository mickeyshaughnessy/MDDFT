/* ----------------------------------------------------------------------

   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */
#include "pair_database.h"
#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "update.h"
#include "integrate.h"
#include "respa.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"
#include <map>
#include <time.h>
#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "database_utility.h"
using namespace SQL_interface;
using std::cout;
using namespace LAMMPS_NS;
using namespace MathConst;
using namespace std;
using namespace Eigen;
//=====================================================================
#define PRINT_DEBUG
//#undef PRINT_DEBUG
#define VERBOSE
//#define MAX_CLUSTERS_BREAK

// TO DO:
//* time parallel
//* use sql to winnow on type
// * find a way to conserve momentum -- I think this can be done by enforcing it as a constraint in this way:
  //Sum_of_forces = eps at any time step
  //if eps > energy tolerance
  //Then add -(eps/N)*weight_i to each force_i that is passed to the LAMMPS integrator, where weight i is some function of the error designed to keep the sum of the constraint forces equal to eps
//=====================================================================

/* ---------------------------------------------------------------------- */
PairDatabase::PairDatabase(LAMMPS *lmp) : Pair(lmp), db_(NULL) { }
/* ---------------------------------------------------------------------- */
PairDatabase::~PairDatabase()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(cut);
  }
  if(db_) delete db_;
}
/* ---------------------------------------------------------------------- */
void PairDatabase::compute(int eflag, int vflag)
{
  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;

  // loop over neighbors of my atoms
  int inum = list->inum;
  int * ilist = list->ilist;
  int * numneigh = list->numneigh;
  int ** firstneigh = list->firstneigh;
  int * jlist = NULL;
  
  int jnum = 0, j = 0, i = 0;
  double xi,yi,zi,dx,dy,dz,rsq;
  int itype, jtype;

  nDistComps_ = 0;
  INDICES indices;
  DATA data;
  MatrixXd C_(n_neighbs_,3);
  for (int ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itype = type[i];
    double * Xi = x[i];
    xi = Xi[0];
    yi = Xi[1];
    zi = Xi[2];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    int k = 0;
    for (int jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      jtype = type[j];
      j &= NEIGHMASK;
      double * Xj = x[j];
      dx = xi - Xj[0];
      dy = yi - Xj[1];
      dz = zi - Xj[2];
      rsq = dx*dx+dy*dy+dz*dz;
      if (rsq < cutsq[itype][jtype]) {
        INDEX index; index.idx = k; index.type = jtype; index.r = rsq;
        indices.push_back(index);
        triple p; p.x1=-dx; p.x2=-dy; p.x3=-dz;
#ifdef PRINT_DEBUG
//      printf("%i type %i dx %g %g %g %g < %g \n",jj,jtype,dx,dy,dz,rsq,cutsq[itype][jtype]);
#endif
        data.push_back(p);
        k++;
      }
    }
    if (indices.size() < (unsigned int) n_neighbs_) {
      error->all(FLERR,"cluster is smaller than expected number of neighbors: increase neighborhood cutoff");
    }
    // sort indices by distance
    sort(indices.begin(),indices.end(),comp_index);
#ifdef PRINT_DEBUG
    //int step = update->ntimestep;
    cout << "\n==============================================================\n";
    cout << i+1 << " sorted cluster: rank, sort_id, type, distance\n";
    for (unsigned int jj=0; jj<indices.size(); jj++) {
      INDEX index = indices[jj];
      if ((int) jj==n_neighbs_) printf("-----------------\n");
      printf("%2d %4d %2d %7.4f\n",jj+1,index.idx,index.type,sqrt(index.r));
    }
#endif
    // pack
    string C_type = type_string(indices,n_neighbs_);
    for (int jj=0;jj<n_neighbs_;jj++) {
      int index = indices[jj].idx;
      triple d = data[index];
      C_(jj,0) = d.x1;
      C_(jj,1) = d.x2;
      C_(jj,2) = d.x3;
    }
    db_->force(C_, itype, C_type, f[i]);
    indices.clear();
    data.clear(); 
  } // central atom loop
#ifdef VERBOSE
  db_->report();
#endif
}
/* ----------------------------------------------------------------------
   global settings 
------------------------------------------------------------------------- */
void PairDatabase::settings(int narg, char **arg)
{
  if (narg != 5) error->all(FLERR,"Illegal pair_style command");
  if (!allocated) allocate();
  double c1 = force->numeric(FLERR,arg[0]);
  n_neighbs_ = force->numeric(FLERR,arg[1]);
  database_ = arg[2];
  rmsdTol_ = force->numeric(FLERR,arg[3]);
  maxLocalClusters_ = force->numeric(FLERR,arg[4]);
#ifdef PRINT_DEBUG
  cout << "> database " << database_ << " with cluster size " << n_neighbs_ << ",\n    tol " << rmsdTol_ << ", max local clusters " << maxLocalClusters_ << "> cut_global " << c1 << "\n";
#endif
  double c2 = c1*c1;
  // reset cutoffs that have been explicitly set
  if (allocated) {
    int n = atom->ntypes;
    for (int i = 1; i <= n; i++) {
      for (int j = i; j <= n; j++) {
        if (setflag[i][j]) {
          cut[i][j] = c1;
          cutsq[i][j] = c2;
        }
      }
    }
  }
}
/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */
void PairDatabase::coeff(int narg, char **arg)
{
  if (!allocated) allocate();
  int ilo,ihi,jlo,jhi;
  force->bounds(arg[0],atom->ntypes,ilo,ihi);
  force->bounds(arg[1],atom->ntypes,jlo,jhi);

  if (narg < 3) return;
  double cut_one = force->numeric(FLERR,arg[2]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      cut[i][j] = cut_one;
      setflag[i][j] = 1;
      count++;
    }
  }
  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}
/* ----------------------------------------------------------------------
   allocate all arrays 
------------------------------------------------------------------------- */
void PairDatabase::allocate()
{
  allocated = 1;
  int n = atom->ntypes;
  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 1; // apply to all types
  memory->create(cutsq,n+1,n+1,"pair:cutsq");
  memory->create(cut,n+1,n+1,"pair:cut");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */
void PairDatabase::init_style()
{
  // load database
  bool useExhaustive = (maxLocalClusters_ < 1);
  maxLocalClusters_=abs(maxLocalClusters_);
  db_ = new ClusterDatabase(database_,n_neighbs_,maxLocalClusters_,rmsdTol_);
  if (useExhaustive) db_->use_exhaustive_search();
  int n_clusters  = db_->size();
  if (n_clusters == 0) { 
    stringstream ss; ss << "database " << database_ << " has no entries";
    error->all(FLERR,(ss.str()).c_str());
  }
  int size = db_->configuration_size();
  if (comm->me == 0 && screen) {
    printf("database %s has %i clusters of size %i\n",database_.c_str(),n_clusters,size);
    if (db_->has_distances()) { printf("using metric search\n"); }
    else                      { printf("using exhaustive search\n"); }
  }
  if (size < n_neighbs_) {
    stringstream ss; ss << "database " << database_ << " has cluster size " << size << "  but " << n_neighbs_ << " neighbors have been requested";
    error->all(FLERR,(ss.str()).c_str());
  }
  // need a full neighbor list from LAMMPS
  int irequest = neighbor->request(this);
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
}

/* ----------------------------------------------------------------------
   neighbor callback to inform pair style of neighbor list to use
------------------------------------------------------------------------- */

void PairDatabase::init_list(int id, NeighList *ptr)
{
  list = ptr;
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairDatabase::init_one(int i, int j)
{
  return cut[i][j];
}
