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

#include "pair_ext.h"
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

// needed for "sort" in later GCC's
#include <algorithm>

#include <iostream>
#include <fstream>
#include <iomanip>

#include "database_utility.h"

using std::cout;
using std::setprecision;

//#define PRINT_DEBUG

using namespace LAMMPS_NS;
using namespace MathConst;
using namespace SQL_interface;

// TO DO:
// index,dist sort
// propagate sort to angles (upper triangle matrix)

/* ---------------------------------------------------------------------- */

PairExt::PairExt(LAMMPS *lmp) : Pair(lmp)
{
  one_coeff = 1;
  nclusters = 0;
}

/* ---------------------------------------------------------------------- */

PairExt::~PairExt()
{
  cluster_file.close();
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(cut);
  }
}

/* ---------------------------------------------------------------------- */

void PairExt::compute(int eflag, int vflag)
{
  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int *tag = atom->tag;
  //int nlocal = atom->nlocal;

  // loop over neighbors of my atoms
  int inum = list->inum;
  int * ilist = list->ilist;
  int * numneigh = list->numneigh;
  int ** firstneigh = list->firstneigh;
  int * jlist = NULL;
  int jnum = 0, j = 0, i = 0;
  double xi,yi,zi,dx,dy,dz,rsq;
  int itype, jtype;
  INDICES indices;
  DATA data;
  double fsum[3] = {0,0,0};
  for (int ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itype = type[i];
    xi = x[i][0];
    yi = x[i][1];
    zi = x[i][2];

    jlist = firstneigh[i];
    jnum = numneigh[i];
    int k = 0;
#ifdef PRINT_DEBUG
cout << "\natom " << i << "\n";
#endif
    for (int jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      jtype = type[j];
      j &= NEIGHMASK;
      double * xj = x[j];
      dx = xj[0]-xi;
      dy = xj[1]-yi;
      dz = xj[2]-zi;
      rsq = dx*dx+dy*dy+dz*dz;
      if (rsq < cutsq[itype][jtype]) {
        INDEX index; index.idx = k; index.type=jtype; index.r = rsq;
#ifdef PRINT_DEBUG
cout << k << " " << jj <<"/"<< jnum << " type " << jtype << " rsq " << rsq << " xi " << xi << " xj " << xj[0] << "\n";
#endif
        indices.push_back(index);
        triple p; p.x1=dx; p.x2=dy; p.x3=dz;
        data.push_back(p);
        k++;
      }
    }
    sort(indices.begin(),indices.end(),comp_index);
    string cluster_type = type_string(indices,n_neighbs_);
    write_cluster(indices,data,itype,cluster_type,f[i]);
    for (int j = 0; j < 3; ++j) fsum[j] += f[i][j];
#ifdef PRINT_DEBUG
    cout << (i+1) << " tag:" << atom->tag[i] << "\n";
    //cout << (i+1) << " " << x[i][0] << "  " << f[i][0] << "\n";
#endif
    indices.clear();
    data.clear();
  }
  printf("> sum f = %g %g %g nclusters = %d/%d\n",fsum[0],fsum[1],fsum[2],inum,nclusters);
}

/* ----------------------------------------------------------------------
   write cluster
// [cid], ctype, type_tag, fx, fy, fz, dx_1, dy_1, dz_1, dx_2, ... 
------------------------------------------------------------------------- */

void PairExt::write_cluster(INDICES indices, DATA data, 
 int central_type, string cluster_type, double * f) 
{
  INDICES::iterator it; 
  cluster_file << central_type  << " " << cluster_type << " ";
  cluster_file << f[0]<< " " << f[1] << " " << f[2] << " ";
  int i = 0;
  for (it=indices.begin(); it < indices.end(); it++,i++) {
    if (i == n_neighbs_) break;
    int index = it->idx;
    triple d = data[index];
    cluster_file << d.x1 << " " << d.x2 << " " << d.x3 << " ";
#ifdef PRINT_DEBUG
cout << i << " index:" << index << " " << it->r << "\n";
#endif
  }
  cluster_file << "\n";
  nclusters += 1;
}

/* ----------------------------------------------------------------------
   global settings 
------------------------------------------------------------------------- */

void PairExt::settings(int narg, char **arg)
{
  if (narg != 3) error->all(FLERR,"Illegal pair_style command");
  if (!allocated) allocate();
  double c1 = force->numeric(FLERR,arg[0]);
  database_ = arg[1];
  n_neighbs_ = force->numeric(FLERR,arg[2]);
  cluster_file.open(database_.c_str(),fstream::out);
  //cluster_file.open(database_.c_str(),fstream::out | fstream::app );
  cluster_file << setprecision (fp_precision_) << std::scientific;
  cout << setprecision (fp_precision_) << std::scientific;
#ifdef PRINT_DEBUG
  cout << ">>>> cut_global " << c1 << "\n";
#endif
  double c2 = c1*c1;
  // reset cutoffs that have been explicitly set
  if (allocated) {
    int n = atom->ntypes;
    for (int i = 1; i <= n; i++) {
      for (int j = i; j <= n; j++) {
#ifdef PRINT_DEBUG
        cout << i << " " << j << " f" << setflag[i][j] << "\n";
#endif
        if (setflag[i][j]) {
          cut[i][j] = c1;
          cutsq[i][j] = c2;
#ifdef PRINT_DEBUG
          cout << i << " " << j << " " << cutsq[i][j] << "\n";
#endif
        }
      }
    }
  }
}


/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairExt::coeff(int narg, char **arg)
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

void PairExt::allocate()
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

void PairExt::init_style()
{
  // need a full neighbor list from LAMMPS
  int irequest = neighbor->request(this);
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
  nclusters = 0;
}



/* ----------------------------------------------------------------------
   neighbor callback to inform pair style of neighbor list to use
------------------------------------------------------------------------- */

void PairExt::init_list(int id, NeighList *ptr)
{
  list = ptr;
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairExt::init_one(int i, int j)
{
  return cut[i][j];
}

