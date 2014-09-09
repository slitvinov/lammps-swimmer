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

#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "pair_sph_bn_rhosum.h"
#include "atom.h"
#include "force.h"
#include "comm.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "memory.h"
#include "error.h"
#include "neighbor.h"
#include "update.h"
#include "domain.h"
#include "sph_kernel_dispatch.h"
#include "sph_kernel.h"
#include "sph_bn_utils.h"
#include <algorithm>    // std::max

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairSPHBNRhoSum::PairSPHBNRhoSum(LAMMPS *lmp) : Pair(lmp)
{
  restartinfo = 0;

  // set comm size needed by this Pair

  comm_forward = 1;
  first = 1;
}

/* ---------------------------------------------------------------------- */

PairSPHBNRhoSum::~PairSPHBNRhoSum() {
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut);

    int n = atom->ntypes;
    for (int i=0; i<=n; ++i) {
      for (int j=0; j<=n; ++j) 	delete ker[i][j];
      delete[] ker[i];
    }
    delete[] ker;

    memory->destroy(T_initial_set);
    memory->destroy(T_target);
  }
}

/* ----------------------------------------------------------------------
 init specific to this pair style
 ------------------------------------------------------------------------- */

void PairSPHBNRhoSum::init_style() {
  // need a full neighbor list
  int irequest = neighbor->request(this);
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
}

/* ---------------------------------------------------------------------- */

void PairSPHBNRhoSum::compute(int eflag, int vflag) {
  int i, j, ii, jj, jnum, itype, jtype;
  double xtmp, ytmp, ztmp, delx, dely, delz;
  double rsq, imass;
  int *jlist;
  double wf;
  // neighbor list variables
  int inum, *ilist, *numneigh, **firstneigh;

  if (eflag || vflag)
    ev_setup(eflag, vflag);
  else
    evflag = vflag_fdotr = 0;

  double **x = atom->x;
  double *rho = atom->rho;
  int *type = atom->type;
  double *mass = atom->mass;

  // check consistency of pair coefficients

  if (first) {
    for (i = 1; i <= atom->ntypes; i++) {
      for (j = 1; i <= atom->ntypes; i++) {
        if (cutsq[i][j] > 0.0) {
          if (!setflag[i][i] || !setflag[j][j]) {
            if (comm->me == 0) {
              printf(
                  "SPH particle types %d and %d interact, but not all of their single particle properties are set.\n",
                  i, j);
            }
          }
        }
      }
    }
    first = 0;
  }

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // recompute density
  // we use a full neighborlist here

  if (nstep != 0) {
    if ((update->ntimestep % nstep) == 0) {

      // initialize density with self-contribution,
      for (ii = 0; ii < inum; ii++) {
        i = ilist[ii];
        itype = type[i];
        imass = mass[itype];
	double rho0i = get_target_field(x[i], domain , T_target,
					nxnodes, nynodes, nznodes);
	double cuti = std::max(get_target_cutoff(imass, nneighbors, rho0i), 
			       cut[itype][itype]); 
	wf = ker[itype][itype]->w(0.0, cuti);
        rho[i] = imass * wf;
      }

      // add density at each atom via kernel function overlap
      for (ii = 0; ii < inum; ii++) {
        i = ilist[ii];
        xtmp = x[i][0];
        ytmp = x[i][1];
        ztmp = x[i][2];
        itype = type[i];
        jlist = firstneigh[i];
        jnum = numneigh[i];

        for (jj = 0; jj < jnum; jj++) {
          j = jlist[jj];
          j &= NEIGHMASK;

          jtype = type[j];
	  double jmass = mass[itype];
          delx = xtmp - x[j][0];
          dely = ytmp - x[j][1];
          delz = ztmp - x[j][2];
          rsq = delx * delx + dely * dely + delz * delz;

          if (rsq < cutsq[itype][jtype]) {
	    double rho0j = get_target_field(x[j], domain , T_target,
					    nxnodes, nynodes, nznodes);
	    double cutj = std::max(get_target_cutoff(jmass, nneighbors, rho0j),
				   cut[itype][itype]); 
	    wf = ker[itype][jtype]->w(sqrt(rsq), cutj);

            rho[i] += mass[jtype] * wf;
          }

        }
      }
    }
  }

  // communicate densities
  comm->forward_comm_pair(this);
}

/* ----------------------------------------------------------------------
 allocate all arrays
 ------------------------------------------------------------------------- */

void PairSPHBNRhoSum::allocate() {
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag, n + 1, n + 1, "pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq, n + 1, n + 1, "pair:cutsq");
  memory->create(cut, n + 1, n + 1, "pair:cut");

  ker = new pSPHKernel*[n+1];
  for (int i=0; i<=n; ++i) {
    ker[i] = new pSPHKernel[n+1];
    for (int j=0; j<=n; ++j)
      ker[i][j] = NULL;
  }
}

/* ----------------------------------------------------------------------
 global settings
 ------------------------------------------------------------------------- */

void PairSPHBNRhoSum::settings(int narg, char **arg) {
  if (narg != 6)
    error->all(FLERR,
        "Illegal number of setting arguments for pair_style sph/bn/rhosum");
  nstep = force->inumeric(FLERR,arg[0]);
  nneighbors = force->inumeric(FLERR,arg[1]);
  nxnodes = force->inumeric(FLERR,arg[2]);
  nynodes = force->inumeric(FLERR,arg[3]);
  nznodes = force->inumeric(FLERR,arg[4]);

  FILE* fpr = fopen(arg[5],"r");
  if (fpr == NULL) {
    char str[128];
    sprintf(str,"Cannot open file %s",arg[5]);
    error->one(FLERR,str);
  }

  if (nxnodes <= 0 || nynodes <= 0 || nznodes <= 0)
    error->all(FLERR,"pair/sph/bn number of nodes must be > 0");


  // allocate 3d grid variables
  total_nnodes = nxnodes*nynodes*nznodes;

  memory->create(T_initial_set,nxnodes,nynodes,nznodes,"ttm:T_initial_set");
  memory->create(T_target,nxnodes,nynodes,nznodes,"ttm:T_target");

  // set target field from input file
  MPI_Comm_rank(world,&me);
  if (me == 0) read_initial_target_field(fpr, nxnodes, nynodes, nznodes, 
					 T_initial_set, T_target, error);
  MPI_Bcast(&T_target[0][0][0],total_nnodes,MPI_DOUBLE,0,world);
}

/* ----------------------------------------------------------------------
 set coeffs for one or more type pairs
 ------------------------------------------------------------------------- */

void PairSPHBNRhoSum::coeff(int narg, char **arg) {
  if (narg != 4)
    error->all(FLERR,"Incorrect number of args for sph/bn/rhosum coefficients");
  if (!allocated)
    allocate();

  int ilo, ihi, jlo, jhi;
  force->bounds(arg[0], atom->ntypes, ilo, ihi);
  force->bounds(arg[1], atom->ntypes, jlo, jhi);

  int i_kernel_name = 2;
  char *kernel_name_one;
  int n_kernel_name = strlen(arg[i_kernel_name]) + 1;
  kernel_name_one = new char[n_kernel_name];
  strcpy(kernel_name_one, arg[i_kernel_name]);
  double cut_one = force->numeric(FLERR,arg[3]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      //printf("setting cut[%d][%d] = %f\n", i, j, cut_one);
      cut[i][j] = cut_one;

      ker[i][j] = sph_kernel_dispatch(kernel_name_one, domain->dimension, error);
      if (i!=j) ker[j][i] = sph_kernel_dispatch(kernel_name_one, domain->dimension, error);

      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0)
    error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
 init for one type pair i,j and corresponding j,i
 ------------------------------------------------------------------------- */

double PairSPHBNRhoSum::init_one(int i, int j) {
  if (setflag[i][j] == 0) {
    error->all(FLERR,"All pair sph/bn/rhosum coeffs are not set");
  }

  cut[j][i] = cut[i][j];

  return cut[i][j];
}

/* ---------------------------------------------------------------------- */

double PairSPHBNRhoSum::single(int i, int j, int itype, int jtype, double rsq,
    double factor_coul, double factor_lj, double &fforce) {
  fforce = 0.0;

  return 0.0;
}

/* ---------------------------------------------------------------------- */

int PairSPHBNRhoSum::pack_forward_comm(int n, int *list, double *buf, 
                                     int pbc_flag, int *pbc) {
  int i, j, m;
  double *rho = atom->rho;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = rho[j];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void PairSPHBNRhoSum::unpack_forward_comm(int n, int first, double *buf) {
  int i, m, last;
  double *rho = atom->rho;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++)
    rho[i] = buf[m++];
}
