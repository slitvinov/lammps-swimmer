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
#include "pair_sph_heatconduction.h"
#include "sph_energy_equation.h"
#include "atom.h"
#include "force.h"
#include "comm.h"
#include "memory.h"
#include "error.h"
#include "neigh_list.h"
#include "domain.h"
#include "sph_kernel_dispatch.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairSPHHeatConduction::PairSPHHeatConduction(LAMMPS *lmp) : Pair(lmp)
{
  restartinfo = 0;
}

/* ---------------------------------------------------------------------- */

PairSPHHeatConduction::~PairSPHHeatConduction() {
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(cut);
    memory->destroy(alpha);

    int n = atom->ntypes;
    for (int i=0; i<=n; ++i) {
      for (int j=0; j<=n; ++j) 	delete ker[i][j];
      delete[] ker[i];
    }
    delete[] ker;
  }
}

/* ---------------------------------------------------------------------- */

void PairSPHHeatConduction::compute(int eflag, int vflag) {
  int i, j, ii, jj, inum, jnum, itype, jtype;
  double xtmp, ytmp, ztmp, delx, dely, delz;

  int *ilist, *jlist, *numneigh, **firstneigh;
  double imass, jmass, h;
  double rsq, wfd, D, deltaE;

  if (eflag || vflag)
    ev_setup(eflag, vflag);
  else
    evflag = vflag_fdotr = 0;

  double **x = atom->x;
  double *e = atom->e;
  double *de = atom->de;
  double *rho = atom->rho;
  double *rmass, *mass;
  if (atom->rmass_flag) {
    rmass = atom->rmass;
    mass  = NULL;
  } else {
    rmass = NULL;
    mass = atom->mass;
  }
  double *cv = atom->cv;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms and do heat diffusion

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itype = type[i];

    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

    jlist = firstneigh[i];
    jnum = numneigh[i];

    if (atom->rmass) {
      imass = rmass[i];
    } else {
      imass = mass[itype];
    }


    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
        // kernel function
        wfd = ker[itype][jtype]->dw_per_r(sqrt(rsq), cut[itype][jtype]);

	if (atom->rmass) {
	  jmass = rmass[j];
	} else {
	  jmass = mass[jtype];
	}

        D = alpha[itype][jtype]; // diffusion coefficient

	double Ti;
	int region_i;
	int region_j;
	if ( (x[i][0]>0.25) && (x[i][0]<0.75) ) {	
	  Ti = sph_energy2t(e[i], 2.0);
	  region_i = 1;
	} else {
	  Ti = sph_energy2t(e[i], 1.0);
	  region_i = 2;
	}

	double Tj;
	if ( (x[j][0]>0.25) && (x[j][0]<0.75) ) {
	  Tj = sph_energy2t(e[j], 2.0);
	  region_j = 1;
	} else {
	  Tj = sph_energy2t(e[j], 1.0);
	  region_j = 2;	  
	}

	//	if (region_j!=region_i) return;
	  
        double deltaE = 2.0*D*(Ti - Tj)*wfd/(rho[i]*rho[j]);
        de[i] += deltaE; //*jmass;
        if (newton_pair || j < nlocal) {
          de[j] -= deltaE; // *imass;
        }

      }
    }
  }
}

/* ----------------------------------------------------------------------
 allocate all arrays
 ------------------------------------------------------------------------- */

void PairSPHHeatConduction::allocate() {
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag, n + 1, n + 1, "pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq, n + 1, n + 1, "pair:cutsq");
  memory->create(cut, n + 1, n + 1, "pair:cut");
  memory->create(alpha, n + 1, n + 1, "pair:alpha");

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

void PairSPHHeatConduction::settings(int narg, char **arg) {
  if (narg != 0)
    error->all(FLERR,
        "Illegal number of setting arguments for pair_style sph/heatconduction");
}

/* ----------------------------------------------------------------------
 set coeffs for one or more type pairs
 ------------------------------------------------------------------------- */

void PairSPHHeatConduction::coeff(int narg, char **arg) {
  if (narg != 5)
    error->all(FLERR,"Incorrect number of args for pair_style sph/heatconduction coefficients");
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

  double alpha_one = force->numeric(FLERR,arg[3]);
  double cut_one   = force->numeric(FLERR,arg[4]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      //printf("setting cut[%d][%d] = %f\n", i, j, cut_one);
      cut[i][j] = cut_one;

      ker[i][j] = sph_kernel_dispatch(kernel_name_one, domain->dimension, error);
      if (i!=j) ker[j][i] = sph_kernel_dispatch(kernel_name_one, domain->dimension, error);

      alpha[i][j] = alpha_one;
      setflag[i][j] = 1;
      count++;
    }
  }
  delete[] kernel_name_one;

  if (count == 0)
    error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
 init for one type pair i,j and corresponding j,i
 ------------------------------------------------------------------------- */

double PairSPHHeatConduction::init_one(int i, int j) {

  if (setflag[i][j] == 0) {
    error->all(FLERR,"All pair sph/heatconduction coeffs are not set");
  }

  cut[j][i] = cut[i][j];
  alpha[j][i] = alpha[i][j];

  return cut[i][j];
}

/* ---------------------------------------------------------------------- */

double PairSPHHeatConduction::single(int i, int j, int itype, int jtype,
    double rsq, double factor_coul, double factor_lj, double &fforce) {
  fforce = 0.0;

  return 0.0;
}
