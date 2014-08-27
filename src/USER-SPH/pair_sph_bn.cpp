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
#include "pair_sph_bn.h"
#include "atom.h"
#include "force.h"
#include "comm.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"
#include "domain.h"
#include "sph_kernel_dispatch.h"
#include "sph_bn_utils.h"

using namespace LAMMPS_NS;


/* ---------------------------------------------------------------------- */

PairSPHBN::PairSPHBN(LAMMPS *lmp) : Pair(lmp)
{
  restartinfo = 0;
  first = 1;
}

/* ---------------------------------------------------------------------- */

PairSPHBN::~PairSPHBN() {
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut);
    memory->destroy(soundspeed);
    memory->destroy(B);
    memory->destroy(viscosity);

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

/* ---------------------------------------------------------------------- */

void PairSPHBN::compute(int eflag, int vflag) {
  int i, j, ii, jj, inum, jnum, itype, jtype;
  double xtmp, ytmp, ztmp, delx, dely, delz, fpair;

  int *ilist, *jlist, *numneigh, **firstneigh;
  double vxtmp, vytmp, vztmp, imass, jmass, fi, fj, fvisc, velx, vely, velz;
  double rsq, tmp, wfd, delVdotDelR, deltaE;

  if (eflag || vflag)
    ev_setup(eflag, vflag);
  else
    evflag = vflag_fdotr = 0;
  
  double **v = atom->vest;
  double **x = atom->x;
  double **f = atom->f;
  double *rho = atom->rho;
  double *mass = atom->mass;
  double *de = atom->de;
  double *drho = atom->drho;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;
  // check consistency of pair coefficients

  if (first) {
    for (i = 1; i <= atom->ntypes; i++) {
      for (j = 1; i <= atom->ntypes; i++) {
        if (cutsq[i][j] > 1.e-32) {
          if (!setflag[i][i] || !setflag[j][j]) {
            if (comm->me == 0) {
              printf(
                  "SPH particle types %d and %d interact with cutoff=%g, but not all of their single particle properties are set.\n",
                  i, j, sqrt(cutsq[i][j]));
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

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];

    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    vxtmp = v[i][0];
    vytmp = v[i][1];
    vztmp = v[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    imass = mass[itype];

    // compute pressure of atom i
    double rho0i = get_target_field(x[i], domain , T_target,
      nxnodes, nynodes, nznodes);
    double pi = bn_eos(rho[i], rho0i, B[itype]);
    fi = pi  / (rho[i] * rho[i]);

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;
      jtype = type[j];
      jmass = mass[jtype];

      if (rsq < cutsq[itype][jtype]) {
	wfd = ker[itype][jtype]->dw_per_r(sqrt(rsq), cut[itype][jtype]);

        // compute pressure  of atom j
	double rho0j = get_target_field(x[j], domain , T_target,
	  nxnodes, nynodes, nznodes);
        double pj = bn_eos(rho[j], rho0j, B[jtype]);
        fj = pj  / (rho[j] * rho[j]);

        velx=vxtmp - v[j][0];
        vely=vytmp - v[j][1];
        velz=vztmp - v[j][2];

        // dot product of velocity delta and distance vector
        delVdotDelR = delx * velx + dely * vely + delz * velz;

        // Morris Viscosity (Morris, 1996)

        fvisc = 2 * viscosity[itype][jtype] / (rho[i] * rho[j]);

        fvisc *= imass * jmass * wfd;

        // total pair force & thermal energy increment
        fpair = -imass * jmass * (fi + fj) * wfd;
        deltaE = -0.5 *(fpair * delVdotDelR + fvisc * (velx*velx + vely*vely + velz*velz));

       // printf("testvar= %f, %f \n", delx, dely);

        f[i][0] += delx * fpair + velx * fvisc;
        f[i][1] += dely * fpair + vely * fvisc;
        f[i][2] += delz * fpair + velz * fvisc;

        // and change in density
        drho[i] += jmass * delVdotDelR * wfd;

        // change in thermal energy
        de[i] += deltaE;

        if (newton_pair || j < nlocal) {
          f[j][0] -= delx * fpair + velx * fvisc;
          f[j][1] -= dely * fpair + vely * fvisc;
          f[j][2] -= delz * fpair + velz * fvisc;
          de[j] += deltaE;
          drho[j] += imass * delVdotDelR * wfd;
        }

        if (evflag)
          ev_tally(i, j, nlocal, newton_pair, 0.0, 0.0, fpair, delx, dely, delz);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
 allocate all arrays
 ------------------------------------------------------------------------- */

void PairSPHBN::allocate() {
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag, n + 1, n + 1, "pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq, n + 1, n + 1, "pair:cutsq");

  memory->create(soundspeed, n + 1, "pair:soundspeed");
  memory->create(B, n + 1, "pair:B");
  memory->create(cut, n + 1, n + 1, "pair:cut");
  memory->create(viscosity, n + 1, n + 1, "pair:viscosity");

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
void PairSPHBN::settings(int narg, char **arg) {
  if (narg != 4) error->all(FLERR,"Illegal pair_style pair/sph/bn command"
			    " (4 arguments required)");

  nxnodes = force->inumeric(FLERR,arg[0]);
  nynodes = force->inumeric(FLERR,arg[1]);
  nznodes = force->inumeric(FLERR,arg[2]);

  FILE* fpr = fopen(arg[3],"r");
  if (fpr == NULL) {
    char str[128];
    sprintf(str,"Cannot open file %s",arg[3]);
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

 void PairSPHBN::coeff(int narg, char **arg) {
  if (narg != 6)
    error->all(FLERR,
        "Incorrect args for pair_style pair/sph/bn coefficients");
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

  double soundspeed_one = force->numeric(FLERR,arg[3]);
  double viscosity_one = force->numeric(FLERR,arg[4]);
  double cut_one = force->numeric(FLERR,arg[5]);
  double B_one = soundspeed_one * soundspeed_one;

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    soundspeed[i] = soundspeed_one;
    B[i] = B_one;
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      viscosity[i][j] = viscosity_one;
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

double PairSPHBN::init_one(int i, int j) {

  if (setflag[i][j] == 0) {
    error->all(FLERR,"Not all pair pair/sph/bn coeffs are not set");
  }

  cut[j][i] = cut[i][j];
  viscosity[j][i] = viscosity[i][j];

  return cut[i][j];
}

/* ---------------------------------------------------------------------- */

double PairSPHBN::single(int i, int j, int itype, int jtype,
    double rsq, double factor_coul, double factor_lj, double &fforce) {
  fforce = 0.0;

  return 0.0;
}


double PairSPHBN::bn_eos (double rho, double rho0, double B) {
    return B  * rho / rho0;
}
