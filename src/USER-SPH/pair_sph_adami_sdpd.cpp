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
#include "pair_sph_adami_sdpd.h"
#include "atom.h"
#include "force.h"
#include "comm.h"
#include "neigh_list.h"
#include "memory.h"
#include "update.h"
#include "random_mars.h"
#include "error.h"
#include "domain.h"
#include "sph_kernel_dispatch.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairSPHAdamiSDPD::PairSPHAdamiSDPD(LAMMPS *lmp) : Pair(lmp)
{
  random = NULL;
  restartinfo = 0;
  first = 1;
}

/* ---------------------------------------------------------------------- */

PairSPHAdamiSDPD::~PairSPHAdamiSDPD() {
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut);
    memory->destroy(rho0);
    memory->destroy(soundspeed);
    memory->destroy(B);
    memory->destroy(viscosity);
    memory->destroy(pb);
    memory->destroy(temperature);

    int n = atom->ntypes;
    for (int i=0; i<=n; ++i) {
      for (int j=0; j<=n; ++j) 	delete ker[i][j];
      delete[] ker[i];
    }
    delete[] ker;
  }
}

/* ---------------------------------------------------------------------- */

void PairSPHAdamiSDPD::compute(int eflag, int vflag) {
  int i, j, ii, jj, inum, jnum, itype, jtype;
  double xtmp, ytmp, ztmp, delx, dely, delz;

  int *ilist, *jlist, *numneigh, **firstneigh;
  double vxtmp, vytmp, vztmp, imass, jmass, fvisc, velx, vely, velz;
  double rsq, wfd, delVdotDelR, deltaE;

  if (eflag || vflag)
    ev_setup(eflag, vflag);
  else
    evflag = vflag_fdotr = 0;
  
  double **v = atom->v;
  double **x = atom->x;
  double **f = atom->f;
  double **fb = atom->fb;
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
    double pi  = B[itype] * (rho[i] / rho0[itype] - 1.0);
    double Vi  = imass/rho[i];
    double Vi2 = Vi * Vi;

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
	double rabs = sqrt(rsq);
	wfd = ker[itype][jtype]->dw_per_r(rabs, cut[itype][jtype]);
	double Vj  = jmass/rho[j];
	double Vj2 = Vj * Vj;

        // compute pressure
	double pj  = B[jtype] * (rho[j] / rho0[jtype] - 1.0);
	double pij_wave = (rho[j]*pi + rho[i]*pj)/(rho[i] + rho[j]);
	double pij_b    = pb[jtype];

        velx=vxtmp - v[j][0];
        vely=vytmp - v[j][1];
        velz=vztmp - v[j][2];

        // dot product of velocity delta and distance vector
        delVdotDelR = delx * velx + dely * vely + delz * velz;

        fvisc = (Vi2 + Vj2) * viscosity[itype][jtype] * wfd;

        // total pair force & thermal energy increment
        double fpair =   - (Vi2 + Vj2) * pij_wave * wfd;
	double fpair_b = - (Vi2 + Vj2) * pij_b    * wfd;

	double frx, fry, frz, krnd;
	double dt = update->dt;
	krnd = 2*sqrt(Vj2+Vi2)
	  *sqrt(force->boltz*viscosity[itype][jtype]*temperature[itype][jtype]*abs(wfd))
	  /(sqrt(dt)*rabs);
	if (domain->dimension == 3) {
	  double W11, W12, W13;
	  double W21, W22, W23;
	  double W31, W32, W33;
	  W11 = random->gaussian();  W12 = random->gaussian(); W13 = random->gaussian();
	  W21 = random->gaussian();  W22 = random->gaussian(); W23 = random->gaussian();
	  W31 = random->gaussian();  W32 = random->gaussian(); W33 = random->gaussian();
	  //	  W11 = 1*1;  W12 = 1*2; W13 = 1*3;
	  //	  W21 = 2*1;  W22 = 2*2; W23 = 2*3;
	  //	  W31 = 3*1;  W32 = 3*2; W33 = 3*3;

	  frx = krnd*(delx*(W11-(W33+W22+W11)/3.0)+delz*(W31+W13)/2.0+dely*(W21+W12)/2.0);
	  fry = krnd*(dely*(W22-(W33+W22+W11)/3.0)+delz*(W32+W23)/2.0+delx*(W21+W12)/2.0);
	  frz = krnd*(delz*(W33-(W33+W22+W11)/3.0)+dely*(W32+W23)/2.0+delx*(W31+W13)/2.0);
	} else {
	  double W11, W12, W21, W22;
	  W11 = random->gaussian();  W12 = random->gaussian();
	  W21 = random->gaussian();  W22 = random->gaussian();
	  frx = krnd*(delx*(W11-(W22+W11)/2.0)+dely*(W21+W12)/2.0);
	  fry = krnd*(dely*(W22-(W22+W11)/2.0)+delx*(W21+W12)/2.0);
	  frz = 0.0;
	}

        deltaE = -0.5 *(fpair * delVdotDelR + fvisc * (velx*velx + vely*vely + velz*velz));

	// printf("testvar= %f, %f \n", delx, dely);
        f[i][0] += delx * fpair + velx * fvisc + frx;
        f[i][1] += dely * fpair + vely * fvisc + fry;
        f[i][2] += delz * fpair + velz * fvisc + frz;

	// change in background pressure
        fb[i][0] += delx * fpair_b;
        fb[i][1] += dely * fpair_b;
        fb[i][2] += delz * fpair_b;

        // and change in density
        drho[i] += jmass * delVdotDelR * wfd;

        // change in thermal energy
        de[i] += deltaE;

        if (newton_pair || j < nlocal) {
          f[j][0] -= delx * fpair + velx * fvisc + frx;
          f[j][1] -= dely * fpair + vely * fvisc + fry;
          f[j][2] -= delz * fpair + velz * fvisc + frz;

	  fb[j][0] -= delx * fpair_b;
	  fb[j][1] -= dely * fpair_b;
	  fb[j][2] -= delz * fpair_b;

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

void PairSPHAdamiSDPD::allocate() {
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag, n + 1, n + 1, "pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq, n + 1, n + 1, "pair:cutsq");
  memory->create(temperature, n + 1, n + 1, "pair:temperature");

  memory->create(rho0, n + 1, "pair:rho0");
  memory->create(soundspeed, n + 1, "pair:soundspeed");
  memory->create(B, n + 1, "pair:B");
  memory->create(cut, n + 1, n + 1, "pair:cut");
  memory->create(viscosity, n + 1, n + 1, "pair:viscosity");
  memory->create(pb, n + 1, "pair:pb");

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

void PairSPHAdamiSDPD::settings(int narg, char **arg) {
  if (narg != 1)
    error->all(FLERR,
        "Illegal number of setting arguments for pair_style sph/adami/sdpd");

  seed = force->inumeric(FLERR,arg[0]);
  
  if (seed <= 0) error->all(FLERR,"Illegal pair_style command (seed <= 0)");
  delete random;

  random = new RanMars(lmp,seed + comm->me);
}

/* ----------------------------------------------------------------------
 set coeffs for one or more type pairs
 ------------------------------------------------------------------------- */

void PairSPHAdamiSDPD::coeff(int narg, char **arg) {
  if (narg != 9)
    error->all(FLERR,
        "Incorrect args for pair_style sph/adami coefficients");
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

  double rho0_one = force->numeric(FLERR,arg[3]);
  double soundspeed_one = force->numeric(FLERR,arg[4]);
  double viscosity_one = force->numeric(FLERR,arg[5]);
  double pb_one = force->numeric(FLERR,arg[6]);
  double temperature_one = force->numeric(FLERR,arg[7]);
  double cut_one = force->numeric(FLERR,arg[8]);

  double B_one = soundspeed_one * soundspeed_one * rho0_one ;

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    rho0[i] = rho0_one;
    soundspeed[i] = soundspeed_one;
    B[i] = B_one;
    pb[i] = pb_one;

    for (int j = MAX(jlo,i); j <= jhi; j++) {
      viscosity[i][j] = viscosity_one;
      temperature[i][j] = temperature_one;
      cut[i][j] = cut_one;

      ker[i][j] = sph_kernel_dispatch(kernel_name_one, domain->dimension, error);
      if (i!=j) ker[j][i] = sph_kernel_dispatch(kernel_name_one, domain->dimension, error);

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

double PairSPHAdamiSDPD::init_one(int i, int j) {

  if (setflag[i][j] == 0) {
    error->all(FLERR,"Not all pair sph/adami/sdpd coeffs are not set");
  }

  cut[j][i] = cut[i][j];
  temperature[j][i] = temperature[i][j];
  viscosity[j][i] = viscosity[i][j];

  return cut[i][j];
}

/* ---------------------------------------------------------------------- */

double PairSPHAdamiSDPD::single(int i, int j, int itype, int jtype,
    double rsq, double factor_coul, double factor_lj, double &fforce) {
  fforce = 0.0;

  return 0.0;
}