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

/* ----------------------------------------------------------------------
   Contributing author: Carsten Svaneborg, science@zqex.dk
------------------------------------------------------------------------- */

#include "math.h"
#include "stdlib.h"
#include "bond_harmonic_swimmer_extended.h"
#include "atom.h"
#include "neighbor.h"
#include "domain.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "error.h"
#include "update.h"

using namespace LAMMPS_NS;

#define EPSILON 1.0e-20

/* ---------------------------------------------------------------------- */

BondHarmonicSwimmerExtended::BondHarmonicSwimmerExtended(LAMMPS *lmp) : Bond(lmp) {
  time_origin = update->ntimestep;
}

/* ---------------------------------------------------------------------- */

BondHarmonicSwimmerExtended::~BondHarmonicSwimmerExtended()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(k);
    memory->destroy(r0);
    memory->destroy(r1);

    memory->destroy(A_alpha);
    memory->destroy(A_beta);

    memory->destroy(omega_alpha);
    memory->destroy(omega_beta);

    memory->destroy(phi);
    memory->destroy(vel_sw);

    memory->destroy(n1);
    memory->destroy(n2);
  }
}

/* ---------------------------------------------------------------------- */

void BondHarmonicSwimmerExtended::compute(int eflag, int vflag)
{
  int i1,i2,n,type;
  tagint tag1, tag2;
  double delx,dely,delz,ebond,fbond;
  double rsq,r,dr,rk;
  double r0_local;

  ebond = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = 0;

  double **x = atom->x;
  double **f = atom->f;
  tagint *tag = atom->tag;

  int **bondlist = neighbor->bondlist;
  int nbondlist = neighbor->nbondlist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;
  double delta = (update->ntimestep - time_origin) * update->dt;

  for (n = 0; n < nbondlist; n++) {
    i1 = bondlist[n][0];
    i2 = bondlist[n][1];

    tag1 = tag[i1];
    tag2 = tag[i2];

    if (tag1>=tag2) {
       char str[128];
       sprintf(str,"tag1>=tag2: something wrong with a bond between %i and %i", i1, i2);
       error->all(FLERR, str);
    }
 
    type = bondlist[n][2];

    delx = x[i1][0] - x[i2][0];
    dely = x[i1][1] - x[i2][1];
    delz = x[i1][2] - x[i2][2];

    rsq = delx*delx + dely*dely + delz*delz;
    r = sqrt(rsq);
   
    if ( (tag1>=n1[type]) && (tag1<=n2[type]) && ((tag2-tag1)==1) ) {
      double dn = static_cast<double>(tag1) - n1[type];
      double omega = omega_beta[type]*dn + omega_alpha[type];
      double A     = A_beta[type]*dn     + A_alpha[type];
      double s_aux = sin(omega*dn + phi[type] - vel_sw[type]*delta);
      r0_local = r0[type] + A*s_aux ;
    } else {
       r0_local = r0[type];
    }

    dr = r - r0_local;
    rk = k[type] * dr;

    // force & energy

    if (r > 0.0) fbond = -2.0*rk/r;
    else fbond = 0.0;

    if (eflag)
      ebond = k[type]*(dr*dr -(r0_local-r1[type])*(r0_local-r1[type]) );

    // apply force to each of 2 atoms

    if (newton_bond || i1 < nlocal) {
      f[i1][0] += delx*fbond;
      f[i1][1] += dely*fbond;
      f[i1][2] += delz*fbond;
    }

    if (newton_bond || i2 < nlocal) {
      f[i2][0] -= delx*fbond;
      f[i2][1] -= dely*fbond;
      f[i2][2] -= delz*fbond;
    }

    if (evflag) ev_tally(i1,i2,nlocal,newton_bond,ebond,fbond,delx,dely,delz);
  }
}

/* ---------------------------------------------------------------------- */

void BondHarmonicSwimmerExtended::allocate()
{
  allocated = 1;
  int n = atom->nbondtypes;

  memory->create(k ,    n+1,"bond:k");
  memory->create(r0,    n+1,"bond:r0");
  memory->create(r1,    n+1,"bond:r1");

  memory->create(A_alpha,    n+1,"bond:A_alpha");
  memory->create(A_beta,    n+1,"bond:A_beta");

  memory->create(omega_alpha,    n+1,"bond:omega_alpha");
  memory->create(omega_beta,     n+1,"bond:omega_beta");

  memory->create(phi,    n+1,"bond:phi");
  memory->create(vel_sw,    n+1,"bond:vel_sw");

  memory->create(n1,    n+1,"bond:n1");
  memory->create(n2,    n+1,"bond:n2");

  memory->create(setflag,n+1,"bond:setflag");

  for (int i = 1; i <= n; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more types
------------------------------------------------------------------------- */

void BondHarmonicSwimmerExtended::coeff(int narg, char **arg)
{
  if (narg != 12) error->all(FLERR,"Incorrect args for bond coefficients");
  if (!allocated) allocate();
  
  if (atom->tag_enable==0) {
    error->all(FLERR,"Bond harmonic/swimmer/extended requires tag_enable=1");
  }

  int ilo,ihi;
  force->bounds(arg[0],atom->nbondtypes,ilo,ihi);

  double Umin = force->numeric(FLERR,arg[1]);   // energy at minimum
  double r0_one = force->numeric(FLERR,arg[2]); // position of minimum
  double r1_one = force->numeric(FLERR,arg[3]);  // position where energy = 0

  // swimmer wave parameters A*sin(omega*N + phi - vel_sw*time)
  double A_alpha_one = force->numeric(FLERR,arg[4]);
  double A_beta_one =  force->numeric(FLERR,arg[5]);

  double omega_alpha_one = force->numeric(FLERR,arg[6]);
  double omega_beta_one =  force->numeric(FLERR,arg[7]);

  double phi_one = force->numeric(FLERR,arg[8]);
  double vel_sw_one = force->numeric(FLERR,arg[9]);

  tagint n1_one = force->numeric(FLERR,arg[10]);
  tagint n2_one = force->numeric(FLERR,arg[11]);

  if (r0_one == r1_one)
    error->all(FLERR,"Bond harmonic/swimmer/extended r0 and r1 must be different");

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    k[i] = Umin/((r0_one-r1_one)*(r0_one-r1_one));
    r0[i] = r0_one;
    r1[i] = r1_one;

    A_alpha[i] = A_alpha_one;
    A_beta[i] =  A_beta_one;

    omega_alpha[i] = omega_alpha_one;
    omega_beta[i] = omega_beta_one;

    phi[i] = phi_one;
    vel_sw[i] = vel_sw_one;

    n1[i] = n1_one;
    n2[i] = n2_one;

    setflag[i] = 1;
    count++;
  }

  if (count == 0) error->all(FLERR,"Incorrect args for bond coefficients");
}

/* ----------------------------------------------------------------------
   return an equilbrium bond length
------------------------------------------------------------------------- */

double BondHarmonicSwimmerExtended::equilibrium_distance(int i)
{
  return r0[i];
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file
------------------------------------------------------------------------- */

void BondHarmonicSwimmerExtended::write_restart(FILE *fp)
{
  fwrite(&k[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&r0[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&r1[1],sizeof(double),atom->nbondtypes,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

void BondHarmonicSwimmerExtended::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
    fread(&k[1],sizeof(double),atom->nbondtypes,fp);
    fread(&r0[1],sizeof(double),atom->nbondtypes,fp);
    fread(&r1[1],sizeof(double),atom->nbondtypes,fp);
  }
  MPI_Bcast(&k[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&r0[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&r1[1],atom->nbondtypes,MPI_DOUBLE,0,world);

  for (int i = 1; i <= atom->nbondtypes; i++) setflag[i] = 1;
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void BondHarmonicSwimmerExtended::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->nbondtypes; i++) {
    double d2 = (r0[i]-r1[i])*(r0[i]-r1[i]);
    fprintf(fp,"%d %g %g %g\n",i,k[i]*d2,r0[i],r1[i]);
  }
}

/* ---------------------------------------------------------------------- */

double BondHarmonicSwimmerExtended::single(int type, double rsq, int i, int j,
				 double &fforce)
{
  double r = sqrt(rsq);
  double dr = r - r0[type];
  double dr2=r0[type]-r1[type];

  fforce =  -2.0*k[type]*dr/r;
  return k[type]*(dr*dr - dr2*dr2);
}