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

#include "stdio.h"
#include "string.h"
#include "fix_meso_novel.h"
#include "math.h"
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
#include "memory.h"
#include "error.h"
#include "pair.h"

#include "assert.h"
#include <iostream>
using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixMesoNoVel::FixMesoNoVel(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg) {

  if (atom->rho_flag != 1)
    error->all(FLERR,
        "fix meso/novel command requires atom_style with density");

  if (narg != 4)
    error->all(FLERR,"Illegal number of arguments for fix meso command");

  max_dr = force->numeric(FLERR,arg[3]);

  time_integrate = 1;
}

/* ---------------------------------------------------------------------- */

int FixMesoNoVel::setmask() {
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixMesoNoVel::init() {
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
}

/* ----------------------------------------------------------------------
 allow for both per-type and per-atom mass
 ------------------------------------------------------------------------- */

void FixMesoNoVel::initial_integrate(int vflag) {
  // update v and x and rho and e of atoms in group

  double **x = atom->x;
  double **f = atom->f;
  double *rho = atom->rho;
  double *drho = atom->drho;
  double *e = atom->e;
  double *de = atom->de;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int rmass_flag = atom->rmass_flag;

  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int i;
  double dtfm;

  if (igroup == atom->firstgroup)
    nlocal = atom->nfirst;

  // find df
  double maxA2 = 0.0;
  double max;
  double mass_atom;
  for (i=0; i<nlocal; i++) {
    if (mask[i] & groupbit) {
      if (rmass_flag) { 
	mass_atom = rmass[i];
      } else {
        mass_atom = mass[type[i]];
      }
      double A2 = (f[i][0]*f[i][0] + f[i][1]*f[i][1] + f[i][2]*f[i][2])/(mass_atom*mass_atom);
      if (A2>maxA2) maxA2 = A2;
    }
  }
  MPI_Allreduce(&maxA2,&max,1,MPI_DOUBLE,MPI_MAX,world);
  double dt2 = max_dr / sqrt(max);

  for (i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      if (rmass_flag) { 
	mass_atom = rmass[i];
      } else {
        mass_atom = mass[type[i]];
      }

      x[i][0] += dt2 * f[i][0]/mass_atom;
      x[i][1] += dt2 * f[i][1]/mass_atom;
      x[i][2] += dt2 * f[i][2]/mass_atom;
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixMesoNoVel::reset_dt() {
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
}
