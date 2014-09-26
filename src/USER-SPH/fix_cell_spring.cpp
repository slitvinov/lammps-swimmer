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

#include "string.h"
#include "fix_cell_spring.h"
#include "atom.h"
#include "update.h"
#include "domain.h"
#include "respa.h"
#include "error.h"
#include "force.h"
#include <algorithm> // std::max()

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixCellSpring::FixCellSpring(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg != 4) error->all(FLERR,"Illegal fix cell/spring command");
  k_fene = force->numeric(FLERR,arg[3]);
}

/* ---------------------------------------------------------------------- */

int FixCellSpring::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixCellSpring::init()
{
  if ( (atom->xc_flag != 1) || (atom->rc_flag != 1) )
    error->all(FLERR,"Atom style should have `xc' and 'rc' attributes to used cell/stick");
}

/* ---------------------------------------------------------------------- */

void FixCellSpring::setup(int vflag)
{
  if (strstr(update->integrate_style,"verlet"))
    post_force(vflag);
  else {
    int nlevels_respa = ((Respa *) update->integrate)->nlevels;
    for (int ilevel = 0; ilevel < nlevels_respa; ilevel++) {
      ((Respa *) update->integrate)->copy_flevel_f(ilevel);
      post_force_respa(vflag,ilevel,0);
      ((Respa *) update->integrate)->copy_f_flevel(ilevel);
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixCellSpring::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixCellSpring::post_force(int vflag)
{
  double **x = atom->x;
  double **xc = atom->xc;
  double **f = atom->f;
  double *rc = atom->rc;

  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      double delx = x[i][0] - xc[i][0];
      double dely = x[i][1] - xc[i][1];
      double delz = x[i][2] - xc[i][2];
      
      double rsq = delx*delx + dely*dely + delz*delz;
      double r0sq = rc[i] * rc[i];
      double rlogarg = std::min(1.0 - rsq/r0sq, 0.01);
      double fbond = -k_fene/rlogarg;

      f[i][0] += delx*fbond;
      f[i][1] += dely*fbond;
      f[i][2] += delz*fbond;
    }
}

/* ---------------------------------------------------------------------- */

void FixCellSpring::post_force_respa(int vflag, int ilevel, int iloop)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixCellSpring::min_post_force(int vflag)
{
  post_force(vflag);
}
