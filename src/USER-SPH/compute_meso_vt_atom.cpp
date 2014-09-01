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
#include "compute_meso_vt_atom.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeMesoVTAtom::ComputeMesoVTAtom(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 3) error->all(FLERR,"Number of arguments for compute meso_vt/atom command != 3");
  if (atom->e_flag != 1) error->all(FLERR,"compute meso_vt/atom command requires atom_style with energy (e.g. meso)");

  peratom_flag = 1;
  size_peratom_cols = 3;

  nmax = 0;
  vtvector = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeMesoVTAtom::~ComputeMesoVTAtom()
{
  memory->destroy(vtvector);
}

/* ---------------------------------------------------------------------- */

void ComputeMesoVTAtom::init()
{

  int count = 0;
  for (int i = 0; i < modify->ncompute; i++)
    if (strcmp(modify->compute[i]->style,"vtvector/atom") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute vtvector/atom");
}

/* ---------------------------------------------------------------------- */

void ComputeMesoVTAtom::compute_peratom()
{
  double dtfm;
  
  invoked_peratom = update->ntimestep;

  // grow vtvector array if necessary

  if (atom->nlocal > nmax) {
    memory->destroy(vtvector);
    nmax = atom->nmax;
    memory->create(vtvector,nmax,3,"vtvector/atom:vtvector");
    array_atom = vtvector;
  }

  double **fb = atom->fb;
  double **v = atom->v;
  int *type = atom->type;
  int *mask = atom->mask;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int rmass_flag = atom->rmass_flag;
  int nlocal = atom->nlocal;
  double dtf = 0.5 * update->dt * force->ftm2v;

    for (int i = 0; i < nlocal; i++) {
      if (rmass_flag) {
        dtfm = dtf / rmass[i];
      } else {
        dtfm = dtf / mass[type[i]];
      }

      if (mask[i] & groupbit) {
              vtvector[i][0] = v[i][0] + dtfm * fb[i][0];
              vtvector[i][1] = v[i][1] + dtfm * fb[i][1];
              vtvector[i][2] = v[i][2] + dtfm * fb[i][2];
      }
      else {
              vtvector[i][0] = 0.0;
              vtvector[i][1] = 0.0;
              vtvector[i][2] = 0.0;
      }
    }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeMesoVTAtom::memory_usage()
{
  double bytes = nmax * sizeof(double);
  return bytes;
}
