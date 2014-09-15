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
#include "compute_target_field.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "error.h"
#include "sph_bn_utils.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeTargetField::ComputeTargetField(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 8) error->all(FLERR,"Wrong number of arguments for compute target/field command");

  nxnodes = force->inumeric(FLERR,arg[3]);
  nynodes = force->inumeric(FLERR,arg[4]);
  nznodes = force->inumeric(FLERR,arg[5]);

  FILE* fpr = fopen(arg[6],"r");
  if (fpr == NULL) {
    char str[128];
    sprintf(str,"Cannot open file %s",arg[6]);
    error->one(FLERR,str);
  }
  ntime_smooth = force->inumeric(FLERR,arg[7]);

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

  peratom_flag = 1;
  size_peratom_cols = 0;

  nmax = 0;
  evector = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeTargetField::~ComputeTargetField()
{
  memory->sfree(evector);
}

/* ---------------------------------------------------------------------- */

void ComputeTargetField::init()
{

  int count = 0;
  for (int i = 0; i < modify->ncompute; i++)
    if (strcmp(modify->compute[i]->style,"evector/atom") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute evector/atom");
}

/* ---------------------------------------------------------------------- */

void ComputeTargetField::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  // grow evector array if necessary

  if (atom->nlocal > nmax) {
    memory->sfree(evector);
    nmax = atom->nmax;
    evector = (double *) memory->smalloc(nmax*sizeof(double),"evector/atom:evector");
    vector_atom = evector;
  }

  double *e = atom->e;
  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      evector[i] = get_target_field(x[i], domain , T_target,
				    nxnodes, nynodes, nznodes,
				    ntime_smooth,     update->ntimestep);
    }
    else {
      evector[i] = 0.0;
    }
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeTargetField::memory_usage()
{
  double bytes = nmax * sizeof(double);
  return bytes;
}
