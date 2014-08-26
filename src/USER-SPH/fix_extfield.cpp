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
   Contributing authors: Paul Crozier (SNL)
                         Carolyn Phillips (University of Michigan)
------------------------------------------------------------------------- */

#include "lmptype.h"
#include "mpi.h"
#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "fix_extfield.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "domain.h"
#include "region.h"
#include "respa.h"
#include "comm.h"
#include "random_mars.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

#define MAXLINE 1024

/* ---------------------------------------------------------------------- */

FixExtField::FixExtField(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg != 7) error->all(FLERR,"Illegal fix_extfield command");

  vector_flag = 1;
  size_vector = 1;
  global_freq = 1;
  extvector = 1;
  nevery = 1;
  restart_peratom = 1;
  restart_global = 1;

  nxnodes = force->inumeric(FLERR,arg[3]);
  nynodes = force->inumeric(FLERR,arg[4]);
  nznodes = force->inumeric(FLERR,arg[5]);

  fpr = fopen(arg[6],"r");
  if (fpr == NULL) {
    char str[128];
    sprintf(str,"Cannot open file %s",arg[6]);
    error->one(FLERR,str);
  }


  if (nxnodes <= 0 || nynodes <= 0 || nznodes <= 0)
    error->all(FLERR,"fix_extfield number of nodes must be > 0");


  // allocate 3d grid variables
  total_nnodes = nxnodes*nynodes*nznodes;

  memory->create(T_initial_set,nxnodes,nynodes,nznodes,"ttm:T_initial_set");
  memory->create(T_electron,nxnodes,nynodes,nznodes,"ttm:T_electron");

  // set initial electron temperatures from user input file

  if (me == 0) read_initial_electron_temperatures();
  MPI_Bcast(&T_electron[0][0][0],total_nnodes,MPI_DOUBLE,0,world);
}

/* ---------------------------------------------------------------------- */

FixExtField::~FixExtField()
{
  memory->destroy(T_initial_set);
  memory->destroy(T_electron);
}

/* ---------------------------------------------------------------------- */

int FixExtField::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixExtField::post_force(int vflag)
{
  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {

      double xscale = (x[i][0] - domain->boxlo[0])/domain->xprd;
      double yscale = (x[i][1] - domain->boxlo[1])/domain->yprd;
      double zscale = (x[i][2] - domain->boxlo[2])/domain->zprd;
      int ixnode = static_cast<int>(xscale*nxnodes);
      int iynode = static_cast<int>(yscale*nynodes);
      int iznode = static_cast<int>(zscale*nznodes);
      while (ixnode > nxnodes-1) ixnode -= nxnodes;
      while (iynode > nynodes-1) iynode -= nynodes;
      while (iznode > nznodes-1) iznode -= nznodes;
      while (ixnode < 0) ixnode += nxnodes;
      while (iynode < 0) iynode += nynodes;
      while (iznode < 0) iznode += nznodes;

      printf("T_electron: %g\n", T_electron[ixnode][iynode][iznode]);
      
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixExtField::post_force_respa(int vflag, int ilevel, int iloop)
{
  if (ilevel == nlevels_respa-1) post_force(vflag);
}


/* ---------------------------------------------------------------------- */

void FixExtField::init()
{
  if (domain->dimension == 2)
    error->all(FLERR,"Cannot use fix_extfield with 2d simulation");
  if (domain->nonperiodic != 0)
    error->all(FLERR,"Cannot use nonperiodic boundares with fix_extfield");
  if (domain->triclinic)
    error->all(FLERR,"Cannot use fix_extfield with triclinic box");

  if (strstr(update->integrate_style,"respa"))
    nlevels_respa = ((Respa *) update->integrate)->nlevels;
}



/* ----------------------------------------------------------------------
   read in initial electron temperatures from a user-specified file
   only called by proc 0
------------------------------------------------------------------------- */

void FixExtField::read_initial_electron_temperatures()
{
  char line[MAXLINE];

  for (int ixnode = 0; ixnode < nxnodes; ixnode++)
    for (int iynode = 0; iynode < nynodes; iynode++)
      for (int iznode = 0; iznode < nznodes; iznode++)
        T_initial_set[ixnode][iynode][iznode] = 0;

  // read initial electron temperature values from file

  int ixnode,iynode,iznode;
  double T_tmp;
  while (1) {
    if (fgets(line,MAXLINE,fpr) == NULL) break;
    sscanf(line,"%d %d %d %lg",&ixnode,&iynode,&iznode,&T_tmp);
    if (T_tmp < 0.0) error->one(FLERR,"fix_extfield electron temperatures must be > 0.0");
    T_electron[ixnode][iynode][iznode] = T_tmp;
    T_initial_set[ixnode][iynode][iznode] = 1;
  }

  for (int ixnode = 0; ixnode < nxnodes; ixnode++)
    for (int iynode = 0; iynode < nynodes; iynode++)
      for (int iznode = 0; iznode < nznodes; iznode++)
        if (T_initial_set[ixnode][iynode][iznode] == 0)
          error->one(FLERR,"Initial temperatures not all set in fix_extfield");

  // close file

  fclose(fpr);
}


/* ----------------------------------------------------------------------
   memory usage of 3d grid
------------------------------------------------------------------------- */

double FixExtField::memory_usage()
{
  double bytes = 0.0;
  bytes += 5*total_nnodes * sizeof(int);
  bytes += 14*total_nnodes * sizeof(double);
  return bytes;
}


/* ----------------------------------------------------------------------
  return the energy of the electronic subsystem or the net_energy transfer
   between the subsystems
------------------------------------------------------------------------- */

double FixExtField::compute_vector(int n)
{
  double e_energy = 0.0;

  double dx = domain->xprd/nxnodes;
  double dy = domain->yprd/nynodes;
  double dz = domain->zprd/nznodes;
  double del_vol = dx*dy*dz;

  for (int ixnode = 0; ixnode < nxnodes; ixnode++)
    for (int iynode = 0; iynode < nynodes; iynode++)
      for (int iznode = 0; iznode < nznodes; iznode++) {
        e_energy +=
          T_electron[ixnode][iynode][iznode];
  }

  if (n == 0) return e_energy;
  return 0.0;
}


