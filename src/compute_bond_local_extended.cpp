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
#include "compute_bond_local_extended.h"
#include "atom.h"
#include "atom_vec.h"
#include "molecule.h"
#include "update.h"
#include "domain.h"
#include "force.h"
#include "bond.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define DELTA 10000

enum{DIST,ENG,FORCE,POWER};

/* ---------------------------------------------------------------------- */

ComputeBondLocalExtended::ComputeBondLocalExtended(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg < 4) error->all(FLERR,"Illegal compute bond/local/extended command");

  if (atom->avec->bonds_allow == 0)
    error->all(FLERR,"Compute bond/local/extended used when bonds are not allowed");

  local_flag = 1;
  nvalues = narg - 3;
  if (nvalues == 1) size_local_cols = 0;
  else size_local_cols = nvalues;

  bstyle = new int[nvalues];

  nvalues = 0;
  for (int iarg = 3; iarg < narg; iarg++) {
    if (strcmp(arg[iarg],"dist") == 0) bstyle[nvalues++] = DIST;
    else if (strcmp(arg[iarg],"eng") == 0) bstyle[nvalues++] = ENG;
    else if (strcmp(arg[iarg],"force") == 0) bstyle[nvalues++] = FORCE;
    else if (strcmp(arg[iarg],"power") == 0) bstyle[nvalues++] = POWER;
    else error->all(FLERR,"Invalid keyword in compute bond/local/extended command");
  }

  // set singleflag if need to call bond->single()

  singleflag = 0;
  for (int i = 0; i < nvalues; i++)
    if (bstyle[i] != DIST) singleflag = 1;

  nmax = 0;
  vector = NULL;
  array = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeBondLocalExtended::~ComputeBondLocalExtended()
{
  memory->destroy(vector);
  memory->destroy(array);
  delete [] bstyle;
}

/* ---------------------------------------------------------------------- */

void ComputeBondLocalExtended::init()
{
  if (force->bond == NULL)
    error->all(FLERR,"No bond style is defined for compute bond/local/extended");

  // do initial memory allocation so that memory_usage() is correct

  ncount = compute_bonds(0);
  if (ncount > nmax) reallocate(ncount);
  size_local_rows = ncount;
}

/* ---------------------------------------------------------------------- */

void ComputeBondLocalExtended::compute_local()
{
  invoked_local = update->ntimestep;

  // count local entries and compute bond info

  ncount = compute_bonds(0);
  if (ncount > nmax) reallocate(ncount);
  size_local_rows = ncount;
  ncount = compute_bonds(1);
}

/* ----------------------------------------------------------------------
   count bonds and compute bond info on this proc
   only count bond once if newton_bond is off
   all atoms in interaction must be in group
   all atoms in interaction must be known to proc
   if bond is deleted (type = 0), do not count
   if bond is turned off (type < 0), still count
   if flag is set, compute requested info about bond
   if bond is turned off (type < 0), energy = 0.0
------------------------------------------------------------------------- */

int ComputeBondLocalExtended::compute_bonds(int flag)
{
  int i,m,n,nb,atom1,atom2,imol,iatom,btype;
  tagint tagprev;
  double delx,dely,delz,rsq;
  double *ptr;
  double dtfm1, dtfm2;

  double **x = atom->x;
  double **fb = atom->fb;
  double **v = atom->v;
  tagint *tag = atom->tag;
  int *num_bond = atom->num_bond;
  tagint **bond_atom = atom->bond_atom;
  int **bond_type = atom->bond_type;
  int *type = atom->type;
  int *mask = atom->mask;

  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int rmass_flag = atom->rmass_flag;
  double dtf = 0.5 * update->dt * force->ftm2v;

  int *molindex = atom->molindex;
  int *molatom = atom->molatom;
  Molecule **onemols = atom->avec->onemols;

  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;
  int molecular = atom->molecular;

  Bond *bond = force->bond;
  double eng,fbond;

  m = n = 0;
  for (atom1 = 0; atom1 < nlocal; atom1++) {
    if (!(mask[atom1] & groupbit)) continue;

    if (molecular == 1) nb = num_bond[atom1];
    else {
      if (molindex[atom1] < 0) continue;
      imol = molindex[atom1];
      iatom = molatom[atom1];
      nb = onemols[imol]->num_bond[iatom];
    }

    for (i = 0; i < nb; i++) {
      if (molecular == 1) {
        btype = bond_type[atom1][i];
        atom2 = atom->map(bond_atom[atom1][i]);
      } else {
        tagprev = tag[atom1] - iatom - 1;
        btype = atom->map(onemols[imol]->bond_type[atom1][i]);
        atom2 = atom->map(onemols[imol]->bond_atom[atom1][i]+tagprev);
      }

      if (atom2 < 0 || !(mask[atom2] & groupbit)) continue;
      if (newton_bond == 0 && tag[atom1] > tag[atom2]) continue;
      if (btype == 0) continue;

      if (flag) {
        delx = x[atom1][0] - x[atom2][0];
        dely = x[atom1][1] - x[atom2][1];
        delz = x[atom1][2] - x[atom2][2];
        domain->minimum_image(delx,dely,delz);
        rsq = delx*delx + dely*dely + delz*delz;

        if (singleflag) {
          if (btype > 0)
            eng = bond->single(btype,rsq,atom1,atom2,fbond);
          else eng = fbond = 0.0;
        }

        if (nvalues == 1) ptr = &vector[m];
        else ptr = array[m];

        for (n = 0; n < nvalues; n++) {
          switch (bstyle[n]) {
          case DIST:
            ptr[n] = sqrt(rsq);
            break;
          case ENG:
            ptr[n] = eng;
            break;
          case FORCE:
            ptr[n] = sqrt(rsq)*fbond;
            break;
          case POWER:
	    double vtmp1[3];
	    double vtmp2[3];
	    if (atom->fb_flag==1) { 
	      if (rmass_flag) {
		dtfm1 = dtf / rmass[atom1];
		dtfm2 = dtf / rmass[atom2];
	      } else {
		dtfm1 = dtf / mass[type[atom1]];
		dtfm2 = dtf / mass[type[atom2]];
	      }
	      vtmp1[0] = v[atom1][0] + dtfm1 * fb[atom1][0];
              vtmp1[1] = v[atom1][1] + dtfm1 * fb[atom1][1];
              vtmp1[2] = v[atom1][2] + dtfm1 * fb[atom1][2];
	      
	      vtmp2[0] = v[atom2][0] + dtfm2 * fb[atom2][0];
	      vtmp2[1] = v[atom2][1] + dtfm2 * fb[atom2][1];
	      vtmp2[2] = v[atom2][2] + dtfm2 * fb[atom2][2];

	    } else {
	      vtmp1[0] = v[atom1][0];
              vtmp1[1] = v[atom1][1];
              vtmp1[2] = v[atom1][2];
	      
	      vtmp2[0] = v[atom2][0];
              vtmp2[1] = v[atom2][1];
              vtmp2[2] = v[atom2][2];
	    }
	    
	    double dvx = vtmp1[0] - vtmp2[0];
	    double dvy = vtmp1[1] - vtmp2[1];
	    double dvz = vtmp1[2] - vtmp2[2];
            ptr[n] = (delz*dvz+dely*dvy+delx*dvx)*fbond;
            break;
          }
        }
      }

      m++;
    }
  }

  return m;
}

/* ---------------------------------------------------------------------- */

void ComputeBondLocalExtended::reallocate(int n)
{
  // grow vector or array and indices array

  while (nmax < n) nmax += DELTA;

  if (nvalues == 1) {
    memory->destroy(vector);
    memory->create(vector,nmax,"bond/local/extended:vector");
    vector_local = vector;
  } else {
    memory->destroy(array);
    memory->create(array,nmax,nvalues,"bond/local/extended:array");
    array_local = array;
  }
}

/* ----------------------------------------------------------------------
   memory usage of local data
------------------------------------------------------------------------- */

double ComputeBondLocalExtended::memory_usage()
{
  double bytes = nmax*nvalues * sizeof(double);
  return bytes;
}
