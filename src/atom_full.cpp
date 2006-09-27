/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   www.cs.sandia.gov/~sjplimp/lammps.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "atom_full.h"
#include "domain.h"
#include "modify.h"
#include "fix.h"

/* ---------------------------------------------------------------------- */

AtomFull::AtomFull(int narg, char **arg) : Atom(narg, arg) {}

/* ---------------------------------------------------------------------- */

void AtomFull::copy(int i, int j)
{
  int k;

  tag[j] = tag[i];
  type[j] = type[i];
  mask[j] = mask[i];
  image[j] = image[i];
  x[j][0] = x[i][0];
  x[j][1] = x[i][1];
  x[j][2] = x[i][2];
  v[j][0] = v[i][0];
  v[j][1] = v[i][1];
  v[j][2] = v[i][2];

  q[j] = q[i];

  molecule[j] = molecule[i];
  num_bond[j] = num_bond[i];
  num_angle[j] = num_angle[i];
  num_dihedral[j] = num_dihedral[i];
  num_improper[j] = num_improper[i];
  nspecial[j][0] = nspecial[i][0];
  nspecial[j][1] = nspecial[i][1];
  nspecial[j][2] = nspecial[i][2];

  for (k = 0; k < num_bond[j]; k++) {
    bond_type[j][k] = bond_type[i][k];
    bond_atom[j][k] = bond_atom[i][k];
  }

  for (k = 0; k < num_angle[j]; k++) {
    angle_type[j][k] = angle_type[i][k];
    angle_atom1[j][k] = angle_atom1[i][k];
    angle_atom2[j][k] = angle_atom2[i][k];
    angle_atom3[j][k] = angle_atom3[i][k];
  }

  for (k = 0; k < num_dihedral[j]; k++) {
    dihedral_type[j][k] = dihedral_type[i][k];
    dihedral_atom1[j][k] = dihedral_atom1[i][k];
    dihedral_atom2[j][k] = dihedral_atom2[i][k];
    dihedral_atom3[j][k] = dihedral_atom3[i][k];
    dihedral_atom4[j][k] = dihedral_atom4[i][k];
  }

  for (k = 0; k < num_improper[j]; k++) {
    improper_type[j][k] = improper_type[i][k];
    improper_atom1[j][k] = improper_atom1[i][k];
    improper_atom2[j][k] = improper_atom2[i][k];
    improper_atom3[j][k] = improper_atom3[i][k];
    improper_atom4[j][k] = improper_atom4[i][k];
  }

  for (k = 0; k < nspecial[j][2]; k++) special[j][k] = special[i][k];

  if (nextra_grow)
    for (int iextra = 0; iextra < nextra_grow; iextra++) 
      modify->fix[extra_grow[iextra]]->copy_arrays(i,j);
}

/* ---------------------------------------------------------------------- */

void AtomFull::pack_comm(int n, int *list, double *buf, int *pbc_flags)
{
  int i,j,m;

  m = 0;
  if (pbc_flags[0] == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0];
      buf[m++] = x[j][1];
      buf[m++] = x[j][2];
    }
  } else {
    double xprd = domain->xprd;
    double yprd = domain->yprd;
    double zprd = domain->zprd;
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0] + pbc_flags[1]*xprd;
      buf[m++] = x[j][1] + pbc_flags[2]*yprd;
      buf[m++] = x[j][2] + pbc_flags[3]*zprd;
    }
  }
}

/* ---------------------------------------------------------------------- */

void AtomFull::unpack_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

void AtomFull::pack_reverse(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    buf[m++] = f[i][0];
    buf[m++] = f[i][1];
    buf[m++] = f[i][2];
  }
}

/* ---------------------------------------------------------------------- */

void AtomFull::unpack_reverse(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    f[j][0] += buf[m++];
    f[j][1] += buf[m++];
    f[j][2] += buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

void AtomFull::pack_border(int n, int *list, double *buf, int *pbc_flags)
{
  int i,j,m;

  m = 0;
  if (pbc_flags[0] == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0];
      buf[m++] = x[j][1];
      buf[m++] = x[j][2];
      buf[m++] = tag[j];
      buf[m++] = type[j];
      buf[m++] = mask[j];
      buf[m++] = q[j];
      buf[m++] = molecule[j];
    }
  } else {
    double xprd = domain->xprd;
    double yprd = domain->yprd;
    double zprd = domain->zprd;
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0] + pbc_flags[1]*xprd;
      buf[m++] = x[j][1] + pbc_flags[2]*yprd;
      buf[m++] = x[j][2] + pbc_flags[3]*zprd;
      buf[m++] = tag[j];
      buf[m++] = type[j];
      buf[m++] = mask[j];
      buf[m++] = q[j];
      buf[m++] = molecule[j];
    }
  }
}

/* ---------------------------------------------------------------------- */

void AtomFull::unpack_border(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    if (i == nmax) grow(0);
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
    tag[i] = static_cast<int> (buf[m++]);
    type[i] = static_cast<int> (buf[m++]);
    mask[i] = static_cast<int> (buf[m++]);
    q[i] = buf[m++];
    molecule[i] = static_cast<int> (buf[m++]);
  }
}

/* ----------------------------------------------------------------------
   pack all atom quantities for shipping to another proc
   xyz must be 1st 3 values, so that comm::exchange can test on them
------------------------------------------------------------------------- */

int AtomFull::pack_exchange(int i, double *buf)
{
  int k;

  int m = 1;
  buf[m++] = x[i][0];
  buf[m++] = x[i][1];
  buf[m++] = x[i][2];
  buf[m++] = tag[i];
  buf[m++] = type[i];
  buf[m++] = mask[i];
  buf[m++] = image[i];
  buf[m++] = v[i][0];
  buf[m++] = v[i][1];
  buf[m++] = v[i][2];

  buf[m++] = q[i];

  buf[m++] = molecule[i];

  buf[m++] = num_bond[i];
  for (k = 0; k < num_bond[i]; k++) {
    buf[m++] = bond_type[i][k];
    buf[m++] = bond_atom[i][k];
  }

  buf[m++] = num_angle[i];
  for (k = 0; k < num_angle[i]; k++) {
    buf[m++] = angle_type[i][k];
    buf[m++] = angle_atom1[i][k];
    buf[m++] = angle_atom2[i][k];
    buf[m++] = angle_atom3[i][k];
  }

  buf[m++] = num_dihedral[i];
  for (k = 0; k < num_dihedral[i]; k++) {
    buf[m++] = dihedral_type[i][k];
    buf[m++] = dihedral_atom1[i][k];
    buf[m++] = dihedral_atom2[i][k];
    buf[m++] = dihedral_atom3[i][k];
    buf[m++] = dihedral_atom4[i][k];
  }

  buf[m++] = num_improper[i];
  for (k = 0; k < num_improper[i]; k++) {
    buf[m++] = improper_type[i][k];
    buf[m++] = improper_atom1[i][k];
    buf[m++] = improper_atom2[i][k];
    buf[m++] = improper_atom3[i][k];
    buf[m++] = improper_atom4[i][k];
  }

  buf[m++] = nspecial[i][0];
  buf[m++] = nspecial[i][1];
  buf[m++] = nspecial[i][2];
  for (k = 0; k < nspecial[i][2]; k++) buf[m++] = special[i][k];

  if (nextra_grow)
    for (int iextra = 0; iextra < nextra_grow; iextra++) 
      m += modify->fix[extra_grow[iextra]]->pack_exchange(i,&buf[m]);

  buf[0] = m;
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomFull::unpack_exchange(double *buf)
{
  int k;
  if (nlocal == nmax) grow(0);

  int m = 1;
  x[nlocal][0] = buf[m++];
  x[nlocal][1] = buf[m++];
  x[nlocal][2] = buf[m++];
  tag[nlocal] = static_cast<int> (buf[m++]);
  type[nlocal] = static_cast<int> (buf[m++]);
  mask[nlocal] = static_cast<int> (buf[m++]);
  image[nlocal] = static_cast<int> (buf[m++]);
  v[nlocal][0] = buf[m++];
  v[nlocal][1] = buf[m++];
  v[nlocal][2] = buf[m++];

  q[nlocal] = buf[m++];

  molecule[nlocal] = static_cast<int> (buf[m++]);

  num_bond[nlocal] = static_cast<int> (buf[m++]);
  for (k = 0; k < num_bond[nlocal]; k++) {
    bond_type[nlocal][k] = static_cast<int> (buf[m++]);
    bond_atom[nlocal][k] = static_cast<int> (buf[m++]);
  }

  num_angle[nlocal] = static_cast<int> (buf[m++]);
  for (k = 0; k < num_angle[nlocal]; k++) {
    angle_type[nlocal][k] = static_cast<int> (buf[m++]);
    angle_atom1[nlocal][k] = static_cast<int> (buf[m++]);
    angle_atom2[nlocal][k] = static_cast<int> (buf[m++]);
    angle_atom3[nlocal][k] = static_cast<int> (buf[m++]);
  }

  num_dihedral[nlocal] = static_cast<int> (buf[m++]);
  for (k = 0; k < num_dihedral[nlocal]; k++) {
    dihedral_type[nlocal][k] = static_cast<int> (buf[m++]);
    dihedral_atom1[nlocal][k] = static_cast<int> (buf[m++]);
    dihedral_atom2[nlocal][k] = static_cast<int> (buf[m++]);
    dihedral_atom3[nlocal][k] = static_cast<int> (buf[m++]);
    dihedral_atom4[nlocal][k] = static_cast<int> (buf[m++]);
  }

  num_improper[nlocal] = static_cast<int> (buf[m++]);
  for (k = 0; k < num_improper[nlocal]; k++) {
    improper_type[nlocal][k] = static_cast<int> (buf[m++]);
    improper_atom1[nlocal][k] = static_cast<int> (buf[m++]);
    improper_atom2[nlocal][k] = static_cast<int> (buf[m++]);
    improper_atom3[nlocal][k] = static_cast<int> (buf[m++]);
    improper_atom4[nlocal][k] = static_cast<int> (buf[m++]);
  }

  nspecial[nlocal][0] = static_cast<int> (buf[m++]);
  nspecial[nlocal][1] = static_cast<int> (buf[m++]);
  nspecial[nlocal][2] = static_cast<int> (buf[m++]);
  for (k = 0; k < nspecial[nlocal][2]; k++) 
    special[nlocal][k] = static_cast<int> (buf[m++]);

  if (nextra_grow)
    for (int iextra = 0; iextra < nextra_grow; iextra++) 
      m += modify->fix[extra_grow[iextra]]->unpack_exchange(nlocal,&buf[m]);

  nlocal++;
  return m;
}
