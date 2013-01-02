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

#include "stdlib.h"
#include "string.h"
#include "kspace.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "pair.h"
#include "memory.h"
#include "atom_masks.h"
#include "error.h"
#include "suffix.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

KSpace::KSpace(LAMMPS *lmp, int narg, char **arg) : Pointers(lmp)
{
  energy = 0.0;
  virial[0] = virial[1] = virial[2] = virial[3] = virial[4] = virial[5] = 0.0;

  ewaldflag = pppmflag = msmflag = dispersionflag = tip4pflag = 0;
  compute_flag = 1;
  group_group_enable = 0;

  order = 5;
  gridflag = 0;
  gewaldflag = 0;
  minorder = 2;
  overlap_allowed = 1;

  order_6 = 5;
  gridflag_6 = 0;
  gewaldflag_6 = 0;
    
  slabflag = 0;
  differentiation_flag = 0;
  slab_volfactor = 1;
  suffix_flag = Suffix::NONE;
  adjust_cutoff_flag = 1;

  accuracy_absolute = -1.0;
  two_charge_force = force->qqr2e *
    (force->qelectron * force->qelectron) /
    (force->angstrom * force->angstrom);

  maxeatom = maxvatom = 0;
  eatom = NULL;
  vatom = NULL;

  datamask = ALL_MASK;
  datamask_ext = ALL_MASK;

  memory->create(gcons,7,7,"kspace:gcons");
  gcons[2][0] = 15.0 / 8.0;
  gcons[2][1] = -5.0 / 4.0;
  gcons[2][2] = 3.0 / 8.0;
  gcons[3][0] = 35.0 / 16.0;
  gcons[3][1] = -35.0 / 16.0;
  gcons[3][2] = 21.0 / 16.0;
  gcons[3][3] = -5.0 / 16.0;
  gcons[4][0] = 315.0 / 128.0;
  gcons[4][1] = -105.0 / 32.0;
  gcons[4][2] = 189.0 / 64.0;
  gcons[4][3] = -45.0 / 32.0;
  gcons[4][4] = 35.0 / 128.0;
  gcons[5][0] = 693.0 / 256.0;
  gcons[5][1] = -1155.0 / 256.0;
  gcons[5][2] = 693.0 / 128.0;
  gcons[5][3] = -495.0 / 128.0;
  gcons[5][4] = 385.0 / 256.0;
  gcons[5][5] = -63.0 / 256.0;
  gcons[6][0] = 3003.0 / 1024.0;
  gcons[6][1] = -3003.0 / 512.0;
  gcons[6][2] = 9009.0 / 1024.0;
  gcons[6][3] = -2145.0 / 256.0;
  gcons[6][4] = 5005.0 / 1024.0;
  gcons[6][5] = -819.0 / 512.0;
  gcons[6][6] = 231.0 / 1024.0;

  memory->create(dgcons,7,6,"kspace:dgcons");
  dgcons[2][0] = -5.0 / 2.0;
  dgcons[2][1] = 3.0 / 2.0;
  dgcons[3][0] = -35.0 / 8.0;
  dgcons[3][1] = 21.0 / 4.0;
  dgcons[3][2] = -15.0 / 8.0;
  dgcons[4][0] = -105.0 / 16.0;
  dgcons[4][1] = 189.0 / 16.0;
  dgcons[4][2] = -135.0 / 16.0;
  dgcons[4][3] = 35.0 / 16.0;
  dgcons[5][0] = -1155.0 / 128.0;
  dgcons[5][1] = 693.0 / 32.0;
  dgcons[5][2] = -1485.0 / 64.0;
  dgcons[5][3] = 385.0 / 32.0;
  dgcons[5][4] = -315.0 / 128.0;
  dgcons[6][0] = -3003.0 / 256.0;
  dgcons[6][1] = 9009.0 / 256.0;
  dgcons[6][2] = -6435.0 / 128.0;
  dgcons[6][3] = 5005.0 / 128.0;
  dgcons[6][4] = -4095.0 / 256.0;
  dgcons[6][5] = 693.0 / 256.0;
}

/* ---------------------------------------------------------------------- */

KSpace::~KSpace()
{
  memory->destroy(eatom);
  memory->destroy(vatom);
  memory->destroy(gcons);
  memory->destroy(dgcons);
}

/* ---------------------------------------------------------------------- */

void KSpace::compute_dummy(int eflag, int vflag)
{
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = evflag_atom = eflag_global = vflag_global =
         eflag_atom = vflag_atom = 0;
}

/* ----------------------------------------------------------------------
   check that pair style is compatible with long-range solver
------------------------------------------------------------------------- */

void KSpace::pair_check()
{
  if (force->pair == NULL)
    error->all(FLERR,"KSpace solver requires a pair style");
  if (ewaldflag && force->pair->ewaldflag == 0)
    error->all(FLERR,"KSpace style is incompatible with Pair style");
  if (pppmflag && force->pair->pppmflag == 0)
    error->all(FLERR,"KSpace style is incompatible with Pair style");
  if (msmflag && force->pair->msmflag == 0)
    error->all(FLERR,"KSpace style is incompatible with Pair style");
  if (dispersionflag && force->pair->dispersionflag == 0)
    error->all(FLERR,"KSpace style is incompatible with Pair style");
  if (tip4pflag && force->pair->tip4pflag == 0)
    error->all(FLERR,"KSpace style is incompatible with Pair style");
}

/* ----------------------------------------------------------------------
   setup for energy, virial computation
   see integrate::ev_set() for values of eflag (0-3) and vflag (0-6)
------------------------------------------------------------------------- */

void KSpace::ev_setup(int eflag, int vflag)
{
  int i,n;

  evflag = 1;

  eflag_either = eflag;
  eflag_global = eflag % 2;
  eflag_atom = eflag / 2;

  vflag_either = vflag;
  vflag_global = vflag % 4;
  vflag_atom = vflag / 4;

  if (eflag_atom || vflag_atom) evflag_atom = 1;
  else evflag_atom = 0;

  // reallocate per-atom arrays if necessary

  if (eflag_atom && atom->nlocal > maxeatom) {
    maxeatom = atom->nmax;
    memory->destroy(eatom);
    memory->create(eatom,maxeatom,"kspace:eatom");
  }
  if (vflag_atom && atom->nlocal > maxvatom) {
    maxvatom = atom->nmax;
    memory->destroy(vatom);
    memory->create(vatom,maxvatom,6,"kspace:vatom");
  }

  // zero accumulators

  if (eflag_global) energy = 0.0;
  if (vflag_global) for (i = 0; i < 6; i++) virial[i] = 0.0;
  if (eflag_atom) {
    n = atom->nlocal;
    for (i = 0; i < n; i++) eatom[i] = 0.0;
  }
  if (vflag_atom) {
    n = atom->nlocal;
    for (i = 0; i < n; i++) {
      vatom[i][0] = 0.0;
      vatom[i][1] = 0.0;
      vatom[i][2] = 0.0;
      vatom[i][3] = 0.0;
      vatom[i][4] = 0.0;
      vatom[i][5] = 0.0;
    }
  }
}

/* ----------------------------------------------------------------------
   estimate the accuracy of the short-range coulomb tables
------------------------------------------------------------------------- */

double KSpace::estimate_table_accuracy(double q2_over_sqrt, double spr)
{
  double table_accuracy = 0.0;
  int nctb = force->pair->ncoultablebits;
  if (nctb) {
    double empirical_precision[17];
    empirical_precision[6] =  6.99E-03;
    empirical_precision[7] =  1.78E-03;
    empirical_precision[8] =  4.72E-04;
    empirical_precision[9] =  1.17E-04;
    empirical_precision[10] = 2.95E-05;
    empirical_precision[11] = 7.41E-06;
    empirical_precision[12] = 1.76E-06;
    empirical_precision[13] = 9.28E-07;
    empirical_precision[14] = 7.46E-07;
    empirical_precision[15] = 7.32E-07;
    empirical_precision[16] = 7.30E-07;
    if (nctb <= 6) table_accuracy = empirical_precision[6];
    else if (nctb <= 16) table_accuracy = empirical_precision[nctb];
    else table_accuracy = empirical_precision[16];
    table_accuracy *= q2_over_sqrt;
    if ((table_accuracy > spr) && (comm->me == 0))
      error->warning(FLERR,"For better accuracy use 'pair_modify table 0'");
  }

  return table_accuracy;
}

/* ----------------------------------------------------------------------
   compute gamma for MSM and pair styles
   see Eq 4 from Parallel Computing 35 (2009) 164�177
------------------------------------------------------------------------- */

double KSpace::gamma(const double &rho)
{
  if (rho <= 1) {
    int split_order = order/2;
    double g = gcons[split_order][0];
    double rho2 = rho*rho;
    double rho_n = rho2;
    for (int n=1; n<=split_order; n++) {
      g += gcons[split_order][n]*rho_n;
      rho_n *= rho2;
    }
    return g;
  }
  else
    return (1.0/rho);
}

/* ----------------------------------------------------------------------
   compute the derivative of gamma for MSM and pair styles
   see Eq 4 from Parallel Computing 35 (2009) 164-177
------------------------------------------------------------------------- */

double KSpace::dgamma(const double &rho)
{
  if (rho <= 1) {
    int split_order = order/2;
    double dg = dgcons[split_order][0]*rho;
    double rho2 = rho*rho;
    double rho_n = rho*rho2;
    for (int n=1; n<split_order; n++) {
      dg += dgcons[split_order][n]*rho_n;
      rho_n *= rho2;
    }
    return dg;
  }
  else
    return (-1.0/rho/rho);
}

/* ----------------------------------------------------------------------
   modify parameters of the KSpace style
------------------------------------------------------------------------- */

void KSpace::modify_params(int narg, char **arg)
{
  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"mesh") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal kspace_modify command");
      nx_pppm = nx_msm_max = atoi(arg[iarg+1]);
      ny_pppm = ny_msm_max = atoi(arg[iarg+2]);
      nz_pppm = nz_msm_max = atoi(arg[iarg+3]);
      if (nx_pppm == 0 && ny_pppm == 0 && nz_pppm == 0) gridflag = 0;
      else gridflag = 1;
      iarg += 4;
    } else if (strcmp(arg[iarg],"mesh/disp") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal kspace_modify command");
      nx_pppm_6 = atoi(arg[iarg+1]);
      ny_pppm_6 = atoi(arg[iarg+2]);
      nz_pppm_6 = atoi(arg[iarg+3]);
      if (nx_pppm_6 == 0 || ny_pppm_6 == 0 || nz_pppm_6 == 0) gridflag_6 = 0;
      else gridflag_6 = 1;
      iarg += 4;
    } else if (strcmp(arg[iarg],"order") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal kspace_modify command");
      order = atoi(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"order/disp") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal kspace_modify command");
      order_6 = atoi(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"minorder") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal kspace_modify command");
      minorder = atoi(arg[iarg+1]);
      if (minorder < 2) error->all(FLERR,"Illegal kspace_modify command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"overlap") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal kspace_modify command");
      if (strcmp(arg[iarg+1],"yes") == 0) overlap_allowed = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) overlap_allowed = 0;
      else error->all(FLERR,"Illegal kspace_modify command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"force") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal kspace_modify command");
      accuracy_absolute = atof(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"gewald") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal kspace_modify command");
      g_ewald = atof(arg[iarg+1]);
      if (g_ewald == 0.0) gewaldflag = 0;
      else gewaldflag = 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"gewald/disp") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal kspace_modify command");
      g_ewald_6 = atof(arg[iarg+1]);
      if (g_ewald_6 == 0.0) gewaldflag_6 = 0;
      else gewaldflag_6 = 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"slab") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal kspace_modify command");
      if (strcmp(arg[iarg+1],"nozforce") == 0) {
        slabflag = 2;
      } else {
        slabflag = 1;
        slab_volfactor = atof(arg[iarg+1]);
        if (slab_volfactor <= 1.0)
          error->all(FLERR,"Bad kspace_modify slab parameter");
        if (slab_volfactor < 2.0 && comm->me == 0)
          error->warning(FLERR,"Kspace_modify slab param < 2.0 may "
                         "cause unphysical behavior");
      }
      iarg += 2;
    } else if (strcmp(arg[iarg],"compute") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal kspace_modify command");
      if (strcmp(arg[iarg+1],"yes") == 0) compute_flag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) compute_flag = 0;
      else error->all(FLERR,"Illegal kspace_modify command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"diff") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal kspace_modify command");
      if (strcmp(arg[iarg+1],"ad") == 0) differentiation_flag = 1;
      else if (strcmp(arg[iarg+1],"ik") == 0) differentiation_flag = 0;
      else error->all(FLERR, "Illegal kspace_modify command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"cutoff/adjust") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal kspace_modify command");
      if (strcmp(arg[iarg+1],"yes") == 0) adjust_cutoff_flag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) adjust_cutoff_flag = 0;
      else error->all(FLERR,"Illegal kspace_modify command");
      iarg += 2;
    } else error->all(FLERR,"Illegal kspace_modify command");
  }
}

/* ---------------------------------------------------------------------- */

void *KSpace::extract(const char *str)
{
  if (strcmp(str,"scale") == 0) return (void *) &scale;
  return NULL;
}
