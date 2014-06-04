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

#include "modify_kokkos.h"
#include "atom_kokkos.h"
#include "update.h"
#include "fix.h"
#include "compute.h"

using namespace LAMMPS_NS;

#define BIG 1.0e20

/* ---------------------------------------------------------------------- */

ModifyKokkos::ModifyKokkos(LAMMPS *lmp) : Modify(lmp) 
{
  atomKK = (AtomKokkos *) atom;
}

/* ----------------------------------------------------------------------
   setup for run, calls setup() of all fixes and computes
   called from Verlet, RESPA, Min
------------------------------------------------------------------------- */

void ModifyKokkos::setup(int vflag)
{
  // compute setup needs to come before fix setup
  // b/c NH fixes need use DOF of temperature computes

  for (int i = 0; i < ncompute; i++) compute[i]->setup();

  if (update->whichflag == 1)
    for (int i = 0; i < nfix; i++) {
      atomKK->sync(fix[i]->execution_space,fix[i]->datamask_read);
      atomKK->modified(fix[i]->execution_space,fix[i]->datamask_modify);
      fix[i]->setup(vflag);
    }
  else if (update->whichflag == 2)
    for (int i = 0; i < nfix; i++) {
      atomKK->sync(fix[i]->execution_space,fix[i]->datamask_read);
      atomKK->modified(fix[i]->execution_space,fix[i]->datamask_modify);
      fix[i]->min_setup(vflag);
    }
}

/* ----------------------------------------------------------------------
   setup pre_exchange call, only for fixes that define pre_exchange
   called from Verlet, RESPA, Min, and WriteRestart with whichflag = 0
------------------------------------------------------------------------- */

void ModifyKokkos::setup_pre_exchange()
{
  if (update->whichflag <= 1)
    for (int i = 0; i < n_pre_exchange; i++) {
      atomKK->sync(fix[list_pre_exchange[i]]->execution_space,
                   fix[list_pre_exchange[i]]->datamask_read);
      atomKK->modified(fix[list_pre_exchange[i]]->execution_space,
                       fix[list_pre_exchange[i]]->datamask_modify);
      fix[list_pre_exchange[i]]->setup_pre_exchange();
    }
  else if (update->whichflag == 2)
    for (int i = 0; i < n_min_pre_exchange; i++) {
      atomKK->sync(fix[list_min_pre_exchange[i]]->execution_space,
                   fix[list_min_pre_exchange[i]]->datamask_read);
      atomKK->modified(fix[list_min_pre_exchange[i]]->execution_space,
                       fix[list_min_pre_exchange[i]]->datamask_modify);
      fix[list_min_pre_exchange[i]]->min_setup_pre_exchange();
    }
}

/* ----------------------------------------------------------------------
   setup pre_neighbor call, only for fixes that define pre_neighbor
   called from Verlet, RESPA
------------------------------------------------------------------------- */

void ModifyKokkos::setup_pre_neighbor()
{
  if (update->whichflag == 1)
    for (int i = 0; i < n_pre_neighbor; i++) {
      atomKK->sync(fix[list_pre_neighbor[i]]->execution_space,
                   fix[list_pre_neighbor[i]]->datamask_read);
      atomKK->modified(fix[list_pre_neighbor[i]]->execution_space,
                       fix[list_pre_neighbor[i]]->datamask_modify);
      fix[list_pre_neighbor[i]]->setup_pre_neighbor();
    }
  else if (update->whichflag == 2)
    for (int i = 0; i < n_min_pre_neighbor; i++) {
      atomKK->sync(fix[list_min_pre_neighbor[i]]->execution_space,
                   fix[list_min_pre_neighbor[i]]->datamask_read);
      atomKK->modified(fix[list_min_pre_neighbor[i]]->execution_space,
                       fix[list_min_pre_neighbor[i]]->datamask_modify);
      fix[list_min_pre_neighbor[i]]->min_setup_pre_neighbor();
    }
}

/* ----------------------------------------------------------------------
   setup pre_force call, only for fixes that define pre_force
   called from Verlet, RESPA, Min
------------------------------------------------------------------------- */

void ModifyKokkos::setup_pre_force(int vflag)
{
  if (update->whichflag == 1)
    for (int i = 0; i < n_pre_force; i++) {
      atomKK->sync(fix[list_pre_force[i]]->execution_space,
                   fix[list_pre_force[i]]->datamask_read);
      atomKK->modified(fix[list_pre_force[i]]->execution_space,
                       fix[list_pre_force[i]]->datamask_modify);
      fix[list_pre_force[i]]->setup_pre_force(vflag);
    }
  else if (update->whichflag == 2)
    for (int i = 0; i < n_min_pre_force; i++) {
      atomKK->sync(fix[list_min_pre_force[i]]->execution_space,
                   fix[list_min_pre_force[i]]->datamask_read);
      atomKK->modified(fix[list_min_pre_force[i]]->execution_space,
                       fix[list_min_pre_force[i]]->datamask_modify);
      fix[list_min_pre_force[i]]->min_setup_pre_force(vflag);
    }
}

/* ----------------------------------------------------------------------
   1st half of integrate call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyKokkos::initial_integrate(int vflag)
{
  for (int i = 0; i < n_initial_integrate; i++) {
    atomKK->sync(fix[list_initial_integrate[i]]->execution_space,
                 fix[list_initial_integrate[i]]->datamask_read);
    atomKK->modified(fix[list_initial_integrate[i]]->execution_space,
                     fix[list_initial_integrate[i]]->datamask_modify);
    fix[list_initial_integrate[i]]->initial_integrate(vflag);
  }
}

/* ----------------------------------------------------------------------
   post_integrate call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyKokkos::post_integrate()
{
  for (int i = 0; i < n_post_integrate; i++) {
    atomKK->sync(fix[list_post_integrate[i]]->execution_space,
                 fix[list_post_integrate[i]]->datamask_read);
    atomKK->modified(fix[list_post_integrate[i]]->execution_space,
                     fix[list_post_integrate[i]]->datamask_modify);
    fix[list_post_integrate[i]]->post_integrate();
  }
}

/* ----------------------------------------------------------------------
   pre_exchange call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyKokkos::pre_exchange()
{
  for (int i = 0; i < n_pre_exchange; i++) {
    atomKK->sync(fix[list_pre_exchange[i]]->execution_space,
                 fix[list_pre_exchange[i]]->datamask_read);
    atomKK->modified(fix[list_pre_exchange[i]]->execution_space,
                     fix[list_pre_exchange[i]]->datamask_modify);
    fix[list_pre_exchange[i]]->pre_exchange();
  }
}

/* ----------------------------------------------------------------------
   pre_neighbor call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyKokkos::pre_neighbor()
{
  for (int i = 0; i < n_pre_neighbor; i++) {
    atomKK->sync(fix[list_pre_neighbor[i]]->execution_space,
                 fix[list_pre_neighbor[i]]->datamask_read);
    atomKK->modified(fix[list_pre_neighbor[i]]->execution_space,
                     fix[list_pre_neighbor[i]]->datamask_modify);
    fix[list_pre_neighbor[i]]->pre_neighbor();
  }
}

/* ----------------------------------------------------------------------
   pre_force call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyKokkos::pre_force(int vflag)
{
  for (int i = 0; i < n_pre_force; i++) {
    atomKK->sync(fix[list_pre_force[i]]->execution_space,
                 fix[list_pre_force[i]]->datamask_read);
    atomKK->modified(fix[list_pre_force[i]]->execution_space,
                     fix[list_pre_force[i]]->datamask_modify);
    fix[list_pre_force[i]]->pre_force(vflag);
  }
}

/* ----------------------------------------------------------------------
   post_force call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyKokkos::post_force(int vflag)
{
  for (int i = 0; i < n_post_force; i++) {
    atomKK->sync(fix[list_post_force[i]]->execution_space,
                 fix[list_post_force[i]]->datamask_read);
    atomKK->modified(fix[list_post_force[i]]->execution_space,
                     fix[list_post_force[i]]->datamask_modify);
    fix[list_post_force[i]]->post_force(vflag);
  }
}

/* ----------------------------------------------------------------------
   2nd half of integrate call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyKokkos::final_integrate()
{
  for (int i = 0; i < n_final_integrate; i++) {
    atomKK->sync(fix[list_final_integrate[i]]->execution_space,
                 fix[list_final_integrate[i]]->datamask_read);
    atomKK->modified(fix[list_final_integrate[i]]->execution_space,
                     fix[list_final_integrate[i]]->datamask_modify);
    fix[list_final_integrate[i]]->final_integrate();
  }
}

/* ----------------------------------------------------------------------
   end-of-timestep call, only for relevant fixes
   only call fix->end_of_step() on timesteps that are multiples of nevery
------------------------------------------------------------------------- */

void ModifyKokkos::end_of_step()
{
  for (int i = 0; i < n_end_of_step; i++)
    if (update->ntimestep % end_of_step_every[i] == 0) {
      atomKK->sync(fix[list_end_of_step[i]]->execution_space,
                   fix[list_end_of_step[i]]->datamask_read);
      atomKK->modified(fix[list_end_of_step[i]]->execution_space,
                       fix[list_end_of_step[i]]->datamask_modify);
      fix[list_end_of_step[i]]->end_of_step();
    }
}

/* ----------------------------------------------------------------------
   thermo energy call, only for relevant fixes
   called by Thermo class
   compute_scalar() is fix call to return energy
------------------------------------------------------------------------- */

double ModifyKokkos::thermo_energy()
{
  double energy = 0.0;
  for (int i = 0; i < n_thermo_energy; i++) {
    atomKK->sync(fix[list_thermo_energy[i]]->execution_space,
                 fix[list_thermo_energy[i]]->datamask_read);
    atomKK->modified(fix[list_thermo_energy[i]]->execution_space,
                     fix[list_thermo_energy[i]]->datamask_modify);
    energy += fix[list_thermo_energy[i]]->compute_scalar();
  }
  return energy;
}

/* ----------------------------------------------------------------------
   post_run call
------------------------------------------------------------------------- */

void ModifyKokkos::post_run()
{
  for (int i = 0; i < nfix; i++) {
    atomKK->sync(fix[i]->execution_space,
                 fix[i]->datamask_read);
    atomKK->modified(fix[i]->execution_space,
                     fix[i]->datamask_modify);
    fix[i]->post_run();
  }
}

/* ----------------------------------------------------------------------
   setup rRESPA pre_force call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyKokkos::setup_pre_force_respa(int vflag, int ilevel)
{
  for (int i = 0; i < n_pre_force; i++) {
    atomKK->sync(fix[list_pre_force[i]]->execution_space,
                 fix[list_pre_force[i]]->datamask_read);
    atomKK->modified(fix[list_pre_force[i]]->execution_space,
                     fix[list_pre_force[i]]->datamask_modify);
    fix[list_pre_force[i]]->setup_pre_force_respa(vflag,ilevel);
  }
}

/* ----------------------------------------------------------------------
   1st half of rRESPA integrate call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyKokkos::initial_integrate_respa(int vflag, int ilevel, int iloop)
{
  for (int i = 0; i < n_initial_integrate_respa; i++) {
    atomKK->sync(fix[list_initial_integrate_respa[i]]->execution_space,
                 fix[list_initial_integrate_respa[i]]->datamask_read);
    atomKK->modified(fix[list_initial_integrate_respa[i]]->execution_space,
                     fix[list_initial_integrate_respa[i]]->datamask_modify);
    fix[list_initial_integrate_respa[i]]->
      initial_integrate_respa(vflag,ilevel,iloop);
  }
}

/* ----------------------------------------------------------------------
   rRESPA post_integrate call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyKokkos::post_integrate_respa(int ilevel, int iloop)
{
  for (int i = 0; i < n_post_integrate_respa; i++) {
    atomKK->sync(fix[list_post_integrate_respa[i]]->execution_space,
                 fix[list_post_integrate_respa[i]]->datamask_read);
    atomKK->modified(fix[list_post_integrate_respa[i]]->execution_space,
                     fix[list_post_integrate_respa[i]]->datamask_modify);
    fix[list_post_integrate_respa[i]]->post_integrate_respa(ilevel,iloop);
  }
}

/* ----------------------------------------------------------------------
   rRESPA pre_force call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyKokkos::pre_force_respa(int vflag, int ilevel, int iloop)
{
  for (int i = 0; i < n_pre_force_respa; i++) {
    atomKK->sync(fix[list_pre_force_respa[i]]->execution_space,
                 fix[list_pre_force_respa[i]]->datamask_read);
    atomKK->modified(fix[list_pre_force_respa[i]]->execution_space,
                     fix[list_pre_force_respa[i]]->datamask_modify);
    fix[list_pre_force_respa[i]]->pre_force_respa(vflag,ilevel,iloop);
  }
}

/* ----------------------------------------------------------------------
   rRESPA post_force call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyKokkos::post_force_respa(int vflag, int ilevel, int iloop)
{
  for (int i = 0; i < n_post_force_respa; i++) {
    atomKK->sync(fix[list_post_force_respa[i]]->execution_space,
                 fix[list_post_force_respa[i]]->datamask_read);
    atomKK->modified(fix[list_post_force_respa[i]]->execution_space,
                     fix[list_post_force_respa[i]]->datamask_modify);
    fix[list_post_force_respa[i]]->post_force_respa(vflag,ilevel,iloop);
  }
}

/* ----------------------------------------------------------------------
   2nd half of rRESPA integrate call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyKokkos::final_integrate_respa(int ilevel, int iloop)
{
  for (int i = 0; i < n_final_integrate_respa; i++) {
    atomKK->sync(fix[list_final_integrate_respa[i]]->execution_space,
                 fix[list_final_integrate_respa[i]]->datamask_read);
    atomKK->modified(fix[list_final_integrate_respa[i]]->execution_space,
                     fix[list_final_integrate_respa[i]]->datamask_modify);
    fix[list_final_integrate_respa[i]]->final_integrate_respa(ilevel,iloop);
  }
}

/* ----------------------------------------------------------------------
   minimizer pre-exchange call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyKokkos::min_pre_exchange()
{
  for (int i = 0; i < n_min_pre_exchange; i++) {
    atomKK->sync(fix[list_min_pre_exchange[i]]->execution_space,
                 fix[list_min_pre_exchange[i]]->datamask_read);
    atomKK->modified(fix[list_min_pre_exchange[i]]->execution_space,
                     fix[list_min_pre_exchange[i]]->datamask_modify);
    fix[list_min_pre_exchange[i]]->min_pre_exchange();
  }
}

/* ----------------------------------------------------------------------
   minimizer pre-neighbor call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyKokkos::min_pre_neighbor()
{
  for (int i = 0; i < n_min_pre_neighbor; i++) {
    atomKK->sync(fix[list_min_pre_neighbor[i]]->execution_space,
                 fix[list_min_pre_neighbor[i]]->datamask_read);
    atomKK->modified(fix[list_min_pre_neighbor[i]]->execution_space,
                     fix[list_min_pre_neighbor[i]]->datamask_modify);
    fix[list_min_pre_neighbor[i]]->min_pre_neighbor();
  }
}

/* ----------------------------------------------------------------------
   minimizer pre-force call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyKokkos::min_pre_force(int vflag)
{
  for (int i = 0; i < n_min_pre_force; i++) {
    atomKK->sync(fix[list_min_pre_force[i]]->execution_space,
                 fix[list_min_pre_force[i]]->datamask_read);
    atomKK->modified(fix[list_min_pre_force[i]]->execution_space,
                     fix[list_min_pre_force[i]]->datamask_modify);
    fix[list_min_pre_force[i]]->min_pre_force(vflag);
  }
}

/* ----------------------------------------------------------------------
   minimizer force adjustment call, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyKokkos::min_post_force(int vflag)
{
  for (int i = 0; i < n_min_post_force; i++) {
    atomKK->sync(fix[list_min_post_force[i]]->execution_space,
                 fix[list_min_post_force[i]]->datamask_read);
    atomKK->modified(fix[list_min_post_force[i]]->execution_space,
                     fix[list_min_post_force[i]]->datamask_modify);
    fix[list_min_post_force[i]]->min_post_force(vflag);
  }
}

/* ----------------------------------------------------------------------
   minimizer energy/force evaluation, only for relevant fixes
   return energy and forces on extra degrees of freedom
------------------------------------------------------------------------- */

double ModifyKokkos::min_energy(double *fextra)
{
  int ifix,index;

  index = 0;
  double eng = 0.0;
  for (int i = 0; i < n_min_energy; i++) {
    ifix = list_min_energy[i];
    atomKK->sync(fix[ifix]->execution_space,fix[ifix]->datamask_read);
    atomKK->modified(fix[ifix]->execution_space,fix[ifix]->datamask_modify);
    eng += fix[ifix]->min_energy(&fextra[index]);
    index += fix[ifix]->min_dof();
  }
  return eng;
}

/* ----------------------------------------------------------------------
   store current state of extra dof, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyKokkos::min_store()
{
  for (int i = 0; i < n_min_energy; i++) {
    atomKK->sync(fix[list_min_energy[i]]->execution_space,
                 fix[list_min_energy[i]]->datamask_read);
    atomKK->modified(fix[list_min_energy[i]]->execution_space,
                     fix[list_min_energy[i]]->datamask_modify);
    fix[list_min_energy[i]]->min_store();
  }
}

/* ----------------------------------------------------------------------
   mange state of extra dof on a stack, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyKokkos::min_clearstore()
{
  for (int i = 0; i < n_min_energy; i++) {
    atomKK->sync(fix[list_min_energy[i]]->execution_space,
                 fix[list_min_energy[i]]->datamask_read);
    atomKK->modified(fix[list_min_energy[i]]->execution_space,
                     fix[list_min_energy[i]]->datamask_modify);
    fix[list_min_energy[i]]->min_clearstore();
  }
}

void ModifyKokkos::min_pushstore()
{
  for (int i = 0; i < n_min_energy; i++) {
    atomKK->sync(fix[list_min_energy[i]]->execution_space,
                 fix[list_min_energy[i]]->datamask_read);
    atomKK->modified(fix[list_min_energy[i]]->execution_space,
                     fix[list_min_energy[i]]->datamask_modify);
    fix[list_min_energy[i]]->min_pushstore();
  }
}

void ModifyKokkos::min_popstore()
{
  for (int i = 0; i < n_min_energy; i++) {
    atomKK->sync(fix[list_min_energy[i]]->execution_space,
                 fix[list_min_energy[i]]->datamask_read);
    atomKK->modified(fix[list_min_energy[i]]->execution_space,
                     fix[list_min_energy[i]]->datamask_modify);
    fix[list_min_energy[i]]->min_popstore();
  }
}

/* ----------------------------------------------------------------------
   displace extra dof along vector hextra, only for relevant fixes
------------------------------------------------------------------------- */

void ModifyKokkos::min_step(double alpha, double *hextra)
{
  int ifix,index;

  index = 0;
  for (int i = 0; i < n_min_energy; i++) {
    ifix = list_min_energy[i];
    atomKK->sync(fix[ifix]->execution_space,fix[ifix]->datamask_read);
    atomKK->modified(fix[ifix]->execution_space,fix[ifix]->datamask_modify);
    fix[ifix]->min_step(alpha,&hextra[index]);
    index += fix[ifix]->min_dof();
  }
}

/* ----------------------------------------------------------------------
   compute max allowed step size along vector hextra, only for relevant fixes
------------------------------------------------------------------------- */

double ModifyKokkos::max_alpha(double *hextra)
{
  int ifix,index;

  double alpha = BIG;
  index = 0;
  for (int i = 0; i < n_min_energy; i++) {
    ifix = list_min_energy[i];
    atomKK->sync(fix[ifix]->execution_space,fix[ifix]->datamask_read);
    atomKK->modified(fix[ifix]->execution_space,fix[ifix]->datamask_modify);
    double alpha_one = fix[ifix]->max_alpha(&hextra[index]);
    alpha = MIN(alpha,alpha_one);
    index += fix[ifix]->min_dof();
  }
  return alpha;
}

/* ----------------------------------------------------------------------
   extract extra dof for minimization, only for relevant fixes
------------------------------------------------------------------------- */

int ModifyKokkos::min_dof()
{
  int ndof = 0;
  for (int i = 0; i < n_min_energy; i++) {
    atomKK->sync(fix[list_min_energy[i]]->execution_space,
                 fix[list_min_energy[i]]->datamask_read);
    atomKK->modified(fix[list_min_energy[i]]->execution_space,
                     fix[list_min_energy[i]]->datamask_modify);
    ndof += fix[list_min_energy[i]]->min_dof();
  }
  return ndof;
}

/* ----------------------------------------------------------------------
   reset reference state of fix, only for relevant fixes
------------------------------------------------------------------------- */

int ModifyKokkos::min_reset_ref()
{
  int itmp,itmpall;
  itmpall = 0;
  for (int i = 0; i < n_min_energy; i++) {
    atomKK->sync(fix[list_min_energy[i]]->execution_space,
                 fix[list_min_energy[i]]->datamask_read);
    atomKK->modified(fix[list_min_energy[i]]->execution_space,
                     fix[list_min_energy[i]]->datamask_modify);
    itmp = fix[list_min_energy[i]]->min_reset_ref();
    if (itmp) itmpall = 1;
  }
  return itmpall;
}