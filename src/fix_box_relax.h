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

#ifndef FIX_BOX_RELAX_H
#define FIX_BOX_RELAX_H

#include "fix.h"

namespace LAMMPS_NS {

class FixBoxRelax : public Fix {
 public:
  FixBoxRelax(class LAMMPS *, int, char **);
  ~FixBoxRelax();
  int setmask();
  void init();
  double min_energy(double *);
  void min_store();
  void min_step(double, double *);
  double max_alpha(double *);
  int min_dof();

 private:
  int p_flag[3];
  int press_couple,allremap;
  int dimension;
  double p_target[3],p_current[3];
  double dilation[3];
  double volinit,xprdinit,yprdinit,zprdinit;
  double smax,pv2e,pflagsum;

  double boxlo0[3],boxhi0[3];     // box bounds at start of line search
  double s0[6];	          	  // scale matrix in Voigt notation
  double ds[6];                   // increment in scale matrix

  char *id_temp,*id_press;
  class Compute *temperature,*pressure;
  int tflag,pflag;

  int nrigid;
  int *rfix;

  void remap();
  void couple();
};

}

#endif
