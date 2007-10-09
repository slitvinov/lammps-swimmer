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

#ifndef FIX_INDENT_H
#define FIX_INDENT_H

#include "fix.h"

namespace LAMMPS_NS {

class FixIndent : public Fix {
 public:
  FixIndent(class LAMMPS *, int, char **);
  int setmask();
  void init();
  void setup();
  void min_setup();
  virtual void post_force(int);
  void post_force_respa(int, int, int);
  void min_post_force(int);
  double compute_scalar();
  double compute_vector(int);

 protected:
  int istyle,scaleflag,radflag,thermo_flag,eflag_enable;
  double k,k3;
  double x0,y0,z0,r0_stop,r0_start;
  int indenter_flag;
  double indenter[4],indenter_all[4];
  int cdim;
  double c1,c2;
  double vx,vy,vz;
  int nlevels_respa;

  void options(int, char **);
};

}

#endif
