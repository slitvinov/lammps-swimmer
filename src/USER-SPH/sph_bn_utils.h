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

#ifndef LMP_SPH_BN_UTILS_H
#define LMP_SPH_BN_UTILS_H

#include "pointers.h"
#include "domain.h"

namespace  LAMMPS_NS {
  void read_initial_target_field(FILE* fpr,
				 int nxnodes, int nynodes, int nznodes,
				 int ***T_initial_set, double ***T_target,
				 Error *error);
  double get_target_field (double* xi, Domain *&domain, double ***T_target,
			   int nxnodes, int nynodes, int nznodes);
  
  double get_target_cutoff (double m, int nn, double rhot);
};

#endif
