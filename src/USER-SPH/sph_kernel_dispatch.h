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

#ifndef LMP_SPH_DISPATCH_H
#define LMP_SPH_DISPATCH_H

namespace  LAMMPS_NS {
  class SPHKernel;
  enum SPHKernelCodeType { Lucy2D };

  SPHKernelCodeType sph_kernel_code(char* kernel_name, int kernel_dimension,
				    Error *error);

  SPHKernel* sph_kernel_decode(SPHKernelCodeType kernel_code,
			       Error *error);
};

#endif
