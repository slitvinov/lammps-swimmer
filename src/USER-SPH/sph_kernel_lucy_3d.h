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

#ifndef LMP_SPH_KERNEL_LUCY_3D_H
#define LMP_SPH_KERNEL_LUCY_3D_H
#include "sph_kernel.h"

namespace LAMMPS_NS {
  class SPHKernelLucy3D : public SPHKernel {
  public:
    SPHKernelLucy3D() {};
    virtual double w  (double r, double h);
    virtual double dw (double r, double h);
    virtual double dw_per_r (double r, double h);
  };
}

#endif