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

#ifndef LMP_SPH_KERNEL_LIQ_2D_H
#define LMP_SPH_KERNEL_LIQ_2D_H
#include "sph_kernel.h"
#include "math.h"

namespace LAMMPS_NS {
  namespace sph_kernel_liq_2d
  {
    const double xs = 0.3;
    const double A = 1.0/(2*xs*xs*xs-6*xs*xs+6*xs-2);
    const double B = -(xs+1)/(xs*xs*xs-3*xs*xs+3*xs-1);
    const double C = 3.0*xs/(xs*xs*xs-3*xs*xs+3*xs-1);
    const double D = -(3*xs-1)/(xs*xs*xs-3*xs*xs+3*xs-1);
    const double E = (2*xs-1)/(2*xs*xs*xs-6*xs*xs+6*xs-2);
    const double F = (xs+1)/2.0;
    const double norm2d = 30.0/(2*xs*xs*xs*M_PI+3*xs*xs*M_PI+3*xs*M_PI+2*M_PI);
  }

  class SPHKernelLIQ2D : public SPHKernel {
  public:
    SPHKernelLIQ2D() {};
    virtual double w  (double r, double h);
    virtual double dw (double r, double h);
    virtual double dw_per_r (double r, double h);
  };
}

#endif
