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

#include "sph_kernel_wendldand4_2d.h"
#include "math.h"
using namespace LAMMPS_NS;

double SPHKernelWendland42D::w (double r, double h) {
  double norm2d = 2.864788975654116/(h*h);
  double s = r/h;
  if (s<1.0) {
    return norm2d*pow(1.0-s, 6)*(35.0*s*s/3.0+6.0*s+1.0);
  }
  return 0.0;
}

double SPHKernelWendland42D::dw (double r, double h) {
  double norm2d = 2.864788975654116/(h*h*h);
  double s = r/h;
  double wfd;
  if (s<1) {
    wfd = -56*pow(1.0-s, 5)*s*(5.0*s+1.0)/3.0;
  } else {
    wfd = 0.0;
  }
  return norm2d*wfd;
}

double SPHKernelWendland42D::dw_per_r (double r, double h) {
  return dw(r, h)/r;
}
