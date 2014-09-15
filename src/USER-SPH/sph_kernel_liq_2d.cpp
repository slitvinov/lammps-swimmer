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

#include "sph_kernel_liq_2d.h"
#include "math.h"
using namespace LAMMPS_NS;
using namespace sph_kernel_liq_2d;

double SPHKernelLIQ2D::w (double r, double h) {
  double s = r/h;
  double w;
  if (s<xs) {
    w = F-s;
  } else if (s<1)  {
    w = E+s*D+s*s*C+s*s*s*B+pow(s,4)*A;
  } else {
    w = 0;
  }
  return norm2d*w/(h*h);
}

double SPHKernelLIQ2D::dw (double r, double h) {
  double s = r/h;
  double wfd;
  if (s<xs) {
    wfd = -1;
  } else if (s<1) {
    wfd = D+2*s*C+3*s*s*B+4*s*s*s*A;
  } else {
    wfd = 0;
  }
  return norm2d*wfd/(h*h*h);
}

double SPHKernelLIQ2D::dw_per_r (double r, double h) {
  return dw(r, h)/r;
}
