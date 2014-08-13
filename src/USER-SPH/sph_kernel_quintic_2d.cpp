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

#include "sph_kernel_quintic_2d.h"
#include "math.h"
using namespace LAMMPS_NS;

double SPHKernelQuintic2D::w (double r, double h) {
  double norm2d = 0.04195297663091802/(h*h);
  double s = 3.0*r/h;
  if (s<1.0) {
    return norm2d*(pow(3 - s, 5) - 6*pow(2 - s, 5) + 15*pow(1 - s, 5));
  } else if (s<2.0) {
    return norm2d*(pow(3 - s, 5) - 6*pow(2 - s, 5));
  } else if (s<3.0) {
    return norm2d*pow(3 - s, 5);
  }
  return 0.0;
}

double SPHKernelQuintic2D::dw (double r, double h) {
  double norm2d = 3.0*0.04195297663091802/(h*h*h);
  double s = 3.0*r/h;
  double wfd;
  if (s<1) {
    wfd = -50*pow(s,4)+120*pow(s,3)-120*s;
  } else if (s<2) {
    wfd = 25*pow(s,4)-180*pow(s,3)+450*pow(s,2)-420*s+75;
  } else if (s<3) {
    wfd = -5*pow(s,4)+60*pow(s,3)-270*pow(s,2)+540*s-405;
  } else {
    wfd = 0.0;
  }
  return norm2d*wfd;
}

double SPHKernelQuintic2D::dw_per_r (double r, double h) {
  return dw(r, h)/r;
}
