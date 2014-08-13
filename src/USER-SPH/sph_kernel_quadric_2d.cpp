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

#include "sph_kernel_quadric_2d.h"
#include "math.h"

using namespace LAMMPS_NS;

double SPHKernelQuadric2D::w (double r, double h) {
  return 1.591549430918954*pow(r-h,4)*pow(r+h,4)/pow(h,10);
}

double SPHKernelQuadric2D::dw (double r, double h) {
  return 12.73239544735162*r*pow(r-h,3)*pow(r+h,3)/pow(h,10);
}

double SPHKernelQuadric2D::dw_per_r (double r, double h) {
  return 12.73239544735162*pow(r-h,3)*pow(r+h,3)/pow(h,10);
}
