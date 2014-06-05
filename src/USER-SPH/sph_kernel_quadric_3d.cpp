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

#include "sph_kernel_quadric_3d.h"
#include "math.h"

using namespace LAMMPS_NS;

double SPHKernelQuadric3D::w (double r, double h) {
  return 2.154187022708661*pow(r-h,4)*pow(r+h,4)/pow(h,11);
}

double SPHKernelQuadric3D::dw (double r, double h) {
  return 17.23349618166929*r*pow(r-h,3)*pow(r+h,3)/pow(h,11);
}

double SPHKernelQuadric3D::dw_per_r (double r, double h) {
  return 17.23349618166929*pow(r-h,3)*pow(r+h,3)/pow(h,11);
}
