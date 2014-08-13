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

#include "sph_kernel_lucy_3d.h"
using namespace LAMMPS_NS;

double SPHKernelLucy3D::w (double r, double h) {
  return -2.088908628081126*(r-h)*(r-h)*(r-h)*(3*r+h)/(h*h*h*h*h*h*h);
}

double SPHKernelLucy3D::dw (double r, double h) {
  return -25.06690353697351*r*(r-h)*(r-h)/(h*h*h*h*h*h*h);
}

double SPHKernelLucy3D::dw_per_r (double r, double h) {
  return -25.06690353697351*(r-h)*(r-h)/(h*h*h*h*h*h*h);
}
