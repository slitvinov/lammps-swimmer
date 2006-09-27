/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   www.cs.sandia.gov/~sjplimp/lammps.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef FIX_INSERT_H
#define FIX_INSERT_H

#include "fix.h"

class RanPark;

class FixInsert : public Fix {
  friend class PairGranHistory;
  friend class PairGranHertzian;
  friend class PairGranNoHistory;

 public:
  FixInsert(int, char **);
  ~FixInsert();
  int setmask();
  void init();
  void pre_exchange();

 private:
  int ninsert,ntype,seed;
  double radius_lo,radius_hi;
  double density_lo,density_hi;
  double volfrac;
  int maxattempt;
  int region_style;
  double rate;
  double vxlo,vxhi,vylo,vyhi,vy,vz;
  double xlo,xhi,ylo,yhi,zlo,zhi;
  double xc,yc,rc;

  int me,nprocs;
  int *recvcounts,*displs;
  double PI;
  int nfreq,nfirst,ninserted,nper;
  double lo_current,hi_current;
  int ifix_history;

  RanPark *random;

  int overlap(int);
  void xyz_random(double, double &, double &, double &);
};

#endif
