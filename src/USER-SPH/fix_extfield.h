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

#ifdef FIX_CLASS

FixStyle(extfield,FixExtField)

#else

#ifndef LMP_FIX_EXT_FIELD_H
#define LMP_FIX_EXT_FIELD_H

#include "fix.h"

namespace LAMMPS_NS {

class FixExtField : public Fix {
 public:
  FixExtField(class LAMMPS *, int, char **);
  ~FixExtField();
  void post_force(int);
  void post_force_respa(int, int, int);
  int setmask();
  void init();
  double memory_usage();
  double compute_vector(int);

 private:
  int me;
  int nlevels_respa;
  FILE *fpr;
  int nxnodes,nynodes,nznodes,total_nnodes;
  int ***T_initial_set;
  double ***T_electron;

  void read_initial_electron_temperatures();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Cannot open file %s

The specified file cannot be opened.  Check that the path and name are
correct. If the file is a compressed file, also check that the gzip
executable can be found and run.

E: Cannot open fix ttm file %s

The output file for the fix ttm command cannot be opened.  Check that
the path and name are correct.

E: Invalid random number seed in fix ttm command

Random number seed must be > 0.

E: Fix ttm electronic_specific_heat must be > 0.0

Self-explanatory.

E: Fix ttm electronic_density must be > 0.0

Self-explanatory.

E: Fix ttm electronic_thermal_conductivity must be >= 0.0

Self-explanatory.

E: Fix ttm gamma_p must be > 0.0

Self-explanatory.

E: Fix ttm gamma_s must be >= 0.0

Self-explanatory.

E: Fix ttm v_0 must be >= 0.0

Self-explanatory.

E: Fix ttm number of nodes must be > 0

Self-explanatory.

E: Cannot use fix ttm with 2d simulation

This is a current restriction of this fix due to the grid it creates.

E: Cannot use nonperiodic boundares with fix ttm

This fix requires a fully periodic simulation box.

E: Cannot use fix ttm with triclinic box

This is a current restriction of this fix due to the grid it creates.

E: Electronic temperature dropped below zero

Something has gone wrong with the fix ttm electron temperature model.

E: Fix ttm electron temperatures must be > 0.0

Self-explanatory.

E: Initial temperatures not all set in fix ttm

Self-explantory.

W: Too many inner timesteps in fix ttm

Self-explanatory.

*/
