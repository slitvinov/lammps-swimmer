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

/* ----------------------------------------------------------------------
   Contributing authors: Leo Silbert (SNL), Gary Grest (SNL)
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "string.h"
#include "pair_gran_no_history.h"
#include "atom.h"
#include "force.h"
#include "neighbor.h"

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

/* ---------------------------------------------------------------------- */

void PairGranNoHistory::compute(int eflag, int vflag)
{
  int i,j,k,numneigh;
  double xtmp,ytmp,ztmp,delx,dely,delz;
  double radi,radj,radsum,rsq,r,rinv;
  double vr1,vr2,vr3,vnnr,vn1,vn2,vn3,vt1,vt2,vt3;
  double wr1,wr2,wr3;
  double vtr1,vtr2,vtr3,vrel;
  double xmeff,damp,ccel,ccelx,ccely,ccelz,tor1,tor2,tor3;
  double fn,fs,ft,fs1,fs2,fs3;
  int *neighs;

  double **f = atom->f;
  double **x = atom->x;
  double **v = atom->v;
  double **phiv = atom->phiv;
  double **phia = atom->phia;
  double *radius = atom->radius;
  double *rmass = atom->rmass;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;

  // loop over neighbors of my atoms

  for (i = 0; i < nlocal; i++) {
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    radi = radius[i];
    neighs = neighbor->firstneigh[i];
    numneigh = neighbor->numneigh[i];

    for (k = 0; k < numneigh; k++) {
      j = neighs[k];

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      radj = radius[j];
      radsum = radi + radj;

      if (rsq < radsum*radsum) {
	r = sqrt(rsq);

	// relative translational velocity

	vr1 = v[i][0] - v[j][0];
	vr2 = v[i][1] - v[j][1];
	vr3 = v[i][2] - v[j][2];

	vr1 *= dt;
	vr2 *= dt;
	vr3 *= dt;

	// normal component

	vnnr = vr1*delx + vr2*dely + vr3*delz;
	vn1 = delx*vnnr / rsq;
	vn2 = dely*vnnr / rsq;
	vn3 = delz*vnnr / rsq;

	// tangential component

	vt1 = vr1 - vn1;
	vt2 = vr2 - vn2;
	vt3 = vr3 - vn3;

	// relative rotational velocity

	wr1 = radi*phiv[i][0] + radj*phiv[j][0];
	wr2 = radi*phiv[i][1] + radj*phiv[j][1];
	wr3 = radi*phiv[i][2] + radj*phiv[j][2];

	wr1 *= dt/r;
	wr2 *= dt/r;
	wr3 *= dt/r;

	// normal damping term
	// this definition of DAMP includes the extra 1/r term

	xmeff = rmass[i]*rmass[j] / (rmass[i]+rmass[j]);
	if (mask[i] & freeze_group_bit) xmeff = rmass[j];
	if (mask[j] & freeze_group_bit) xmeff = rmass[i];
	damp = xmeff*gamman_dl*vnnr/rsq;
	ccel = xkk*(radsum-r)/r - damp;

	// relative velocities

	vtr1 = vt1 - (delz*wr2-dely*wr3);
	vtr2 = vt2 - (delx*wr3-delz*wr1);
	vtr3 = vt3 - (dely*wr1-delx*wr2);
	vrel = vtr1*vtr1 + vtr2*vtr2 + vtr3*vtr3;
	vrel = sqrt(vrel);

	// force normalization

	fn = xmu * fabs(ccel*r);
	fs = xmeff*gammas_dl*vrel;
	if (vrel != 0.0) ft = MIN(fn,fs) / vrel;
	else ft = 0.0;

	// shear friction forces

	fs1 = -ft*vtr1;
	fs2 = -ft*vtr2;
	fs3 = -ft*vtr3;

	// forces & torques

	ccelx = delx*ccel + fs1;
	ccely = dely*ccel + fs2;
	ccelz = delz*ccel + fs3;
	f[i][0] += ccelx;
	f[i][1] += ccely;
	f[i][2] += ccelz;

	rinv = 1/r;
	tor1 = rinv * (dely*fs3 - delz*fs2);
	tor2 = rinv * (delz*fs1 - delx*fs3);
	tor3 = rinv * (delx*fs2 - dely*fs1);
	phia[i][0] -= radi*tor1;
	phia[i][1] -= radi*tor2;
	phia[i][2] -= radi*tor3;

	if (newton_pair || j < nlocal) {
	  f[j][0] -= ccelx;
	  f[j][1] -= ccely;
	  f[j][2] -= ccelz;
	  phia[j][0] -= radj*tor1;
	  phia[j][1] -= radj*tor2;
	  phia[j][2] -= radj*tor3;
	}
      }
    }
  }
}
