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

/* ----------------------------------------------------------------------
   Contributing author: Mike Brown (SNL)
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "pair_resquared.h"
#include "math_extra.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "integrate.h"
#include "memory.h"
#include "error.h"

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

using namespace LAMMPS_NS;

enum{SPHERE_SPHERE,SPHERE_ELLIPSE,ELLIPSE_SPHERE,ELLIPSE_ELLIPSE};

/* ---------------------------------------------------------------------- */

PairRESquared::PairRESquared(LAMMPS *lmp) : Pair(lmp), 
					    b_alpha(45.0/56.0),
                                            cr60(pow(60.0,1.0/3.0))
{
  single_enable = 0;
  cr60 = pow(60.0,1.0/3.0);
  b_alpha = 45.0/56.0;
  solv_f_a = 3.0/(16.0*atan(1.0)*-36.0);
  solv_f_r = 3.0/(16.0*atan(1.0)*2025.0);
}

/* ----------------------------------------------------------------------
   free all arrays
------------------------------------------------------------------------- */

PairRESquared::~PairRESquared()
{
  if (allocated) {
    memory->destroy_2d_int_array(setflag);
    memory->destroy_2d_double_array(cutsq);

    memory->destroy_2d_int_array(form);
    memory->destroy_2d_double_array(epsilon);
    memory->destroy_2d_double_array(sigma);
    memory->destroy_2d_double_array(shape2);
    memory->destroy_2d_double_array(well);
    memory->destroy_2d_double_array(cut);
    memory->destroy_2d_double_array(lj1);
    memory->destroy_2d_double_array(lj2);
    memory->destroy_2d_double_array(lj3);
    memory->destroy_2d_double_array(lj4);
    memory->destroy_2d_double_array(offset);
    delete [] lshape;
    delete [] setwell;
  }
}

/* ---------------------------------------------------------------------- */

void PairRESquared::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double evdwl,one_eng,rsq,r2inv,r6inv,forcelj,factor_lj;
  double fforce[3],ttor[3],rtor[3],r12[3];
  int *ilist,*jlist,*numneigh,**firstneigh;
  RE2Vars wi,wj;

  evdwl = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  double **x = atom->x;
  double **f = atom->f;
  double **tor = atom->torque;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  
  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itype = type[i];

    // not a LJ sphere

    if (lshape[itype] != 0.0) precompute_i(i,wi);

    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];

      if (j < nall) factor_lj = 1.0;
      else {
        factor_lj = special_lj[j/nall];
        j %= nall;
      }

      // r12 = center to center vector

      r12[0] = x[j][0]-x[i][0];
      r12[1] = x[j][1]-x[i][1];
      r12[2] = x[j][2]-x[i][2];
      rsq = MathExtra::dot3(r12,r12);
      jtype = type[j];

      // compute if less than cutoff

      if (rsq < cutsq[itype][jtype]) {
        switch (form[itype][jtype]) {

         case SPHERE_SPHERE:
          r2inv = 1.0/rsq;
          r6inv = r2inv*r2inv*r2inv;
          forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);
          forcelj *= -r2inv;
          if (eflag) one_eng =
              r6inv*(r6inv*lj3[itype][jtype]-lj4[itype][jtype]) -
              offset[itype][jtype];
          fforce[0] = r12[0]*forcelj;
          fforce[1] = r12[1]*forcelj;
          fforce[2] = r12[2]*forcelj;
          break;

         case SPHERE_ELLIPSE:
          precompute_i(j,wj);
          if (newton_pair || j < nlocal) {
            one_eng = resquared_lj(j,i,wj,r12,rsq,fforce,rtor,true);
            tor[j][0] += rtor[0]*factor_lj;
            tor[j][1] += rtor[1]*factor_lj;
            tor[j][2] += rtor[2]*factor_lj;
          } else
            one_eng = resquared_lj(j,i,wj,r12,rsq,fforce,rtor,false);
          break;

         case ELLIPSE_SPHERE:
          one_eng = resquared_lj(i,j,wi,r12,rsq,fforce,ttor,true);
          tor[i][0] += ttor[0]*factor_lj;
          tor[i][1] += ttor[1]*factor_lj;
          tor[i][2] += ttor[2]*factor_lj;
          break;

         default:
          precompute_i(j,wj);
          one_eng = resquared_analytic(i,j,wi,wj,r12,rsq,fforce,ttor,rtor);
          tor[i][0] += ttor[0]*factor_lj;
          tor[i][1] += ttor[1]*factor_lj;
          tor[i][2] += ttor[2]*factor_lj;
          if (newton_pair || j < nlocal) {
            tor[j][0] += rtor[0]*factor_lj;
            tor[j][1] += rtor[1]*factor_lj;
            tor[j][2] += rtor[2]*factor_lj;
          }
         break;
        }

        fforce[0] *= factor_lj;
        fforce[1] *= factor_lj;
        fforce[2] *= factor_lj;
        f[i][0] += fforce[0];
        f[i][1] += fforce[1];
        f[i][2] += fforce[2];

        if (newton_pair || j < nlocal) {
          f[j][0] -= fforce[0];
          f[j][1] -= fforce[1];
          f[j][2] -= fforce[2];
        }

        if (eflag) evdwl = factor_lj*one_eng;

	if (evflag) ev_tally_xyz(i,j,nlocal,newton_pair,
				 evdwl,0.0,fforce[0],fforce[1],fforce[2],
				 -r12[0],-r12[1],-r12[2]);
      }
    }
  }

  if (vflag_fdotr) virial_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairRESquared::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  setflag = memory->create_2d_int_array(n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  cutsq = memory->create_2d_double_array(n+1,n+1,"pair:cutsq");

  form = memory->create_2d_int_array(n+1,n+1,"pair:form");
  epsilon = memory->create_2d_double_array(n+1,n+1,"pair:epsilon");
  sigma = memory->create_2d_double_array(n+1,n+1,"pair:sigma");
  shape2 = memory->create_2d_double_array(n+1,3,"pair:shape2");
  well = memory->create_2d_double_array(n+1,3,"pair:well");
  cut = memory->create_2d_double_array(n+1,n+1,"pair:cut");
  lj1 = memory->create_2d_double_array(n+1,n+1,"pair:lj1");
  lj2 = memory->create_2d_double_array(n+1,n+1,"pair:lj2");
  lj3 = memory->create_2d_double_array(n+1,n+1,"pair:lj3");
  lj4 = memory->create_2d_double_array(n+1,n+1,"pair:lj4");
  offset = memory->create_2d_double_array(n+1,n+1,"pair:offset");

  lshape = new double[n+1];
  setwell = new int[n+1];
  for (int i = 1; i <= n; i++) setwell[i] = 0;
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairRESquared::settings(int narg, char **arg)
{
  if (narg != 1) error->all("Illegal pair_style command");

  cut_global = atof(arg[0]);
  
  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i+1; j <= atom->ntypes; j++)
        if (setflag[i][j]) cut[i][j] = cut_global;
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairRESquared::coeff(int narg, char **arg)
{
  if (narg < 10 || narg > 11)
    error->all("Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(arg[0],atom->ntypes,ilo,ihi);
  force->bounds(arg[1],atom->ntypes,jlo,jhi);

  double epsilon_one = atof(arg[2]);
  double sigma_one = atof(arg[3]);
  double eia_one = atof(arg[4]);
  double eib_one = atof(arg[5]);
  double eic_one = atof(arg[6]);
  double eja_one = atof(arg[7]);
  double ejb_one = atof(arg[8]);
  double ejc_one = atof(arg[9]);
  
  double cut_one = cut_global;
  if (narg == 11) cut_one = atof(arg[10]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      epsilon[i][j] = epsilon_one;
      sigma[i][j] = sigma_one;
      cut[i][j] = cut_one;
      if (eia_one != 0.0 || eib_one != 0.0 || eic_one != 0.0) {
        well[i][0] = eia_one;
        well[i][1] = eib_one;
        well[i][2] = eic_one;
        if (eia_one == 1.0 && eib_one == 1.0 && eic_one == 1.0) setwell[i] = 2;
        else setwell[i] = 1;
      }
      if (eja_one != 0.0 || ejb_one != 0.0 || ejc_one != 0.0) {
        well[j][0] = eja_one;
        well[j][1] = ejb_one;
        well[j][2] = ejc_one;
        if (eja_one == 1.0 && ejb_one == 1.0 && ejc_one == 1.0) setwell[j] = 2;
        else setwell[j] = 1;
      }
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all("Incorrect args for pair coefficients");
}


/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairRESquared::init_style()
{
  if (!atom->quat_flag || !atom->torque_flag)
    error->all("Pair resquared requires atom attributes quat, torque");

  int irequest = neighbor->request(this);

  // per-type shape precalculations

  for (int i = 1; i <= atom->ntypes; i++) {
    if (setwell[i]) {
      double *one = atom->shape[i];
      shape2[i][0] = one[0]*one[0];
      shape2[i][1] = one[1]*one[1];
      shape2[i][2] = one[2]*one[2];
      lshape[i] = one[0]*one[1]*one[2];
    }
  }
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairRESquared::init_one(int i, int j)
{
  double **shape = atom->shape;
  
  if (setwell[i] == 0 || setwell[j] == 0)
    error->all("Pair resquared epsilon a,b,c coeffs are not all set");

  int ishape = 0;
  if (shape[i][0] != 0 && shape[i][1] != 0 && shape[i][2] != 0) ishape = 1;
  int jshape = 0;
  if (shape[j][0] != 0 && shape[j][1] != 0 && shape[j][2] != 0) jshape = 1;
  
  if (ishape == 0 && jshape == 0) {
    form[i][j] = SPHERE_SPHERE;
    form[j][i] = SPHERE_SPHERE;
  } else if (ishape == 0) {
    form[i][j] = SPHERE_ELLIPSE;
    form[j][i] = ELLIPSE_SPHERE;
  } else if (jshape == 0) {
    form[i][j] = ELLIPSE_SPHERE;
    form[j][i] = SPHERE_ELLIPSE;
  } else {
    form[i][j] = ELLIPSE_ELLIPSE;
    form[j][i] = ELLIPSE_ELLIPSE;
  }

  // allow mixing only for LJ spheres

  if (setflag[i][j] == 0) {
    if (setflag[j][i] == 0) {
      if (ishape == 0 && jshape == 0) {
        epsilon[i][j] = mix_energy(epsilon[i][i],epsilon[j][j],
                                   sigma[i][i],sigma[j][j]);
        sigma[i][j] = mix_distance(sigma[i][i],sigma[j][j]);
        cut[i][j] = mix_distance(cut[i][i],cut[j][j]);
      } else
        error->all("Pair resquared epsilon and sigma coeffs are not all set");
    }
    epsilon[i][j] = epsilon[j][i];
    sigma[i][j] = sigma[j][i];
    cut[i][j] = cut[j][i];
  }
  
  lj1[i][j] = 48.0 * epsilon[i][j] * pow(sigma[i][j],12.0);
  lj2[i][j] = 24.0 * epsilon[i][j] * pow(sigma[i][j],6.0);
  lj3[i][j] = 4.0 * epsilon[i][j] * pow(sigma[i][j],12.0);
  lj4[i][j] = 4.0 * epsilon[i][j] * pow(sigma[i][j],6.0);
     
  if (offset_flag) {
    double ratio = sigma[i][j] / cut[i][j];
    offset[i][j] = 4.0 * epsilon[i][j] * (pow(ratio,12.0) - pow(ratio,6.0));
  } else offset[i][j] = 0.0;

  epsilon[j][i] = epsilon[i][j];
  sigma[j][i] = sigma[i][j];
  cut[j][i] = cut[i][j];
  lj1[j][i] = lj1[i][j];
  lj2[j][i] = lj2[i][j];
  lj3[j][i] = lj3[i][j];
  lj4[j][i] = lj4[i][j];
  offset[j][i] = offset[i][j];

  return cut[i][j];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairRESquared::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++) {
    fwrite(&setwell[i],sizeof(int),1,fp);
    if (setwell[i]) fwrite(&well[i][0],sizeof(double),3,fp);
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&epsilon[i][j],sizeof(double),1,fp);
        fwrite(&sigma[i][j],sizeof(double),1,fp);
        fwrite(&cut[i][j],sizeof(double),1,fp);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairRESquared::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++) {
    if (me == 0) fread(&setwell[i],sizeof(int),1,fp);
    MPI_Bcast(&setwell[i],1,MPI_INT,0,world);
    if (setwell[i]) {
      if (me == 0) fread(&well[i][0],sizeof(double),3,fp);
      MPI_Bcast(&well[i][0],3,MPI_DOUBLE,0,world);
    }
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) fread(&setflag[i][j],sizeof(int),1,fp);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0) {
          fread(&epsilon[i][j],sizeof(double),1,fp);
          fread(&sigma[i][j],sizeof(double),1,fp);
          fread(&cut[i][j],sizeof(double),1,fp);
        }
        MPI_Bcast(&epsilon[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&sigma[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairRESquared::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairRESquared::read_restart_settings(FILE *fp)
{
  int me = comm->me;
  if (me == 0) {
    fread(&cut_global,sizeof(double),1,fp);
    fread(&mix_flag,sizeof(int),1,fp);
  }
  MPI_Bcast(&cut_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
}

/* ----------------------------------------------------------------------
   Precompute per-particle temporaries for RE-squared calculation
------------------------------------------------------------------------- */

void PairRESquared::precompute_i(const int i,RE2Vars &ws) {
  double aTs[3][3];       // A1'*S1^2

  MathExtra::quat_to_mat_trans(atom->quat[i],ws.A);
  MathExtra::transpose_times_diag3(ws.A,well[atom->type[i]],ws.aTe);
  MathExtra::transpose_times_diag3(ws.A,shape2[atom->type[i]],aTs);
  MathExtra::diag_times3(shape2[atom->type[i]],ws.A,ws.sa);
  MathExtra::times3(aTs,ws.A,ws.gamma);
  MathExtra::rotation_generator_x(ws.A,ws.lA[0]);
  MathExtra::rotation_generator_y(ws.A,ws.lA[1]);
  MathExtra::rotation_generator_z(ws.A,ws.lA[2]);
  for (int i=0; i<3; i++) {
    MathExtra::times3(aTs,ws.lA[i],ws.lAtwo[i]);
    MathExtra::transpose_times3(ws.lA[i],ws.sa,ws.lAsa[i]);
    MathExtra::plus3(ws.lAsa[i],ws.lAtwo[i],ws.lAsa[i]);
  }
}
                                      
/* ----------------------------------------------------------------------
   Compute the derivative of the determinant of m, using m and the
   derivative of m (m2)
------------------------------------------------------------------------- */

double PairRESquared::det_prime(const double m[3][3], const double m2[3][3])
{
  double ans;
  ans = m2[0][0]*m[1][1]*m[2][2] - m2[0][0]*m[1][2]*m[2][1] -
        m[1][0]*m2[0][1]*m[2][2] + m[1][0]*m2[0][2]*m[2][1] +
        m[2][0]*m2[0][1]*m[1][2] - m[2][0]*m2[0][2]*m[1][1] +
        m[0][0]*m2[1][1]*m[2][2] - m[0][0]*m2[1][2]*m[2][1] -
        m2[1][0]*m[0][1]*m[2][2] + m2[1][0]*m[0][2]*m[2][1] +
        m[2][0]*m[0][1]*m2[1][2] - m[2][0]*m[0][2]*m2[1][1] +
        m[0][0]*m[1][1]*m2[2][2] - m[0][0]*m[1][2]*m2[2][1] -
        m[1][0]*m[0][1]*m2[2][2] + m[1][0]*m[0][2]*m2[2][1] +
        m2[2][0]*m[0][1]*m[1][2] - m2[2][0]*m[0][2]*m[1][1];
  return ans;
}

/* ----------------------------------------------------------------------
   Compute the energy, force, torque for a pair (INTEGRATED-INTEGRATED)
------------------------------------------------------------------------- */

double PairRESquared::resquared_analytic(const int i, const int j,
                                         const RE2Vars &wi, const RE2Vars &wj,
                                         const double *r, const double rsq,
                                         double *fforce, double *ttor,
                                         double *rtor)
{
  int *type = atom->type;
  double **shape = atom->shape;
    
  // pair computations for energy, force, torque

  double z1[3],z2[3];        // A1*rhat  # don't need to store
  double v1[3],v2[3];        // inv(S1^2)*z1 # don't need to store
  double sigma1,sigma2;      // 1/sqrt(z1'*v1)
  double sigma1p2,sigma2p2;  // sigma1^2
  double rnorm;              // L2 norm of r
  double rhat[3];            // r/rnorm
  double s[3];               // inv(gamma1+gamma2)*rhat
  double sigma12;            // 1/sqrt(0.5*s'*rhat)
  double H12[3][3];          // gamma1/sigma1+gamma2/sigma2
  double dH;                 // det(H12)
  double lambda;             // dS1/sigma1p2+dS2/sigma2p2
  double nu;                 // sqrt(dH/(sigma1+sigma2))
  double w[3];               // inv(A1'*E1*A1+A2'*E2*A2)*rhat
  double h12;                // rnorm-sigma12;
  double eta;                // lambda/nu
  double chi;                // 2*rhat'*w
  double sprod;              // dS1*dS2
  double sigh;               // sigma/h12
  double tprod;              // eta*chi*sigh
  double Ua,Ur;              // attractive/repulsive parts of potential
    
  // pair computations for force, torque

  double sec;                          // sigma*eta*chi
  double sigma1p3, sigma2p3;           // sigma1^3
  double vsigma1[3], vsigma2[3];       // sigma1^3*v1;
  double sigma12p3;                    // sigma12^3
  double gsigma1[3][3], gsigma2[3][3]; // -gamma1/sigma1^2
  double tsig1sig2;                    // eta/(2*(sigma1+sigma2))
  double tdH;                          // eta/(2*dH)
  double teta1,teta2;                  // 2*eta/lambda*dS1/sigma1p3
  double fourw[3];                     // 4*w;
  double spr[3];                       // 0.5*sigma12^3*s
  double hsec;                         // h12+[3,b_alpha]*sec
  double dspu;                         // 1/h12 - 1/hsec + temp
  double pbsu;                         // 3*sigma/hsec
  double dspr;                         // 7/h12-1/hsec+temp
  double pbsr;                         // b_alpha*sigma/hsec;
  double u[3];                         // (-rhat(i)*rhat+eye(:,i))/rnorm
  double u1[3],u2[3];                  // A1*u
  double dsigma1,dsigma2;              // u1'*vsigma1 (force) p'*vsigma1 (tor)
  double dH12[3][3];                   // dsigma1*gsigma1 + dsigma2*gsigma2
  double ddH;                          // derivative of det(H12)
  double deta,dchi,dh12;               // derivatives of eta,chi,h12
  double dUr,dUa;                      // derivatives of Ua,Ur
      
  // pair computations for torque

  double fwae[3];        // -fourw'*aTe
  double p[3];           // lA*rhat

  rnorm = sqrt(rsq);
  rhat[0] = r[0]/rnorm;
  rhat[1] = r[1]/rnorm;
  rhat[2] = r[2]/rnorm;

  // energy

  double temp[3][3];
  MathExtra::plus3(wi.gamma,wj.gamma,temp);
  MathExtra::mldivide3(temp,rhat,s,error);
  sigma12 = 1.0/sqrt(0.5*MathExtra::dot3(s,rhat));
  MathExtra::times_column3(wi.A,rhat,z1);
  MathExtra::times_column3(wj.A,rhat,z2);
  v1[0] = z1[0]/shape2[type[i]][0];
  v1[1] = z1[1]/shape2[type[i]][1];
  v1[2] = z1[2]/shape2[type[i]][2];
  v2[0] = z2[0]/shape2[type[j]][0];
  v2[1] = z2[1]/shape2[type[j]][1];
  v2[2] = z2[2]/shape2[type[j]][2];
  sigma1 = 1.0/sqrt(MathExtra::dot3(z1,v1));
  sigma2 = 1.0/sqrt(MathExtra::dot3(z2,v2));
  H12[0][0] = wi.gamma[0][0]/sigma1+wj.gamma[0][0]/sigma2;
  H12[0][1] = wi.gamma[0][1]/sigma1+wj.gamma[0][1]/sigma2;
  H12[0][2] = wi.gamma[0][2]/sigma1+wj.gamma[0][2]/sigma2;
  H12[1][0] = wi.gamma[1][0]/sigma1+wj.gamma[1][0]/sigma2;
  H12[1][1] = wi.gamma[1][1]/sigma1+wj.gamma[1][1]/sigma2;
  H12[1][2] = wi.gamma[1][2]/sigma1+wj.gamma[1][2]/sigma2;
  H12[2][0] = wi.gamma[2][0]/sigma1+wj.gamma[2][0]/sigma2;
  H12[2][1] = wi.gamma[2][1]/sigma1+wj.gamma[2][1]/sigma2;
  H12[2][2] = wi.gamma[2][2]/sigma1+wj.gamma[2][2]/sigma2;
  dH=MathExtra::det3(H12);
  sigma1p2 = sigma1*sigma1;
  sigma2p2 = sigma2*sigma2;
  lambda = lshape[type[i]]/sigma1p2 + lshape[type[j]]/sigma2p2;
  nu = sqrt(dH/(sigma1+sigma2));
  MathExtra::times3(wi.aTe,wi.A,temp);
  double temp2[3][3];
  MathExtra::times3(wj.aTe,wj.A,temp2);
  MathExtra::plus3(temp,temp2,temp);
  MathExtra::mldivide3(temp,rhat,w,error);
  h12 = rnorm-sigma12;
  eta = lambda/nu;
  chi = 2.0*MathExtra::dot3(rhat,w);
  sprod = lshape[type[i]] * lshape[type[j]];
  sigh = sigma[type[i]][type[j]]/h12;
  tprod = eta*chi*sigh;

  double stemp = h12/2.0;
  Ua = (shape[type[i]][0]+stemp)*(shape[type[i]][1]+stemp)*
       (shape[type[i]][2]+stemp)*(shape[type[j]][0]+stemp)*
       (shape[type[j]][1]+stemp)*(shape[type[j]][2]+stemp);
  Ua = (1.0+3.0*tprod)*sprod/Ua;
  Ua = epsilon[type[i]][type[j]]*Ua/-36.0;

  stemp = h12/cr60;
  Ur = (shape[type[i]][0]+stemp)*(shape[type[i]][1]+stemp)*
       (shape[type[i]][2]+stemp)*(shape[type[j]][0]+stemp)*
       (shape[type[j]][1]+stemp)*(shape[type[j]][2]+stemp);
  Ur = (1.0+b_alpha*tprod)*sprod/Ur;
  Ur = epsilon[type[i]][type[j]]*Ur*pow(sigh,6.0)/2025.0;

  // force

  sec = sigma[type[i]][type[j]]*eta*chi;
  sigma12p3 = pow(sigma12,3.0);
  sigma1p3 = sigma1p2*sigma1;
  sigma2p3 = sigma2p2*sigma2;
  vsigma1[0] = -sigma1p3*v1[0];
  vsigma1[1] = -sigma1p3*v1[1];
  vsigma1[2] = -sigma1p3*v1[2];
  vsigma2[0] = -sigma2p3*v2[0];
  vsigma2[1] = -sigma2p3*v2[1];
  vsigma2[2] = -sigma2p3*v2[2];
  gsigma1[0][0] = -wi.gamma[0][0]/sigma1p2;
  gsigma1[0][1] = -wi.gamma[0][1]/sigma1p2;
  gsigma1[0][2] = -wi.gamma[0][2]/sigma1p2;
  gsigma1[1][0] = -wi.gamma[1][0]/sigma1p2;
  gsigma1[1][1] = -wi.gamma[1][1]/sigma1p2;
  gsigma1[1][2] = -wi.gamma[1][2]/sigma1p2;
  gsigma1[2][0] = -wi.gamma[2][0]/sigma1p2;
  gsigma1[2][1] = -wi.gamma[2][1]/sigma1p2;
  gsigma1[2][2] = -wi.gamma[2][2]/sigma1p2;
  gsigma2[0][0] = -wj.gamma[0][0]/sigma2p2;
  gsigma2[0][1] = -wj.gamma[0][1]/sigma2p2;
  gsigma2[0][2] = -wj.gamma[0][2]/sigma2p2;
  gsigma2[1][0] = -wj.gamma[1][0]/sigma2p2;
  gsigma2[1][1] = -wj.gamma[1][1]/sigma2p2;
  gsigma2[1][2] = -wj.gamma[1][2]/sigma2p2;
  gsigma2[2][0] = -wj.gamma[2][0]/sigma2p2;
  gsigma2[2][1] = -wj.gamma[2][1]/sigma2p2;
  gsigma2[2][2] = -wj.gamma[2][2]/sigma2p2;
  tsig1sig2 = eta/(2.0*(sigma1+sigma2));
  tdH = eta/(2.0*dH);
  teta1 = 2.0*eta/lambda;
  teta2 = teta1*lshape[type[j]]/sigma2p3;
  teta1 = teta1*lshape[type[i]]/sigma1p3;
  fourw[0] = 4.0*w[0];
  fourw[1] = 4.0*w[1];
  fourw[2] = 4.0*w[2];
  spr[0] = 0.5*sigma12p3*s[0];
  spr[1] = 0.5*sigma12p3*s[1];
  spr[2] = 0.5*sigma12p3*s[2];

  stemp = 1.0/(shape[type[i]][0]*2.0+h12)+
          1.0/(shape[type[i]][1]*2.0+h12)+
          1.0/(shape[type[i]][2]*2.0+h12)+
          1.0/(shape[type[j]][0]*2.0+h12)+
          1.0/(shape[type[j]][1]*2.0+h12)+
          1.0/(shape[type[j]][2]*2.0+h12);
  hsec = h12+3.0*sec;
  dspu = 1.0/h12-1.0/hsec+stemp;
  pbsu = 3.0*sigma[type[i]][type[j]]/hsec;
  
  stemp = 1.0/(shape[type[i]][0]*cr60+h12)+
          1.0/(shape[type[i]][1]*cr60+h12)+
          1.0/(shape[type[i]][2]*cr60+h12)+
          1.0/(shape[type[j]][0]*cr60+h12)+
          1.0/(shape[type[j]][1]*cr60+h12)+
          1.0/(shape[type[j]][2]*cr60+h12);
  hsec = h12+b_alpha*sec;
  dspr = 7.0/h12-1.0/hsec+stemp;
  pbsr = b_alpha*sigma[type[i]][type[j]]/hsec;
  
  for (int i=0; i<3; i++) {
    u[0] = -rhat[i]*rhat[0];
    u[1] = -rhat[i]*rhat[1];
    u[2] = -rhat[i]*rhat[2];
    u[i] += 1.0;
    u[0] /= rnorm;
    u[1] /= rnorm;
    u[2] /= rnorm;
    MathExtra::times_column3(wi.A,u,u1);
    MathExtra::times_column3(wj.A,u,u2);
    dsigma1=MathExtra::dot3(u1,vsigma1);
    dsigma2=MathExtra::dot3(u2,vsigma2);
    dH12[0][0] = dsigma1*gsigma1[0][0]+dsigma2*gsigma2[0][0];
    dH12[0][1] = dsigma1*gsigma1[0][1]+dsigma2*gsigma2[0][1];
    dH12[0][2] = dsigma1*gsigma1[0][2]+dsigma2*gsigma2[0][2];
    dH12[1][0] = dsigma1*gsigma1[1][0]+dsigma2*gsigma2[1][0];
    dH12[1][1] = dsigma1*gsigma1[1][1]+dsigma2*gsigma2[1][1];
    dH12[1][2] = dsigma1*gsigma1[1][2]+dsigma2*gsigma2[1][2];
    dH12[2][0] = dsigma1*gsigma1[2][0]+dsigma2*gsigma2[2][0];
    dH12[2][1] = dsigma1*gsigma1[2][1]+dsigma2*gsigma2[2][1];
    dH12[2][2] = dsigma1*gsigma1[2][2]+dsigma2*gsigma2[2][2];
    ddH = det_prime(H12,dH12);
    deta = (dsigma1+dsigma2)*tsig1sig2;
    deta -= ddH*tdH;
    deta -= dsigma1*teta1+dsigma2*teta2;
    dchi = MathExtra::dot3(u,fourw);
    dh12 = rhat[i]+MathExtra::dot3(u,spr);
    dUa = pbsu*(eta*dchi+deta*chi)-dh12*dspu;
    dUr = pbsr*(eta*dchi+deta*chi)-dh12*dspr;
    fforce[i]=dUr*Ur+dUa*Ua;
  }
    
  // torque on i

  MathExtra::row_times3(fourw,wi.aTe,fwae);

  for (int i=0; i<3; i++) {
    MathExtra::times_column3(wi.lA[i],rhat,p);
    dsigma1 = MathExtra::dot3(p,vsigma1);
    dH12[0][0] = wi.lAsa[i][0][0]/sigma1+dsigma1*gsigma1[0][0];
    dH12[0][1] = wi.lAsa[i][0][1]/sigma1+dsigma1*gsigma1[0][1];
    dH12[0][2] = wi.lAsa[i][0][2]/sigma1+dsigma1*gsigma1[0][2];
    dH12[1][0] = wi.lAsa[i][1][0]/sigma1+dsigma1*gsigma1[1][0];
    dH12[1][1] = wi.lAsa[i][1][1]/sigma1+dsigma1*gsigma1[1][1];
    dH12[1][2] = wi.lAsa[i][1][2]/sigma1+dsigma1*gsigma1[1][2];
    dH12[2][0] = wi.lAsa[i][2][0]/sigma1+dsigma1*gsigma1[2][0];
    dH12[2][1] = wi.lAsa[i][2][1]/sigma1+dsigma1*gsigma1[2][1];
    dH12[2][2] = wi.lAsa[i][2][2]/sigma1+dsigma1*gsigma1[2][2];
    ddH = det_prime(H12,dH12);
    deta = tsig1sig2*dsigma1-tdH*ddH;
    deta -= teta1*dsigma1;
    double tempv[3];
    MathExtra::times_column3(wi.lA[i],w,tempv);
    dchi = -MathExtra::dot3(fwae,tempv);
    MathExtra::times_column3(wi.lAtwo[i],spr,tempv);
    dh12 = -MathExtra::dot3(s,tempv);

    dUa = pbsu*(eta*dchi + deta*chi)-dh12*dspu;
    dUr = pbsr*(eta*dchi + deta*chi)-dh12*dspr;
    ttor[i] = -(dUa*Ua+dUr*Ur);
  }
  
  // torque on j

  if (!(force->newton_pair || j < atom->nlocal))
    return Ua+Ur;

  MathExtra::row_times3(fourw,wj.aTe,fwae);

  for (int i=0; i<3; i++) {
    MathExtra::times_column3(wj.lA[i],rhat,p);
    dsigma2 = MathExtra::dot3(p,vsigma2);
    dH12[0][0] = wj.lAsa[i][0][0]/sigma2+dsigma2*gsigma2[0][0];
    dH12[0][1] = wj.lAsa[i][0][1]/sigma2+dsigma2*gsigma2[0][1];
    dH12[0][2] = wj.lAsa[i][0][2]/sigma2+dsigma2*gsigma2[0][2];
    dH12[1][0] = wj.lAsa[i][1][0]/sigma2+dsigma2*gsigma2[1][0];
    dH12[1][1] = wj.lAsa[i][1][1]/sigma2+dsigma2*gsigma2[1][1];
    dH12[1][2] = wj.lAsa[i][1][2]/sigma2+dsigma2*gsigma2[1][2];
    dH12[2][0] = wj.lAsa[i][2][0]/sigma2+dsigma2*gsigma2[2][0];
    dH12[2][1] = wj.lAsa[i][2][1]/sigma2+dsigma2*gsigma2[2][1];
    dH12[2][2] = wj.lAsa[i][2][2]/sigma2+dsigma2*gsigma2[2][2];
    ddH = det_prime(H12,dH12);
    deta = tsig1sig2*dsigma2-tdH*ddH;
    deta -= teta2*dsigma2;
    double tempv[3];
    MathExtra::times_column3(wj.lA[i],w,tempv);
    dchi = -MathExtra::dot3(fwae,tempv);
    MathExtra::times_column3(wj.lAtwo[i],spr,tempv);
    dh12 = -MathExtra::dot3(s,tempv);

    dUa = pbsu*(eta*dchi + deta*chi)-dh12*dspu;
    dUr = pbsr*(eta*dchi + deta*chi)-dh12*dspr;
    rtor[i] = -(dUa*Ua+dUr*Ur);
  }

  return Ua+Ur;
}

/* ----------------------------------------------------------------------
   Compute the energy, force, torque for a pair (INTEGRATED-LJ)
------------------------------------------------------------------------- */

double PairRESquared::resquared_lj(const int i, const int j,
                                   const RE2Vars &wi, const double *r, 
                                   const double rsq, double *fforce, 
                                   double *ttor, bool calc_torque)
{
  int *type = atom->type;
  double **shape = atom->shape;
    
  // pair computations for energy, force, torque

  double rnorm;              // L2 norm of r
  double rhat[3];            // r/rnorm
  double s[3];               // inv(gamma1)*rhat
  double sigma12;            // 1/sqrt(0.5*s'*rhat)
  double w[3];               // inv(A1'*E1*A1+I)*rhat
  double h12;                // rnorm-sigma12;
  double chi;                // 2*rhat'*w
  double sigh;               // sigma/h12
  double tprod;              // chi*sigh
  double Ua,Ur;              // attractive/repulsive parts of potential
    
  // pair computations for force, torque

  double sec;                          // sigma*chi
  double sigma12p3;                    // sigma12^3
  double fourw[3];                     // 4*w;
  double spr[3];                       // 0.5*sigma12^3*s
  double hsec;                         // h12+[3,b_alpha]*sec
  double dspu;                         // 1/h12 - 1/hsec + temp
  double pbsu;                         // 3*sigma/hsec
  double dspr;                         // 7/h12-1/hsec+temp
  double pbsr;                         // b_alpha*sigma/hsec;
  double u[3];                         // (-rhat(i)*rhat+eye(:,i))/rnorm
  double dchi,dh12;                    // derivatives of chi,h12
  double dUr,dUa;                      // derivatives of Ua,Ur
  double h12p3;                        // h12^3
      
  // pair computations for torque

  double fwae[3];        // -fourw'*aTe
  double p[3];           // lA*rhat

  // distance of closest approach correction

  double aTs[3][3];       // A1'*S1^2
  double gamma[3][3];     // A1'*S1^2*A
  double lAtwo[3][3][3];  // A1'*S1^2*wi.lA
  double scorrect[3];
  double half_sigma=sigma[type[i]][type[j]] / 2.0;
  scorrect[0] = shape[type[i]][0]+half_sigma;
  scorrect[1] = shape[type[i]][1]+half_sigma;
  scorrect[2] = shape[type[i]][2]+half_sigma;
  scorrect[0] = scorrect[0] * scorrect[0] / 2.0;
  scorrect[1] = scorrect[1] * scorrect[1] / 2.0;
  scorrect[2] = scorrect[2] * scorrect[2] / 2.0;
  MathExtra::transpose_times_diag3(wi.A,scorrect,aTs);
  MathExtra::times3(aTs,wi.A,gamma);
  for (int ii=0; ii<3; ii++)
    MathExtra::times3(aTs,wi.lA[ii],lAtwo[ii]);

  rnorm=sqrt(rsq);
  rhat[0] = r[0]/rnorm;
  rhat[1] = r[1]/rnorm;
  rhat[2] = r[2]/rnorm;

  // energy

  MathExtra::mldivide3(gamma,rhat,s,error);
  sigma12 = 1.0/sqrt(0.5*MathExtra::dot3(s,rhat));
  double temp[3][3];
  MathExtra::times3(wi.aTe,wi.A,temp);
  temp[0][0] += 1.0;
  temp[1][1] += 1.0;
  temp[2][2] += 1.0;
  MathExtra::mldivide3(temp,rhat,w,error);
  h12 = rnorm-sigma12;
  chi = 2.0*MathExtra::dot3(rhat,w);
  sigh = sigma[type[i]][type[j]]/h12;
  tprod = chi*sigh;

  h12p3 = pow(h12,3.0);
  double sigmap3 = pow(sigma[type[i]][type[j]],3.0);
  double stemp = h12/2.0;
  Ua = (shape[type[i]][0]+stemp)*(shape[type[i]][1]+stemp)*
       (shape[type[i]][2]+stemp)*h12p3/8.0;
  Ua = (1.0+3.0*tprod)*lshape[type[i]]/Ua;
  Ua = epsilon[type[i]][type[j]]*Ua*sigmap3*solv_f_a;

  stemp = h12/cr60;
  Ur = (shape[type[i]][0]+stemp)*(shape[type[i]][1]+stemp)*
       (shape[type[i]][2]+stemp)*h12p3/60.0;
  Ur = (1.0+b_alpha*tprod)*lshape[type[i]]/Ur;
  Ur = epsilon[type[i]][type[j]]*Ur*sigmap3*pow(sigh,6.0)*solv_f_r;

  // force

  sec = sigma[type[i]][type[j]]*chi;
  sigma12p3 = pow(sigma12,3.0);
  fourw[0] = 4.0*w[0];
  fourw[1] = 4.0*w[1];
  fourw[2] = 4.0*w[2];
  spr[0] = 0.5*sigma12p3*s[0];
  spr[1] = 0.5*sigma12p3*s[1];
  spr[2] = 0.5*sigma12p3*s[2];

  stemp = 1.0/(shape[type[i]][0]*2.0+h12)+
          1.0/(shape[type[i]][1]*2.0+h12)+
          1.0/(shape[type[i]][2]*2.0+h12)+
          3.0/h12;
  hsec = h12+3.0*sec;
  dspu = 1.0/h12-1.0/hsec+stemp;
  pbsu = 3.0*sigma[type[i]][type[j]]/hsec;
  
  stemp = 1.0/(shape[type[i]][0]*cr60+h12)+
          1.0/(shape[type[i]][1]*cr60+h12)+
          1.0/(shape[type[i]][2]*cr60+h12)+
          3.0/h12;
  hsec = h12+b_alpha*sec;
  dspr = 7.0/h12-1.0/hsec+stemp;
  pbsr = b_alpha*sigma[type[i]][type[j]]/hsec;
  
  for (int i=0; i<3; i++) {
    u[0] = -rhat[i]*rhat[0];
    u[1] = -rhat[i]*rhat[1];
    u[2] = -rhat[i]*rhat[2];
    u[i] += 1.0;
    u[0] /= rnorm;
    u[1] /= rnorm;
    u[2] /= rnorm;
    dchi = MathExtra::dot3(u,fourw);
    dh12 = rhat[i]+MathExtra::dot3(u,spr);
    dUa = pbsu*dchi-dh12*dspu;
    dUr = pbsr*dchi-dh12*dspr;
    fforce[i]=dUr*Ur+dUa*Ua;
  }
    
  // torque on i

  if (calc_torque) {
    MathExtra::row_times3(fourw,wi.aTe,fwae);

    for (int i=0; i<3; i++) {
      MathExtra::times_column3(wi.lA[i],rhat,p);
      double tempv[3];
      MathExtra::times_column3(wi.lA[i],w,tempv);
      dchi = -MathExtra::dot3(fwae,tempv);
      MathExtra::times_column3(lAtwo[i],spr,tempv);
      dh12 = -MathExtra::dot3(s,tempv);

      dUa = pbsu*dchi-dh12*dspu;
      dUr = pbsr*dchi-dh12*dspr;
      ttor[i] = -(dUa*Ua+dUr*Ur);
    }
  }

  return Ua+Ur;
}
