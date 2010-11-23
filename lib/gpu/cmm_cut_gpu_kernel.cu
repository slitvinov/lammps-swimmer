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
   Contributing authors: Mike Brown (ORNL), brownw@ornl.gov
------------------------------------------------------------------------- */

#ifndef CMM_GPU_KERNEL
#define CMM_GPU_KERNEL

#define MAX_SHARED_TYPES 8

#ifdef _DOUBLE_DOUBLE
#define numtyp double
#define numtyp2 double2
#define numtyp4 double4
#define acctyp double
#define acctyp4 double4
#endif

#ifdef _SINGLE_DOUBLE
#define numtyp float
#define numtyp2 float2
#define numtyp4 float4
#define acctyp double
#define acctyp4 double4
#endif

#ifndef numtyp
#define numtyp float
#define numtyp2 float2
#define numtyp4 float4
#define acctyp float
#define acctyp4 float4
#endif

#ifdef NV_KERNEL

#include "geryon/ucl_nv_kernel.h"
texture<float4> pos_tex;

#ifdef _DOUBLE_DOUBLE
__inline double4 fetch_pos(const int& i, const double4 *pos)
{
  return pos[i];
}
#else
__inline float4 fetch_pos(const int& i, const float4 *pos)
{
  return tex1Dfetch(pos_tex, i);
}
#endif

#else

#pragma OPENCL EXTENSION cl_khr_fp64: enable
#define GLOBAL_ID_X get_global_id(0)
#define THREAD_ID_X get_local_id(0)
#define BLOCK_ID_X get_group_id(0)
#define BLOCK_SIZE_X get_local_size(0)
#define __syncthreads() barrier(CLK_LOCAL_MEM_FENCE)
#define __inline inline

#define fetch_pos(i,y) x_[i]

#endif

__kernel void kernel_pair(__global numtyp4 *x_, __global numtyp4 *lj1,
                          __global numtyp4* lj3, const int lj_types, 
                          __global numtyp *sp_lj_in, __global int *dev_nbor, 
                          __global acctyp4 *ans, __global acctyp *engv, 
                          const int eflag, const int vflag, const int inum, 
                          const int nall, const int nbor_pitch) {
  // ii indexes the two interacting particles in gi
  int ii=GLOBAL_ID_X;
  __local numtyp sp_lj[4];
  sp_lj[0]=sp_lj_in[0];
  sp_lj[1]=sp_lj_in[1];
  sp_lj[2]=sp_lj_in[2];
  sp_lj[3]=sp_lj_in[3];

  if (ii<inum) {
  
    acctyp energy=(numtyp)0;
    acctyp4 f;
    f.x=(numtyp)0;
    f.y=(numtyp)0;
    f.z=(numtyp)0;
    acctyp virial[6];
    for (int i=0; i<6; i++)
      virial[i]=(numtyp)0;
  
    __global int *nbor=dev_nbor+ii;
    int i=*nbor;
    nbor+=nbor_pitch;
    int numj=*nbor;
    nbor+=nbor_pitch;
    __global int *list_end=nbor+mul24(numj,nbor_pitch);
  
    numtyp4 ix=fetch_pos(i,x_); //x_[i];
    int itype=ix.w;

    numtyp factor_lj;
    for ( ; nbor<list_end; nbor+=nbor_pitch) {
  
      int j=*nbor;
      if (j < nall) 
        factor_lj = (numtyp)1.0;
      else {
        factor_lj = sp_lj[j/nall];
        j %= nall;
      }
      numtyp4 jx=fetch_pos(j,x_); //x_[j];
      int jtype=jx.w;

      // Compute r12
      numtyp delx = ix.x-jx.x;
      numtyp dely = ix.y-jx.y;
      numtyp delz = ix.z-jx.z;
      numtyp r2inv = delx*delx+dely*dely+delz*delz;
        
      int mtype=itype*lj_types+jtype;
      if (r2inv<lj1[mtype].x) {
        r2inv=(numtyp)1.0/r2inv;
        numtyp inv1,inv2;
        
        if (lj1[mtype].y == 2) {
          inv1=r2inv*r2inv;
          inv2=inv1*inv1;
        } else if (lj1[mtype].y == 1) {
          inv2=r2inv*sqrt(r2inv);
          inv1=inv2*inv2;
        } else {
          inv1=r2inv*r2inv*r2inv;
          inv2=inv1;
        }
        numtyp force = factor_lj*r2inv*inv1*(lj1[mtype].z*inv2-lj1[mtype].w);
      
        f.x+=delx*force;
        f.y+=dely*force;
        f.z+=delz*force;
        if (eflag>0)
          energy += factor_lj*inv1*(lj3[mtype].x*inv2-lj3[mtype].y)-
                    lj3[mtype].z;
        if (vflag>0) {
          virial[0] += delx*delx*force;
          virial[1] += dely*dely*force;
          virial[2] += delz*delz*force;
          virial[3] += delx*dely*force;
          virial[4] += delx*delz*force;
          virial[5] += dely*delz*force;
        }
      }

    } // for nbor

    // Store answers
    __global acctyp *ap1=engv+ii;
    if (eflag>0) {
      *ap1=energy;
      ap1+=inum;
    }
    if (vflag>0) {
      for (int i=0; i<6; i++) {
        *ap1=virial[i];
        ap1+=inum;
      }
    }
    ans[ii]=f;
  } // if ii
}

__kernel void kernel_pair_fast(__global numtyp4 *x_, __global numtyp4 *lj1_in,
                               __global numtyp4* lj3_in, 
                               __global numtyp* sp_lj_in,__global int *dev_nbor, 
                               __global acctyp4 *ans, __global acctyp *engv, 
                               const int eflag, const int vflag, const int inum, 
                               const int nall, const int nbor_pitch) {
  // ii indexes the two interacting particles in gi
  int ii=THREAD_ID_X;
  __local numtyp4 lj1[MAX_SHARED_TYPES*MAX_SHARED_TYPES];
  __local numtyp4 lj3[MAX_SHARED_TYPES*MAX_SHARED_TYPES];
  __local numtyp sp_lj[4];
  if (ii<4)
    sp_lj[ii]=sp_lj_in[ii];
  if (ii<MAX_SHARED_TYPES*MAX_SHARED_TYPES) {
    lj1[ii]=lj1_in[ii];
    if (eflag>0)
      lj3[ii]=lj3_in[ii];
  }
  ii+=mul24((int)BLOCK_ID_X,(int)BLOCK_SIZE_X);
  __syncthreads();
  
  if (ii<inum) {
  
    acctyp energy=(numtyp)0;
    acctyp4 f;
    f.x=(numtyp)0;
    f.y=(numtyp)0;
    f.z=(numtyp)0;
    acctyp virial[6];
    for (int i=0; i<6; i++)
      virial[i]=(numtyp)0;
  
    __global int *nbor=dev_nbor+ii;
    int i=*nbor;
    nbor+=nbor_pitch;
    int numj=*nbor;
    nbor+=nbor_pitch;
    __global int *list_end=nbor+mul24(numj,nbor_pitch);
  
    numtyp4 ix=fetch_pos(i,x_); //x_[i];
    int iw=ix.w;
    int itype=mul24((int)MAX_SHARED_TYPES,iw);

    numtyp factor_lj;
    for ( ; nbor<list_end; nbor+=nbor_pitch) {
  
      int j=*nbor;
      if (j < nall) 
        factor_lj = (numtyp)1.0;
      else {
        factor_lj = sp_lj[j/nall];
        j %= nall;
      }
      numtyp4 jx=fetch_pos(j,x_); //x_[j];
      int mtype=itype+jx.w;

      // Compute r12
      numtyp delx = ix.x-jx.x;
      numtyp dely = ix.y-jx.y;
      numtyp delz = ix.z-jx.z;
      numtyp r2inv = delx*delx+dely*dely+delz*delz;
        
      if (r2inv<lj1[mtype].x) {
        r2inv=(numtyp)1.0/r2inv;
        numtyp inv1,inv2;
        
        if (lj1[mtype].y == (numtyp)2) {
          inv1=r2inv*r2inv;
          inv2=inv1*inv1;
        } else if (lj1[mtype].y == (numtyp)1) {
          inv2=r2inv*sqrt(r2inv);
          inv1=inv2*inv2;
        } else {
          inv1=r2inv*r2inv*r2inv;
          inv2=inv1;
        }
        numtyp force = factor_lj*r2inv*inv1*(lj1[mtype].z*inv2-lj1[mtype].w);
      
        f.x+=delx*force;
        f.y+=dely*force;
        f.z+=delz*force;
        if (eflag>0)
          energy += factor_lj*inv1*(lj3[mtype].x*inv2-lj3[mtype].y)-
                    lj3[mtype].z;
        if (vflag>0) {
          virial[0] += delx*delx*force;
          virial[1] += dely*dely*force;
          virial[2] += delz*delz*force;
          virial[3] += delx*dely*force;
          virial[4] += delx*delz*force;
          virial[5] += dely*delz*force;
        }
      }

    } // for nbor

    // Store answers
    __global acctyp *ap1=engv+ii;
    if (eflag>0) {
      *ap1=energy;
      ap1+=inum;
    }
    if (vflag>0) {
      for (int i=0; i<6; i++) {
        *ap1=virial[i];
        ap1+=inum;
      }
    }
    ans[ii]=f;
  } // if ii*/
}

#endif

