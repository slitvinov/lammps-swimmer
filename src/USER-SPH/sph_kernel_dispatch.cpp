#include "string.h"
#include "error.h"
#include "sph_kernel_dispatch.h"
#include "sph_kernel_lucy_2d.h"
#include "sph_kernel_lucy_3d.h"
#include "sph_kernel_quadric_2d.h"
#include "sph_kernel_quadric_3d.h"
#include "sph_kernel_quintic_2d.h"
#include "sph_kernel_quintic_3d.h"
#include "sph_kernel_wendland4_2d.h"

namespace LAMMPS_NS {
  SPHKernel* sph_kernel_dispatch(char* kernel_name, int kernel_dimension,
				 Error *error) {
    if ( (strcmp(kernel_name, "lucy")==0) && (kernel_dimension==2) ) {
      return new SPHKernelLucy2D();
    }
    else if ( (strcmp(kernel_name, "lucy")==0) && (kernel_dimension==3) ) {
      return new SPHKernelLucy3D();
    }
    else if ( (strcmp(kernel_name, "quadric")==0) && (kernel_dimension==2) ) {
      return new SPHKernelQuadric2D();
    }
    else if ( (strcmp(kernel_name, "quadric")==0) && (kernel_dimension==3) ) {
      return new SPHKernelQuadric3D();
    }
    else if ( (strcmp(kernel_name, "quintic")==0) && (kernel_dimension==2) ) {
      return new SPHKernelQuintic2D();
    }
    else if ( (strcmp(kernel_name, "wendland4")==0) && (kernel_dimension==2) ) {
      return new SPHKernelWendland42D();
    }
    else if ( (strcmp(kernel_name, "quintic")==0) && (kernel_dimension==3) ) {
      return new SPHKernelQuintic3D();
    }
    //    else if ( (strcmp(kernel_name, "wendland4")==0) && (kernel_dimension==3) ) {
    //      return new SPHKernelQuintic3D();
    //    }
    else {
      char str[128];
      sprintf(str, "Unknown kernel type and/or dimension: %s and %i", kernel_name , kernel_dimension);
      error->all(FLERR, str);
    }
  }

}
