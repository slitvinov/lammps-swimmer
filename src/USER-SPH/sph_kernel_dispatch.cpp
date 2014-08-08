#include "string.h"
#include "error.h"
#include "sph_kernel_dispatch.h"
#include "sph_kernel_lucy_2d.h"
#include "sph_kernel_lucy_3d.h"

namespace LAMMPS_NS {
  SPHKernelCodeType sph_kernel_code(char* kernel_name, int kernel_dimension,
				    Error *error) {
    if ( (strcmp(kernel_name, "lucy")==0) && (kernel_dimension==2) ) {
      return Lucy2D;
    }
    else if ( (strcmp(kernel_name, "lucy")==0) && (kernel_dimension==3) ) {
      return Lucy3D;
    }
    else {
      char str[128];
      sprintf(str, "Unknown kernel type and/or dimension: %s and %i", kernel_name , kernel_dimension);
      error->all(FLERR, str);
    }
  };

  SPHKernel* sph_kernel_decode(SPHKernelCodeType kernel_code,
			       Error *error) {
    if (kernel_code == Lucy2D) {
      return new SPHKernelLucy2D();
    }
    else if (kernel_code == Lucy3D) {
      return new SPHKernelLucy3D();
    }
    else {
      error->all(FLERR, "Unknown kernel code");
    }
  }
}
