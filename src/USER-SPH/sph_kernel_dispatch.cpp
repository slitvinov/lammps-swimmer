#include "string.h"
#include "error.h"
#include "sph_kernel_lucy_2d.h"

namespace LAMMPS_NS {
  SPHKernel* sph_kernel_dispatch(char* kernel_name, int kernel_dimension,
				 Error *error) {
    if ( (strcmp(kernel_name, "lucy")==0) && (kernel_dimension==2) ) {
      return new SPHKernelLucy2D();    
    }
    else {
      char str[128];
      sprintf(str, "Unknown kernel type and/or dimension: %s and %i", kernel_name , kernel_dimension);
      error->all(FLERR, str);
      return NULL;
    }
  }
}
