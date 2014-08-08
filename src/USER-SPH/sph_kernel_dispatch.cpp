#include "string.h"
#include "sph_kernel_lucy_2d.h"

namespace LAMMPS_NS {
  SPHKernel* sph_kernel_dispatch(char* kernel_name, int kernel_dimension) {
    if ( (strcpy(kernel_name, "lucy")==0) && (kernel_dimension==2) ) {
      return new SPHKernelLucy2D();    
    }
    else {
      return NULL;
    }
  }
}
