#include "string.h"
#include "error.h"
#include "sph_bn_utils.h"
#include <algorithm>    // std::max

#define MAXLINE 1024

namespace LAMMPS_NS {
  /* ----------------------------------------------------------------------
     read target field from a user-specified file
     only called by proc 0
     ------------------------------------------------------------------------- */
  
  void read_initial_target_field(FILE* fpr,
				 int nxnodes, int nynodes, int nznodes,
				 int ***T_initial_set, double ***T_target,
				 Error *error) {
    char line[MAXLINE];
    
    for (int ixnode = 0; ixnode < nxnodes; ixnode++)
      for (int iynode = 0; iynode < nynodes; iynode++)
	for (int iznode = 0; iznode < nznodes; iznode++)
	  T_initial_set[ixnode][iynode][iznode] = 0;
    
    // read target values from file
    
    int ixnode,iynode,iznode;
    double T_tmp;
    while (1) {
      if (fgets(line,MAXLINE,fpr) == NULL) break;
      sscanf(line,"%d %d %d %lg",&ixnode,&iynode,&iznode,&T_tmp);
      if (T_tmp < 0.0) error->one(FLERR,"pair/sph/bn target filed must be > 0.0");
      T_target[ixnode][iynode][iznode] = T_tmp;
      T_initial_set[ixnode][iynode][iznode] = 1;
    }
    
    for (int ixnode = 0; ixnode < nxnodes; ixnode++)
      for (int iynode = 0; iynode < nynodes; iynode++)
	for (int iznode = 0; iznode < nznodes; iznode++)
	  if (T_initial_set[ixnode][iynode][iznode] == 0)
	    error->one(FLERR,"Initial filed not all set in pair/sph/bn");
    
    // close file
    fclose(fpr);
  };

  double get_target_field (double* xi, Domain *&domain, double ***T_target,
			   int nxnodes, int nynodes, int nznodes, 
			   int ntime_smooth, bigint ntimestep) {
    double xscale = (xi[0] - domain->boxlo[0])/domain->xprd;
    double yscale = (xi[1] - domain->boxlo[1])/domain->yprd;
    double zscale = (xi[2] - domain->boxlo[2])/domain->zprd;
    int ixnode = static_cast<int>(xscale*nxnodes);
    int iynode = static_cast<int>(yscale*nynodes);
    int iznode = static_cast<int>(zscale*nznodes);
    while (ixnode > nxnodes-1) ixnode -= nxnodes;
    while (iynode > nynodes-1) iynode -= nynodes;
    while (iznode > nznodes-1) iznode -= nznodes;
    while (ixnode < 0) ixnode += nxnodes;
    while (iynode < 0) iynode += nynodes;
    while (iznode < 0) iznode += nznodes;

    double Tfin = T_target[ixnode][iynode][iznode];
    if (ntimestep>ntime_smooth) {
      return Tfin;
    }
    double Tr  = (1.0*(ntime_smooth-ntimestep) + Tfin*ntimestep)
      /static_cast<double>(ntime_smooth);
    return Tr;
  }

  double get_target_cutoff (double m, int nn, double rhot) {
    double rhot_cut = std::max(rhot, 1e-9);
    return sqrt(m*nn/rhot_cut)/sqrt(3.141592653589793);
  }
  
}
