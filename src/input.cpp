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

#include "mpi.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "input.h"
#include "system.h"
#include "universe.h"
#include "atom.h"
#include "comm.h"
#include "group.h"
#include "domain.h"
#include "output.h"
#include "thermo.h"
#include "force.h"
#include "pair.h"
#include "min.h"
#include "modify.h"
#include "bond.h"
#include "angle.h"
#include "dihedral.h"
#include "improper.h"
#include "kspace.h"
#include "update.h"
#include "neighbor.h"
#include "special.h"
#include "variable.h"
#include "error.h"
#include "memory.h"

#define AtomInclude
#define CommandInclude
#include "style.h"
#undef AtomInclude
#undef CommandInclude

#define MAXLINE 1024
#define DELTA 4

/* ---------------------------------------------------------------------- */

Input::Input(int argc, char **argv)
{
  MPI_Comm_rank(world,&me);

  line = new char[MAXLINE];
  copy = new char[MAXLINE];
  work = new char[MAXLINE];
  narg = maxarg = 0;
  arg = NULL;

  echo_screen = 0;
  echo_log = 1;

  label_active = 0;
  labelstr = NULL;
  jump_skip = 0;

  if (me == 0) {
    nfile = maxfile = 1;
    infiles = (FILE **) memory->smalloc(sizeof(FILE *),"input:infiles");
    infiles[0] = infile;
  } else infiles = NULL;

  variable = new Variable;

  // process command-line args
  // check for args "-var" and "-echo"
  // caller has already checked that sufficient arguments exist

  int iarg = 0;
  while (iarg < argc) {
    if (strcmp(argv[iarg],"-var") == 0) {
      variable->set(argv[iarg+1],argv[iarg+2]);
      iarg += 2;
    } else if (strcmp(argv[iarg],"-echo") == 0) {
      narg = 1;
      char **tmp = arg;        // trick echo() into using argv instead of arg
      arg = &argv[iarg+1];
      echo();
      arg = tmp;
      iarg += 3;
    } else iarg++;
  }
}

/* ---------------------------------------------------------------------- */

Input::~Input()
{
  // don't free command and arg strings
  // they just point to other allocated memory

  delete variable;
  delete [] line;
  delete [] copy;
  delete [] work;
  if (labelstr) delete [] labelstr;
  if (arg) memory->sfree(arg);
  if (infiles) memory->sfree(infiles);
}

/* ----------------------------------------------------------------------
   process all input from infile
   infile = stdin or file if command-line arg "-in" was used
------------------------------------------------------------------------- */

void Input::file()
{
  int n,flag;

  while (1) {
    
    // read one line from input script
    // if line ends in continuation char '&', concatenate next line(s)
    // n = str length of line
    
    if (me == 0) {
      if (fgets(line,MAXLINE,infile) == NULL) n = 0;
      else n = strlen(line) + 1;
      while (n >= 3 && line[n-3] == '&') {
	if (fgets(&line[n-3],MAXLINE-n+3,infile) == NULL) n = 0;
	else n = strlen(line) + 1;
      }
    }

    // bcast the line
    // if n = 0, end-of-file
    // error if label_active is set, since label wasn't encountered
    // if original input file, code is done
    // else go back to previous input file

    MPI_Bcast(&n,1,MPI_INT,0,world);
    if (n == 0) {
      if (label_active) error->all("Label wasn't found in input script");
      if (me == 0) {
	if (infile != stdin) fclose(infile);
	nfile--;
      }
      MPI_Bcast(&nfile,1,MPI_INT,0,world);
      if (nfile == 0) break;
      if (me == 0) infile = infiles[nfile-1];
      continue;
    }

    MPI_Bcast(line,n,MPI_CHAR,0,world);

    // if n = MAXLINE, line is too long

    if (n == MAXLINE) {
      char str[MAXLINE+32];
      sprintf(str,"Input line too long: %s",line);
      error->all(str);
    }

    // echo the command unless scanning for label

    if (me == 0 && label_active == 0) {
      if (echo_screen && screen) fprintf(screen,"%s",line); 
      if (echo_log && logfile) fprintf(logfile,"%s",line);
    }

    // parse the line
    // if no command, skip to next line in input script

    parse();
    if (command == NULL) continue;

    // if scanning for label, skip command unless it's a label command

    if (label_active && strcmp(command,"label") != 0) continue;

    // execute the command

    if (execute_command()) {
      char str[MAXLINE];
      sprintf(str,"Unknown command: %s",line);
      error->all(str);
    }
  }
}

/* ----------------------------------------------------------------------
   process all input from filename
------------------------------------------------------------------------- */

void Input::file(char *filename)
{
  // error if another nested file still open
  // if single open file is not stdin, close it
  // open new filename and set infile, infiles[0]

  if (me == 0) {
    if (nfile > 1)
      error->one("Another input script is already being processed");
    if (infile != stdin) fclose(infile);
    infile = fopen(filename,"r");
    if (infile == NULL) {
      char str[128];
      sprintf(str,"Cannot open input script %s",filename);
      error->one(str);
    }
    infiles[0] = infile;
  } else infile = NULL;

  file();
}

/* ----------------------------------------------------------------------
   parse the command in single and execute it
   return command name to caller
------------------------------------------------------------------------- */

char *Input::one(char *single)
{
  strcpy(line,single);

  // echo the command unless scanning for label
  
  if (me == 0 && label_active == 0) {
    if (echo_screen && screen) fprintf(screen,"%s",line); 
    if (echo_log && logfile) fprintf(logfile,"%s",line);
  }

  // parse the line
  // if no command, just return NULL

  parse();
  if (command == NULL) return NULL;

  // if scanning for label, skip command unless it's a label command

  if (label_active && strcmp(command,"label") != 0) return NULL;

  // execute the command and return its name

  if (execute_command()) {
    char str[MAXLINE];
    sprintf(str,"Unknown command: %s",line);
    error->all(str);
  }

  return command;
}

/* ----------------------------------------------------------------------
   parse copy of command line
   strip comment = all chars from # on
   replace all $ via variable substitution
   command = first word
   narg = # of args
   arg[] = individual args
   treat multiple args between double quotes as one arg
------------------------------------------------------------------------- */

void Input::parse()
{
  // make a copy to work on

  strcpy(copy,line);

  // strip any # comment by resetting string terminator
  // do not strip # inside double quotes

  int level = 0;
  char *ptr = copy;
  while (*ptr) {
    if (*ptr == '#' && level == 0) {
      *ptr = '\0';
      break;
    }
    if (*ptr == '"') {
      if (level == 0) level = 1;
      else level = 0;
    }
    ptr++;
  }

  // perform $ variable substitution (print changes)

  substitute(copy,1);

  // command = 1st arg

  command = strtok(copy," \t\n\r\f");
  if (command == NULL) return;

  // point arg[] at each subsequent arg
  // treat text between double quotes as one arg
  // insert string terminators in copy to delimit args

  narg = 0;
  while (1) {
    if (narg == maxarg) {
      maxarg += DELTA;
      arg = (char **) memory->srealloc(arg,maxarg*sizeof(char *),"input:arg");
    }
    arg[narg] = strtok(NULL," \t\n\r\f");
    if (arg[narg] && arg[narg][0] == '\"') {
      arg[narg] = &arg[narg][1];
      if (strchr(arg[narg],'\"')) error->all("Quotes in a single arg");
      arg[narg][strlen(arg[narg])] = ' ';
      ptr = strtok(NULL,"\"");
      if (ptr == NULL) error->all("Unbalanced quotes in input line");
    }
    if (arg[narg]) narg++;
    else break;
  }
}

/* ----------------------------------------------------------------------
   substitute for $ variables in str
   print updated string if flag is set and not searching for label
------------------------------------------------------------------------- */

void Input::substitute(char *str, int flag)
{
  // use work[] as scratch space to expand str
  // do not replace $ inside double quotes as flagged by level
  // var = pts at variable name, ended by NULL
  //   if $ is followed by '{', trailing '}' becomes NULL
  //   else $x becomes x followed by NULL
  // beyond = pts at text following variable

  char *var,*value,*beyond;
  int level = 0;
  char *ptr = str;

  while (*ptr) {
    if (*ptr == '$' && level == 0) {
      if (*(ptr+1) == '{') {
	var = ptr+2;
	int i = 0;
	while (var[i] != '\0' && var[i] != '}') i++;
	if (var[i] == '\0') error->one("Invalid variable name");
	var[i] = '\0';
	beyond = ptr + strlen(var) + 3;
      } else {
	var = ptr;
	var[0] = var[1];
	var[1] = '\0';
	beyond = ptr + strlen(var) + 1;
      }
      value = variable->retrieve(var);
      if (value == NULL) error->one("Substitution for undefined variable");

      *ptr = '\0';
      strcpy(work,str);
      strcat(work,value);
      strcat(work,beyond);
      strcpy(str,work);
      ptr += strlen(value);
      if (flag && me == 0 && label_active == 0) {
	if (echo_screen && screen) fprintf(screen,"%s",str); 
	if (echo_log && logfile) fprintf(logfile,"%s",str);
      }
      continue;
    }
    if (*ptr == '"') {
      if (level == 0) level = 1;
      else level = 0;
    }
    ptr++;
  }
}

/* ----------------------------------------------------------------------
   process a single parsed command
   return 0 if successful, -1 if did not recognize command
------------------------------------------------------------------------- */

int Input::execute_command()
{
  int flag = 1;

  if (!strcmp(command,"angle_coeff")) angle_coeff();
  else if (!strcmp(command,"angle_style")) angle_style();
  else if (!strcmp(command,"atom_modify")) atom_modify();
  else if (!strcmp(command,"atom_style")) atom_style();
  else if (!strcmp(command,"bond_coeff")) bond_coeff();
  else if (!strcmp(command,"bond_style")) bond_style();
  else if (!strcmp(command,"boundary")) boundary();
  else if (!strcmp(command,"clear")) clear();
  else if (!strcmp(command,"dielectric")) dielectric();
  else if (!strcmp(command,"dihedral_coeff")) dihedral_coeff();
  else if (!strcmp(command,"dihedral_style")) dihedral_style();
  else if (!strcmp(command,"dimension")) dimension();
  else if (!strcmp(command,"dipole")) dipole();
  else if (!strcmp(command,"dump")) dump();
  else if (!strcmp(command,"dump_modify")) dump_modify();
  else if (!strcmp(command,"echo")) echo();
  else if (!strcmp(command,"fix")) fix();
  else if (!strcmp(command,"fix_modify")) fix_modify();
  else if (!strcmp(command,"group")) group_command();
  else if (!strcmp(command,"improper_coeff")) improper_coeff();
  else if (!strcmp(command,"improper_style")) improper_style();
  else if (!strcmp(command,"include")) include();
  else if (!strcmp(command,"jump")) jump();
  else if (!strcmp(command,"kspace_modify")) kspace_modify();
  else if (!strcmp(command,"kspace_style")) kspace_style();
  else if (!strcmp(command,"label")) label();
  else if (!strcmp(command,"lattice")) lattice();
  else if (!strcmp(command,"log")) log();
  else if (!strcmp(command,"mass")) mass();
  else if (!strcmp(command,"min_modify")) min_modify();
  else if (!strcmp(command,"min_style")) min_style();
  else if (!strcmp(command,"neigh_modify")) neigh_modify();
  else if (!strcmp(command,"neighbor")) neighbor_command();
  else if (!strcmp(command,"newton")) newton();
  else if (!strcmp(command,"next")) next_command();
  else if (!strcmp(command,"orient")) orient();
  else if (!strcmp(command,"origin")) origin();
  else if (!strcmp(command,"pair_coeff")) pair_coeff();
  else if (!strcmp(command,"pair_modify")) pair_modify();
  else if (!strcmp(command,"pair_style")) pair_style();
  else if (!strcmp(command,"pair_write")) pair_write();
  else if (!strcmp(command,"print")) print();
  else if (!strcmp(command,"processors")) processors();
  else if (!strcmp(command,"region")) region();
  else if (!strcmp(command,"reset_timestep")) reset_timestep();
  else if (!strcmp(command,"restart")) restart();
  else if (!strcmp(command,"run_style")) run_style();
  else if (!strcmp(command,"special_bonds")) special_bonds();
  else if (!strcmp(command,"temp_modify")) temp_modify();
  else if (!strcmp(command,"temperature")) temperature();
  else if (!strcmp(command,"thermo")) thermo();
  else if (!strcmp(command,"thermo_modify")) thermo_modify();
  else if (!strcmp(command,"thermo_style")) thermo_style();
  else if (!strcmp(command,"timestep")) timestep();
  else if (!strcmp(command,"undump")) undump();
  else if (!strcmp(command,"unfix")) unfix();
  else if (!strcmp(command,"units")) units();
  else if (!strcmp(command,"variable")) variable_command();
  else flag = 0;

  if (flag) return 0;

  if (0) return 0;      // dummy line to enable else-if macro expansion

#define CommandClass
#define CommandStyle(key,Class)         \
  else if (strcmp(command,#key) == 0) { \
    Class key;                          \
    key.command(narg,arg);              \
    return 0;                           \
  }
#include "style.h"
#undef CommandClass

  return -1;
}

/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   one function for each input command executed here
------------------------------------------------------------------------- */

void Input::angle_coeff()
{
  if (domain->box_exist == 0)
    error->all("Angle_coeff command before simulation box is defined");
  if (force->angle == NULL) 
    error->all("Angle_coeff command before angle_style is defined");
  if (atom->angles_allow == 0) 
    error->all("Angle_coeff command when no angles allowed");
  force->angle->coeff(0,narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::angle_style()
{
  if (narg < 1) error->all("Illegal angle_style command");
  force->create_angle(arg[0]);
  if (force->angle) force->angle->settings(narg-1,&arg[1]);
}

/* ---------------------------------------------------------------------- */

void Input::atom_modify()
{
  if (domain->box_exist) 
    error->all("Atom_modify command after simulation box is defined");
  if (!atom) error->all("Atom_modify command before atom_style command");
  atom->modify_params(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::atom_style()
{
  if (domain->box_exist) 
    error->all("Atom_style command after simulation box is defined");
  if (narg < 1) error->all("Illegal atom_style command");

  Atom *old = atom;

  if (strcmp(arg[0],"none") == 0) error->all("Invalid atom style");

#define AtomClass
#define AtomStyle(key,Class) \
  else if (strcmp(arg[0],#key) == 0) atom = new Class(narg,arg);
#include "style.h"
#undef AtomClass

  else error->all("Invalid atom style");

  if (old) {
    atom->settings(old);
    delete old;
  }
}

/* ---------------------------------------------------------------------- */

void Input::bond_coeff()
{
  if (domain->box_exist == 0)
    error->all("Bond_coeff command before simulation box is defined");
  if (force->bond == NULL) 
    error->all("Bond_coeff command before bond_style is defined");
  if (atom->bonds_allow == 0) 
    error->all("Bond_coeff command when no bonds allowed");
  force->bond->coeff(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::bond_style()
{
  if (narg < 1) error->all("Illegal bond_style command");
  force->create_bond(arg[0]);
  if (force->bond) force->bond->settings(narg-1,&arg[1]);
}

/* ---------------------------------------------------------------------- */

void Input::boundary()
{
  if (domain->box_exist)
    error->all("Boundary command after simulation box is defined");
  domain->set_boundary(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::clear()
{
  if (narg > 0) error->all("Illegal clear command");
  sys->destroy();
  sys->create();
}

/* ---------------------------------------------------------------------- */

void Input::dielectric()
{
  if (narg != 1) error->all("Illegal dielectric command");
  force->dielectric = atof(arg[0]);
}

/* ---------------------------------------------------------------------- */

void Input::dihedral_coeff()
{
  if (domain->box_exist == 0)
    error->all("Dihedral_coeff command before simulation box is defined");
  if (force->dihedral == NULL) 
    error->all("Dihedral_coeff command before dihedral_style is defined");
  if (atom->dihedrals_allow == 0) 
    error->all("Dihedral_coeff command when no dihedrals allowed");
  force->dihedral->coeff(0,narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::dihedral_style()
{
  if (narg < 1) error->all("Illegal dihedral_style command");
  force->create_dihedral(arg[0]);
  if (force->dihedral) force->dihedral->settings(narg-1,&arg[1]);
}

/* ---------------------------------------------------------------------- */

void Input::dimension()
{
  if (domain->box_exist) 
    error->all("Dimension command after simulation box is defined");
  if (narg != 1) error->all("Illegal dimension command");
  force->dimension = atoi(arg[0]);
}

/* ---------------------------------------------------------------------- */

void Input::dipole()
{
  if (domain->box_exist == 0)
    error->all("Dipole command before simulation box is defined");
  if (narg != 2) error->all("Illegal dipole command");
  atom->set_dipole(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::dump()
{
  output->add_dump(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::dump_modify()
{
  output->modify_dump(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::echo()
{
  if (narg != 1) error->all("Illegal echo command");

  if (strcmp(arg[0],"none") == 0) {
    echo_screen = 0;
    echo_log = 0;
  } else if (strcmp(arg[0],"screen") == 0) {
    echo_screen = 1;
    echo_log = 0;
  } else if (strcmp(arg[0],"log") == 0) {
    echo_screen = 0;
    echo_log = 1;
  } else if (strcmp(arg[0],"both") == 0) {
    echo_screen = 1;
    echo_log = 1;
  } else error->all("Illegal echo command");
}

/* ---------------------------------------------------------------------- */

void Input::fix()
{
  modify->add_fix(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::fix_modify()
{
  modify->modify_fix(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::group_command()
{
  group->assign(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::improper_coeff()
{
  if (domain->box_exist == 0)
    error->all("Improper_coeff command before simulation box is defined");
  if (force->improper == NULL) 
    error->all("Improper_coeff command before improper_style is defined");
  if (atom->impropers_allow == 0) 
    error->all("Improper_coeff command when no impropers allowed");
  force->improper->coeff(0,narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::improper_style()
{
  if (narg < 1) error->all("Illegal improper_style command");
  force->create_improper(arg[0]);
  if (force->improper) force->improper->settings(narg-1,&arg[1]);
}

/* ---------------------------------------------------------------------- */

void Input::include()
{
  if (narg != 1) error->all("Illegal include command");

  if (me == 0) {
    if (nfile == maxfile) {
      maxfile++;
      infiles = (FILE **) 
        memory->srealloc(infiles,maxfile*sizeof(FILE *),"input:infiles");
    }
    infile = fopen(arg[0],"r");
    if (infile == NULL) {
      char str[128];
      sprintf(str,"Cannot open input script %s",arg[0]);
      error->one(str);
    }
    infiles[nfile++] = infile;
  }
}

/* ---------------------------------------------------------------------- */

void Input::jump()
{
  if (narg < 1 || narg > 2) error->all("Illegal jump command");

  if (jump_skip) {
    jump_skip = 0;
    return;
  }

  if (me == 0) {
    if (infile != stdin) fclose(infile);
    infile = fopen(arg[0],"r");
    if (infile == NULL) {
      char str[128];
      sprintf(str,"Cannot open input script %s",arg[0]);
      error->one(str);
    }
    infiles[nfile-1] = infile;
  }

  if (narg == 2) {
    label_active = 1;
    if (labelstr) delete [] labelstr;
    int n = strlen(arg[1]) + 1;
    labelstr = new char[n];
    strcpy(labelstr,arg[1]);
  }
}

/* ---------------------------------------------------------------------- */

void Input::kspace_modify()
{
  if (force->kspace == NULL) error->all("KSpace style has not yet been set");
  force->kspace->modify_params(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::kspace_style()
{
  force->create_kspace(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::label()
{
  if (narg != 1) error->all("Illegal label command");
  if (label_active && strcmp(labelstr,arg[0]) == 0) label_active = 0;
}

/* ---------------------------------------------------------------------- */

void Input::lattice()
{
  domain->set_lattice(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::log()
{
  if (narg != 1) error->all("Illegal log command");

  if (me == 0) {
    if (logfile) fclose(logfile);
    if (strcmp(arg[0],"none") == 0) logfile = NULL;
    else {
      char fname[128];
      if (universe->nworlds == 1) strcpy(fname,arg[0]);
      else sprintf(fname,"%s.%d",arg[0],universe->iworld);
      logfile = fopen(fname,"w");
      if (logfile == NULL) {
	char str[128];
	sprintf(str,"Cannot open logfile %s",fname);
	error->one(str);
      }
    }
    if (universe->nworlds == 1) universe->ulogfile = logfile;
  }
}

/* ---------------------------------------------------------------------- */

void Input::mass()
{
  if (domain->box_exist == 0)
    error->all("Mass command before simulation box is defined");
  if (narg != 2) error->all("Illegal mass command");
  atom->set_mass(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::min_modify()
{
  update->minimize->modify_params(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::min_style()
{
  if (domain->box_exist == 0)
    error->all("Min_style command before simulation box is defined");
  update->create_minimize(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::neigh_modify()
{
  neighbor->modify_params(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::neighbor_command()
{
  if (narg != 2) error->all("Illegal neighbor command");

  neighbor->skin = atof(arg[0]);
  if (strcmp(arg[1],"nsq") == 0) neighbor->style = 0;
  else if (strcmp(arg[1],"bin") == 0) neighbor->style = 1;
  else error->all("Illegal neighbor command");
}

/* ---------------------------------------------------------------------- */

void Input::newton()
{
  int newton_pair,newton_bond;

  if (narg == 1) {
    if (strcmp(arg[0],"off") == 0) newton_pair = newton_bond = 0;
    else if (strcmp(arg[0],"on") == 0) newton_pair = newton_bond = 1;
    else error->all("Illegal newton command");
  } else if (narg == 2) {
    if (strcmp(arg[0],"off") == 0) newton_pair = 0;
    else if (strcmp(arg[0],"on") == 0) newton_pair= 1;
    else error->all("Illegal newton command");
    if (strcmp(arg[1],"off") == 0) newton_bond = 0;
    else if (strcmp(arg[1],"on") == 0) newton_bond = 1;
    else error->all("Illegal newton command");
  } else error->all("Illegal newton command");

  force->newton_pair = newton_pair;

  if (newton_bond == 0) {
    if (domain->box_exist && force->newton_bond == 1) 
      error->all("Newton bond change after simulation box is defined");
    force->newton_bond = 0;
  } else {
    if (domain->box_exist && force->newton_bond == 0) 
      error->all("Newton bond change after simulation box is defined");
    force->newton_bond = 1;
  }

  if (newton_pair || newton_bond) force->newton = 1;
  else force->newton = 0;
}

/* ---------------------------------------------------------------------- */

void Input::next_command()
{
  if (variable->next(narg,arg)) jump_skip = 1;
}

/* ---------------------------------------------------------------------- */

void Input::orient()
{
  if (narg != 4) error->all("Illegal orient command");

  if (strcmp(arg[0],"x") == 0) {
    domain->orient_x[0] = atoi(arg[1]);
    domain->orient_x[1] = atoi(arg[2]);
    domain->orient_x[2] = atoi(arg[3]);
  } else if (strcmp(arg[0],"y") == 0) {
    domain->orient_y[0] = atoi(arg[1]);
    domain->orient_y[1] = atoi(arg[2]);
    domain->orient_y[2] = atoi(arg[3]);
  } else if (strcmp(arg[0],"z") == 0) {
    domain->orient_z[0] = atoi(arg[1]);
    domain->orient_z[1] = atoi(arg[2]);
    domain->orient_z[2] = atoi(arg[3]);
  } else error->all("Illegal orient command");
}

/* ---------------------------------------------------------------------- */

void Input::origin()
{
  if (narg != 3) error->all("Illegal origin command");

  domain->origin_x = atof(arg[0]);
  domain->origin_y = atof(arg[1]);
  domain->origin_z = atof(arg[2]);

  if (domain->origin_x < 0.0 || domain->origin_x > 1.0 ||
      domain->origin_y < 0.0 || domain->origin_y > 1.0 ||
      domain->origin_x < 0.0 || domain->origin_z > 1.0)
    error->all("Illegal origin command");
}

/* ---------------------------------------------------------------------- */

void Input::pair_coeff()
{
  if (domain->box_exist == 0)
    error->all("Pair_coeff command before simulation box is defined");
  if (force->pair == NULL) 
    error->all("Pair_coeff command before pair_style is defined");
  force->pair->coeff(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::pair_modify()
{
  if (force->pair == NULL) 
    error->all("Pair_modify command before pair_style is defined");
  force->pair->modify_params(narg,arg);
}

/* ----------------------------------------------------------------------
   if old pair style exists and new style is same, just change settings
   else create new pair class
------------------------------------------------------------------------- */

void Input::pair_style()
{
  if (narg < 1) error->all("Illegal pair_style command");
  if (force->pair && strcmp(arg[0],force->pair_style) == 0) {
    force->pair->settings(narg-1,&arg[1]);
    return;
  }
  force->create_pair(arg[0]);
  if (force->pair) force->pair->settings(narg-1,&arg[1]);
}

/* ---------------------------------------------------------------------- */

void Input::pair_write()
{
  if (force->pair == NULL) 
    error->all("Pair_write command before pair_style is defined");
  force->pair->write_file(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::print()
{
  if (narg < 1) error->all("Illegal print command");

  char *str = new char[MAXLINE];
  str[0] = '\0';
  for (int iarg = 0; iarg < narg; iarg++) {
    strcat(str,arg[iarg]);
    strcat(str," ");
  }
  
  if (me == 0) {
    if (screen) fprintf(screen,"%s\n",str);
    if (logfile) fprintf(logfile,"%s\n",str);
  }

  delete [] str;
}

/* ---------------------------------------------------------------------- */

void Input::processors()
{
  if (domain->box_exist)
    error->all("Processors command after simulation box is defined");
  if (narg != 3) error->all("Illegal processors command");
  comm->user_procgrid[0] = atoi(arg[0]);
  comm->user_procgrid[1] = atoi(arg[1]);
  comm->user_procgrid[2] = atoi(arg[2]);
}

/* ---------------------------------------------------------------------- */

void Input::region()
{
  domain->add_region(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::reset_timestep()
{
  if (narg != 1) error->all("Illegal reset_timestep command");
  update->ntimestep = atoi(arg[0]);
  if (update->ntimestep < 0) error->all("Timestep must be >= 0");
}

/* ---------------------------------------------------------------------- */

void Input::restart()
{
  output->create_restart(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::run_style()
{
  if (domain->box_exist == 0)
    error->all("Run_style command before simulation box is defined");
  update->create_integrate(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::special_bonds()
{
  // store 1-3,1-4 values before change

  double lj2 = force->special_lj[2];
  double lj3 = force->special_lj[3];
  double coul2 = force->special_coul[2];
  double coul3 = force->special_coul[3];

  force->set_special(narg,arg);

  // if simulation box defined and 1-3,1-4 values changed, redo special list

  if (domain->box_exist && atom->molecular) {
    if (lj2 != force->special_lj[2] || lj3 != force->special_lj[3] ||
	coul2 != force->special_coul[2] || coul3 != force->special_coul[3]) {
      Special special;
      special.build();
    }
  }
}

/* ---------------------------------------------------------------------- */

void Input::temp_modify()
{
  force->modify_temp(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::temperature()
{
  force->add_temp(narg,arg,0);
}

/* ---------------------------------------------------------------------- */

void Input::thermo()
{
  if (narg != 1) error->all("Illegal thermo command");
  output->thermo_every = atoi(arg[0]);
}

/* ---------------------------------------------------------------------- */

void Input::thermo_modify()
{
  output->thermo->modify_params(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::thermo_style()
{
  output->create_thermo(narg,arg);
}

/* ---------------------------------------------------------------------- */

void Input::timestep()
{
  if (narg != 1) error->all("Illegal timestep command");
  update->dt = atof(arg[0]);
}

/* ---------------------------------------------------------------------- */

void Input::undump()
{
  if (narg != 1) error->all("Illegal undump command");
  output->delete_dump(arg[0]);
}

/* ---------------------------------------------------------------------- */

void Input::unfix()
{
  if (narg != 1) error->all("Illegal unfix command");
  modify->delete_fix(arg[0]);
}

/* ---------------------------------------------------------------------- */

void Input::units()
{
  if (domain->box_exist) 
    error->all("Units command after simulation box is defined");
  if (narg != 1) error->all("Illegal units command");
  update->set_units(arg[0]);
}

/* ---------------------------------------------------------------------- */

void Input::variable_command()
{
  variable->set(narg,arg);
}
