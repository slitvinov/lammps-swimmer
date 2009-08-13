# Install/unInstall package classes in LAMMPS
# edit Makefile.package to include/exclude REAX library

if ($1 == 1) then

  sed -i -e 's/[^ \t]*reax //' ../Makefile.package
  sed -i -e 's|^PKGINC =[ \t]*|&-I../../lib/reax |' ../Makefile.package
  sed -i -e 's|^PKGPATH =[ \t]*|&-L../../lib/reax |' ../Makefile.package
  sed -i -e 's|^PKGLIB =[ \t]*|&-lreax |' ../Makefile.package

  cp style_reax.h ..

  cp pair_reax.cpp ..

  cp pair_reax.h ..
  cp pair_reax_fortran.h ..

  cp fix_reax_bonds.h ..
  cp fix_reax_bonds.cpp ..

else if ($1 == 0) then

  sed -i -e 's/[^ \t]*reax //' ../Makefile.package

  rm ../style_reax.h
  touch ../style_reax.h

  rm ../pair_reax.cpp

  rm ../pair_reax.h
  rm ../pair_reax_fortran.h

  rm ../fix_reax_bonds.h
  rm ../fix_reax_bonds.cpp

endif
