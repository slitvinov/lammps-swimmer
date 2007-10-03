# Install/unInstall package classes in LAMMPS

if ($1 == 1) then

  cp -p style_user_ackland.h ..

  cp -p compute_ackland_atom.cpp ..

  cp -p compute_ackland_atom.h ..

else if ($1 == 0) then

  rm ../style_user_ackland.h
  touch ../style_user_ackland.h

  rm ../compute_ackland_atom.cpp

  rm ../compute_ackland_atom.h

endif
