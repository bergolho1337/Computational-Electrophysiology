#!/bin/bash

if [ ! -f meshGenerator ]; then
	make
fi

SK_FOLDER="Skeletons"
MSH_FOLDER="Meshes"

./meshGenerator $SK_FOLDER/cable2-05cm.vtk $MSH_FOLDER/cable2-05cm-alien.msh --alien
./meshGenerator $SK_FOLDER/cable2-05cm.vtk $MSH_FOLDER/cable2-05cm-dog.msh --dog
./meshGenerator $SK_FOLDER/cable2-05cm.vtk $MSH_FOLDER/cable2-05cm-orc.msh --orc
./meshGenerator $SK_FOLDER/cable2-05cm.vtk $MSH_FOLDER/cable2-05cm-pig.msh --pig

./meshGenerator $SK_FOLDER/cable2-1cm.vtk $MSH_FOLDER/cable2-1cm-alien.msh --alien
./meshGenerator $SK_FOLDER/cable2-1cm.vtk $MSH_FOLDER/cable2-1cm-dog.msh --dog
./meshGenerator $SK_FOLDER/cable2-1cm.vtk $MSH_FOLDER/cable2-1cm-orc.msh --orc
./meshGenerator $SK_FOLDER/cable2-1cm.vtk $MSH_FOLDER/cable2-1cm-pig.msh --pig

./meshGenerator $SK_FOLDER/cable2-2cm.vtk $MSH_FOLDER/cable2-2cm-alien.msh --alien
./meshGenerator $SK_FOLDER/cable2-2cm.vtk $MSH_FOLDER/cable2-2cm-dog.msh --dog
./meshGenerator $SK_FOLDER/cable2-2cm.vtk $MSH_FOLDER/cable2-2cm-orc.msh --orc
./meshGenerator $SK_FOLDER/cable2-2cm.vtk $MSH_FOLDER/cable2-2cm-pig.msh --pig

./meshGenerator $SK_FOLDER/cable2-5cm.vtk $MSH_FOLDER/cable2-5cm-alien.msh --alien
./meshGenerator $SK_FOLDER/cable2-5cm.vtk $MSH_FOLDER/cable2-5cm-dog.msh --dog
./meshGenerator $SK_FOLDER/cable2-5cm.vtk $MSH_FOLDER/cable2-5cm-orc.msh --orc
./meshGenerator $SK_FOLDER/cable2-5cm.vtk $MSH_FOLDER/cable2-5cm-pig.msh --pig

./meshGenerator $SK_FOLDER/biff2-05cm.vtk $MSH_FOLDER/biff2-05cm-alien.msh --alien
./meshGenerator $SK_FOLDER/biff2-05cm.vtk $MSH_FOLDER/biff2-05cm-dog.msh --dog
./meshGenerator $SK_FOLDER/biff2-05cm.vtk $MSH_FOLDER/biff2-05cm-orc.msh --orc
./meshGenerator $SK_FOLDER/biff2-05cm.vtk $MSH_FOLDER/biff2-05cm-pig.msh --pig

./meshGenerator $SK_FOLDER/biff2-1cm.vtk $MSH_FOLDER/biff2-1cm-alien.msh --alien
./meshGenerator $SK_FOLDER/biff2-1cm.vtk $MSH_FOLDER/biff2-1cm-dog.msh --dog
./meshGenerator $SK_FOLDER/biff2-1cm.vtk $MSH_FOLDER/biff2-1cm-orc.msh --orc
./meshGenerator $SK_FOLDER/biff2-1cm.vtk $MSH_FOLDER/biff2-1cm-pig.msh --pig

./meshGenerator $SK_FOLDER/biff2-2cm.vtk $MSH_FOLDER/biff2-2cm-alien.msh --alien
./meshGenerator $SK_FOLDER/biff2-2cm.vtk $MSH_FOLDER/biff2-2cm-dog.msh --dog
./meshGenerator $SK_FOLDER/biff2-2cm.vtk $MSH_FOLDER/biff2-2cm-orc.msh --orc
./meshGenerator $SK_FOLDER/biff2-2cm.vtk $MSH_FOLDER/biff2-2cm-pig.msh --pig

./meshGenerator $SK_FOLDER/biff2-5cm.vtk $MSH_FOLDER/biff2-5cm-alien.msh --alien
./meshGenerator $SK_FOLDER/biff2-5cm.vtk $MSH_FOLDER/biff2-5cm-dog.msh --dog
./meshGenerator $SK_FOLDER/biff2-5cm.vtk $MSH_FOLDER/biff2-5cm-orc.msh --orc
./meshGenerator $SK_FOLDER/biff2-5cm.vtk $MSH_FOLDER/biff2-5cm-pig.msh --pig
