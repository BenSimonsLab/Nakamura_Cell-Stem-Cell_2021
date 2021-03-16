# This script compiles two versions of the simulation taking different data files
# (/data/init_gfra1.csv and /data/init_ngn3.csv, respectively) as initial condition.

# To run the simulation, consult comments in the simulation code or see run.sh for
# examples.

gfortran -cpp -Dgfr transplantation.f90 -ffree-line-length-none -fbackslash -Ofast -o transplantation_gfra1

gfortran -cpp -Dngn transplantation.f90 -ffree-line-length-none -fbackslash -Ofast -o transplantation_ngn3
