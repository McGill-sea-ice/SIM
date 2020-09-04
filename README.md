# SIM
McGill sea ice model

CONFIGURATIONS:
  In order to compile in Fortran 4.8.5 (master version is 4.4.7)
  PROGRAM test_datetime >> PROGRAM datetime_test in "datetime/tests/test_datetime.f90",
  function "to_s_1d_logical_" had twice "CHARACTER(len=500) :: out",  comment line 163 in "fruit/fruit_util.f90"
  comment "use datetime, only: operator(==), operator(/=)" in  "src/stepper.f90"
  netCDF and HDF5 libraries in the Sconscript settings need to compiled with gfortran 4.8.5, 

OPERATIONS:
compile with "scons"
run experience with "./zoupa input_norestart >outEXPERIENCENo" or "./zoupa nano input_norestart >outEXPERIENCENo" 


  
KNOWN BUGS:
Sometimes when adding changes the compilation doesn't include them. Before compiling try : rm -r build/ zoupa libs/*
