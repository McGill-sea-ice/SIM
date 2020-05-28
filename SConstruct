from glob import glob

# This file sets up the compilation environment for the sea ice model and its libraries.

# Refer questions to 
# David Huard <david.huard@gmail.com>
# McGill U.

# Sept. 18, 2008


# ----------------------------------
#         SCONS Variables
# ----------------------------------

# LIBPATH : The path searched for libraries.
# FORTRANMODDIR : The path where .mod files are stored, ie compiled modules.
# FORTRANMODDIRPREFIX : The compiler prefix that prefixes FORTRANMODDIR (-J for gfortran)
# F90PATH : The path searched for include files.
# F90FLAGS : Flags for the F90 compiler.


# ----------------------------------
#              Notes
# ----------------------------------
#
# Flags
# -----
#   -fbackslash : Allow escape characters \n, \t, etc in strings. 
#                 In gcc-4.3, -fno-backslash is the default.
#   -Wall : Print all warnings.

# ----------------------------------
#  Standard Compilation Environment
# ----------------------------------

include = Dir('#/include')
libs = Dir('#/libs')
FC = F90 = 'gfortran'

env = Environment(LIBPATH=[libs],
                  FORTRANMODDIR=[include,], 
                  FORTRANMODDIRPREFIX='-J', 
                  F90PATH=[include,], 
                  F90FLAGS=['-fPIC','-fbackslash', '-O3', '-ffast-math'], 
                  FORTRAN=FC,
                  LINKFLAGS=[],
                  F90=FC)

# Command line options to modify the build environment.
# Use scons debug=1 to compile with debugging flags. 


DEBUG = ARGUMENTS.get('debug', 0)
EXE = ARGUMENTS.get('exe', 'zoupa')


if int(DEBUG):  
    env.Append(F90FLAGS = '-g')
    env.Append(F90FLAGS = '-static')
    env.Append(F90FLAGS = '-fbacktrace')
    env.Append(LINKFLAGS = '-fbacktrace')


# ---------------------------------------
#           NetCDF4 Librairies
# ---------------------------------------

env.Append(LIBPATH='/aos/shared/lib/el7/netcdf/default/gnu/lib')
env.Append(F90PATH='/aos/shared/lib/el7/netcdf/default/gnu/include')
env.Append(LIBPATH='/aos/shared/lib/el7/hdf5/default/gnu/lib')
env.Append(F90PATH='/aos/shared/lib/el7/hdf5/default/gnu/include')
#env.Append(LIBPATH='/usr/lib64')

# ----------------------------------------
#             Create Executables
# ----------------------------------------

# Export the compilation environment so SConscript files in subdirectories can access it. 

Export('env')
Export('EXE')
# The datetime module
datetime = SConscript('datetime/SConscript', variant_dir='build/datetime')

# Tests for the datetime module
datetime_tests = SConscript('datetime/tests/SConscript', variant_dir='build/datetime_tests')

# The Fortran Unit Test library
fruit = SConscript('fruit/SConscript', variant_dir='build/fruit')

# Tests for the Fortran Unit Test Library
fruit_tests = SConscript('fruit/tests/SConscript', variant_dir='build/fruit_tests')

# The Sea Ice Model library and executable
sim = SConscript('src/SConscript', variant_dir='build/sim')

# Tests for the Sea Ice Model
# test_sim = SConscript('src/tests/SConscript', variant_dir='build/sim_tests')

