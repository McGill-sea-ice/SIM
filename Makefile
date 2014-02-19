## INSTRUCTION: make EXE='executable_filename' ice

# Makefile for my sea-ice model
# compilation by 
FC	=	gfortran -O3 -ffast-math
#FC	=	ifc -w# fortran 95

# -w above means that the warning messages are ignored...there are a few left
# but they don't affect the compilation

#LIBDIR = -L/home/hdx1/gavin/lib
LIBS  = #-llapack -lblas #-lnag
#FFLAGS 	= -s -O -Pv -lc    # used preprocessing
#FFLAGS  = -g  #  used for debug
#FFLAGS  = -pg  #  used for profiling run the code and then type gprof
#FFlAGS  = -O0  
OBJ	= src/ice.o src/var_analysis.o src/lsor.o src/leapyear.o src/julianday.o src/oceanTadv.o src/ddot.o src/ini_get.o src/advection.o src/vect_v_3dag.o src/prep_res_norm.o src/fluxy.o src/T_atm.o src/thermodynamic.o src/lagrangiantracer.o src/prep_fgmres.o src/res_norm_conv.o src/bc_get.o src/vect_u_3dag.o src/r1pr2p.o src/matvec.o src/uvsolve_sor.o src/KEcalc.o src/ocn_Tclim.o src/heatfluxes.o src/precond_lsor.o src/sea_ice_post.o src/stepper.o src/par_get.o src/update_tracer.o src/res_norm.o src/oceanTthermo.o src/dcopy.o src/ocn_current.o src/r1ppr2pp.o src/wind_forcing.o src/viscouscoefficient.o src/precond_sor.o  src/uvsolve_NR_B.o src/fgmresD.o src/chdate.o src/transformer.o src/atmosphere.o src/julian2date.o src/pressure.o src/daxpy.o src/shortwave.o src/fluxx.o src/tridag.o src/uvchecker.o  src/res_normNL.o

 

ice: $(OBJ) $(COMMON) 
	$(FC)  $(FFLAGS) -o $(EXE) $(OBJ)  $(LIBS)

lib: ice
	ar rv libice.a $(OBJ)

tidy:	
	rm src/*.o src/*~ 


# $(LIBDIR)  $(LIBS)

 
#$@
