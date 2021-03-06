FC = mpif90

# Common compiler flags
NCFLAG = -I $(NETCDF_ROOT)/include 
#MPIFLAG = -Dusempi3
FFLAGS = -xHost -ftz -fp-model precise $(MPIFLAG) $(NCFLAG)
LIBS = -L $(NETCDF_ROOT)/lib -lnetcdf -lnetcdff
PPFLAG90 = -fpp
PPFLAG77 = -fpp
PPFLAG90F = -fpp
REAL8FLAG = -r8
INT8FLAG = -i8
DEBUGFLAG = -check all -debug all -traceback -fpe0

# Gfortran compiler options
ifeq ($(GFORTRAN),yes)
MPIFC = gfortran
MPIF77 = gfortran
FC = mpif90
NCFLAG += -Dusempif -Dusenc3
FFLAGS = -O2 -mtune=native -march=native $(MPIFLAG) $(NCFLAG)
PPFLAG90 = -x f95-cpp-input
PPFLAG77 = -x f77-cpp-input
PPFLAG90F =
REAL8FLAG = -fdefault-real-8
INT8FLAG = -fdefault-int-8
DEBUGFLAG = -g -Wall -Wextra -fbounds-check -fbacktrace
endif

# Options for building with VAMPIRTrace
ifeq ($(VT),yes)
FC = vtfort -vt:fc mpif90 -vt:inst manual
FFLAGS += -Dvampir -DVTRACE
else
FFLAGS += -Dsimple_timer
endif

#Decomposition method
ifeq ($(DECOMP),uniform)
FFLAGS += -Duniform_decomp
endif

# Testing - I/O and fpmodel
ifeq ($(TEST),yes)
FFLAGS += -Doutsync $(DEBUGFLAG)
endif

# Build with 64 ints/reals
ifeq ($(I8R8),yes)
FFLAGS += $(REAL8FLAG) $(INT8FLAG) -Di8r8
endif

OBJS = gettin.o \
globpe.o  indata.o infile.o jimcc.o outcdf.o setxyz.o sflux.o soilsnow.o latltoij.o \
zenith.o cc_mpi.o diag_m.o sumdd_m.o utilities.o onthefly.o stacklimit.o \
xyzinfo_m.o vecsuv_m.o map_m.o latlong_m.o indices_m.o bigxy4_m.o arrays_m.o \
histave_m.o morepbl_m.o nsibd_m.o parmhdff_m.o pbl_m.o \
permsurf_m.o prec_m.o savuvt_m.o savuv1_m.o sigs_m.o soil_m.o soilsnow_m.o \
workglob_m.o work2_m.o work3_m.o work3b_m.o \
cable_ccam2.o cable_common.o cable_data.o cable_define_types.o cable_roughness.o \
cable_soilsnow.o river.o netcdf_m.o mpif_m.o

globpea: $(OBJS)
	$(FC) -o globpea $(FFLAGS) $(OBJS) $(LIBS)

clean:
	rm *.o *.i *.mod globpea

.SUFFIXES:.f90 .F90

netcdf_m.o: netcdf_m.f90
	$(FC) -c $(PPFLAG90) $(NCFLAG) $<
mpif_m.o: mpif_m.f90
	$(FC) -c $(PPFLAG90) $(MPIFLAG) $<
stacklimit.o: stacklimit.c
	cc -c stacklimit.c
version.h: FORCE
	rm -f brokenver tmpver
	echo "      character(len=*), parameter :: version ='CCAM r'" > brokenver
	echo "      character(len=*), parameter :: version ='CCAM r`svnversion .`'" > tmpver
	grep exported tmpver || grep Unversioned tmpver || cmp tmpver brokenver || cmp tmpver version.h || mv tmpver version.h
FORCE:

.f90.o:
	$(FC) -c $(FFLAGS) $(PPFLAG90) $<
.F90.o:
	$(FC) -c $(FFLAGS) $(PPFLAG90F) $<	
.f.o:
	$(FC) -c $(FFLAGS) $(PPFLAG77) $<

# Remove mod rule from Modula 2 so GNU make doesn't get confused
%.o : %.mod

# Dependencies
cable_common.o : cable_define_types.o
cable_ccam2.o : arrays_m.o cable_common.o cable_define_types.o cable_roughness.o cable_soilsnow.o cc_mpi.o infile.o latlong_m.o morepbl_m.o nsibd_m.o pbl_m.o permsurf_m.o prec_m.o sigs_m.o soil_m.o soilsnow_m.o work2_m.o work3_m.o zenith.o const_phys.h darcdf.h dates.h newmpar.h parm.h parmgeom.h soilv.h
cable_roughness.o : cable_common.o cable_data.o cable_define_types.o
cable_soilsnow.o : cable_common.o cable_data.o cable_define_types.o
cc_mpi.o : arrays_m.o indices_m.o latlong_m.o map_m.o mpif_m.o sigs_m.o sumdd_m.o vecsuv_m.o xyzinfo_m.o newmpar.h parm.h 
diag_m.o : cc_mpi.o sigs_m.o sumdd_m.o xyzinfo_m.o newmpar.h parm.h
gettin.o : arrays_m.o savuvt_m.o newmpar.h 
globpe.o : arrays_m.o bigxy4_m.o cable_ccam2.o cc_mpi.o diag_m.o histave_m.o indata.o indices_m.o infile.o latlong_m.o map_m.o morepbl_m.o nsibd_m.o outcdf.o parmhdff_m.o pbl_m.o permsurf_m.o prec_m.o river.o savuvt_m.o savuv1_m.o sigs_m.o soil_m.o soilsnow_m.o vecsuv_m.o workglob_m.o work2_m.o work3_m.o xyzinfo_m.o const_phys.h darcdf.h dates.h filnames.h kuocom.h newmpar.h parm.h parmdyn.h parmgeom.h parmhor.h parmsurf.h soilv.h stime.h version.h
indata.o : arrays_m.o bigxy4_m.o cable_ccam2.o cc_mpi.o diag_m.o indices_m.o infile.o latlong_m.o map_m.o morepbl_m.o nsibd_m.o onthefly.o pbl_m.o permsurf_m.o river.o sigs_m.o soil_m.o soilsnow_m.o vecsuv_m.o xyzinfo_m.o const_phys.h darcdf.h dates.h filnames.h newmpar.h parm.h parmdyn.h parmgeom.h soilv.h stime.h
indices_m.o : newmpar.h
infile.o : cc_mpi.o netcdf_m.o dates.h newmpar.h parm.h parmgeom.h
latltoij.o : utilities.o const_phys.h newmpar.h parm.h parmdyn.h
onthefly.o : cc_mpi.o cable_define_types.o diag_m.o infile.o latlong_m.o morepbl_m.o nsibd_m.o river.o savuvt_m.o savuv1_m.o sigs_m.o soil_m.o utilities.o vecsuv_m.o workglob_m.o work2_m.o const_phys.h darcdf.h kuocom.h newmpar.h parm.h parmdyn.h parmgeom.h soilv.h stime.h
outcdf.o : arrays_m.o cable_ccam2.o cable_define_types.o cc_mpi.o histave_m.o infile.o latlong_m.o map_m.o morepbl_m.o nsibd_m.o parmhdff_m.o pbl_m.o prec_m.o river.o savuvt_m.o savuv1_m.o sigs_m.o soil_m.o soilsnow_m.o work2_m.o const_phys.h dates.h filnames.h kuocom.h newmpar.h parm.h parmdyn.h parmgeom.h parmhor.h parmsurf.h soilv.h version.h
river.o : arrays_m.o cable_ccam2.o cc_mpi.o indices_m.o map_m.o nsibd_m.o soil_m.o soilsnow_m.o xyzinfo_m.o const_phys.h newmpar.h parm.h soilv.h
setxyz.o : cc_mpi.o indices_m.o latlong_m.o map_m.o utilities.o workglob_m.o const_phys.h newmpar.h parm.h
sflux.o : arrays_m.o cable_ccam2.o cc_mpi.o diag_m.o latlong_m.o map_m.o morepbl_m.o nsibd_m.o pbl_m.o permsurf_m.o prec_m.o river.o savuvt_m.o sigs_m.o soil_m.o soilsnow_m.o vecsuv_m.o work2_m.o work3_m.o xyzinfo_m.o const_phys.h dates.h newmpar.h parm.h parmgeom.h parmsurf.h soilv.h
soilsnow.o : arrays_m.o cc_mpi.o diag_m.o morepbl_m.o nsibd_m.o permsurf_m.o sigs_m.o soil_m.o soilsnow_m.o work2_m.o work3_m.o work3b_m.o const_phys.h newmpar.h parm.h soilv.h
utilities.o : const_phys.h 
