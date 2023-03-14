###FC = mpinfort
FC = mpinfort -compiler /opt/nec/ve/nfort/3.0.4/bin/nfort
FFLAGS  = -report-all -O3 -fpp $(EXTRA) #-mparallel
FFLAGS  += -fdiag-vector=2 -fdiag-inline=2 -fdiag-parallel=2
FFLAGS_OMP1= #-fopenmp #around FFT
FFLAGS_OMP2= #-fopenmp #others
FFLAGS_OMP3= #-fopenmp #bndry&advnc

ifneq ("x$(FFLAGS_OMP1)_$(FFLAGS_OMP2)_$(FFLAGS_OMP3)_","x___")
   FFLAGS_LOMP= -fopenmp
endif

PROG = 'gkvp.exe'

SDIR = src
SRC = ../$(SDIR)/
MYL = ../lib/

MATH = gkvp_math_portable

FFT = gkvp_f0.56_fft_fftw_tune2r_0813
### Usage of FFTW
ifeq ($(FFT),gkvp_f0.56_fft_fftw_tune2r_0813)
  NLC_HOME=/opt/nec/ve/nlc/2.1.0
  INC = -I$(NLC_HOME)/include
  LIB = -L$(NLC_HOME)/lib  -laslfftw3 -lasl_sequential -ftrace
  #LIB = -L$(NLC_HOME)/lib  -laslfftw3 -lasl_openmp
endif

FILEIO=gkvp_fileio_fortran
#FILEIO=gkvp_fileio_netcdf
### Usage of NetCDF (module load netcdf-parallelIO-fortran-sx)
ifeq ($(FILEIO),gkvp_fileio_netcdf)
  FC = mpinfort
  #INC += -I$(NFORT_INCLUDE_PATH)
  #LIB += -L$(NFORT_LIBRARY_PATH) -lnetcdff -lnetcdf -lhdf5_hl -lhdf5
  LIB += -lnetcdff -lnetcdf -lhdf5_hl -lhdf5
endif


gkvp:	$(SRC)gkvp_header.f90\
	$(SRC)gkvp_mpienv.f90\
	$(MYL)$(MATH).f90\
	$(SRC)gkvp_clock.f90\
	$(SRC)$(FILEIO).f90\
	$(SRC)gkvp_intgrl.f90\
	$(SRC)gkvp_tips.f90\
	$(SRC)gkvp_vmecbzx.f90\
	$(SRC)gkvp_igs.f90\
	$(SRC)gkvp_ring.f90\
	$(SRC)gkvp_f0.56_bndry_tune_nec1.f90\
	$(SRC)gkvp_f0.56_colli_tune_nifs.f90\
	$(SRC)$(FFT).f90\
	$(SRC)gkvp_fld.f90\
	$(SRC)gkvp_colliimp.f90\
	$(SRC)gkvp_freq.f90\
	$(SRC)gkvp_f0.56_zfilter_tune_nec1.f90\
	$(SRC)gkvp_f0.56_exb_tune2r_0813.f90\
	$(SRC)gkvp_trans.f90\
	$(SRC)gkvp_shearflow.f90\
	$(SRC)gkvp_geom.f90\
	$(SRC)gkvp_advnc.f90\
	$(SRC)gkvp_dtc.f90\
	$(SRC)gkvp_out.f90\
	$(SRC)gkvp_set.f90\
	$(SRC)gkvp_main.f90

	$(FC) $(FFLAGS) -c $(SRC)gkvp_header.f90
	$(FC) $(FFLAGS) -c $(SRC)gkvp_mpienv.f90
	$(FC) $(FFLAGS) -c $(MYL)$(MATH).f90
	$(FC) $(FFLAGS) -c $(SRC)gkvp_clock.f90
	$(FC) $(FFLAGS) -c $(SRC)$(FILEIO).f90 $(INC)
	$(FC) $(FFLAGS) -c $(SRC)gkvp_intgrl.f90
	$(FC) $(FFLAGS) -c $(SRC)gkvp_tips.f90
	$(FC) $(FFLAGS) -c $(SRC)gkvp_vmecbzx.f90
	$(FC) $(FFLAGS) -c $(SRC)gkvp_igs.f90
	$(FC) $(FFRAGS) -c $(SRC)gkvp_ring.f90
	$(FC) $(FFLAGS) $(FFLAGS_OMP3) -c $(SRC)gkvp_f0.56_bndry_tune_nec1.f90
	$(FC) $(FFLAGS) $(FFLAGS_OMP2) -c $(SRC)gkvp_f0.56_colli_tune_nifs.f90
	$(FC) $(FFLAGS) $(FFLAGS_OMP1) -c $(SRC)$(FFT).f90 $(INC)
	$(FC) $(FFLAGS) -c $(SRC)gkvp_fld.f90
	$(FC) $(FFLAGS) -c $(SRC)gkvp_colliimp.f90
	$(FC) $(FFLAGS) -c $(SRC)gkvp_freq.f90
	$(FC) $(FFLAGS) $(FFLAGS_OMP2) -c $(SRC)gkvp_f0.56_zfilter_tune_nec1.f90
	$(FC) $(FFLAGS) $(FFLAGS_OMP1) -c $(SRC)gkvp_f0.56_exb_tune2r_0813.f90
	$(FC) $(FFLAGS) $(FFLAGS_OMP1) -c $(SRC)gkvp_trans.f90
	$(FC) $(FFLAGS) -c $(SRC)gkvp_shearflow.f90
	$(FC) $(FFLAGS) -c $(SRC)gkvp_geom.f90
	$(FC) $(FFLAGS) $(FFLAGS_OMP3) -c $(SRC)gkvp_advnc.f90
	$(FC) $(FFLAGS) -c $(SRC)gkvp_dtc.f90
	$(FC) $(FFLAGS) $(FFLAGS_OMP2) -c $(SRC)gkvp_out.f90
	$(FC) $(FFLAGS) -c $(SRC)gkvp_set.f90
	$(FC) $(FFLAGS) -c $(SRC)gkvp_main.f90

	$(FC) $(FFLAGS)	$(FFLAGS_LOMP) \
			gkvp_header.o\
			gkvp_mpienv.o\
			$(MATH).o\
			gkvp_clock.o\
			$(FILEIO).o\
			gkvp_intgrl.o\
			gkvp_tips.o\
			gkvp_vmecbzx.o\
			gkvp_igs.o\
			gkvp_ring.o\
			gkvp_f0.56_bndry_tune_nec1.o\
			gkvp_f0.56_colli_tune_nifs.o\
			$(FFT).o\
			gkvp_fld.o\
			gkvp_colliimp.o\
			gkvp_freq.o\
			gkvp_f0.56_zfilter_tune_nec1.o\
			gkvp_f0.56_exb_tune2r_0813.o\
			gkvp_trans.o\
			gkvp_shearflow.o\
			gkvp_geom.o\
			gkvp_advnc.o\
			gkvp_dtc.o\
			gkvp_out.o\
			gkvp_set.o\
                        gkvp_main.o\
			-o $(PROG) $(LIB)

	\cp *.L *.o *.mod ../$(SDIR)/
	\rm -f *.L *.o *.mod

clean:
	rm -f  ../$(SDIR)/*.LL ../$(SDIR)/*.L ../$(SDIR)/*.o ../$(SDIR)/*.mod ../$(SDIR)/*.lst 
	rm -f ./*.exe ./sub.q.* ./gkvp_namelist.*

clear:
	rm -f ./*.o ./*.mod ./*.L ./*.LL
 
