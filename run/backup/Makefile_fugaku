### Fujitsu Fortran Compiler ###
FC = mpifrtpx
FFLAGS = -Kfast,parallel # Optimization
FFLAGS += -X9 # Fortran95
FFLAGS += -Koptmsg=2 -Nlst=t # Optimization report
FFLAGS += -fw # Suppress message
FFLAGS += -Kopenmp #-Nfjomplib # OpenMP
FFLAGS += -mcmodel=large # Static memory larger than 2GB
#FFLAGS += -Haefosux -NRtrap #-O0 # Debug
OPTRPT = 'lst'
#FFLAGS += -Nfjprof # Fujitsu profiler fapp
#FFLAGS += -Ksimd_nouse_multiple_structures # Specific option for compiler tcs1.2.26 to avoid slowing down GKV
#FFLAGS += -Knosch_pre_ra # Specific option for compiler tcs1.2.26 to avoid slowing down GKV


PROG = 'gkvp.exe'

SRC = ../src/
MYL = ../lib/

MATH = gkvp_math_portable

FFT = gkvp_fft_fftw
### Usage of FFTW
ifeq ($(FFT),gkvp_fft_fftw)
  ### FFTW-SVE
  FFTW_DIR=/home/apps/r/OSS_CN/fftw-3.3.8/
  INC = -I$(FFTW_DIR)/include
  LIB = -L$(FFTW_DIR)/lib64 -lfftw3 -lm -SSL2
  #### FFTW-SPACK (. /home/apps/oss/spack/share/spack/setup-env.sh; spack load fftw) ###
  #FFTW_DIR=`spack location -i fftw`
  #INC = -I$(FFTW_DIR)/include
  #LIB = -L$(FFTW_DIR)/lib -lfftw3 -lm -SSL2
endif

FILEIO=gkvp_fileio_fortran
#FILEIO=gkvp_fileio_netcdf
### Usage of NetCDF (. /home/apps/oss/spack/share/spack/setup-env.sh; spack load netcdf-fortran%fj) ###
###   Operation of NetCDF has not yet been checked on Fugaku, Jan 26 2021
ifeq ($(FILEIO),gkvp_fileio_netcdf)
  NETCDF_DIR=`spack location -i netcdf-fortran%fj`
  INC += -I$(NETCDF_DIR)/include
  LIB += -L$(NETCDF_DIR)/lib -lnetcdff -lnetcdf -lhdf5_hl -lhdf5
endif

OBJS =	gkvp_header.o\
	gkvp_mpienv.o\
	$(MATH).o\
	gkvp_clock.o\
	$(FILEIO).o\
	gkvp_intgrl.o\
	gkvp_tips.o\
	gkvp_vmecbzx.o\
	gkvp_igs.o\
	gkvp_ring.o\
	gkvp_bndry.o\
	gkvp_colli.o\
	$(FFT).o\
	gkvp_fld.o\
	gkvp_colliimp.o\
	gkvp_freq.o\
	gkvp_zfilter.o\
	gkvp_geom.o\
	gkvp_exb.o\
	gkvp_trans.o\
	gkvp_advnc.o\
	gkvp_shearflow.o\
	gkvp_dtc.o\
	gkvp_out.o\
	gkvp_set.o\
	gkvp_main.o

main:
	(cp Makefile $(SRC); cd $(SRC); make gkvp)

gkvp: $(OBJS)
	$(FC) $(FFLAGS)	$(OBJS) -o $(PROG) $(LIB)
	mv $(PROG) ../run/

#------------------------------>
gkvp_advnc.o : gkvp_advnc.f90    gkvp_geom.o gkvp_tips.o gkvp_zfilter.o gkvp_clock.o gkvp_bndry.o gkvp_colliimp.o gkvp_colli.o gkvp_exb.o gkvp_fld.o gkvp_mpienv.o gkvp_header.o
	$(FC) $(FFLAGS) -c $<
gkvp_bndry.o : gkvp_bndry.f90    gkvp_clock.o gkvp_mpienv.o gkvp_header.o
	$(FC) $(FFLAGS) -c $<
gkvp_clock.o : gkvp_clock.f90    gkvp_mpienv.o gkvp_header.o
	$(FC) $(FFLAGS) -c $<
gkvp_colli.o : gkvp_colli.f90    gkvp_bndry.o gkvp_clock.o gkvp_mpienv.o gkvp_header.o
	$(FC) $(FFLAGS) -c $<
gkvp_colliimp.o : gkvp_colliimp.f90    gkvp_fld.o $(MATH).o gkvp_clock.o gkvp_mpienv.o gkvp_header.o
	$(FC) $(FFLAGS) -c $<
gkvp_dtc.o : gkvp_dtc.f90    gkvp_colliimp.o gkvp_exb.o gkvp_mpienv.o gkvp_header.o
	$(FC) $(FFLAGS) -c $<
gkvp_exb.o : gkvp_exb.f90    gkvp_clock.o $(FFT).o gkvp_mpienv.o gkvp_header.o
	$(FC) $(FFLAGS) -c $<
$(FFT).o : $(FFT).f90    gkvp_clock.o gkvp_mpienv.o gkvp_header.o
	$(FC) $(FFLAGS) -c $< $(INC)
$(FILEIO).o : $(FILEIO).f90    gkvp_mpienv.o gkvp_header.o
	$(FC) $(FFLAGS) -c $< $(INC)
gkvp_fld.o : gkvp_fld.f90    gkvp_clock.o gkvp_intgrl.o gkvp_mpienv.o gkvp_header.o
	$(FC) $(FFLAGS) -c $<
gkvp_freq.o : gkvp_freq.f90    gkvp_mpienv.o gkvp_header.o
	$(FC) $(FFLAGS) -c $<
gkvp_geom.o : gkvp_geom.f90    gkvp_ring.o gkvp_igs.o gkvp_vmecbzx.o gkvp_intgrl.o $(MATH).o gkvp_mpienv.o gkvp_header.o
	$(FC) $(FFLAGS) -c $<
gkvp_header.o : gkvp_header.f90
	$(FC) $(FFLAGS) -c $<
gkvp_igs.o : gkvp_igs.f90    gkvp_mpienv.o gkvp_header.o
	$(FC) $(FFLAGS) -c $<
gkvp_intgrl.o : gkvp_intgrl.f90    gkvp_mpienv.o gkvp_header.o
	$(FC) $(FFLAGS) -c $<
gkvp_main.o : gkvp_main.f90    gkvp_shearflow.o gkvp_tips.o gkvp_freq.o $(FFT).o gkvp_colliimp.o gkvp_advnc.o gkvp_fld.o gkvp_dtc.o gkvp_out.o gkvp_clock.o gkvp_set.o gkvp_mpienv.o gkvp_header.o
	$(FC) $(FFLAGS) -c $<
gkvp_mpienv.o : gkvp_mpienv.f90    gkvp_header.o
	$(FC) $(FFLAGS) -c $<
gkvp_out.o : gkvp_out.f90    $(FILEIO).o gkvp_tips.o gkvp_dtc.o gkvp_colliimp.o gkvp_advnc.o gkvp_freq.o gkvp_trans.o gkvp_fld.o gkvp_intgrl.o gkvp_mpienv.o gkvp_header.o
	$(FC) $(FFLAGS) -c $<
gkvp_ring.o : gkvp_ring.f90    $(MATH).o
	$(FC) $(FFLAGS) -c $<
gkvp_set.o : gkvp_set.f90    gkvp_geom.o $(FILEIO).o gkvp_tips.o gkvp_colliimp.o gkvp_colli.o gkvp_dtc.o gkvp_advnc.o gkvp_bndry.o gkvp_fld.o $(MATH).o gkvp_mpienv.o gkvp_header.o
	$(FC) $(FFLAGS) -c $<
gkvp_shearflow.o : gkvp_shearflow.f90    gkvp_tips.o gkvp_fld.o gkvp_mpienv.o gkvp_header.o
	$(FC) $(FFLAGS) -c $<
gkvp_tips.o : gkvp_tips.f90    gkvp_mpienv.o gkvp_header.o
	$(FC) $(FFLAGS) -c $<
gkvp_trans.o : gkvp_trans.f90    $(FILEIO).o gkvp_exb.o gkvp_clock.o gkvp_intgrl.o gkvp_mpienv.o gkvp_header.o
	$(FC) $(FFLAGS) -c $<
gkvp_vmecbzx.o : gkvp_vmecbzx.f90    gkvp_mpienv.o gkvp_header.o
	$(FC) $(FFLAGS) -c $<
gkvp_vmecin.o : gkvp_vmecin.f90    gkvp_mpienv.o gkvp_header.o
	$(FC) $(FFLAGS) -c $<
gkvp_zfilter.o : gkvp_zfilter.f90    gkvp_clock.o gkvp_mpienv.o gkvp_header.o
	$(FC) $(FFLAGS) -c $<
$(MATH).o : $(MYL)$(MATH).f90 $(MYL)Bessel0_Zeros.f90    gkvp_header.o
	$(FC) $(FFLAGS) -c $<
#------------------------------<

clean:
	rm -f ../src/Makefile ../src/*.o ../src/*.mod ../src/*.$(OPTRPT) ./*.exe ./sub.q.*.o* \
	      ./*.o ./*.mod ./*.$(OPTRPT) ./*namelist.* ./sub.q.*

clear:
	rm -f ./*.o ./*.mod ./*.$(OPTRPT) ./*namelist.* ./sub.q.*
 
