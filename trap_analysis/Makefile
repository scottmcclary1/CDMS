#
# Makefile for the trap analysis program.
#

  MAKEFILE  = Makefile

#===============================================================================
# The following definitions are for the Intel compiler, releases 13 through 15.
# 
#-------------------------------------------------------------------------------
# Intel Fortran 95 compiler
 F95    = ifort
 OMPF95 = ifort -qopenmp

#-------------------------------------------------------------------------------
# Flags other than for optimization.
 F95FLAGS = -nologo -gen-interfaces -debug-parameters all \
   -warn declarations -warn interfaces -module work -traceback \
   -check uninit

#F95FLAGS = -nologo -gen-interfaces -debug-parameters all \
#  -warn declarations -warn interfaces -module work -traceback -check bounds \
#  -check uninit

#F95FLAGS = -nologo -gen-interfaces -debug-parameters all \
#  -warn declarations -warn interfaces -module work -traceback -check bounds \
#  -check uninit -heap-arrays:2

#F95FLAGS = -nologo -g -O0 -gen-interfaces -openmp -debug-parameters all \
#  -warn declarations -warn interfaces -module work -traceback -check bounds \
#  -check uninit -threads

#-------------------------------------------------------------------------------
# Optimization flags
#F95OPT = -g -O0
 F95OPT =  -O2 -ipo -no-prec-div -xHost
#===============================================================================

#-------------------------------------------------------------------------------
 FCFLAGS    = $(F95OPT) $(TARCH) $(F95FLAGS) $(ATP)

#OBJC = mkl_dfti.o omp_lib.o FFT_of_Selected_Section.o Find_Frequency_Shift.o Find_Peak_Integral.o \
#      Frequency_Data.o Harmonic_Analysis.o Mag_Cutoff_317k_Datapoints.o \
#      Mag_Cutoff_69k_Datapoints.o Noise_Analysis.o \
#      Simulations.o Trap_Analysis.o Trap_Charge_Calibration.o Voltage_Finder.o

 OBJC      = FFT_of_Selected_Section.o Find_Frequency_Shift.o Find_Peak_Integral.o \
             Frequency_Data.o Harmonic_Analysis.o Mag_Cutoff_317k_Datapoints.o \
             Mag_Cutoff_69k_Datapoints.o Noise_Analysis.o \
             Simulations.o Trap_Charge_Calibration.o Voltage_Finder.o

#OBJC      = omp_lib.o FFT_of_Selected_Section.o Find_Frequency_Shift.o Find_Peak_Integral.o \
#            Frequency_Data.o Harmonic_Analysis.o Mag_Cutoff_317k_Datapoints.o \
#            Mag_Cutoff_69k_Datapoints.o Noise_Analysis.o \
#            Simulations.o Trap_Charge_Calibration.o Voltage_Finder.o

#MODOBJ    = 
#MODOBJ    = mkl_dfti.o
 MODOBJ    = Trap_Mod.o mkl_dfti.o

#OBJL = omp_lib.o FFT_of_Selected_Section.o Find_Frequency_Shift.o Find_Peak_Integral.o \
#      Frequency_Data.o Harmonic_Analysis.o Mag_Cutoff_317k_Datapoints.o \
#      Mag_Cutoff_69k_Datapoints.o Noise_Analysis.o \
#      Simulations.o Trap_Analysis.o Trap_Charge_Calibration.o Voltage_Finder.o

#OBJL = FFT_of_Selected_Section.o Find_Frequency_Shift.o Find_Peak_Integral.o \
#      Frequency_Data.o Harmonic_Analysis.o Mag_Cutoff_317k_Datapoints.o \
#      Mag_Cutoff_69k_Datapoints.o Noise_Analysis.o \
#      Simulations.o Trap_Charge_Calibration.o Voltage_Finder.o


#-----------------------------------------------------------------------------

trap: dummy
	make -f $(MAKEFILE) target "TARGET=trap" "FC=$(F95)" \
              "FEATURES=" \
              "MODOBJ=$(MODOBJ)" \
              "OBJ=$(OBJC) NCC_Peakfinder_VI.o Trap_Analysis.o" \
              "EXTRAOBJ=" \
              "LIBS= -mkl=sequential"

trap_omp: dummy
	make -f $(MAKEFILE) target "TARGET=trap_omp" "FC=$(OMPF95)" \
              "FEATURES=" \
              "MODOBJ=$(MODOBJ) omp_lib.o" \
              "OBJ=$(OBJC) NCC_Peakfinder_VI.o Trap_Analysis.o" \
              "EXTRAOBJ=" \
              "LIBS= -mkl=sequential"

dummy:

clean:
	rm -f *.o work/*.f90 work/*.mod

realclean:
	rm -f *.o work/*.f90 work/*.mod trap trap_omp

#===============================================================================
#     Each make of a program sets several macros, and then recursively calls
#     make. The target of this recursive call is "target":

target: $(MODOBJ) $(OBJ) $(EXTRAOBJ)
	$(FC) $(FLDFLAGS) -o $(TARGET) $(OBJ) $(EXTRAOBJ) $(LIBS)
#$(FC) $(FLDFLAGS) -o $(TARGET) $(OBJ) $(MODOBJ) $(EXTRAOBJ) $(LIBS)


#===============================================================================
#     Object files


#-----------------------------------------------------------------------------
FFT_of_Selected_Section.o: FFT_of_Selected_Section.F90
	$(FC) -c $(FCFLAGS) FFT_of_Selected_Section.F90

Find_Frequency_Shift.o: Find_Frequency_Shift.F90
	$(FC) -c $(FCFLAGS) Find_Frequency_Shift.F90

Find_Peak_Integral.o: Find_Peak_Integral.F90
	$(FC) -c $(FCFLAGS) Find_Peak_Integral.F90

Frequency_Data.o: Frequency_Data.F90
	$(FC) -c $(FCFLAGS) Frequency_Data.F90

Harmonic_Analysis.o: Harmonic_Analysis.F90
	$(FC) -c $(FCFLAGS) Harmonic_Analysis.F90

Mag_Cutoff_317k_Datapoints.o: Mag_Cutoff_317k_Datapoints.F90
	$(FC) -c $(FCFLAGS) Mag_Cutoff_317k_Datapoints.F90

Mag_Cutoff_69k_Datapoints.o: Mag_Cutoff_69k_Datapoints.F90
	$(FC) -c $(FCFLAGS) Mag_Cutoff_69k_Datapoints.F90

NCC_Peakfinder_VI.o: NCC_Peakfinder_VI.F90
	$(FC) -c $(FCFLAGS) NCC_Peakfinder_VI.F90

Noise_Analysis.o: Noise_Analysis.F90
	$(FC) -c $(FCFLAGS) Noise_Analysis.F90

Simulations.o: Simulations.F90
	$(FC) -c $(FCFLAGS) Simulations.F90

Trap_Analysis.o: Trap_Analysis.f90
	$(FC) -c $(FCFLAGS) Trap_Analysis.f90

Trap_Analysis_03.o: Trap_Analysis_03.F90
	$(FC) -c $(FCFLAGS) Trap_Analysis_03.F90

Trap_Charge_Calibration.o: Trap_Charge_Calibration.f90
	$(FC) -c $(FCFLAGS) Trap_Charge_Calibration.f90

Voltage_Finder.o: Voltage_Finder.F90
	$(FC) -c $(FCFLAGS) Voltage_Finder.F90

Trap_Mod.o: Trap_Mod.F90
	$(FC) -c $(FCFLAGS) Trap_Mod.F90

mkl_dfti.o: mkl_dfti.f90
	$(FC) -c $(FCFLAGS) mkl_dfti.f90

omp_lib.o: omp_lib.f
	$(FC) -c $(FCFLAGS) omp_lib.f

#omp_lib_ser.o: omp_lib_ser.f90
#	$(FC) -c $(FCFLAGS) omp_lib_ser.f90
