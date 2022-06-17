
#--------------------------------------------------------------------------
#------- USER DEFINED -----------------------------------------------------
#--------------------------------------------------------------------------

NAME = xPOD3D

F90FILES = modfile.f90 modfile_FFT.f90 cgns_routines.f90 \
           export.f90 averaging.f90 probe.f90 profile.f90 compute_area.f90 \
           POD_Fourier.f90 POD_snapshot.f90 POD_harmonic.f90 \
           general.f90 der.f90 TKE.f90 lumley.f90 main.f90

OFILES = $(F90FILES:.f90=.o)

export hostname
HOST = $(shell hostname)

ifeq "$(HOST)" "Nazgul"
	CYCLE = 3
	LIBS_CGNS = -I/usr/local/include /usr/local/lib/libcgns.a
endif

ifeq "$(HOST)" "gabnog"
	CYCLE = 4
	LIBS_CGNS = -I/usr/local/include /usr/local/lib/libcgns.a
endif

ifeq "$(HOST)" "ice"
	CYCLE = 1
	LIBS_CGNS = -I../../libs/include -L../../libs/ -lcgns
endif

#CYCLE = 3

#LIBS_CFF = -L $(CFDROOT)/$(SYSTEM)/$(SYSTEM_CPU)/lib/ -lcfd-intel10 -L $(CFDROOT)/$(SYSTEM)/$(SYSTEM_CPU)/lib/ -ladf-intel10 -L /opt/intel/lib/intel64 -lirc
#LIBS_FFT = -I./newfft/ -L./newfft -lnewfft

LIBS = $(LIBS_CFF) $(LIBS_CGNS) $(LIBS_FFT)

ifeq "$(CYCLE)" "0"
FC = ifort
OPTLEVEL = -O0
FLAGS1 = -g -traceback -ftrapuv -debug all -warn all -implicitnone
FLAGS2 = -check bounds -check uninit -check format -check output_conversion -check power -check args
FLAGS3 = -fp-speculation=off -heap-arrays -openmp 
#FLAGS4 = -convert big_endian
FLAGS = $(OPTLEVEL) $(FLAGS1) $(FLAGS2) $(FLAGS3) #$(FLAGS4)
endif

ifeq "$(CYCLE)" "1"
FC = ifort
OPTLEVEL = -O2
FLAGS1 = -heap-arrays -openmp -static
FLAGS = $(OPTLEVEL) $(FLAGS1)
endif

ifeq "$(CYCLE)" "2"
FC = ifort
OPTLEVEL = -O2
FLAGS1 = -heap-arrays -openmp
FLAGS2 = #-static 
FLAGS = $(OPTLEVEL) $(FLAGS1) $(FLAGS2)
endif

ifeq "$(CYCLE)" "3"
FC = gfortran 
OPTLEVEL = -O2
FLAGS1 = -fopenmp #-static #-fautomatic
FLAGS2 = #-fconvert=big-endian
FLAGS3 = -ffree-line-length-none
FLAGS = $(OPTLEVEL) $(FLAGS1) $(FLAGS2) $(FLAGS3)
endif

ifeq "$(CYCLE)" "4"
FC = gfortran 
OPTLEVEL = -O0
FLAGS1 = -fcheck=all -fimplicit-none -Wall -fbacktrace -static
FLAGS2 = #-fopenmp
FLAGS3 = -ffree-line-length-none
FLAGS4 = #-fconvert=big-endian
FLAGS = $(OPTLEVEL) $(FLAGS1) $(FLAGS2) $(FLAGS3) $(FLAGS4)
endif

#--------------------------------------------------------------------------
#------- COMPILATION RULES ------------------------------------------------
#--------------------------------------------------------------------------

%.o : %.f90
	$(FC) $(FLAGS) -c $*.f90 $(LIBS)

$(NAME) : $(OFILES)
	$(FC) $(FLAGS) $(OFILES) -o $(NAME) $(LIBS)

print : 
	echo $(HOST)

remove : 
	clear; rm -f ?akefile~ *.o *.mod *mod.f90 *.inp *.inp~ *.out *.out~ *.f90~ *.dat~ *.out~ fort* $(NAME)

clean :
	clear; rm -f ?akefile~ *.o *.mod *mod.f90 *.inp *.inp~ *.out *.out~ *.f90~ *.dat~ *.out~ fort*

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------

