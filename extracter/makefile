##########################################################
##       Makefile to compile extracter.f90	##
##       Written by Oscar Barragán Oct-16-2014	        ##
##########################################################

FC=gfortran 
MPIFC=mpif90
SRC=../sources_f90/
FLAGS= -cpp -finit-local-zero -fno-automatic -ffloat-store
SOURCES=extracter.f90\
        trans_coords.f90\
        create_names.f90\
        input_file.f90
OBJS=$(SOURCES:.f90=.o)
EXECUTABLE=extracter

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJS)
	${MPIFC} ${FLAGS} $(OBJS) -o extracter

extracter.o: $(SRC)extracter.f90
	${MPIFC} ${FLAGS} -c $(SRC)extracter.f90

%.o: $(SRC)%.f90
	${MPIFC} ${FLAGS} -c $<


clean:
	rm -rf *.o $(EXECUTABLE)
