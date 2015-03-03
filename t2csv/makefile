##########################################################
##       Makefile to compile t2csv.f90		##
##       Written by Oscar Barragán Oct-23-2014	        ##
##########################################################

FC=gfortran 
MPIFC=mpif90
SRC=../sources_f90/
FLAGS= -cpp -finit-local-zero -fno-automatic -ffloat-store
SOURCES=t2csv.f90\
        trans_coords.f90\
        create_names.f90\
        input_file.f90
OBJS=$(SOURCES:.f90=.o)
EXECUTABLE=t2csv

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJS)
	${MPIFC} ${FLAGS} $(OBJS) -o t2csv

t2csv.o: $(SRC)t2csv.f90
	${MPIFC} ${FLAGS} -c $(SRC)t2csv.f90

%.o: $(SRC)%.f90
	${MPIFC} ${FLAGS} -c $<


clean:
	rm -rf *.o $(EXECUTABLE)
