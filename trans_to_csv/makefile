##########################################################
##       Makefile to compile trans_to_csv.f90		##
##       Written by Oscar Barragán Oct-23-2014	        ##
##########################################################

FC=gfortran 
MPIFC=mpif90
SRC=../sources_f90/
FLAGS= -cpp -finit-local-zero -fno-automatic -ffloat-store
SOURCES=trans_to_csv.f90\
        trans_coords.f90\
        create_names.f90\
        input_file.f90
OBJS=$(SOURCES:.f90=.o)
EXECUTABLE=trans_to_csv

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJS)
	${MPIFC} ${FLAGS} $(OBJS) -o trans_to_csv

trans_to_csv.o: $(SRC)trans_to_csv.f90
	${MPIFC} ${FLAGS} -c $(SRC)trans_to_csv.f90

%.o: $(SRC)%.f90
	${MPIFC} ${FLAGS} -c $<


clean:
	rm -rf *.o $(EXECUTABLE)
