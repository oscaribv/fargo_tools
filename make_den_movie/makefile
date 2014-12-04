##########################################################
##       Makefile to compile make_den_movie.f90		##
##       Written by Oscar Barragán Oct-16-2014	        ##
##########################################################

FC=gfortran 
MPIFC=mpif90
SRC=../sources_f90/
FLAGS= -cpp -finit-local-zero -fno-automatic -ffloat-store
SOURCES=readbinaryall.f90\
        trans_coords.f90\
        create_names.f90\
        input_file.f90
OBJS=$(SOURCES:.f90=.o)
EXECUTABLE=readbinaryall

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJS)
	${MPIFC} ${FLAGS} $(OBJS) -o readbinaryall

readbinaryall.o: $(SRC)readbinaryall.f90
	${MPIFC} ${FLAGS} -c $(SRC)readbinaryall.f90

%.o: $(SRC)%.f90
	${MPIFC} ${FLAGS} -c $<


clean:
	rm -rf *.o $(EXECUTABLE)
