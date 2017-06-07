#
# KD_TREE2
#      A Fortran 95 module for finding close neighbors in data sets of
#      points in k-dimensional Euclidean space. 
#



# 
# For Intel Fortran 8.1
#
#FLAGSMACHINE= -xK -arch pn4 -rounding_mode=chopped  # Specific to your machine architecture
#FLAGSALWAYS = -static-libcxa -u -warn all  # for both optimization and debugging
FLAGSALWAYS = -u -warn all  # for both optimization and debugging
FLAGSMACHINE=
#FLAGSOPT=-O3 -ipo -fpe0 ${FLAGSMACHINE} -fno-alias 
FLAGSOPT=-O2 -g -traceback -fpe0 ${FLAGSMACHINE} -fno-alias
FLAGSDEBUG= -check all -traceback  #-g #-traceback # -C
F90=ifort

INCLUDE=-I${NETCDF_ROOT}/Intel/include
LDFLAGS=-L${NETCDF_ROOT}/Intel/lib
LIBS=-lnetcdff -lnetcdf


#

# for gfortran
#
#FLAGSMACHINE=#-march=athlon-xp -m3dnow -msse -mno-ieee-fp -mfpmath=sse # Athlon XP for example
#FLAGSMACHINE=
#FLAGSALWAYS = -u -Wall
#FLAGSOPT= -O3
#FLAGSDEBUG= -g
#F90=gfortran


# 
# choose debugging or optimization
FLAGS= ${FLAGSMACHINE} ${FLAGSALWAYS} ${FLAGSOPT} ${INCLUDE}  #  change the last to ${FLAGSOPT} or ${FLAGSDEBUG}


MY_DIR=`basename ${PWD}`

all:	create_runoff_nn create_runoff_weights process_runoff create_model_coast create_model_wet create_runoff_weights_spread

create_model_wet: create_model_wet.o runoff_modules.o
	${F90} ${FLAGS} -o create_model_wet create_model_wet.o  runoff_modules.o ${LDFLAGS} ${LIBS}

create_model_coast: create_model_coast.o runoff_modules.o
	${F90} ${FLAGS} -o create_model_coast create_model_coast.o runoff_modules.o ${LDFLAGS} ${LIBS}

process_runoff: process_runoff.o
	${F90} ${FLAGS} -o process_runoff process_runoff.o  ${LDFLAGS} ${LIBS}

create_runoff_weights_spread: create_runoff_weights_spread.o kdtree2.o
	${F90} ${FLAGS} -o create_runoff_weights_spread create_runoff_weights_spread.o kdtree2.o  ${LDFLAGS} ${LIBS}

create_runoff_weights: create_runoff_weights.o kdtree2.o
	${F90} ${FLAGS} -o create_runoff_weights create_runoff_weights.o kdtree2.o  ${LDFLAGS} ${LIBS}

create_runoff_nn: create_runoff_nn.o kdtree2.o runoff_modules.o
	${F90} ${FLAGS} -o create_runoff_nn create_runoff_nn.o kdtree2.o  runoff_modules.o ${LDFLAGS} ${LIBS}


create_model_coast.o: runoff_modules.o

create_model_wet.o: runoff_modules.o

create_runoff_nn.o: kdtree2.o runoff_modules.o

create_runoff_weights.o: kdtree2.o

create_runoff_weights_spread.o: kdtree2.o

%.o :: %.f90
	${F90} ${FLAGS} -c $<

tar:
	cd ..; tar zcvf ${MY_DIR}.tgz ${MY_DIR}

clean:
	/bin/rm -f *.o *.mod runoff
