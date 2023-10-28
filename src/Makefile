PETSC_DIR=/home/petsc/
PETSC_ARCH=arch-linux-c-debug

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
CFLAGS = -Wall -Wwrite-strings -Wno-unknown-pragmas -Wno-lto-type-mismatch -Wno-stringop-overflow -fstack-protector -fvisibility=hidden -g3 -O0
# CC = ${PETSC_DIR}/${PETSC_ARCH}/bin/mpicc
CC = /usr/bin/mpicc
# LINKLIBRARIES = -Wl,-rpath,${PETSC_DIR}/${PETSC_ARCH}/lib -L${PETSC_DIR}/${PETSC_ARCH}/lib -Wl,-rpath,${PETSC_DIR}/${PETSC_ARCH}/lib -L${PETSC_DIR}/${PETSC_ARCH}/lib -Wl,-rpath,/usr/lib/gcc/x86_64-linux-gnu/11 -L/usr/lib/gcc/x86_64-linux-gnu/11 -lpetsc -llapack -lblas -lm -lstdc++ -ldl -lmpifort -lmpi -lgfortran -lm -lgfortran -lm -lgcc_s -lquadmath -lstdc++ -ldl
LINKLIBRARIES = -Wl,-rpath,/home/petsc//arch-linux-c-debug//lib -L/home/petsc//arch-linux-c-debug//lib -Wl,-rpath,/usr/lib/x86_64-linux-gnu/openmpi/lib/fortran/gfortran -L/usr/lib/x86_64-linux-gnu/openmpi/lib/fortran/gfortran -Wl,-rpath,/usr/lib/gcc/x86_64-linux-gnu/11 -L/usr/lib/gcc/x86_64-linux-gnu/11 -lpetsc -llapack -lblas -lm -lstdc++ -ldl -lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh -lmpi -lopen-rte -lopen-pal -lhwloc -levent_core -levent_pthreads -lgfortran -lm -lz -lgfortran -lm -lgfortran -lgcc_s -lquadmath -lstdc++ -ldl

my_funcs.so: my_funcs.cpp
	g++ -O3 -shared -fPIC -o my_funcs.so my_funcs.cpp

my_petsc.o: my_petsc.c my_petsc.h
	${CC} -o my_petsc.o \
		-c ${CFLAGS} -I${PETSC_DIR}/include -I${PETSC_DIR}/${PETSC_ARCH}/include    \
		my_petsc.c

isotherm.o: isotherm.c isotherm.h
	${CC} -o isotherm.o \
		-c ${CFLAGS} -I${PETSC_DIR}/include -I${PETSC_DIR}/${PETSC_ARCH}/include    \
		isotherm.c

run_petsc.exe: run_petsc.c my_petsc.o my_funcs.so isotherm.o
	${CC} -o run_petsc.o \
		-c ${CFLAGS} -I${PETSC_DIR}/include -I${PETSC_DIR}/${PETSC_ARCH}/include    \
		run_petsc.c 
	${CC} ${CFLAGS} \
		-o run_petsc.exe run_petsc.o my_petsc.o my_funcs.so isotherm.o  ${LINKLIBRARIES}
	rm run_petsc.o

run_s_meshes.exe: run_s_meshes.c my_petsc.o my_funcs.so isotherm.o
	${CC} -o run_s_meshes.o \
		-c ${CFLAGS} -I${PETSC_DIR}/include -I${PETSC_DIR}/${PETSC_ARCH}/include    \
		run_s_meshes.c 
	${CC} ${CFLAGS} \
		-o run_s_meshes.exe run_s_meshes.o my_petsc.o my_funcs.so isotherm.o ${LINKLIBRARIES}

run_s_timesteps.exe: run_s_timesteps.c my_petsc.o my_funcs.so isotherm.o
	${CC} -o run_s_timesteps.o \
		-c ${CFLAGS} -I${PETSC_DIR}/include -I${PETSC_DIR}/${PETSC_ARCH}/include    \
		run_s_timesteps.c 
	${CC} ${CFLAGS} \
		-o run_s_timesteps.exe run_s_timesteps.o my_petsc.o my_funcs.so isotherm.o ${LINKLIBRARIES}