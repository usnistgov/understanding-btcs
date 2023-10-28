out_dir := out
PETSC_DIR=/home/petsc/
PETSC_ARCH=arch-linux-c-debug

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
CFLAGS = -Wall -Wwrite-strings -Wno-unknown-pragmas -Wno-lto-type-mismatch -Wno-stringop-overflow -fstack-protector -fvisibility=hidden -g3 -O0
# CC = ${PETSC_DIR}/lib/bin/mpic++


run_petsc.exe: run_petsc.o
	${PETSC_DIR}/${PETSC_ARCH}/bin/mpicc -o run_petsc.o \
		-c ${CFLAGS} -I${PETSC_DIR}/include -I${PETSC_DIR}/${PETSC_ARCH}/include    \
		run_petsc.c
	${PETSC_DIR}/${PETSC_ARCH}/bin/mpicc ${CFLAGS} \
		-o run_petsc.exe run_petsc.o src/my_funcs.so -Wl,-rpath,${PETSC_DIR}/${PETSC_ARCH}/lib \
		-L${PETSC_DIR}/${PETSC_ARCH}/lib -Wl,-rpath,${PETSC_DIR}/${PETSC_ARCH}/lib -L${PETSC_DIR}/${PETSC_ARCH}/lib -Wl,-rpath,/usr/lib/gcc/x86_64-linux-gnu/11 \
		-L/usr/lib/gcc/x86_64-linux-gnu/11 -lpetsc -llapack -lblas -lm -lstdc++ -ldl -lmpifort -lmpi -lgfortran -lm -lgfortran -lm -lgcc_s -lquadmath -lstdc++ -ldl

run_petsc.so: run_petsc.o
	${PETSC_DIR}/${PETSC_ARCH}/bin/mpicc -o run_petsc.o \
		-c ${CFLAGS} -I${PETSC_DIR}/include -I${PETSC_DIR}/${PETSC_ARCH}/include    \
		run_petsc.c
	${PETSC_DIR}/${PETSC_ARCH}/bin/mpicc ${CFLAGS} \
		-shared -fPIC -o run_petsc.so run_petsc.o src/my_funcs.so -Wl,-rpath,${PETSC_DIR}/${PETSC_ARCH}/lib \
		-L${PETSC_DIR}/${PETSC_ARCH}/lib -Wl,-rpath,${PETSC_DIR}/${PETSC_ARCH}/lib -L${PETSC_DIR}/${PETSC_ARCH}/lib -Wl,-rpath,/usr/lib/gcc/x86_64-linux-gnu/11 \
		-L/usr/lib/gcc/x86_64-linux-gnu/11 -lpetsc -llapack -lblas -lm -lstdc++ -ldl -lmpifort -lmpi -lgfortran -lm -lgfortran -lm -lgcc_s -lquadmath -lstdc++ -ldl