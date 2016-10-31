#-------------------------------------------
# Makefile for any C program with ROOT lib
# Zhihong Ye, 01/09/2009
#-------------------------------------------

PROGRAM = makebin

CFLAGS = $(shell $(ROOTSYS)/bin/root-config --cflags)
LIBS   = $(shell $(ROOTSYS)/bin/root-config --libs)
GLIBS  = $(shell $(ROOTSYS)/bin/root-config --glibs)

CFLAGS += -Wall
INCLUDE := -I${ROOTSYS}/include
RTLIBS := -L${ROOTSYS}/lib -lCore
all: $(PROGRAM)

FFLAGS = -C -g -m64 -fbounds-check -finit-local-zero -fno-automatic


$(PROGRAM): $(PROGRAM).o 
	clang++ -o $(PROGRAM) $(PROGRAM).o $(CFLAGS) $(RTLIBS) $(LIBS) $(GLIBS)


$(PROGRAM).o: $(PROGRAM).C
	clang++ $(CFLAGS) -c -o $(PROGRAM).o $(PROGRAM).C  ${INCLUDE}

clean:
#	rm -f *.o ${GFORTRAN_OBJ}
	rm -f *.o
