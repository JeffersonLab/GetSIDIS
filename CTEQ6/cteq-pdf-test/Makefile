#!/bin/bash

#   The directory where the include files are.
#   By default it was installed into your home directory.
prefix = /v/lustre2/expphy/work/halla/solid/yez/EIC_Gluon/generator/CTEQ6/cteq-pdf-1.0.4

#   include and lib directories
includedir = $(prefix)/include/cteq/
libdir = $(prefix)/lib

#   use you favorite fortran and C compilers 
#   (note, the package was tested with GNU compilers)
F77 = gfortran
CC = gcc


#   targets (you shouldn't modify this part)
all: test1 test2 test3 test4 test5 test6


test1: Cteq6Pdf-2008.f test1.c
	$(CC) -c -O3 -Wall -I$(includedir) test1.c
	$(F77) -O3 Cteq6Pdf-2008.f test1.o -L$(libdir) -o test1 -lcteqpdf

test2: test2.c
	$(CC) -O3 -Wall -I$(includedir) test2.c -o test2 -L$(libdir) -lcteqpdf

test3: test3.f
	$(F77) -O3 test3.f -L$(libdir) -o test3 -lcteqpdf

test4: test4.f
	$(F77) -O3 test4.f -L$(libdir) -o test4 -lcteqpdf

test5: test5.f Cteq6Pdf-2008.f
	$(F77) -O3 test5.f Cteq6Pdf-2008.f -L$(libdir) -o test5 -lcteqpdf

test6: test6.c
	$(CC) -O3 -Wall -I$(includedir) test6.c -o test6 -L$(libdir) -lcteqpdf


clean:
	rm -rf *.o test1 test2 test3 test4 test5 test6


