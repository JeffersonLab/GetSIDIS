#=================================================================================
#
#  File: Makefile
#  Author: Thomas Ullrich (thomas.ullrich@bnl.gov)
#  Last modified: June 16, 2009
#
#  Compile and link eps09query, a command line tool to EPS09 
#
#=================================================================================
CCX = g++
CCXOPT = -O -W -Wall

EPPS16.o: EPPS16.cpp EPPS16.h
	$(CCX) $(CCXOPT) -c EPPS16.cpp

clean:
	rm -rf *.o 

