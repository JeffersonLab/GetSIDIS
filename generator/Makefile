# A Generic Makefile for compiling ROOT programs
# R. Michaels, rom@jlab.org, Aug 2001  See also README !!
# Version of this release
# 
ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs)
ROOTINC       = $(shell root-config --incdir)
LHAPDFLIBS    = $(shell lhapdf-config --cflags --ldflags)
CXXFLAGS      = -Wall -g -frtti -fexceptions -fPIC \
                   -DLINUXVERS -O 
#-std=c++11 
#-Wstatic-float-init

O=GetSIDIS

# Linux with g++
INCLUDES      = -I$(ROOTSYS)/include -I/opt/local/include -I$(ROOTINC) -I$(LHAPDF)/include
#CXX           = clang++
CXX           = g++
#LD            = clang++
LD            = g++
LDFLAGS       = 

LIBS          = $(ROOTLIBS)  -L/usr/lib/  -L$(LHAPDF)/lib  -lLHAPDF
GLIBS         = $(ROOTGLIBS) 

ALL_LIBS =  $(GLIBS) $(LIBS) $(LHAPDFLIBS)

# The following sources comprise the package of scaler classes by R. Michaels.
SRC = $(O).C

HEAD = $(SRC:.C=.h)
DEPS = $(SRC:.C=.d)
SCALER_OBJS = $(SRC:.C=.o)

# Test code executibles
PROGS = $(O)

$(O):  $(O).o $(O).C
	rm -f $@
	$(CXX) $(CXXFLAGS) $(INCLUDES)  -o $@  $(O).o $(ALL_LIBS)

         
clean:
	rm -f *.o core *~ *.d *.tar $(PROGS)

realclean:  clean
	rm -f *.d

###

.SUFFIXES:
.SUFFIXES: .c .cc .cpp .C .o .d

###Adding the EPS09 code
%.o:	%.C
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $<

%.d:	%.C
	@echo Creating dependencies for $<
	@$(SHELL) -ec '$(CXX) -MM $(CXXFLAGS) $(INCLUDES) -c $< \
		| sed '\''s/\($*\)\.o[ :]*/\1.o $@ : /g'\'' > $@; \
		[ -s $@ ] || rm -f $@'

-include $(DEPS)

