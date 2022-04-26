# A Generic Makefile for compiling ROOT programs
# R. Michaels, rom@jlab.org, Aug 2001  See also README !!
# Version of this release
# 
ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs)
ROOTINC       = $(shell root-config --incdir)
CXX           = g++
CXXFLAGS      = -Wall -g -frtti -fno-exceptions -fPIC \
                   -DLINUXVERS -I$(ROOTINC) -I$(LHAPDF)/include -O

#O=makebin_both
#O=makebin
#O=skim
O=sole_inter

# Linux with g++
INCLUDES      = -I$(ROOTSYS)/include 
CXX           = g++
LD            = g++
LDFLAGS       = 

LIBS          = $(ROOTLIBS) -L$(LHAPDF)/lib -lLHAPDF
GLIBS         = $(ROOTGLIBS) -L/usr/X11R6/lib -lXpm -lX11

ALL_LIBS =  $(GLIBS) $(LIBS)

# The following sources comprise the package of scaler classes by R. Michaels.
SRC = $(O).C

HEAD = $(SRC:.C=.h)
DEPS = $(SRC:.C=.d)
SCALER_OBJS = $(SRC:.C=.o)

# Test code executibles
PROGS = $(O)

$(O): $(O).o $(O).C
	rm -f $@
	$(CXX) $(CXXFLAGS) -o $@ $(O).o $(ALL_LIBS)

         
clean:
	rm -f *.o core *~ *.d *.tar $(PROGS)

realclean:  clean
	rm -f *.d

###

.SUFFIXES:
.SUFFIXES: .c .cc .cpp .C .o .d

%.o:	%.C
	$(CXX) $(CXXFLAGS) -c $<

%.d:	%.C
	@echo Creating dependencies for $<
	@$(SHELL) -ec '$(CXX) -MM $(CXXFLAGS) -c $< \
		| sed '\''s/\($*\)\.o[ :]*/\1.o $@ : /g'\'' > $@; \
		[ -s $@ ] || rm -f $@'

-include $(DEPS)








