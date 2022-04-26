# A Generic Makefile for compiling ROOT programs
# R. Michaels, rom@jlab.org, Aug 2001  See also README !!
# Version of this release
# 
ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs)
CXX           = g++
LHAPDF = /home/xqian/lhapdf
CXXFLAGS      = -Wall -frtti -fno-exceptions -fPIC \
                   -DLINUXVERS -I$(ROOTSYS)/include -I$(LHAPDF)/include -O



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

$(O): $(O).o FF.o DF.o vegas_C.o rand.o CrossSection.o $(O).C
	rm -f $@
	$(CXX) $(CXXFLAGS)  -o $@ $(O).o FF.o DF.o vegas_C.o rand.o CrossSection.o  $(ALL_LIBS)

clean:
	rm -f *.o core *~ *.d *.tar $(PROGS)

realclean:  clean
	rm -f *.d

###

.SUFFIXES:
.SUFFIXES: .c .cc .cpp .C .o .d

%.o:	%.C
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $<

%.d:	%.C
	@echo Creating dependencies for $<
	@$(SHELL) -ec '$(CXX) -MM $(CXXFLAGS) -c $< \
		| sed '\''s/\($*\)\.o[ :]*/\1.o $@ : /g'\'' > $@; \
		[ -s $@ ] || rm -f $@'

-include $(DEPS)








