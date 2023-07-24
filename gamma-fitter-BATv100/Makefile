# Makefile for Gamma Analysis

# Compiler and flags
CXX = g++
CXXFLAGS  = -std=c++11 -g -O2 -Wall -fPIC -Wno-deprecated  -fopenmp  
LD        = /usr/bin/ld -m elf_x86_64 
LDFLAGS   = -g -O2   -fopenmp  

# ----------------------------------------------------------------------
# The following definitions rely on the script bat-config being
# available in $PATH. If BAT is not installed in the standard system
# directories or you need to specify custom paths, uncomment the lines below:

# Paths
BAT_INCLUDE = /home/calgaro/bat_install/include # or /opt/bat/include in the container
BAT_LIB = /home/calgaro/bat_install/lib # or /opt/bat/lib in the container
CXXFLAGS += -I$(BAT_INCLUDE)
LIBS += -L$(BAT_LIB) -lBAT
CXXFLAGS += $(shell bat-config --cflags)
LIBS := $(shell bat-config --libs)
CXXFLAGS += -I/lfs/l1/legend/software/root/v06.14/include # or -I/opt/root/include in the container

# Source files
CXXSRCS = GammaLineFit.cxx
PRGSRCS = runGammaAnalysis.cxx

# Object files
CXXOBJS = $(CXXSRCS:.cxx=.o)
PRGOBJS = $(PRGSRCS:.cxx=.o)

# Executable name
EXECUTABLE = runGammaAnalysis

# ----------------------------------------------------------------------
# don't change lines below unless you know what you're doing
#

CXXOBJS = $(addsuffix .o,$(basename $(CXXSRCS)))
MYPROGS = $(basename $(PRGSRCS))
PRGOBJS = $(addsuffix .o,$(basename $(PRGSRCS)))

GARBAGE = $(CXXOBJS) $(PRGOBJS) link.d $(MYPROGS)

# targets
all : $(MYPROGS)

.PHONY : all clean print

# Your existing Makefile rules and variables

# Modify the link.d rule
link.d: $(addsuffix .h,$(basename $(CXXSRCS))) $(CXXSRCS) $(PRGSRCS)
	$(CXX) -M $(CXXFLAGS) $(filter-out %.h,$^) > link.d
	@$(foreach prog,$(MYPROGS), echo $(prog) : $(prog).o >> link.d;)

-include link.d


$(CXXOBJS) $(PRGOBJS) : %.o : %.cxx
	$(CXX) $(CXXFLAGS) -c $< -o $@


$(MYPROGS) : $(CXXOBJS) $(PRGOBJS)
	$(CXX) $(LDFLAGS) $^ -L$(BAT_LIB) $(LIBS) -o $@

clean :
	rm -f $(GARBAGE)

print :
	@echo compiler  : $(CXX)
	@echo c++ srcs  : $(CXXSRCS) $(PRGSRCS)
	@echo c++ objs  : $(CXXOBJS) $(PRGOBJS)
	@echo c++ flags : $(CXXFLAGS)
	@echo ld flags  : $(LDFLAGS)
	@echo libs      : $(LIBS)
