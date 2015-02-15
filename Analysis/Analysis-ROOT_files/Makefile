
# ROOT:
# -flags
 ROOTCFLAGS      := $(shell root-config --cflags)                                
 ROOTGLIBS       := $(shell root-config --glibs)                                 
 GSLLDFLAGS      := $(shell gsl-config --libs)
 ROOT             = $(ROOTCFLAGS) $(ROOTGLIBS) -lEG
 # -lFoam 

ROOTINC          = $(ROOTSYS)/include





#======================================================================
#####  GNU c++ compiler
CXXFLAGS    = -g -ansi -Wall  -I. -I$(ROOTINC) 
LDFLAGS     = $(GSLLDFLAGS)
CXX         =  g++
CLD         =  g++
#-----------------------------


HEADERS = TConfigReader.h TPolicyReader.h
OBJECTS = TConfigReader.o TPolicyReader.o


# Compile all c++ files
%.o : %.cxx
	$(CXX) $(CXXFLAGS) -c $<
	

#Create ROOT dictionary
Dict.cxx: $(HEADERS) LinkDef.h
	@echo "Generating dictionary $@..."
	@rootcint -f $@ -c -p $^	
	

analysis.x: analysis.o $(OBJECTS) Dict.o
	$(CLD) -o analysis.x $^ $(ROOT) $(LDFLAGS)	
	
	
.PHONY: run clean cleanest 

run: analysis.x
	(time ./analysis.x)

clean : 
	-rm *.x
	-rm *.o

cleanest : clean



