#----------------------------------------------------------------------
# NOTES: 
#-ROOTSYS variable should be properly set into directory of ROOT
#-ROOTINC variable shows include directory form ROOT instalation, usually ($ROOTSYS)/include
#----------------------------------------------------------------------
#	make      	- compile TARGET       
#	make run	- compile and run
#	make clean	- clean project
#----------------------------------------------------------------------
# Parallel compilation
#MAKEFLAGS += -j 4  #fixed number of cores
MAKEFLAGS+="-j -l $(shell grep -c ^processor /proc/cpuinfo) "


#----------------------------------------------------------------------
# ROOT:
# -flags
 ROOTCFLAGS      := $(shell root-config --cflags)                                
 ROOTGLIBS       := $(shell root-config --glibs)                                 
 GSLLDFLAGS      := $(shell gsl-config --libs)
 ROOT             = $(ROOTCFLAGS) $(ROOTGLIBS) -lEG -lFoam 

ROOTINC          = $(ROOTSYS)/include

#----------------------------------------------------------------------
CC := g++ 
SRCDIR := ./src
BUILDDIR := ./build
TARGET := ./bin/main.x
 

SRCEXT := cxx
SOURCES := $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))
CFLAGS := -g -ansi -Wall 
LIB := $(ROOT) $(GSLLDFLAGS)
INC := -I include -I$(ROOTINC) -I.
DOC := ./doc

#----------------------------------------------------------------------
#additional checks of changes introduced by user in cuts
ACCEPTANCE := ACCEPTANCE_CUTS/ExplorationCuts.cxx ACCEPTANCE_CUTS/GenerationCuts.cxx
#----------------------------------------------------------------------



$(TARGET): $(OBJECTS)
	@echo " Linking..."
	@echo " $(CC) $^ -o $(TARGET) $(LIB)"; $(CC) $^ -o $(TARGET) $(LIB) 

$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT)  $(ACCEPTANCE)
	@mkdir -p $(BUILDDIR)
	@echo " $(CC) $(CFLAGS) $(INC) -c -o $@ $<"; $(CC) $(CFLAGS) $(INC) $(LIB) -c -o $@ $<

clean:
	@echo " Cleaning..."; 
	@echo " $(RM) -r $(BUILDDIR) $(TARGET)"; $(RM) -r $(BUILDDIR) $(TARGET) Doxyfile*

run: $(TARGET)
	(time $(TARGET))
	
all: $(TARGET)
	@echo 'Created executable: $(TARGET)'  

.PHONY: clean
.PHONY: all
.PHONY: run

#----------------------------------------------------------------------
#Generate documentation
Generate-doc:	
	-doxygen -g 
	-echo 'PROJECT_NAME = "GenEx"' >> Doxyfile
	-echo 'RECURSIVE = YES' >> Doxyfile
	-echo 'INPUT = src doc include'
	-doxygen Doxyfile
	-firefox ./html/index.html

#----------------------------------------------------------------------
#Create configuration file from log
sort_config_log:
	(./sortConfigLog.py)
#----------------------------------------------------------------------
#Run tests
test:
	(./test.py)

#----------------------------------------------------------------------
#Technical checks using valgrind
#----------------------------------------------------------------------
#valgrin suppresion file
VALGRINDSUPP = $(ROOTSYS)/etc/valgrind-root.supp
#----------------------------------------------------------------------
check-memleak: $(TARGET)
	valgrind --leak-check=yes --show-reachable=yes --suppressions=$(VALGRINDSUPP) --log-file=valgrin.log $(TARGET)
#----------------------------------------------------------------------
check-memleakfull: $(TARGET)
	 valgrind --leak-check=full --show-reachable=yes --suppressions=$(VALGRINDSUPP) --log-file=valgrin.log $(TARGET)
#----------------------------------------------------------------------
check-prof: $(TARGET) 
	 valgrind --tool=cachegrind --suppressions=$(VALGRINDSUPP) --log-file=valgrin.log $(TARGET)
	 cg_annotate cachegrind.out.*
#----------------------------------------------------------------------
check-call: $(TARGET)
	valgrind --tool=callgrind --suppressions=$(VALGRINDSUPP) --log-file=valgrin.log $(TARGET)
	callgrind_annotate callgrind.out.*
#----------------------------------------------------------------------
check-alloc: $(TARGET)
	valgrind  --tool=exp-dhat --suppressions=$(VALGRINDSUPP) --log-file=valgrin.log $(TARGET)
#----------------------------------------------------------------------
check-heap: $(TARGET)
	valgrind  --tool=massif --suppressions=$(VALGRINDSUPP)--log-file=valgrin.log $(TARGET)
	ms_print massif.out.*
#----------------------------------------------------------------------
line_count:
	(find . \( -name '*.h' -o -name '*.cxx' \) | xargs wc -l)
#======================================================================
