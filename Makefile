MARCH := $(shell root-config --arch)
CXX := $(shell root-config --cxx)

OBJ=./obj/
SRC=./src/
INC=./include/

EXTRALIBS= -ldl -lssl -lMinuit -lTreePlayer -lProof -lProofPlayer -lRooFit -lRooFitCore
ifeq ($(UNAME), Linux)
EXTRALIBS += -lTMVA -lXMLIO -lMLP -lRFIO #-lNetx
SOFLAGS = -shared
endif
ifeq ($(UNAME), Darwin)
EXTRALIBS += -lTMVA -lXMLIO -lMLP
#SOFLAGS = --dynamic --shared --dynamiclib -undefined dynamic_lookup
SOFLAGS = --shared -undefined dynamic_lookup
endif

LIBS    = $(shell root-config --libs --glibs) $(EXTRALIBS)

ROOTCFLAGS=$(shell root-config --cflags)
ROOTCFLAGSSTUPIDE=-I. #rootcint put #include "include/TrCluster.hh"
CPPFLAGS= $(ROOTCFLAGSSTUPIDE) $(ROOTCFLAGS) -I$(INC) -fstack-protector -fstack-protector-all # avoids problems like arrays used ovr the size writing elsewhere...
#CXXFLAGS= -O3 -Wno-write-strings -fPIC $(CPPFLAGS) 
CXXFLAGS= -g -Wno-write-strings -fPIC $(CPPFLAGS)

LDFLAGS = $(shell root-config --ldflags) #-Wl,-stack_size,0x10000000 # 256MB of stack (default is 8MB...), http://stackoverflow.com/questions/2275550/change-stack-size-for-a-c-application-in-linux-during-compilation-with-gnu-com, http://linuxtoosx.blogspot.it/2010/10/stack-overflow-increasing-stack-limit.html

ifdef ONMACOSX
        CXXFLAGSFORCINT=$(subst -stdlib=libc++ -std=c++11,,$(CXXFLAGS))
else
        CXXFLAGSFORCINT=$(CFLAGS)
endif

default: DataAnalysis DataAnalysisGhosts DataAnalysisMulti

DataAnalysis: $(OBJ)DataAnalysis.o $(OBJ)Dictionary.o
#	echo $(CXXFLAGS)
	$(CXX) -o $@ $(LDFLAGS) $(OBJ)Dictionary.o $(OBJ)DataAnalysis.o $(LIBS)

DataAnalysisGhosts: $(OBJ)DataAnalysisGhosts.o $(OBJ)Dictionary.o
#	echo $(CXXFLAGS)
	$(CXX) -o $@ $(LDFLAGS) $(OBJ)Dictionary.o $(OBJ)DataAnalysisGhosts.o $(LIBS)

DataAnalysisMulti: $(OBJ)DataAnalysisMulti.o $(OBJ)Dictionary.o
#	echo $(CXXFLAGS)
	$(CXX) -o $@ $(LDFLAGS) $(OBJ)Dictionary.o $(OBJ)DataAnalysisMulti.o $(LIBS)

$(OBJ)%.o: $(SRC)%.C
	@echo Compiling  $< ...
	@if ! [ -d $(OBJ) ] ; then mkdir -p $(OBJ); fi
	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(OBJ)%.o: $(SRC)%.cxx
	@echo Compiling  $< ...
	@if ! [ -d $(OBJ) ] ; then mkdir -p $(OBJ); fi
	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(OBJ)%.o: $(SRC)%.cpp
	@echo Compiling  $< ...
	@if ! [ -d $(OBJ) ] ; then mkdir -p $(OBJ); fi
	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(OBJ)Dictionary.o: $(OBJ)Dictionary.cxx
	@echo Compiling  $< ...
	@if ! [ -d $(OBJ) ] ; then mkdir -p $(OBJ); fi
	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(OBJ)Dictionary.cxx: $(INC)TrCluster.hh $(INC)TrCluster_linkdef.h
	@echo Creating  $@ ...
	@if ! [ -d $(OBJ) ] ; then mkdir -p $(OBJ); fi
	$(ROOTSYS)/bin/rootcint -f $@ -c -p $(CXXFLAGSFORCINT) $^

clean:
	rm -fv DataAnalysis
	rm -fv DataAnalysisGhosts
	rm -fv DataAnalysisMulti
	rm -fv $(OBJ)*
