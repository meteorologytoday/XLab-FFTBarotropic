CPP=g++
CPPFLAGS=-Wall -std=c++11
MKDIR=mkdir -p

SRCPATH=src
LIBPATH=lib
BINPATH=bin
INPUTPATH=input
OUTPUTPATH=output
BUILD_DIRS=$(LIBPATH) $(BINPATH) $(INPUTPATH) $(OUTPUTPATH)

EX_LIBS=fftw3f
IN_LIBS=fieldio

CPP_INC_LIBS=$(foreach lib,$(EX_LIBS),-l$(lib)) $(foreach lib,$(IN_LIBS),-l$(lib))

.DEFAULT_GOAL := all

# VPATH specifies the search path
VPATH=$(SRCPATH)

lib%.so: %.cpp
	$(CPP) $(CPPFLAGS) -shared -fPIC -o $(LIBPATH)/$@ $<

%.out: %.cpp
	$(CPP) $(CPPFLAGS) -L$(LIBPATH) -lfieldio -o $(BINPATH)/$@ $<

main.out: main.cpp
	$(CPP) $(CPPFLAGS) -L$(LIBPATH) $(CPP_INC_LIBS) -o $(BINPATH)/$@ $<

$(BUILD_DIRS):
	$(MKDIR) $@

.PHONY: dirs
dirs: $(BUILD_DIRS)

.PHONY: clean
clean:
	for dir in $(BUILD_DIRS); do \
		rm -rf $$dir; \
	done

.PHONY: main
main: main.out

.PHONY: libs
libs: libfieldio.so libfftwfop.so

.PHONY: makefield
makefield: makefield.out


.PHONY: all
all: | dirs libs main makefield
