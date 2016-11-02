CPP=g++
CPPFLAGS=-Wall -std=c++11
MKDIR=mkdir -p

SRCPATH=src
LIBPATH=lib
BINPATH=bin
INPUTPATH=input
OUTPUTPATH=output
BUILD_DIRS=$(LIBPATH) $(BINPATH) $(INPUTPATH) $(OUTPUTPATH)

.DEFAULT_GOAL := all

# VPATH specifies the search path
VPATH=$(SRCPATH)

lib%.so: %.cpp
	$(CPP) $(CPPFLAGS) -shared -fPIC -o $(LIBPATH)/$@ $<

main.out: main.cpp
	$(CPP) $(CPPFLAGS) -L$(LIBPATH) -lfftw3f -lfieldio -o $(BINPATH)/$@ $<

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
libs: libfieldio.so

.PHONY: all
all: | dirs libs main
