# Get all make options
CMDVARS := $(strip $(foreach V,$(.VARIABLES),$(if $(findstring command,$(origin $V)),$V)))

# Test if invalid make options exist
ifneq ($(and $(CMDVARS),$(filter-out $(PARAMS),$(CMDVARS))),)
	$(error Make option $(filter-out $(PARAMS),$(CMDVARS)) is not supported)
endif

CC = gcc
CFLAGS = -g -O0 -Wall -Wextra

# Check CMake version
CMAKE_VERS := $(shell cmake --version 2> /dev/null |awk '{print $$3}')
ifeq ($(CMAKE_VERS),)
	CMAKE_VERS := 0.0.0
endif
CMAKE_VERS_RQD := 3.1.3

# Check if CMake version meets requirement
ifeq ($(shell printf '%s\n' $(CMAKE_VERS) $(CMAKE_VERS_RQD) | sort -V | head -n 1), $(CMAKE_VERS_RQD))
	CMAKE_EXIST := true
	CMAKE := cmake
else
	CMAKE_EXIST := false
	ifeq ($(shell uname), Darwin)
		CMAKE_VERS_INSTALL := cmake-3.7.2-Darwin-x86_64
		CMAKE := $(PWD)/$(CMAKE_VERS_INSTALL)/CMake.app/Contents/bin/cmake
	else
		CMAKE_VERS_INSTALL := cmake-3.7.2-Linux-x86_64
		CMAKE := $(PWD)/$(CMAKE_VERS_INSTALL)/bin/cmake
	endif
endif

# Define CVODE paths
CVODE_PATH = ./cvode/instdir
CVODE_LIB = $(CVODE_PATH)/lib

# Define source directory
SRCDIR = ./src
LIBS = -lm
INCLUDES = \
	-I$(SRCDIR)/include\
	-I$(CVODE_PATH)/include

# Define libraries
LFLAGS = -lsundials_cvode -L$(CVODE_LIB) -lsundials_nvecserial

SRCS_ = \
	main.c \
	custom_io.c \
	cvode_funcs.c \
	initialize.c \
	ode.c \
	read_hydro.c \
	read_input.c \
	soil.c \
	swc.c

HEADERS_ = \
	custom_io.h \
	cycles.h \
	cycles_const.h \
	cycles_func.h \
	cycles_struct.h

EXECUTABLE = swc

SRCS = $(patsubst %,$(SRCDIR)/%,$(SRCS_))
HEADERS = $(patsubst %,$(SRCDIR)/include/%,$(HEADERS_))
OBJS = $(SRCS:.c=.o)

.PHONY: all clean cvode cmake

cycles: $(OBJS)
	@$(CC) $(CFLAGS) $(SFLAGS) $(INCLUDES) -o $(EXECUTABLE) $(OBJS) $(MODULE_OBJS) $(LFLAGS) $(LIBS)

all: cvode cycles

cmake:
ifeq ($(CMAKE_EXIST), false)
	@echo "CVODE installation requires CMake v$(CMAKE_VERS_RQD) or above."
	@echo "Download CMake $(CMAKE_VERS_INSTALL) from cmake.org"
	@curl https://cmake.org/files/v3.7/$(CMAKE_VERS_INSTALL).tar.gz -o $(CMAKE_VERS_INSTALL).tar.gz &> /dev/null
	@echo
	@echo "Extract $(CMAKE_VERS_INSTALL).tar.gz"
	@tar xzf $(CMAKE_VERS_INSTALL).tar.gz
endif

cvode: cmake
	@echo "Install CVODE library"
	@cd cvode && mkdir -p instdir && mkdir -p builddir
	@cd $(CVODE_PATH) && $(CMAKE) -DCMAKE_INSTALL_PREFIX=../instdir -DCMAKE_INSTALL_LIBDIR=lib -DBUILD_SHARED_LIBS=OFF -DEXAMPLES_ENABLE_C=OFF -DEXAMPLES_INSTALL=OFF ../
	@cd $(CVODE_PATH) && make && make install
	@echo "CVODE library installed."
ifeq ($(CMAKE_EXIST), false)
	@echo "Remove CMake files"
	@$(RM) -r $(CMAKE_VERS_INSTALL).tar.gz $(CMAKE_VERS_INSTALL)
endif

%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) $(SFLAGS) $(INCLUDES) -c $< -o $@

clean:		## Clean executables and objects
	@$(RM) $(SRCDIR)/*.o $(SRCDIR)/*/*.o *~
