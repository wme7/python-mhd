


ConfigHeader = \
"""
#ifndef __MARA_CONFIG_HEADER__
#define __MARA_CONFIG_HEADER__

#define __MARA_BASE_VERSION "0.9"
#define __MARA_HG_CHANGESET "0000000"

#define __MARA_INSTALL_DIR "%(install_dir)s"
#define __MARA_USE_MPI %(mpi)d
#define __MARA_USE_HDF5 %(use_hdf5)d
#define __MARA_USE_HDF5_PAR %(use_hdf5_par)d
#define __MARA_USE_FFTW %(use_fftw)d
#define __MARA_USE_GLFW %(use_glfw)d

#ifdef __INTEL_COMPILER
#define Mara_isinf_cxx(x) isinf(x)
#define Mara_isnan_cxx(x) isnan(x)
#else
#define Mara_isinf_cxx(x) std::isinf(x)
#define Mara_isnan_cxx(x) std::isnan(x)
#endif // __INTEL_COMPILER

#endif // __MARA_CONFIG_HEADER__
"""


SystemMakefile = \
"""
# ------------------------------------------------------------------------------
#
#                        Mara Astrophysical gasdynamics code
#
#                          System-specific build macros
#
# ------------------------------------------------------------------------------

HDF5_HOME = %(hdf5)s
FFTW_HOME = %(fftw)s

USE_GLFW = %(use_glfw)s
USE_FFTW = %(use_fftw)s

AR = ar rcu

CC = %(cc)s
CXX = %(cxx)s

CFLAGS = %(cflags)s
CLIBS = %(clibs)s

HDF5_L = -L$(HDF5_HOME)/lib %(hdf5libs)s
HDF5_I = -I$(HDF5_HOME)/include

FFTW_L = -L$(FFTW_HOME)/lib %(fftwlibs)s
FFTW_I = -I$(FFTW_HOME)/include



ifeq ($(shell uname), Linux)
ARCH_LUA  = linux
ARCH_GLFW = x11
GL_L      = -lglfw -lXrandr -lX11 -lGLU -lGL
endif

ifeq ($(shell uname), Darwin)
ARCH_LUA  = macosx
ARCH_GLFW = cocoa
GL_L      = -framework OpenGL -framework Cocoa -lglfw
endif

# ------------------------------------------------------------------------------
"""

