CC= g++
INCLUDE= -Iinclude
CFLAGS= $(INCLUDE) -O -DNOBLAS -DNDEBUG
TESTCFLAGS= $(INCLUDE)
BLASLIBS= 
PREFIX=/usr/local

# The above values should be the only things you will need to change.
# There are some examples for specific systems below, which may help
# you determine what you want to use for your setup.
#
# CC defines what compiler to use.
# 
# INCLUDE specifies the appropriate directories for external .h files
# of BLAS and/or LAPACK libraries.  If you choose to use them,
# specify the directories with -I here.
#
# CFLAGS include both compiler optimizations (like -O) and some
# defines (-D...) which tell the TMV library how to act.
#
# TESTCFLAGS specifies the flags to use when compiling the test suite.
# Typically, this should not include -DNDEBUG.  You might also want a 
# different level of optimization to speed up compilation.
#
# BLASLIBS specifies the linkage requirements for your BLAS/LAPACK
# libraries.  Remember to include any appropriate -L flags to specify
# the directory for these libraries.  See the examples below for some
# specific BLAS and LAPACK library linkages.
#
# PREFIX specifies the directory to install the TMV library files
# when typing make install.
# The header files will be copied to $(PREFIX)/include.
# The library files will be copied to $(PREFIX)/lib.
# Note: PREFIX=. is invalid.
#
# The default (also the first example below) should work for all systems 
# if you don't want to work too hard at getting it to compile.
# It just probably won't be quite as fast as compiling with optimized
# BLAS and/or LAPACK libraries.



#
# For g++ with no BLAS or LAPACK library, with debugging turned off
# This is the default set given above when you first obtain this file.
#

#CC= g++
#INCLUDE= -Iinclude
#CFLAGS= $(INCLUDE) -O -DNOBLAS -DNDEBUG
#BLASLIBS=
#PREFIX=/usr/local



#
# For Intel C Compiler with Intel Math Kernel Library compiling code
# for float, double _and_ long double matrices and vectors:
# (Default is just float and double.)
#

#CC= icpc
#INCLUDE= -Iinclude -I/usr/local/include
#CFLAGS= $(INCLUDE) -O3 -DNDEBUG -DMKL -DINST_LONGDOUBLE
#BLASLIBS= -lmkl_lapack -lmkl_ia32 -lguide 



#
# For Fink g++ 4.0 with ATLAS and CLAPACK (both one directory up), 
# and leaving out the float routines.
# (So using Matrix<float> would give linking errors with this one.)
#

#CC= /sw/bin/g++-4
#INCLUDE= -Iinclude -I../CLAPACK -I../ATLAS/include -I/sw/include
#CFLAGS= $(INCLUDE) -O3 -DATLAS -DCLAPACK -DNDEBUG -DNO_INST_FLOAT
#BLASLIBS= -L../CLAPACK -llapack -lcblaswr -L/sw/lib -lf2c -L../ATLAS/lib -lcblas -latlas



#
# For Portland Group C Compiler with AMD Core Math Library for BLAS
# but the native TMV code for the LAPACK functionality.
# (This isn't necessarily preferable - it's just an example.)
#

#CC= /usr/local/pgi/linux86-64/6.1/bin/pgCC
#INCLUDE= -Iinclude
#CFLAGS= $(INCLUDE) -fast -Mcache_align -DACML -DNDEBUG -DNOLAP
#BLASLIBS= -L/usr/local/pgi/linux86-64/6.1/lib -lacml -lpgftnrtl



#
# There are a few more options which you might want to change, but
# I think the defaults will be fine for almost everyone.
#

RM= rm -f
CP= cp
AR= ar
ARFLAGS= -cr
RANLIB= ranlib
MKDIR= mkdir -p
# RM, CP, AR, ARFLAGS, RANLIB, and MKDIR are all pretty standard,
# but change them if necessary.
LIBDIR= lib
SRCDIR= src
TESTDIR= test
EXDIR= examples
LINKFILE= tmv-link
# LIBDIR indicates where to put the compiled libraries.
# SRCDIR, TESTDIR and EXDIR indicate where the various source files are.
# LINKFILE is the name of the file in the share directory that stores the 
# required linkage flags to compile programs that use the TMV library.

# Make.1 has all the real meat of the makefile.
include Make.1
