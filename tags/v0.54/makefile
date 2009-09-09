
CC= icc
CFLAGS= -Kc++ -O3 -DMKL -I/usr/local/include 
#CFLAGS= -Kc++ -O3 -DLONGDOUBLE
#ODIR= icc-manjula-ld
ODIR= icc-mkl-debug
#BLASLIBS= -lmkl_lapack -lmkl_ia32 -lguide -lpthread
#CFLAGS= -Kc++ -O3 -DNDEBUG
#ODIR= icc-manjula-noblas
#BLASLIBS=
BINDIR= .

#CC= g++-3.4
#CFLAGS= -O3 -DNDEBUG 
#ODIR= gcc-noblas
#BLASLIBS= 

#CC= g++
#CFLAGS= -O -DNOBLAS
#ODIR= gcc-debug
#BLASLIBS=
#BINDIR= .

# The above values should be the only things you will need to change.
# There are some examples for specific systems below, which may help
# you determine what you want to use for your setup.
#
# CC defines what compiler to use.
# 
# CFLAGS include both compiler optimizations (like -O) and some
# defines (-D...) which tell the TMV library how to act.
#
# Remember, if you are using BLAS and/or LAPACK, you may need to specify
# the appropriate directories for the include files with -I here.
#
# ODIR is where the library file (libtmv.a) will be placed.
# My preference is to not put it in the main directory, so I can 
# have several different versions to link with.  
# Hence the ODIR names in the examples below are descriptive of the 
# compiler options used.
#
# BLASLIBS specifies the linkage requirements for your BLAS/LAPACK
# libraries.  Remember to include any appropriate -L flags to specify
# the directory for these libraries.  See the examples below for some
# specific BLAS and LAPACK library linkages.
#
# BINDIR specifies the location for the test suite's executable.
# Probably it's just easiest to put this in the main TMV directory 
# (ie. "."), but you could put it somewhere else if you prefer.
#
# The default (also the first example below) should work for all systems 
# if you don't want to work too hard at getting it to compile.
# It just probably won't be quite as fast as compiling with optimized
# BLAS and/or LAPACK libraries.
# 



#
# For g++ with no BLAS or LAPACK library, with debugging turned on.
# This is the default set given above when you first obtain this file.
#

#CC= g++
#CFLAGS= -O -DNOBLAS
#ODIR= gcc-debug
#BLASLIBS=



#
# For Intel C Compiler with Intel Math Kernel Library compiling code
# for float, double _and_ long double matrices and vectors:
# (Default is just float and double.)
#

#CC= icc
#CFLAGS= -Kc++ -O3 -DNDEBUG -DMKL -I/usr/local/include -DLONGDOUBLE
#ODIR= icc-mkl
#BLASLIBS= -lmkl_lapack -lmkl_ia32 -lguide -lpthread



#
# For Fink g++ 4.0 with ATLAS and CLAPACK, leaving out the float routines.
# (So using Matrix<float> would give linking errors with this one.)
#

#CC= /sw/bin/g++-4
#CFLAGS= -O3 -DATLAS -DCLAPACK -DNDEBUG -I/sw/include -DNOFLOAT
#ODIR= gcc4-clap
#BLASLIBS= -L../CLAPACK -llapack_OSX -lcblaswr -L/sw/lib -lf2c -L../ATLAS/lib/G4 -lcblas -latlas



#
# For Portland Group C Compiler with AMD Core Math Library for BLAS
# but the native TMV code for the LAPACK functionality.
# Note the -DNOSTLSORT.  For my installation of pgCC, the linker couldn't 
# find instantiations of the standard library sort function.  This flag
# tells TMV to use its own sort function (median-of-three quicksort) instead.
#

#CC= /usr/pgi/linux86-64/5.2/bin/pgCC
#CFLAGS= -fastsse -Mcache_align -DACML -DNDEBUG -DNOSTLSORT -DNOLAP
#ODIR= pgcc-acml
#BLASLIBS= -L/usr/pgi/linux86-64/5.2/lib -lacml -lpgftnrtl



# I have several programs which test the speed of some of the algorithms.
# This specifies which one gets made by make speed.
SPEED= TMV_Speed.cpp
#SPEED= TMV_Speed1.cpp
#SPEED= TMV_SymSpeed.cpp
#SPEED= TMV_SymLUSpeed.cpp

# 
# I don't think you should need to modify anything below this point.
#


CLIBS= -L$(ODIR) -ltmv $(BLASLIBS) -lm
SMALLCLIBS= -L$(ODIR) -lsmalltmv $(BLASLIBS) -lm


TMV_VH_FILES= TMV_Base.h TMV_VIt.h TMV_Vector.h TMV_VectorArith.h 
TMV_MH_FILES= TMV_BaseMatrix.h TMV_Matrix.h TMV_MatrixArith.h TMV_Divider.h TMV_LUDiv.h TMV_QRDiv.h TMV_QRPDiv.h TMV_SVDiv.h
TMV_H_FILES= $(TMV_VH_FILES) $(TMV_MH_FILES) $(TMV_PH_FILES) $(TMV_DIV_FILES)

TMV_O_FILES= $(ODIR)/TMV_Vector.o $(ODIR)/TMV_VectorArith.o $(ODIR)/TMV_Matrix.o $(ODIR)/TMV_MatrixArith_A.o $(ODIR)/TMV_MatrixArith_B.o $(ODIR)/TMV_MatrixArith_C.o $(ODIR)/TMV_MatrixArith_D.o $(ODIR)/TMV_MatrixArith_E.o $(ODIR)/TMV_Givens.o $(ODIR)/TMV_Householder.o $(ODIR)/TMV_LUDiv.o $(ODIR)/TMV_LUDiv_A.o $(ODIR)/TMV_LUDiv_B.o $(ODIR)/TMV_LUDiv_C.o $(ODIR)/TMV_QRDiv.o $(ODIR)/TMV_QRDiv_A.o $(ODIR)/TMV_QRDiv_B.o $(ODIR)/TMV_QRDiv_C.o $(ODIR)/TMV_QRDiv_D.o $(ODIR)/TMV_QRDiv_E.o $(ODIR)/TMV_QRDiv_F.o $(ODIR)/TMV_QRPDiv.o $(ODIR)/TMV_QRPDiv_A.o $(ODIR)/TMV_QRPDiv_B.o $(ODIR)/TMV_SVDiv.o $(ODIR)/TMV_SVDiv_A.o $(ODIR)/TMV_SVDiv_B.o $(ODIR)/TMV_SVDiv_C.o

TMV_DIAGH_FILES= TMV_DiagMatrix.h TMV_DiagDiv.h TMV_DiagMatrixArith.h
TMV_DIAGO_FILES= $(ODIR)/TMV_DiagMatrix.o $(ODIR)/TMV_DiagMatrixArith.o $(ODIR)/TMV_DiagDiv.o

TMV_TRIH_FILES= TMV_TriMatrix.h TMV_TriMatrixArith.h TMV_TriDiv.h TMV_TriMatrixArith.h
TMV_TRIO_FILES= $(ODIR)/TMV_TriMatrix.o $(ODIR)/TMV_TriDiv.o $(ODIR)/TMV_TriDiv_A.o $(ODIR)/TMV_TriDiv_B.o $(ODIR)/TMV_TriDiv_C.o $(ODIR)/TMV_TriDiv_D.o $(ODIR)/TMV_TriMatrixArith_A.o $(ODIR)/TMV_TriMatrixArith_B.o $(ODIR)/TMV_TriMatrixArith_C.o $(ODIR)/TMV_TriMatrixArith_D.o $(ODIR)/TMV_TriMatrixArith_E.o $(ODIR)/TMV_TriMatrixArith_F.o 

TMV_BANDH_FILES= TMV_Band.h TMV_BandMatrix.h TMV_BandLUDiv.h TMV_BandQRDiv.h TMV_BandSVDiv.h TMV_BandMatrixArith.h 
TMV_BANDO_FILES= $(ODIR)/TMV_BandMatrix.o $(ODIR)/TMV_BandLUDiv.o $(ODIR)/TMV_BandLUDiv_A.o $(ODIR)/TMV_BandLUDiv_B.o $(ODIR)/TMV_BandLUDiv_C.o $(ODIR)/TMV_BandQRDiv.o $(ODIR)/TMV_BandQRDiv_A.o $(ODIR)/TMV_BandQRDiv_B.o $(ODIR)/TMV_BandQRDiv_C.o $(ODIR)/TMV_BandSVDiv.o $(ODIR)/TMV_BandSVDiv_A.o $(ODIR)/TMV_BandSVDiv_B.o $(ODIR)/TMV_BandMatrixArith_A.o $(ODIR)/TMV_BandMatrixArith_B.o $(ODIR)/TMV_BandMatrixArith_C.o $(ODIR)/TMV_BandMatrixArith_D.o $(ODIR)/TMV_BandMatrixArith_E.o

TMV_SYMH_FILES= TMV_Sym.h TMV_SymMatrix.h TMV_SymMatrixArith.h TMV_SymCHDiv.h TMV_SymSVDiv.h TMV_SymLUDiv.h
TMV_SYMO_FILES= $(ODIR)/TMV_SymMatrix.o $(ODIR)/TMV_SymMatrixArith_A.o $(ODIR)/TMV_SymMatrixArith_B.o $(ODIR)/TMV_SymMatrixArith_C.o $(ODIR)/TMV_SymMatrixArith_D.o $(ODIR)/TMV_SymMatrixArith_E.o $(ODIR)/TMV_SymMatrixArith_F.o $(ODIR)/TMV_SymMatrixArith_G.o $(ODIR)/TMV_SymHouseholder.o $(ODIR)/TMV_SymLUDiv.o $(ODIR)/TMV_SymLUDiv_A.o $(ODIR)/TMV_SymLUDiv_B.o $(ODIR)/TMV_SymLUDiv_C.o $(ODIR)/TMV_SymCHDiv.o $(ODIR)/TMV_SymCHDiv_A.o $(ODIR)/TMV_SymCHDiv_B.o $(ODIR)/TMV_SymSVDiv.o $(ODIR)/TMV_SymSVDiv_A.o $(ODIR)/TMV_SymSVDiv_B.o 

TMV_BASICO_FILES= $(TMV_O_FILES) $(TMV_DIAGO_FILES) $(TMV_TRIO_FILES)
TMV_BASICH_FILES= $(TMV_H_FILES) $(TMV_DIAGH_FILES) $(TMV_TRIH_FILES)
TMV_ALLO_FILES= $(TMV_BASICO_FILES) $(TMV_BANDO_FILES) $(TMV_SYMO_FILES) 
TMV_ALLH_FILES= $(TMV_BASICH_FILES) $(TMV_BANDH_FILES) $(TMV_SYMH_FILES)

TMV_TESTDIAG_FILES=$(ODIR)/TMV_TestDiag.o $(ODIR)/TMV_TestDiagDiv.o 
TMV_TESTTRI_FILES=$(ODIR)/TMV_TestTri.o $(ODIR)/TMV_TestTriArith_A.o $(ODIR)/TMV_TestTriArith_B.o $(ODIR)/TMV_TestTriDiv.o 
TMV_TESTBAND_FILES=$(ODIR)/TMV_TestBand.o $(ODIR)/TMV_TestBandArith_A.o $(ODIR)/TMV_TestBandArith_B.o $(ODIR)/TMV_TestBandDiv.o 
TMV_TESTSYM_FILES= $(ODIR)/TMV_TestSym.o $(ODIR)/TMV_TestSymArith_A.o $(ODIR)/TMV_TestSymArith_B.o $(ODIR)/TMV_TestSymDiv.o
TMV_SMALLTESTO_FILES= $(ODIR)/TMV_TestVector.o $(ODIR)/TMV_TestMatrix.o $(ODIR)/TMV_TestMatrixDiv.o
TMV_TESTO_FILES= $(TMV_SMALLTESTO_FILES) $(TMV_TESTDIAG_FILES) $(TMV_TESTTRI_FILES) $(TMV_TESTBAND_FILES) $(TMV_TESTSYM_FILES) 

test: $(BINDIR)/tmvtest

lib: $(ODIR)/libtmv.a

speed: $(BINDIR)/tmvspeed

small_lib: $(ODIR)/libsmalltmv.a

small_test: $(BINDIR)/tmvsmalltest
 
$(ODIR)/libsmalltmv.a: $(TMV_BASICO_FILES)
	ar -crv $(ODIR)/libsmalltmv.a $(TMV_BASICO_FILES)
	ranlib $(ODIR)/libsmalltmv.a

$(ODIR)/libtmv.a: $(TMV_ALLO_FILES)
	ar -crv $(ODIR)/libtmv.a $(TMV_ALLO_FILES)
	ranlib $(ODIR)/libtmv.a

$(BINDIR)/tmvsmalltest: $(ODIR)/TMV_SmallTest.o $(TMV_SMALLTESTO_FILES) $(ODIR)/libsmalltmv.a 
	$(CC) $(CFLAGS) $(ODIR)/TMV_SmallTest.o $(TMV_SMALLTESTO_FILES) $(SMALLCLIBS) -o $(BINDIR)/tmvsmalltest

$(BINDIR)/tmvtest: $(ODIR)/TMV_Test.o $(TMV_TESTO_FILES) $(ODIR)/libtmv.a 
	$(CC) $(CFLAGS) $(ODIR)/TMV_Test.o $(TMV_TESTO_FILES) $(CLIBS) -o $(BINDIR)/tmvtest

$(BINDIR)/tmvspeed: $(ODIR)/libtmv.a $(ODIR)/TMV_Speed.o
	$(CC) $(CFLAGS) $(ODIR)/TMV_Speed.o $(CLIBS) -o $(BINDIR)/tmvspeed

$(ODIR)/TMV_Speed.o: $(SPEED) $(TMV_H_FILES) 
	$(CC) $(CFLAGS) -c $(SPEED) -o $(ODIR)/TMV_Speed.o

$(ODIR)/TMV_Test.o: TMV_Test.cpp $(TMV_ALLH_FILES) 
	$(CC) $(CFLAGS) -c TMV_Test.cpp -o $(ODIR)/TMV_Test.o

$(ODIR)/TMV_SmallTest.o: TMV_SmallTest.cpp $(TMV_ALLH_FILES) 
	$(CC) $(CFLAGS) -c TMV_SmallTest.cpp -o $(ODIR)/TMV_SmallTest.o

$(ODIR)/TMV_TestVector.o: TMV_TestVector.cpp $(TMV_VH_FILES) 
	$(CC) $(CFLAGS) -c TMV_TestVector.cpp -o $(ODIR)/TMV_TestVector.o

$(ODIR)/TMV_TestMatrix.o: TMV_TestMatrix.cpp $(TMV_MH_FILES)
	$(CC) $(CFLAGS) -c TMV_TestMatrix.cpp -o $(ODIR)/TMV_TestMatrix.o

$(ODIR)/TMV_TestMatrixDiv.o: TMV_TestMatrixDiv.cpp $(TMV_MH_FILES)
	$(CC) $(CFLAGS) -c TMV_TestMatrixDiv.cpp -o $(ODIR)/TMV_TestMatrixDiv.o

$(ODIR)/TMV_TestDiag.o: TMV_TestDiag.cpp $(TMV_DIAGH_FILES)
	$(CC) $(CFLAGS) -c TMV_TestDiag.cpp -o $(ODIR)/TMV_TestDiag.o

$(ODIR)/TMV_TestDiagDiv.o: TMV_TestDiagDiv.cpp $(TMV_DIAGH_FILES)
	$(CC) $(CFLAGS) -c TMV_TestDiagDiv.cpp -o $(ODIR)/TMV_TestDiagDiv.o

$(ODIR)/TMV_TestTriDiv.o: TMV_TestTriDiv.cpp $(TMV_TRIH_FILES)
	$(CC) $(CFLAGS) -c TMV_TestTriDiv.cpp -o $(ODIR)/TMV_TestTriDiv.o

$(ODIR)/TMV_TestTri.o: TMV_TestTri.cpp $(TMV_TRIH_FILES)
	$(CC) $(CFLAGS) -c TMV_TestTri.cpp -o $(ODIR)/TMV_TestTri.o

$(ODIR)/TMV_TestTriArith_A.o: TMV_TestTriArith_A.cpp $(TMV_TRIH_FILES)
	$(CC) $(CFLAGS) -c TMV_TestTriArith_A.cpp -o $(ODIR)/TMV_TestTriArith_A.o

$(ODIR)/TMV_TestTriArith_B.o: TMV_TestTriArith_B.cpp $(TMV_TRIH_FILES)
	$(CC) $(CFLAGS) -c TMV_TestTriArith_B.cpp -o $(ODIR)/TMV_TestTriArith_B.o

$(ODIR)/TMV_TestBand.o: TMV_TestBand.cpp $(TMV_BANDH_FILES) 
	$(CC) $(CFLAGS) -c TMV_TestBand.cpp -o $(ODIR)/TMV_TestBand.o

$(ODIR)/TMV_TestBandArith_A.o: TMV_TestBandArith_A.cpp $(TMV_BANDH_FILES)
	$(CC) $(CFLAGS) -c TMV_TestBandArith_A.cpp -o $(ODIR)/TMV_TestBandArith_A.o

$(ODIR)/TMV_TestBandArith_B.o: TMV_TestBandArith_B.cpp $(TMV_BANDH_FILES)
	$(CC) $(CFLAGS) -c TMV_TestBandArith_B.cpp -o $(ODIR)/TMV_TestBandArith_B.o

$(ODIR)/TMV_TestBandDiv.o: TMV_TestBandDiv.cpp $(TMV_BANDH_FILES)
	$(CC) $(CFLAGS) -c TMV_TestBandDiv.cpp -o $(ODIR)/TMV_TestBandDiv.o

$(ODIR)/TMV_TestSym.o: TMV_TestSym.cpp $(TMV_SYMH_FILES) 
	$(CC) $(CFLAGS) -c TMV_TestSym.cpp -o $(ODIR)/TMV_TestSym.o

$(ODIR)/TMV_TestSymArith_A.o: TMV_TestSymArith_A.cpp $(TMV_SYMH_FILES)
	$(CC) $(CFLAGS) -c TMV_TestSymArith_A.cpp -o $(ODIR)/TMV_TestSymArith_A.o

$(ODIR)/TMV_TestSymArith_B.o: TMV_TestSymArith_B.cpp $(TMV_SYMH_FILES)
	$(CC) $(CFLAGS) -c TMV_TestSymArith_B.cpp -o $(ODIR)/TMV_TestSymArith_B.o

$(ODIR)/TMV_TestSymDiv.o: TMV_TestSymDiv.cpp $(TMV_SYMH_FILES)
	$(CC) $(CFLAGS) -c TMV_TestSymDiv.cpp -o $(ODIR)/TMV_TestSymDiv.o

$(ODIR)/TMV_Vector.o: TMV_Vector.cpp TMV_Vector.h TMV_Vector.inst $(TMV_VH_FILES)
	$(CC) $(CFLAGS) -c TMV_Vector.cpp -o $(ODIR)/TMV_Vector.o

$(ODIR)/TMV_VectorArith.o: TMV_VectorArith.cpp TMV_VectorArith.h TMV_VectorArith.inst $(TMV_VH_FILES) 
	$(CC) $(CFLAGS) -c TMV_VectorArith.cpp -o $(ODIR)/TMV_VectorArith.o

$(ODIR)/TMV_Matrix.o: TMV_Matrix.cpp TMV_Matrix.h TMV_Matrix.inst $(TMV_MH_FILES)
	$(CC) $(CFLAGS) -c TMV_Matrix.cpp -o $(ODIR)/TMV_Matrix.o

$(ODIR)/TMV_MatrixArith_A.o: TMV_MatrixArith_A.cpp TMV_MatrixArith.h TMV_MatrixArith_A.inst TMV_Matrix.h TMV_BaseMatrix.h TMV_Divider.h
	$(CC) $(CFLAGS) -c TMV_MatrixArith_A.cpp -o $(ODIR)/TMV_MatrixArith_A.o

$(ODIR)/TMV_MatrixArith_B.o: TMV_MatrixArith_B.cpp TMV_MatrixArith.h TMV_MatrixArith_B.inst  TMV_Matrix.h TMV_BaseMatrix.h TMV_Divider.h
	$(CC) $(CFLAGS) -c TMV_MatrixArith_B.cpp -o $(ODIR)/TMV_MatrixArith_B.o

$(ODIR)/TMV_MatrixArith_C.o: TMV_MatrixArith_C.cpp TMV_MatrixArith.h TMV_MatrixArith_C.inst  TMV_Matrix.h TMV_BaseMatrix.h TMV_Divider.h
	$(CC) $(CFLAGS) -c TMV_MatrixArith_C.cpp -o $(ODIR)/TMV_MatrixArith_C.o

$(ODIR)/TMV_MatrixArith_D.o: TMV_MatrixArith_D.cpp TMV_MatrixArith.h TMV_MatrixArith_D.inst  TMV_Matrix.h TMV_BaseMatrix.h TMV_Divider.h 
	$(CC) $(CFLAGS) -c TMV_MatrixArith_D.cpp -o $(ODIR)/TMV_MatrixArith_D.o

$(ODIR)/TMV_MatrixArith_E.o: TMV_MatrixArith_E.cpp TMV_MatrixArith.h TMV_MatrixArith_E.inst  TMV_Matrix.h TMV_BaseMatrix.h TMV_Divider.h
	$(CC) $(CFLAGS) -c TMV_MatrixArith_E.cpp -o $(ODIR)/TMV_MatrixArith_E.o

$(ODIR)/TMV_Givens.o: TMV_Givens.cpp TMV_Givens.h TMV_Givens.inst $(TMV_H_FILES)
	$(CC) $(CFLAGS) -c TMV_Givens.cpp -o $(ODIR)/TMV_Givens.o

$(ODIR)/TMV_Householder.o: TMV_Householder.cpp TMV_Householder.h TMV_Householder.inst $(TMV_H_FILES)
	$(CC) $(CFLAGS) -c TMV_Householder.cpp -o $(ODIR)/TMV_Householder.o

$(ODIR)/TMV_SVDiv.o: TMV_SVDiv.cpp TMV_SVDiv.h TMV_SVDiv.inst $(TMV_DIAGH_FILES)
	$(CC) $(CFLAGS) -c TMV_SVDiv.cpp -o $(ODIR)/TMV_SVDiv.o

$(ODIR)/TMV_SVDiv_A.o: TMV_SVDiv_A.cpp TMV_SVDiv.h TMV_SVDiv_A.inst TMV_Householder.h TMV_Givens.h $(TMV_DIAGH_FILES)
	$(CC) $(CFLAGS) -c TMV_SVDiv_A.cpp -o $(ODIR)/TMV_SVDiv_A.o

$(ODIR)/TMV_SVDiv_B.o: TMV_SVDiv_B.cpp TMV_SVDiv.h TMV_SVDiv_B.inst $(TMV_DIAGH_FILES)
	$(CC) $(CFLAGS) -c TMV_SVDiv_B.cpp -o $(ODIR)/TMV_SVDiv_B.o

$(ODIR)/TMV_SVDiv_C.o: TMV_SVDiv_C.cpp TMV_SVDiv.h TMV_SVDiv_C.inst $(TMV_DIAGH_FILES)
	$(CC) $(CFLAGS) -c TMV_SVDiv_C.cpp -o $(ODIR)/TMV_SVDiv_C.o

$(ODIR)/TMV_LUDiv.o: TMV_LUDiv.cpp TMV_LUDiv.h TMV_LUDiv.inst $(TMV_DIAGH_FILES) $(TMV_TRIHFILES)
	$(CC) $(CFLAGS) -c TMV_LUDiv.cpp -o $(ODIR)/TMV_LUDiv.o

$(ODIR)/TMV_LUDiv_A.o: TMV_LUDiv_A.cpp TMV_LUDiv.h TMV_LUDiv_A.inst $(TMV_TRIHFILES)
	$(CC) $(CFLAGS) -c TMV_LUDiv_A.cpp -o $(ODIR)/TMV_LUDiv_A.o

$(ODIR)/TMV_LUDiv_B.o: TMV_LUDiv_B.cpp TMV_LUDiv.h TMV_LUDiv_B.inst $(TMV_TRIHFILES)
	$(CC) $(CFLAGS) -c TMV_LUDiv_B.cpp -o $(ODIR)/TMV_LUDiv_B.o

$(ODIR)/TMV_LUDiv_C.o: TMV_LUDiv_C.cpp TMV_LUDiv.h TMV_LUDiv_C.inst $(TMV_TRIHFILES)
	$(CC) $(CFLAGS) -c TMV_LUDiv_C.cpp -o $(ODIR)/TMV_LUDiv_C.o

$(ODIR)/TMV_QRDiv.o: TMV_QRDiv.cpp TMV_QRDiv.h TMV_QRDiv.inst $(TMV_DIAGH_FILES) $(TMV_TRIHFILES)
	$(CC) $(CFLAGS) -c TMV_QRDiv.cpp -o $(ODIR)/TMV_QRDiv.o

$(ODIR)/TMV_QRDiv_A.o: TMV_QRDiv_A.cpp TMV_QRDiv.h TMV_QRDiv_A.inst TMV_Householder.h $(TMV_TRIHFILES)
	$(CC) $(CFLAGS) -c TMV_QRDiv_A.cpp -o $(ODIR)/TMV_QRDiv_A.o

$(ODIR)/TMV_QRDiv_B.o: TMV_QRDiv_B.cpp TMV_QRDiv.h TMV_QRDiv_B.inst TMV_Householder.h $(TMV_TRIHFILES)
	$(CC) $(CFLAGS) -c TMV_QRDiv_B.cpp -o $(ODIR)/TMV_QRDiv_B.o

$(ODIR)/TMV_QRDiv_C.o: TMV_QRDiv_C.cpp TMV_QRDiv.h TMV_QRDiv_C.inst $(TMV_TRIHFILES)
	$(CC) $(CFLAGS) -c TMV_QRDiv_C.cpp -o $(ODIR)/TMV_QRDiv_C.o

$(ODIR)/TMV_QRDiv_D.o: TMV_QRDiv_D.cpp TMV_QRDiv.h TMV_QRDiv_D.inst TMV_Householder.h $(TMV_DIAGH_FILES) $(TMV_TRIHFILES)
	$(CC) $(CFLAGS) -c TMV_QRDiv_D.cpp -o $(ODIR)/TMV_QRDiv_D.o

$(ODIR)/TMV_QRDiv_E.o: TMV_QRDiv_E.cpp TMV_QRDiv.h TMV_QRDiv_E.inst TMV_Householder.h $(TMV_DIAGH_FILES) $(TMV_TRIHFILES)
	$(CC) $(CFLAGS) -c TMV_QRDiv_E.cpp -o $(ODIR)/TMV_QRDiv_E.o

$(ODIR)/TMV_QRDiv_F.o: TMV_QRDiv_F.cpp TMV_QRDiv.h TMV_QRDiv_F.inst TMV_Householder.h $(TMV_DIAGH_FILES) $(TMV_TRIHFILES)
	$(CC) $(CFLAGS) -c TMV_QRDiv_F.cpp -o $(ODIR)/TMV_QRDiv_F.o

$(ODIR)/TMV_QRPDiv.o: TMV_QRPDiv.cpp TMV_QRPDiv.h TMV_QRPDiv.inst  TMV_QRDiv.h $(TMV_DIAGH_FILES) $(TMV_TRIHFILES)
	$(CC) $(CFLAGS) -c TMV_QRPDiv.cpp -o $(ODIR)/TMV_QRPDiv.o

$(ODIR)/TMV_QRPDiv_A.o: TMV_QRPDiv_A.cpp TMV_QRPDiv.h TMV_QRPDiv_A.inst TMV_QRDiv.h TMV_Householder.h $(TMV_DIAGH_FILES) $(TMV_TRIHFILES)
	$(CC) $(CFLAGS) -c TMV_QRPDiv_A.cpp -o $(ODIR)/TMV_QRPDiv_A.o

$(ODIR)/TMV_QRPDiv_B.o: TMV_QRPDiv_B.cpp TMV_QRPDiv.h TMV_QRPDiv_B.inst TMV_QRDiv.h $(TMV_DIAGH_FILES) $(TMV_TRIHFILES)
	$(CC) $(CFLAGS) -c TMV_QRPDiv_B.cpp -o $(ODIR)/TMV_QRPDiv_B.o

$(ODIR)/TMV_DiagMatrix.o: TMV_DiagMatrix.cpp TMV_DiagMatrix.h TMV_DiagMatrix.inst $(TMV_DH_FILES)
	$(CC) $(CFLAGS) -c TMV_DiagMatrix.cpp -o $(ODIR)/TMV_DiagMatrix.o

$(ODIR)/TMV_DiagMatrixArith.o: TMV_DiagMatrixArith.cpp TMV_DiagMatrixArith.inst $(TMV_DIAGH_FILES) 
	$(CC) $(CFLAGS) -c TMV_DiagMatrixArith.cpp -o $(ODIR)/TMV_DiagMatrixArith.o

$(ODIR)/TMV_DiagDiv.o: TMV_DiagDiv.cpp TMV_DiagDiv.inst $(TMV_DIAGH_FILES)
	$(CC) $(CFLAGS) -c TMV_DiagDiv.cpp -o $(ODIR)/TMV_DiagDiv.o

$(ODIR)/TMV_TriMatrix.o: TMV_TriMatrix.cpp TMV_TriMatrix.inst $(TMV_TRIH_FILES)
	$(CC) $(CFLAGS) -c TMV_TriMatrix.cpp -o $(ODIR)/TMV_TriMatrix.o

$(ODIR)/TMV_TriMatrixArith_A.o: TMV_TriMatrixArith_A.cpp TMV_TriMatrixArith_A.inst $(TMV_TRIH_FILES) 
	$(CC) $(CFLAGS) -c TMV_TriMatrixArith_A.cpp -o $(ODIR)/TMV_TriMatrixArith_A.o

$(ODIR)/TMV_TriMatrixArith_B.o: TMV_TriMatrixArith_B.cpp TMV_TriMatrixArith_B.inst $(TMV_TRIH_FILES) 
	$(CC) $(CFLAGS) -c TMV_TriMatrixArith_B.cpp -o $(ODIR)/TMV_TriMatrixArith_B.o

$(ODIR)/TMV_TriMatrixArith_C.o: TMV_TriMatrixArith_C.cpp TMV_TriMatrixArith_C.inst $(TMV_TRIH_FILES) 
	$(CC) $(CFLAGS) -c TMV_TriMatrixArith_C.cpp -o $(ODIR)/TMV_TriMatrixArith_C.o

$(ODIR)/TMV_TriMatrixArith_D.o: TMV_TriMatrixArith_D.cpp TMV_TriMatrixArith_D.inst $(TMV_TRIH_FILES)
	$(CC) $(CFLAGS) -c TMV_TriMatrixArith_D.cpp -o $(ODIR)/TMV_TriMatrixArith_D.o

$(ODIR)/TMV_TriMatrixArith_E.o: TMV_TriMatrixArith_E.cpp TMV_TriMatrixArith_E.inst $(TMV_TRIH_FILES)
	$(CC) $(CFLAGS) -c TMV_TriMatrixArith_E.cpp -o $(ODIR)/TMV_TriMatrixArith_E.o

$(ODIR)/TMV_TriMatrixArith_F.o: TMV_TriMatrixArith_F.cpp TMV_TriMatrixArith_F.inst $(TMV_TRIH_FILES)
	$(CC) $(CFLAGS) -c TMV_TriMatrixArith_F.cpp -o $(ODIR)/TMV_TriMatrixArith_F.o

$(ODIR)/TMV_TriDiv.o: TMV_TriDiv.cpp TMV_TriDiv.inst $(TMV_TRIH_FILES) 
	$(CC) $(CFLAGS) -c TMV_TriDiv.cpp -o $(ODIR)/TMV_TriDiv.o

$(ODIR)/TMV_TriDiv_A.o: TMV_TriDiv_A.cpp TMV_TriDiv_A.inst $(TMV_TRIH_FILES) 
	$(CC) $(CFLAGS) -c TMV_TriDiv_A.cpp -o $(ODIR)/TMV_TriDiv_A.o

$(ODIR)/TMV_TriDiv_B.o: TMV_TriDiv_B.cpp TMV_TriDiv_B.inst $(TMV_TRIH_FILES) 
	$(CC) $(CFLAGS) -c TMV_TriDiv_B.cpp -o $(ODIR)/TMV_TriDiv_B.o

$(ODIR)/TMV_TriDiv_C.o: TMV_TriDiv_C.cpp TMV_TriDiv_C.inst $(TMV_TRIH_FILES) 
	$(CC) $(CFLAGS) -c TMV_TriDiv_C.cpp -o $(ODIR)/TMV_TriDiv_C.o

$(ODIR)/TMV_TriDiv_D.o: TMV_TriDiv_D.cpp TMV_TriDiv_D.inst $(TMV_TRIH_FILES) 
	$(CC) $(CFLAGS) -c TMV_TriDiv_D.cpp -o $(ODIR)/TMV_TriDiv_D.o

$(ODIR)/TMV_BandMatrix.o: TMV_BandMatrix.cpp TMV_BandMatrix.inst $(TMV_BANDH_FILES)
	$(CC) $(CFLAGS) -c TMV_BandMatrix.cpp -o $(ODIR)/TMV_BandMatrix.o

$(ODIR)/TMV_BandMatrixArith_A.o: TMV_BandMatrixArith_A.cpp TMV_BandMatrixArith_A.inst $(TMV_BANDH_FILES) 
	$(CC) $(CFLAGS) -c TMV_BandMatrixArith_A.cpp -o $(ODIR)/TMV_BandMatrixArith_A.o

$(ODIR)/TMV_BandMatrixArith_B.o: TMV_BandMatrixArith_B.cpp TMV_BandMatrixArith_B.inst $(TMV_BANDH_FILES) 
	$(CC) $(CFLAGS) -c TMV_BandMatrixArith_B.cpp -o $(ODIR)/TMV_BandMatrixArith_B.o

$(ODIR)/TMV_BandMatrixArith_C.o: TMV_BandMatrixArith_C.cpp TMV_BandMatrixArith_C.inst $(TMV_BANDH_FILES) 
	$(CC) $(CFLAGS) -c TMV_BandMatrixArith_C.cpp -o $(ODIR)/TMV_BandMatrixArith_C.o

$(ODIR)/TMV_BandMatrixArith_D.o: TMV_BandMatrixArith_D.cpp TMV_BandMatrixArith_D.inst $(TMV_BANDH_FILES) 
	$(CC) $(CFLAGS) -c TMV_BandMatrixArith_D.cpp -o $(ODIR)/TMV_BandMatrixArith_D.o

$(ODIR)/TMV_BandMatrixArith_E.o: TMV_BandMatrixArith_E.cpp TMV_BandMatrixArith_E.inst $(TMV_BANDH_FILES) 
	$(CC) $(CFLAGS) -c TMV_BandMatrixArith_E.cpp -o $(ODIR)/TMV_BandMatrixArith_E.o

$(ODIR)/TMV_BandLUDiv.o: TMV_BandLUDiv.cpp TMV_BandLUDiv.h TMV_BandLUDiv.inst $(TMV_BANDH_FILES) 
	$(CC) $(CFLAGS) -c TMV_BandLUDiv.cpp -o $(ODIR)/TMV_BandLUDiv.o

$(ODIR)/TMV_BandLUDiv_A.o: TMV_BandLUDiv_A.cpp TMV_BandLUDiv.h TMV_BandLUDiv_A.inst $(TMV_BANDH_FILES) 
	$(CC) $(CFLAGS) -c TMV_BandLUDiv_A.cpp -o $(ODIR)/TMV_BandLUDiv_A.o

$(ODIR)/TMV_BandLUDiv_B.o: TMV_BandLUDiv_B.cpp TMV_BandLUDiv.h TMV_BandLUDiv_B.inst $(TMV_BANDH_FILES) 
	$(CC) $(CFLAGS) -c TMV_BandLUDiv_B.cpp -o $(ODIR)/TMV_BandLUDiv_B.o

$(ODIR)/TMV_BandLUDiv_C.o: TMV_BandLUDiv_C.cpp TMV_BandLUDiv.h TMV_BandLUDiv_C.inst $(TMV_BANDH_FILES) 
	$(CC) $(CFLAGS) -c TMV_BandLUDiv_C.cpp -o $(ODIR)/TMV_BandLUDiv_C.o

$(ODIR)/TMV_BandQRDiv.o: TMV_BandQRDiv.cpp TMV_BandQRDiv.h TMV_BandQRDiv.inst TMV_Householder.h $(TMV_BANDH_FILES)
	$(CC) $(CFLAGS) -c TMV_BandQRDiv.cpp -o $(ODIR)/TMV_BandQRDiv.o

$(ODIR)/TMV_BandQRDiv_A.o: TMV_BandQRDiv_A.cpp TMV_BandQRDiv.h TMV_BandQRDiv_A.inst TMV_Householder.h $(TMV_BANDH_FILES)
	$(CC) $(CFLAGS) -c TMV_BandQRDiv_A.cpp -o $(ODIR)/TMV_BandQRDiv_A.o

$(ODIR)/TMV_BandQRDiv_B.o: TMV_BandQRDiv_B.cpp TMV_BandQRDiv.h TMV_BandQRDiv_B.inst TMV_Householder.h $(TMV_BANDH_FILES)
	$(CC) $(CFLAGS) -c TMV_BandQRDiv_B.cpp -o $(ODIR)/TMV_BandQRDiv_B.o

$(ODIR)/TMV_BandQRDiv_C.o: TMV_BandQRDiv_C.cpp TMV_BandQRDiv.h TMV_BandQRDiv_C.inst $(TMV_BANDH_FILES)
	$(CC) $(CFLAGS) -c TMV_BandQRDiv_C.cpp -o $(ODIR)/TMV_BandQRDiv_C.o

$(ODIR)/TMV_BandSVDiv.o: TMV_BandSVDiv.cpp TMV_BandSVDiv.h TMV_BandSVDiv.inst $(TMV_BANDH_FILES)
	$(CC) $(CFLAGS) -c TMV_BandSVDiv.cpp -o $(ODIR)/TMV_BandSVDiv.o

$(ODIR)/TMV_BandSVDiv_A.o: TMV_BandSVDiv_A.cpp TMV_BandSVDiv.h TMV_BandSVDiv_A.inst TMV_Householder.h $(TMV_BANDH_FILES)
	$(CC) $(CFLAGS) -c TMV_BandSVDiv_A.cpp -o $(ODIR)/TMV_BandSVDiv_A.o

$(ODIR)/TMV_BandSVDiv_B.o: TMV_BandSVDiv_B.cpp TMV_BandSVDiv.h TMV_BandSVDiv_B.inst $(TMV_BANDH_FILES)
	$(CC) $(CFLAGS) -c TMV_BandSVDiv_B.cpp -o $(ODIR)/TMV_BandSVDiv_B.o

$(ODIR)/TMV_SymMatrix.o: TMV_SymMatrix.cpp TMV_SymMatrix.inst $(TMV_SYMH_FILES)
	$(CC) $(CFLAGS) -c TMV_SymMatrix.cpp -o $(ODIR)/TMV_SymMatrix.o

$(ODIR)/TMV_SymMatrixArith_A.o: TMV_SymMatrixArith_A.cpp TMV_SymMatrixArith_A.inst $(TMV_SYMH_FILES) 
	$(CC) $(CFLAGS) -c TMV_SymMatrixArith_A.cpp -o $(ODIR)/TMV_SymMatrixArith_A.o

$(ODIR)/TMV_SymMatrixArith_B.o: TMV_SymMatrixArith_B.cpp TMV_SymMatrixArith_B.inst $(TMV_SYMH_FILES) 
	$(CC) $(CFLAGS) -c TMV_SymMatrixArith_B.cpp -o $(ODIR)/TMV_SymMatrixArith_B.o

$(ODIR)/TMV_SymMatrixArith_C.o: TMV_SymMatrixArith_C.cpp TMV_SymMatrixArith_C.inst $(TMV_SYMH_FILES) 
	$(CC) $(CFLAGS) -c TMV_SymMatrixArith_C.cpp -o $(ODIR)/TMV_SymMatrixArith_C.o

$(ODIR)/TMV_SymMatrixArith_D.o: TMV_SymMatrixArith_D.cpp TMV_SymMatrixArith_D.inst $(TMV_SYMH_FILES) 
	$(CC) $(CFLAGS) -c TMV_SymMatrixArith_D.cpp -o $(ODIR)/TMV_SymMatrixArith_D.o

$(ODIR)/TMV_SymMatrixArith_E.o: TMV_SymMatrixArith_E.cpp TMV_SymMatrixArith_E.inst $(TMV_SYMH_FILES) 
	$(CC) $(CFLAGS) -c TMV_SymMatrixArith_E.cpp -o $(ODIR)/TMV_SymMatrixArith_E.o

$(ODIR)/TMV_SymMatrixArith_F.o: TMV_SymMatrixArith_F.cpp TMV_SymMatrixArith_F.inst $(TMV_SYMH_FILES) 
	$(CC) $(CFLAGS) -c TMV_SymMatrixArith_F.cpp -o $(ODIR)/TMV_SymMatrixArith_F.o

$(ODIR)/TMV_SymMatrixArith_G.o: TMV_SymMatrixArith_G.cpp TMV_SymMatrixArith_G.inst $(TMV_SYMH_FILES) 
	$(CC) $(CFLAGS) -c TMV_SymMatrixArith_G.cpp -o $(ODIR)/TMV_SymMatrixArith_G.o

$(ODIR)/TMV_SymHouseholder.o: TMV_SymHouseholder.cpp TMV_SymHouseholder.h TMV_SymHouseholder.inst $(TMV_H_FILES)
	$(CC) $(CFLAGS) -c TMV_SymHouseholder.cpp -o $(ODIR)/TMV_SymHouseholder.o

$(ODIR)/TMV_SymLUDiv.o: TMV_SymLUDiv.cpp TMV_SymLUDiv.h TMV_SymLUDiv.inst $(TMV_SYMH_FILES) 
	$(CC) $(CFLAGS) -c TMV_SymLUDiv.cpp -o $(ODIR)/TMV_SymLUDiv.o

$(ODIR)/TMV_SymLUDiv_A.o: TMV_SymLUDiv_A.cpp TMV_SymLUDiv.h TMV_SymLUDiv_A.inst $(TMV_SYMH_FILES) 
	$(CC) $(CFLAGS) -c TMV_SymLUDiv_A.cpp -o $(ODIR)/TMV_SymLUDiv_A.o

$(ODIR)/TMV_SymLUDiv_B.o: TMV_SymLUDiv_B.cpp TMV_SymLUDiv.h TMV_SymLUDiv_B.inst $(TMV_SYMH_FILES) 
	$(CC) $(CFLAGS) -c TMV_SymLUDiv_B.cpp -o $(ODIR)/TMV_SymLUDiv_B.o

$(ODIR)/TMV_SymLUDiv_C.o: TMV_SymLUDiv_C.cpp TMV_SymLUDiv.h TMV_SymLUDiv_C.inst $(TMV_SYMH_FILES) 
	$(CC) $(CFLAGS) -c TMV_SymLUDiv_C.cpp -o $(ODIR)/TMV_SymLUDiv_C.o

$(ODIR)/TMV_SymCHDiv.o: TMV_SymCHDiv.cpp TMV_SymCHDiv.h TMV_SymCHDiv.inst $(TMV_SYMH_FILES)
	$(CC) $(CFLAGS) -c TMV_SymCHDiv.cpp -o $(ODIR)/TMV_SymCHDiv.o

$(ODIR)/TMV_SymCHDiv_A.o: TMV_SymCHDiv_A.cpp TMV_SymCHDiv.h TMV_SymCHDiv_A.inst $(TMV_SYMH_FILES)
	$(CC) $(CFLAGS) -c TMV_SymCHDiv_A.cpp -o $(ODIR)/TMV_SymCHDiv_A.o

$(ODIR)/TMV_SymCHDiv_B.o: TMV_SymCHDiv_B.cpp TMV_SymCHDiv.h TMV_SymCHDiv_B.inst $(TMV_SYMH_FILES)
	$(CC) $(CFLAGS) -c TMV_SymCHDiv_B.cpp -o $(ODIR)/TMV_SymCHDiv_B.o

$(ODIR)/TMV_SymSVDiv.o: TMV_SymSVDiv.cpp TMV_SymSVDiv.h TMV_SymSVDiv.inst $(TMV_SYMH_FILES)
	$(CC) $(CFLAGS) -c TMV_SymSVDiv.cpp -o $(ODIR)/TMV_SymSVDiv.o

$(ODIR)/TMV_SymSVDiv_A.o: TMV_SymSVDiv_A.cpp TMV_SymSVDiv.h TMV_SymSVDiv_A.inst TMV_SymHouseholder.h TMV_Givens.h $(TMV_SYMH_FILES)
	$(CC) $(CFLAGS) -c TMV_SymSVDiv_A.cpp -o $(ODIR)/TMV_SymSVDiv_A.o

$(ODIR)/TMV_SymSVDiv_B.o: TMV_SymSVDiv_B.cpp TMV_SymSVDiv.h TMV_SymSVDiv_B.inst $(TMV_SYMH_FILES)
	$(CC) $(CFLAGS) -c TMV_SymSVDiv_B.cpp -o $(ODIR)/TMV_SymSVDiv_B.o

