#CFLAGS1= -O3 -Wall -Werror -ansi -DNOBLAS -DNOFLOAT -DTMVFLDEBUG -DTMV_BLOCKSIZE=3
#ODIR= gcc-3.4-test
#CFLAGS1= -O3 -Wall -Werror -ansi -DNDEBUG -DMKL
#ODIR= gcc-3.4
#CFLAGS1= -O3 -Wall -Werror -ansi -DMKL
#ODIR= gcc-3.4-debug
#CC= /usr/bin/g++-3.4

#CFLAGS1= -wn3 -Kc++ -tpp7 -wd1476,1505 -DNOBLAS -DNOFLOAT -DTMVFLDEBUG -DTMV_BLOCKSIZE=3
#ODIR= icc-test
#CFLAGS1= -wn3 -Kc++ -O3 -tpp7 -wd1476,1505 -DNDEBUG -DNOBLAS
#ODIR= icc-noblas
#CFLAGS1= -wn3 -Kc++ -O3 -tpp7 -wd1476,1505 -DNOBLAS -DNOFLOAT
#ODIR= icc-noblas-debug
#CFLAGS1= -wn3 -Kc++ -O3 -tpp7 -wd1476,1505 -DNDEBUG -DNOBLAS -pg -DTMV_BLOCKSIZE=64
#ODIR= icc-noblas-pg
#CFLAGS1= -wn3 -Kc++ -O3 -tpp7 -wd1476,1505 -DNDEBUG -DMKL -DNOLAP #-DTMV_BLOCKSIZE=64
#ODIR= icc-nolap
CFLAGS1= -wn3 -Kc++ -O3 -tpp7 -wd1476,1505 -DNOLAP -DMKL
ODIR= icc-nolap-debug
#CFLAGS1= -wn3 -Kc++ -O3 -tpp7 -wd1476,1505 -DNDEBUG -DMKL -DNOLAP -pg -DTMV_BLOCKSIZE=32
#ODIR= icc-nolap-pg
#CFLAGS1= -wn3 -Kc++ -O3 -tpp7 -wd1476,1505 -DNDEBUG -DMKL -pg
#ODIR= icc-pg
#CFLAGS1= -wn3 -Kc++ -O3 -tpp7 -wd1476,1505 -DMKL 
#ODIR= icc-debug
#CFLAGS1= -wn3 -Kc++ -O3 -tpp7 -wd1476,1505 -DNDEBUG -DMKL
#ODIR= icc
CC= /usr/local/bin/icc

BLASLIBS= -lmkl_lapack -lmkl_ia32 -lguide -lpthread -lm

CLIBS= -L$(ODIR) -ltmv $(BLASLIBS)
CFLAGS= $(CFLAGS1) -I/usr/local/include -I/usr/local/pkg/intel/mkl721/include -static

#BINDIR= $(ODIR)
BINDIR= .

TMV_VH_FILES= TMV_Base.h TMV_VIt.h TMV_Vector.h TMV_VectorArith.h 
TMV_MH_FILES= TMV_BaseMatrix.h TMV_Matrix.h TMV_MatrixArith.h TMV_Divider.h TMV_LUDiv.h TMV_QRDiv.h TMV_SVDiv.h
TMV_H_FILES= $(TMV_VH_FILES) $(TMV_MH_FILES) $(TMV_PH_FILES) $(TMV_DIV_FILES)

TMV_O_FILES= $(ODIR)/TMV_Vector.o $(ODIR)/TMV_VectorArith.o $(ODIR)/TMV_Matrix.o $(ODIR)/TMV_MatrixArith_A.o $(ODIR)/TMV_MatrixArith_B.o $(ODIR)/TMV_MatrixArith_C.o $(ODIR)/TMV_MatrixArith_D.o $(ODIR)/TMV_MatrixArith_E.o $(ODIR)/TMV_Givens.o $(ODIR)/TMV_Householder.o $(ODIR)/TMV_LUDiv.o $(ODIR)/TMV_QRDiv.o $(ODIR)/TMV_SVDiv.o 

TMV_DIAGH_FILES= TMV_DiagMatrix.h TMV_DiagDiv.h TMV_DiagMatrixArith.h
TMV_DIAGO_FILES= $(ODIR)/TMV_DiagMatrix.o $(ODIR)/TMV_DiagDiv.o $(ODIR)/TMV_DiagMatrixArith.o 

TMV_TRIH_FILES= TMV_TriMatrix.h TMV_TriMatrixArith.h TMV_TriDiv.h TMV_TriMatrixArith.h
TMV_TRIO_FILES= $(ODIR)/TMV_TriMatrix.o $(ODIR)/TMV_TriDiv.o $(ODIR)/TMV_TriMatrixArith_A.o $(ODIR)/TMV_TriMatrixArith_B.o $(ODIR)/TMV_TriMatrixArith_C.o $(ODIR)/TMV_TriMatrixArith_D.o 

TMV_BANDH_FILES= TMV_Band.h TMV_BandMatrix.h TMV_BandLUDiv.h TMV_BandQRDiv.h TMV_BandSVDiv.h TMV_BandMatrixArith.h 
TMV_BANDO_FILES= $(ODIR)/TMV_BandMatrix.o $(ODIR)/TMV_BandLUDiv.o $(ODIR)/TMV_BandQRDiv.o $(ODIR)/TMV_BandSVDiv.o $(ODIR)/TMV_BandMatrixArith_A.o $(ODIR)/TMV_BandMatrixArith_B.o $(ODIR)/TMV_BandMatrixArith_C.o $(ODIR)/TMV_BandMatrixArith_D.o 

TMV_SYMH_FILES= TMV_Sym.h TMV_SymMatrix.h TMV_SymMatrixArith.h TMV_SymCHDiv.h TMV_SymSVDiv.h TMV_SymLUDiv.h 
TMV_SYMO_FILES= $(ODIR)/TMV_SymMatrix.o $(ODIR)/TMV_SymMatrixArith_A.o $(ODIR)/TMV_SymMatrixArith_B.o $(ODIR)/TMV_SymMatrixArith_C.o $(ODIR)/TMV_SymMatrixArith_D.o $(ODIR)/TMV_SymMatrixArith_E.o $(ODIR)/TMV_SymCHDiv.o $(ODIR)/TMV_SymSVDiv.o $(ODIR)/TMV_SymLUDiv.o 

TMV_ALLO_FILES= $(TMV_O_FILES) $(TMV_DIAGO_FILES) $(TMV_TRIO_FILES) $(TMV_BANDO_FILES) $(TMV_SYMO_FILES) 
TMV_ALLH_FILES= $(TMV_H_FILES) $(TMV_DIAGH_FILES) $(TMV_TRIH_FILES) $(TMV_BANDH_FILES) $(TMV_SYMH_FILES)

TMV_TESTDIAG_FILES=$(ODIR)/TMV_TestDiag.o $(ODIR)/TMV_TestDiagDiv.o 
TMV_TESTTRI_FILES=$(ODIR)/TMV_TestTri.o $(ODIR)/TMV_TestTriArith_A.o $(ODIR)/TMV_TestTriArith_B.o $(ODIR)/TMV_TestTriDiv.o 
TMV_TESTBAND_FILES=$(ODIR)/TMV_TestBand.o $(ODIR)/TMV_TestBandArith_A.o $(ODIR)/TMV_TestBandArith_B.o $(ODIR)/TMV_TestBandDiv.o 
TMV_TESTSYM_FILES= $(ODIR)/TMV_TestSym.o $(ODIR)/TMV_TestSymArith_A.o $(ODIR)/TMV_TestSymArith_B.o $(ODIR)/TMV_TestSymDiv.o
TMV_TESTO_FILES= $(ODIR)/TMV_Test.o $(ODIR)/TMV_TestVector.o $(ODIR)/TMV_TestMatrix.o $(ODIR)/TMV_TestMatrixDiv.o $(TMV_TESTDIAG_FILES) $(TMV_TESTTRI_FILES) $(TMV_TESTBAND_FILES) $(TMV_TESTSYM_FILES) 

#all: $(BINDIR)/tmvtest
#all: $(BINDIR)/tmvspeed
all: $(BINDIR)/tmvtest $(BINDIR)/tmvspeed
#all: $(ODIR)/libtmv.a
 
$(ODIR)/libtmv.a: $(TMV_ALLO_FILES)
	ar -crv $(ODIR)/libtmv.a $(TMV_ALLO_FILES)
	ranlib $(ODIR)/libtmv.a

$(BINDIR)/tmvtest: $(ODIR)/libtmv.a $(TMV_TESTO_FILES)
	$(CC) $(CFLAGS) $(TMV_TESTO_FILES) $(CLIBS) -o $(BINDIR)/tmvtest

$(BINDIR)/tmvspeed: $(ODIR)/libtmv.a $(ODIR)/TMV_Speed.o
	$(CC) $(CFLAGS) $(ODIR)/TMV_Speed.o $(CLIBS) -o $(BINDIR)/tmvspeed

$(ODIR)/TMV_Speed.o: TMV_Speed.cpp $(TMV_H_FILES) 
	$(CC) $(CFLAGS) -c TMV_Speed.cpp -o $(ODIR)/TMV_Speed.o

$(ODIR)/TMV_Test.o: TMV_Test.cpp $(TMV_ALLH_FILES) 
	$(CC) $(CFLAGS) -c TMV_Test.cpp -o $(ODIR)/TMV_Test.o

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

$(ODIR)/TMV_VectorArith.o: TMV_VectorArith.cpp TMV_VectorArith.h TMV_VectorArith.inst $(TMV_VH_FILES) TMV_VectorArith_Inline.h
	$(CC) $(CFLAGS) -c TMV_VectorArith.cpp -o $(ODIR)/TMV_VectorArith.o

$(ODIR)/TMV_Matrix.o: TMV_Matrix.cpp TMV_Matrix.h TMV_Matrix.inst $(TMV_MH_FILES)
	$(CC) $(CFLAGS) -c TMV_Matrix.cpp -o $(ODIR)/TMV_Matrix.o

$(ODIR)/TMV_MatrixArith_A.o: TMV_MatrixArith_A.cpp TMV_MatrixArith.h TMV_MatrixArith_A.inst TMV_Matrix.h TMV_BaseMatrix.h TMV_Divider.h  TMV_VectorArith_Inline.h
	$(CC) $(CFLAGS) -c TMV_MatrixArith_A.cpp -o $(ODIR)/TMV_MatrixArith_A.o

$(ODIR)/TMV_MatrixArith_B.o: TMV_MatrixArith_B.cpp TMV_MatrixArith.h TMV_MatrixArith_B.inst  TMV_Matrix.h TMV_BaseMatrix.h TMV_Divider.h TMV_VectorArith_Inline.h
	$(CC) $(CFLAGS) -c TMV_MatrixArith_B.cpp -o $(ODIR)/TMV_MatrixArith_B.o

$(ODIR)/TMV_MatrixArith_C.o: TMV_MatrixArith_C.cpp TMV_MatrixArith.h TMV_MatrixArith_C.inst  TMV_Matrix.h TMV_BaseMatrix.h TMV_Divider.h TMV_VectorArith_Inline.h
	$(CC) $(CFLAGS) -c TMV_MatrixArith_C.cpp -o $(ODIR)/TMV_MatrixArith_C.o

$(ODIR)/TMV_MatrixArith_D.o: TMV_MatrixArith_D.cpp TMV_MatrixArith.h TMV_MatrixArith_D.inst  TMV_Matrix.h TMV_BaseMatrix.h TMV_Divider.h TMV_VectorArith_Inline.h
	$(CC) $(CFLAGS) -c TMV_MatrixArith_D.cpp -o $(ODIR)/TMV_MatrixArith_D.o

$(ODIR)/TMV_MatrixArith_E.o: TMV_MatrixArith_E.cpp TMV_MatrixArith.h TMV_MatrixArith_E.inst  TMV_Matrix.h TMV_BaseMatrix.h TMV_Divider.h TMV_VectorArith_Inline.h
	$(CC) $(CFLAGS) -c TMV_MatrixArith_E.cpp -o $(ODIR)/TMV_MatrixArith_E.o

$(ODIR)/TMV_Givens.o: TMV_Givens.cpp TMV_Givens.h TMV_Givens.inst
	$(CC) $(CFLAGS) -c TMV_Givens.cpp -o $(ODIR)/TMV_Givens.o

$(ODIR)/TMV_Householder.o: TMV_Householder.cpp TMV_Householder.h TMV_Householder.inst
	$(CC) $(CFLAGS) -c TMV_Householder.cpp -o $(ODIR)/TMV_Householder.o

$(ODIR)/TMV_SVDiv.o: TMV_SVDiv.cpp TMV_SVDiv.h TMV_SVDiv.inst TMV_Householder.h TMV_Givens.h $(TMV_DIAGH_FILES)
	$(CC) $(CFLAGS) -c TMV_SVDiv.cpp -o $(ODIR)/TMV_SVDiv.o

$(ODIR)/TMV_LUDiv.o: TMV_LUDiv.cpp TMV_LUDiv.h TMV_LUDiv.inst $(TMV_DIAGH_FILES) $(TMV_TRIHFILES)
	$(CC) $(CFLAGS) -c TMV_LUDiv.cpp -o $(ODIR)/TMV_LUDiv.o

$(ODIR)/TMV_QRDiv.o: TMV_QRDiv.cpp TMV_QRDiv.h TMV_QRDiv.inst TMV_Householder.h TMV_Givens.h $(TMV_DIAGH_FILES) $(TMV_TRIHFILES)
	$(CC) $(CFLAGS) -c TMV_QRDiv.cpp -o $(ODIR)/TMV_QRDiv.o

$(ODIR)/TMV_DiagMatrix.o: TMV_DiagMatrix.cpp TMV_DiagMatrix.h TMV_DiagMatrix.inst $(TMV_DH_FILES)
	$(CC) $(CFLAGS) -c TMV_DiagMatrix.cpp -o $(ODIR)/TMV_DiagMatrix.o

$(ODIR)/TMV_DiagMatrixArith.o: TMV_DiagMatrixArith.cpp TMV_DiagMatrixArith.inst $(TMV_DIAGH_FILES) TMV_VectorArith_Inline.h
	$(CC) $(CFLAGS) -c TMV_DiagMatrixArith.cpp -o $(ODIR)/TMV_DiagMatrixArith.o

$(ODIR)/TMV_DiagDiv.o: TMV_DiagDiv.cpp TMV_DiagDiv.inst $(TMV_DIAGH_FILES)
	$(CC) $(CFLAGS) -c TMV_DiagDiv.cpp -o $(ODIR)/TMV_DiagDiv.o

$(ODIR)/TMV_TriMatrix.o: TMV_TriMatrix.cpp TMV_TriMatrix.inst $(TMV_TRIH_FILES)
	$(CC) $(CFLAGS) -c TMV_TriMatrix.cpp -o $(ODIR)/TMV_TriMatrix.o

$(ODIR)/TMV_TriMatrixArith_A.o: TMV_TriMatrixArith_A.cpp TMV_TriMatrixArith_A.inst $(TMV_TRIH_FILES) TMV_VectorArith_Inline.h
	$(CC) $(CFLAGS) -c TMV_TriMatrixArith_A.cpp -o $(ODIR)/TMV_TriMatrixArith_A.o

$(ODIR)/TMV_TriMatrixArith_B.o: TMV_TriMatrixArith_B.cpp TMV_TriMatrixArith_B.inst $(TMV_TRIH_FILES) TMV_VectorArith_Inline.h
	$(CC) $(CFLAGS) -c TMV_TriMatrixArith_B.cpp -o $(ODIR)/TMV_TriMatrixArith_B.o

$(ODIR)/TMV_TriMatrixArith_C.o: TMV_TriMatrixArith_C.cpp TMV_TriMatrixArith_C.inst $(TMV_TRIH_FILES) TMV_VectorArith_Inline.h
	$(CC) $(CFLAGS) -c TMV_TriMatrixArith_C.cpp -o $(ODIR)/TMV_TriMatrixArith_C.o

$(ODIR)/TMV_TriMatrixArith_D.o: TMV_TriMatrixArith_D.cpp TMV_TriMatrixArith_D.inst $(TMV_TRIH_FILES) TMV_VectorArith_Inline.h
	$(CC) $(CFLAGS) -c TMV_TriMatrixArith_D.cpp -o $(ODIR)/TMV_TriMatrixArith_D.o

$(ODIR)/TMV_TriDiv.o: TMV_TriDiv.cpp TMV_TriDiv.inst $(TMV_TRIH_FILES) TMV_VectorArith_Inline.h
	$(CC) $(CFLAGS) -c TMV_TriDiv.cpp -o $(ODIR)/TMV_TriDiv.o

$(ODIR)/TMV_BandMatrix.o: TMV_BandMatrix.cpp TMV_BandMatrix.inst $(TMV_BANDH_FILES)
	$(CC) $(CFLAGS) -c TMV_BandMatrix.cpp -o $(ODIR)/TMV_BandMatrix.o

$(ODIR)/TMV_BandMatrixArith_A.o: TMV_BandMatrixArith_A.cpp TMV_BandMatrixArith_A.inst $(TMV_BANDH_FILES) TMV_VectorArith_Inline.h
	$(CC) $(CFLAGS) -c TMV_BandMatrixArith_A.cpp -o $(ODIR)/TMV_BandMatrixArith_A.o

$(ODIR)/TMV_BandMatrixArith_B.o: TMV_BandMatrixArith_B.cpp TMV_BandMatrixArith_B.inst $(TMV_BANDH_FILES) TMV_VectorArith_Inline.h
	$(CC) $(CFLAGS) -c TMV_BandMatrixArith_B.cpp -o $(ODIR)/TMV_BandMatrixArith_B.o

$(ODIR)/TMV_BandMatrixArith_C.o: TMV_BandMatrixArith_C.cpp TMV_BandMatrixArith_C.inst $(TMV_BANDH_FILES) TMV_VectorArith_Inline.h
	$(CC) $(CFLAGS) -c TMV_BandMatrixArith_C.cpp -o $(ODIR)/TMV_BandMatrixArith_C.o

$(ODIR)/TMV_BandMatrixArith_D.o: TMV_BandMatrixArith_D.cpp TMV_BandMatrixArith_D.inst $(TMV_BANDH_FILES) TMV_VectorArith_Inline.h
	$(CC) $(CFLAGS) -c TMV_BandMatrixArith_D.cpp -o $(ODIR)/TMV_BandMatrixArith_D.o

$(ODIR)/TMV_BandSVDiv.o: TMV_BandSVDiv.cpp TMV_BandSVDiv.h TMV_BandSVDiv.inst TMV_Householder.h TMV_Givens.h $(TMV_BANDH_FILES)
	$(CC) $(CFLAGS) -c TMV_BandSVDiv.cpp -o $(ODIR)/TMV_BandSVDiv.o

$(ODIR)/TMV_BandLUDiv.o: TMV_BandLUDiv.cpp TMV_BandLUDiv.h TMV_BandLUDiv.inst $(TMV_BANDH_FILES) TMV_VectorArith_Inline.h
	$(CC) $(CFLAGS) -c TMV_BandLUDiv.cpp -o $(ODIR)/TMV_BandLUDiv.o

$(ODIR)/TMV_BandQRDiv.o: TMV_BandQRDiv.cpp TMV_BandQRDiv.h TMV_BandQRDiv.inst TMV_Householder.h TMV_Givens.h $(TMV_BANDH_FILES)
	$(CC) $(CFLAGS) -c TMV_BandQRDiv.cpp -o $(ODIR)/TMV_BandQRDiv.o

$(ODIR)/TMV_SymMatrix.o: TMV_SymMatrix.cpp TMV_SymMatrix.inst $(TMV_SYMH_FILES)
	$(CC) $(CFLAGS) -c TMV_SymMatrix.cpp -o $(ODIR)/TMV_SymMatrix.o

$(ODIR)/TMV_SymMatrixArith_A.o: TMV_SymMatrixArith_A.cpp TMV_SymMatrixArith_A.inst $(TMV_SYMH_FILES) TMV_VectorArith_Inline.h
	$(CC) $(CFLAGS) -c TMV_SymMatrixArith_A.cpp -o $(ODIR)/TMV_SymMatrixArith_A.o

$(ODIR)/TMV_SymMatrixArith_B.o: TMV_SymMatrixArith_B.cpp TMV_SymMatrixArith_B.inst $(TMV_SYMH_FILES) TMV_VectorArith_Inline.h
	$(CC) $(CFLAGS) -c TMV_SymMatrixArith_B.cpp -o $(ODIR)/TMV_SymMatrixArith_B.o

$(ODIR)/TMV_SymMatrixArith_C.o: TMV_SymMatrixArith_C.cpp TMV_SymMatrixArith_C.inst $(TMV_SYMH_FILES) TMV_VectorArith_Inline.h
	$(CC) $(CFLAGS) -c TMV_SymMatrixArith_C.cpp -o $(ODIR)/TMV_SymMatrixArith_C.o

$(ODIR)/TMV_SymMatrixArith_D.o: TMV_SymMatrixArith_D.cpp TMV_SymMatrixArith_D.inst $(TMV_SYMH_FILES) TMV_VectorArith_Inline.h
	$(CC) $(CFLAGS) -c TMV_SymMatrixArith_D.cpp -o $(ODIR)/TMV_SymMatrixArith_D.o

$(ODIR)/TMV_SymMatrixArith_E.o: TMV_SymMatrixArith_E.cpp TMV_SymMatrixArith_E.inst $(TMV_SYMH_FILES) TMV_VectorArith_Inline.h
	$(CC) $(CFLAGS) -c TMV_SymMatrixArith_E.cpp -o $(ODIR)/TMV_SymMatrixArith_E.o

$(ODIR)/TMV_SymSVDiv.o: TMV_SymSVDiv.cpp TMV_SymSVDiv.h TMV_SymSVDiv.inst TMV_Householder.h TMV_Givens.h $(TMV_SYMH_FILES)
	$(CC) $(CFLAGS) -c TMV_SymSVDiv.cpp -o $(ODIR)/TMV_SymSVDiv.o

$(ODIR)/TMV_SymLUDiv.o: TMV_SymLUDiv.cpp TMV_SymLUDiv.h TMV_SymLUDiv.inst $(TMV_SYMH_FILES) TMV_VectorArith_Inline.h
	$(CC) $(CFLAGS) -c TMV_SymLUDiv.cpp -o $(ODIR)/TMV_SymLUDiv.o

$(ODIR)/TMV_SymCHDiv.o: TMV_SymCHDiv.cpp TMV_SymCHDiv.h TMV_SymCHDiv.inst TMV_Householder.h TMV_Givens.h $(TMV_SYMH_FILES)
	$(CC) $(CFLAGS) -c TMV_SymCHDiv.cpp -o $(ODIR)/TMV_SymCHDiv.o


