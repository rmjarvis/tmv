CFLAGS1= -O3 -Wall -Werror -ansi -DNDEBUG
CC= /usr/bin/g++
ODIR= gcc-3.4
#CFLAGS1= -wn3 -Kc++ -tpp7 -wd1476,1505 -DNOFLOAT -DTMVFLDEBUG
#CFLAGS1= -wn3 -Kc++ -O3 -tpp7 -wd1476,1505 -DNDEBUG #-DNOFLOAT
#CFLAGS1= -wn3 -Kc++ -O3 -tpp7 -wd1476,1505 -DNDEBUG -DLONGDOUBLE
#CFLAGS1= -wn3 -Kc++ -O3 -tpp7 -wd1476,1505 -pg -DNDEBUG 
#CC= /usr/local/bin/icc
#ODIR= icc

#BLASFLAGS= -DNOBLAS
BLASLIBS= -lmkl_lapack -lmkl_ia32 -lguide -lpthread -lm

CLIBS= -L$(ODIR) -ltmv $(BLASLIBS)
CFLAGS= $(CFLAGS1) $(BLASFLAGS)

#BINDIR= $(ODIR)
BINDIR= .

TMV_VH_FILES= TMV_Base.h TMV_Blas.h TMV_VIt.h TMV_Vector.h TMV_VectorArith.h 
TMV_MH_FILES= TMV_BaseMatrix.h TMV_Matrix.h TMV_MatrixArith.h TMV_Divider.h TMV_LUDiv.h TMV_QRDiv.h TMV_SVDiv.h
TMV_PH_FILES= TMV_Permutation.h TMV_PermutationArith.h
TMV_H_FILES= $(TMV_VH_FILES) $(TMV_MH_FILES) $(TMV_PH_FILES) $(TMV_DIV_FILES)

TMV_O_FILES= $(ODIR)/TMV_Vector.o $(ODIR)/TMV_VectorArith.o $(ODIR)/TMV_Matrix.o $(ODIR)/TMV_Permutation.o $(ODIR)/TMV_MatrixArith_A.o $(ODIR)/TMV_MatrixArith_B.o $(ODIR)/TMV_MatrixArith_C.o $(ODIR)/TMV_MatrixArith_D.o $(ODIR)/TMV_MatrixArith_E.o $(ODIR)/TMV_Givens.o $(ODIR)/TMV_Householder.o $(ODIR)/TMV_LUDiv.o $(ODIR)/TMV_QRDiv.o $(ODIR)/TMV_SVDiv.o 

TMV_DIAGH_FILES= TMV_DiagMatrix.h TMV_DiagDivider.h TMV_DiagLUDiv.h TMV_DiagSVDiv.h TMV_DiagMatrixArith.h
TMV_DIAGO_FILES= $(ODIR)/TMV_DiagMatrix.o $(ODIR)/TMV_DiagLUDiv.o $(ODIR)/TMV_DiagSVDiv.o $(ODIR)/TMV_DiagMatrixArith.o 

TMV_TRIH_FILES= TMV_TriMatrix.h TMV_TriMatrixArith.h TMV_TriDivider.h TMV_TriLUDiv.h TMV_TriSVDiv.h TMV_TriMatrixArith.h
TMV_TRIO_FILES= $(ODIR)/TMV_TriMatrix.o $(ODIR)/TMV_TriLUDiv.o $(ODIR)/TMV_TriMatrixArith_A.o $(ODIR)/TMV_TriMatrixArith_B.o $(ODIR)/TMV_TriMatrixArith_C.o $(ODIR)/TMV_TriMatrixArith_D.o 

TMV_BANDH_FILES= TMV_Band.h TMV_BandMatrix.h TMV_BandLUDiv.h TMV_BandQRDiv.h TMV_BandSVDiv.h TMV_BandMatrixArith.h 
TMV_BANDO_FILES= $(ODIR)/TMV_BandMatrix.o $(ODIR)/TMV_BandLUDiv.o $(ODIR)/TMV_BandQRDiv.o $(ODIR)/TMV_BandSVDiv.o $(ODIR)/TMV_BandMatrixArith_A.o $(ODIR)/TMV_BandMatrixArith_B.o $(ODIR)/TMV_BandMatrixArith_C.o $(ODIR)/TMV_BandMatrixArith_D.o 

TMV_ALLO_FILES= $(TMV_O_FILES) $(TMV_DIAGO_FILES) $(TMV_TRIO_FILES) $(TMV_BANDO_FILES) 
TMV_ALLH_FILES= $(TMV_H_FILES) $(TMV_DIAGH_FILES) $(TMV_TRIH_FILES) $(TMV_BANDH_FILES) 

TMV_TESTO_FILES= $(ODIR)/TMV_Test.o $(ODIR)/TMV_TestVector.o $(ODIR)/TMV_TestMatrix.o $(ODIR)/TMV_TestPermutation.o $(ODIR)/TMV_TestMatrixDiv.o $(ODIR)/TMV_TestDiag.o $(ODIR)/TMV_TestDiagDiv.o $(ODIR)/TMV_TestTri.o $(ODIR)/TMV_TestTriArith_A.o $(ODIR)/TMV_TestTriArith_B.o $(ODIR)/TMV_TestTriDiv.o $(ODIR)/TMV_TestBand.o $(ODIR)/TMV_TestBandArith_A.o $(ODIR)/TMV_TestBandArith_B.o $(ODIR)/TMV_TestBandDiv.o 

all: $(BINDIR)/tmvtest
#all: $(BINDIR)/tmvtest $(BINDIR)/tmvspeed
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

$(ODIR)/TMV_TestPermutation.o: TMV_TestPermutation.cpp $(TMV_PH_FILES)
	$(CC) $(CFLAGS) -c TMV_TestPermutation.cpp -o $(ODIR)/TMV_TestPermutation.o

$(ODIR)/TMV_TestMatrixDiv.o: TMV_TestMatrixDiv.cpp $(TMV_MH_FILES)
	$(CC) $(CFLAGS) -c TMV_TestMatrixDiv.cpp -o $(ODIR)/TMV_TestMatrixDiv.o

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

$(ODIR)/TMV_TestDiag.o: TMV_TestDiag.cpp $(TMV_DIAGH_FILES)
	$(CC) $(CFLAGS) -c TMV_TestDiag.cpp -o $(ODIR)/TMV_TestDiag.o

$(ODIR)/TMV_TestDiagDiv.o: TMV_TestDiagDiv.cpp $(TMV_DIAGH_FILES)
	$(CC) $(CFLAGS) -c TMV_TestDiagDiv.cpp -o $(ODIR)/TMV_TestDiagDiv.o

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

$(ODIR)/TMV_Permutation.o: TMV_Permutation.cpp TMV_Permutation.h TMV_Permutation.inst $(TMV_PH_FILES)
	$(CC) $(CFLAGS) -c TMV_Permutation.cpp -o $(ODIR)/TMV_Permutation.o

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

$(ODIR)/TMV_DiagSVDiv.o: TMV_DiagSVDiv.cpp TMV_DiagSVDiv.inst $(TMV_DIAGH_FILES)
	$(CC) $(CFLAGS) -c TMV_DiagSVDiv.cpp -o $(ODIR)/TMV_DiagSVDiv.o

$(ODIR)/TMV_DiagLUDiv.o: TMV_DiagLUDiv.cpp TMV_DiagLUDiv.inst $(TMV_DIAGH_FILES)
	$(CC) $(CFLAGS) -c TMV_DiagLUDiv.cpp -o $(ODIR)/TMV_DiagLUDiv.o

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

$(ODIR)/TMV_TriLUDiv.o: TMV_TriLUDiv.cpp TMV_TriLUDiv.inst $(TMV_TRIH_FILES) TMV_VectorArith_Inline.h
	$(CC) $(CFLAGS) -c TMV_TriLUDiv.cpp -o $(ODIR)/TMV_TriLUDiv.o

$(ODIR)/TMV_BandSVDiv.o: TMV_BandSVDiv.cpp TMV_BandSVDiv.h TMV_BandSVDiv.inst TMV_Householder.h TMV_Givens.h $(TMV_BANDH_FILES)
	$(CC) $(CFLAGS) -c TMV_BandSVDiv.cpp -o $(ODIR)/TMV_BandSVDiv.o

$(ODIR)/TMV_BandLUDiv.o: TMV_BandLUDiv.cpp TMV_BandLUDiv.h TMV_BandLUDiv.inst $(TMV_BANDH_FILES) TMV_VectorArith_Inline.h
	$(CC) $(CFLAGS) -c TMV_BandLUDiv.cpp -o $(ODIR)/TMV_BandLUDiv.o

$(ODIR)/TMV_BandQRDiv.o: TMV_BandQRDiv.cpp TMV_BandQRDiv.h TMV_BandQRDiv.inst TMV_Householder.h TMV_Givens.h $(TMV_BANDH_FILES)
	$(CC) $(CFLAGS) -c TMV_BandQRDiv.cpp -o $(ODIR)/TMV_BandQRDiv.o

$(ODIR)/TMV_BandMatrixArith_A.o: TMV_BandMatrixArith_A.cpp TMV_BandMatrixArith_A.inst $(TMV_BANDH_FILES) TMV_VectorArith_Inline.h
	$(CC) $(CFLAGS) -c TMV_BandMatrixArith_A.cpp -o $(ODIR)/TMV_BandMatrixArith_A.o

$(ODIR)/TMV_BandMatrixArith_B.o: TMV_BandMatrixArith_B.cpp TMV_BandMatrixArith_B.inst $(TMV_BANDH_FILES) TMV_VectorArith_Inline.h
	$(CC) $(CFLAGS) -c TMV_BandMatrixArith_B.cpp -o $(ODIR)/TMV_BandMatrixArith_B.o

$(ODIR)/TMV_BandMatrixArith_C.o: TMV_BandMatrixArith_C.cpp TMV_BandMatrixArith_C.inst $(TMV_BANDH_FILES) TMV_VectorArith_Inline.h
	$(CC) $(CFLAGS) -c TMV_BandMatrixArith_C.cpp -o $(ODIR)/TMV_BandMatrixArith_C.o

$(ODIR)/TMV_BandMatrixArith_D.o: TMV_BandMatrixArith_D.cpp TMV_BandMatrixArith_D.inst $(TMV_BANDH_FILES) TMV_VectorArith_Inline.h
	$(CC) $(CFLAGS) -c TMV_BandMatrixArith_D.cpp -o $(ODIR)/TMV_BandMatrixArith_D.o

$(ODIR)/TMV_BandMatrix.o: TMV_BandMatrix.cpp TMV_BandMatrix.inst $(TMV_BANDH_FILES)
	$(CC) $(CFLAGS) -c TMV_BandMatrix.cpp -o $(ODIR)/TMV_BandMatrix.o


