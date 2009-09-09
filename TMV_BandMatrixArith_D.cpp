
#include "TMV_VectorArith_Inline.h"
#include "TMV.h"
#include "TMV_Band.h"

//#define XDEBUG

namespace tmv {

#ifdef TMV_BLOCKSIZE
  const size_t MM_BLOCKSIZE = TMV_BLOCKSIZE;
#else
  const size_t MM_BLOCKSIZE = 64;
#endif

  //
  // MultMM
  //

  template <class T, class Ta, class Tb> inline void RowMultMM(const T alpha,
      const GenBandMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const T beta, const MatrixView<T>& C)
  {
    TMVAssert(A.colsize() == C.colsize());
    TMVAssert(A.rowsize() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(alpha!= T(0));
    TMVAssert(A.rowsize() > 0);
    TMVAssert(C.rowsize() > 0);
    TMVAssert(C.colsize() > 0);
    TMVAssert(C.ct()==NonConj);

    size_t j1=0;
    size_t k=A.nlo();
    size_t j2=A.nhi()+1;

    for(size_t i=0;i<C.colsize(); ++i) {
      // C.row(i) = beta * C.row(i) + alpha * A.row(i,j1,j2) * B.Rows(j1,j2);
      MultMV(alpha,B.Rows(j1,j2).Transpose(),A.row(i,j1,j2),
	  beta,C.row(i));
      if (k>0) --k; else ++j1;
      if (j2<A.rowsize()) ++j2;
      else if (j1==A.rowsize()) {
	C.Rows(i+1,C.colsize()) *= beta;
	break;
      }
    }
  }

  template <class T, class Ta, class Tb> inline void OPMultMM(const T alpha,
      const GenBandMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const T beta, const MatrixView<T>& C)
  {
    TMVAssert(A.colsize() == C.colsize());
    TMVAssert(A.rowsize() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(alpha!= T(0));
    TMVAssert(A.rowsize() > 0);
    TMVAssert(C.rowsize() > 0);
    TMVAssert(C.colsize() > 0);
    TMVAssert(C.ct()==NonConj);

    size_t i1=0;
    size_t k=A.nhi();
    size_t i2=A.nlo()+1;

    C *= beta;
    for(size_t j=0;j<A.rowsize();++j) {
      C.Rows(i1,i2) += alpha * A.col(j,i1,i2) ^ B.row(j);
      if (k>0) --k; else ++i1;
      if (i2<A.colsize()) ++i2;
      else if (i1==A.colsize()) break;
    }
  }

  template <class T, class Ta, class Tb> inline void ColMultMM(const T alpha,
      const GenBandMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const T beta, const MatrixView<T>& C)
  {
    TMVAssert(A.colsize() == C.colsize());
    TMVAssert(A.rowsize() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(alpha!= T(0));
    TMVAssert(A.rowsize() > 0);
    TMVAssert(C.rowsize() > 0);
    TMVAssert(C.colsize() > 0);
    if (C.ct() != NonConj) {
      cerr<<"C = "<<Type(C)<<"  "<<C<<endl;
      cerr<<"C.ct() == "<<Text(C.ct())<<endl;
    }
    TMVAssert(C.ct()==NonConj);

    for(size_t j=0;j<B.rowsize();j++)
      // C.col(j) = beta * C.col(j) + alpha * A * B.col(j);
      MultMV(alpha,A,B.col(j),beta,C.col(j));
  }

  template <class T, class Ta, class Tb> inline void DoMultMM(const T alpha,
      const GenBandMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const T beta, const MatrixView<T>& C)
  {
    TMVAssert(A.colsize() == C.colsize());
    TMVAssert(A.rowsize() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(alpha!= T(0));
    TMVAssert(A.rowsize() > 0);
    TMVAssert(C.rowsize() > 0);
    TMVAssert(C.colsize() > 0);
    TMVAssert(C.ct()==NonConj);

    if (A.isrm() && C.isrm()) RowMultMM(alpha,A,B,beta,C);
    else if (A.iscm() && B.isrm()) OPMultMM(alpha,A,B,beta,C);
    else if (B.iscm() && C.iscm()) ColMultMM(alpha,A,B,beta,C);
    else if (C.colsize() < C.rowsize()) RowMultMM(alpha,A,B,beta,C);
    else ColMultMM(alpha,A,B,beta,C);
  }

  template <class T, class Ta, class Tb> void FullTempMultMM(const T alpha,
      const GenBandMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const T beta, const MatrixView<T>& C)
  {
    if (C.isrm()) {
      Matrix<T,RowMajor> tempC = C;
      DoMultMM(alpha,A,B,beta,tempC.View());
      C = tempC;
    } else {
      Matrix<T,ColMajor> tempC = C;
      DoMultMM(alpha,A,B,beta,tempC.View());
      C = tempC;
    }
  }

  template <class T, class Ta, class Tb> void BlockTempMultMM(const T alpha,
      const GenBandMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const T beta, const MatrixView<T>& C)
  {
    for (size_t j=0;j<C.rowsize();) {
      size_t j2 = min(C.rowsize(),j+MM_BLOCKSIZE);
      if (C.isrm()) {
	Matrix<T,RowMajor> tempB = B.Cols(j,j2);
	DoMultMM(alpha,A,tempB,beta,C.Cols(j,j2));
      } else  {
	Matrix<T,ColMajor> tempB = B.Cols(j,j2);
	DoMultMM(alpha,A,tempB,beta,C.Cols(j,j2));
      }
      j=j2;
    }
  }

  template <class T, class Ta, class Tb> void MultMM(const T alpha,
      const GenBandMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const T beta, const MatrixView<T>& C)
    // C = alpha * A * B + beta * C
  {
    TMVAssert(A.colsize() == C.colsize());
    TMVAssert(A.rowsize() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    if (C.colsize() > 0 && C.rowsize() > 0) {
      if (C.isconj()) MultMM(CONJ(alpha),A.Conjugate(),
	  B.Conjugate(),CONJ(beta),C.Conjugate());
      else if (A.rowsize() == 0 || alpha == T(0)) 
        C *= beta;
      else if (C.SameStorageAs(A)) 
	FullTempMultMM(alpha,A,B,beta,C);
      else if (C.SameStorageAs(B)) 
	if (C.stepi() == B.stepi() && C.stepj() == B.stepj())
	  BlockTempMultMM(alpha,A,B,beta,C);
	else 
	  FullTempMultMM(alpha,A,B,beta,C);
      else DoMultMM(alpha, A, B, beta, C);
    }
  }


  //
  // MultMM (Band * Band)
  //

  template <class T, class Ta, class Tb> inline void RowMultMM(const T alpha,
      const GenBandMatrix<Ta>& A, const GenBandMatrix<Tb>& B,
      const T beta, const BandMatrixView<T>& C)
  {
    TMVAssert(A.colsize() == C.colsize());
    TMVAssert(A.rowsize() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(C.nhi() == A.nhi()+B.nhi());
    TMVAssert(C.nlo() == A.nlo()+B.nlo());
    TMVAssert(alpha!= T(0));
    TMVAssert(A.rowsize() > 0);
    TMVAssert(C.rowsize() > 0);
    TMVAssert(C.colsize() > 0);
    TMVAssert(C.ct()==NonConj);

    size_t m1=0;
    size_t n=A.nlo();
    size_t m2=A.nhi()+1;

    size_t j1=0;
    size_t j2=C.nhi()+1;
    size_t k=C.nlo();

    size_t k2=B.rowsize()-B.nhi();

    int subnlo = min(A.nhi(),B.nlo());
    int subnhi = B.nhi();
    for(size_t i=0;i<C.colsize(); ++i) {

      // C.row(i,j1,j2) = beta * C.row(i,j1,j2) + alpha * A.row(i,m1,m2) * 
      //     B.SubBandMatrix(m1,m2,j1,j2,subnlo,subnhi);
      MultMV(alpha,B.SubBandMatrix(m1,m2,j1,j2,subnlo,subnhi).Transpose(),
	  A.row(i,m1,m2),beta,C.row(i,j1,j2));

      if (k==0) { ++m1; ++j1; }
      else if (n==0) { --k; ++m1; ++subnhi; if(int(m2)>B.nlo()) --subnlo; }
      else { --n; --k; if(subnlo<B.nlo()) ++subnlo; }

      if (j2<C.rowsize()) ++j2;
      else if (j1==C.rowsize()) break;
      else if (m1>=k2) --subnhi; 

      if (m2<A.rowsize()) ++m2;
      else if (m1==A.rowsize()) {
	if (++i < C.colsize()) 
	  C.SubBandMatrix(i,C.colsize(),j1,C.rowsize(),0,j2-j1-1) *= beta;
	break;
      }
    }
  }

  template <class T, class Ta, class Tb> inline void OPMultMM(const T alpha,
      const GenBandMatrix<Ta>& A, const GenBandMatrix<Tb>& B,
      const T beta, const BandMatrixView<T>& C)
  {
    TMVAssert(A.colsize() == C.colsize());
    TMVAssert(A.rowsize() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(C.nhi() == A.nhi()+B.nhi());
    TMVAssert(C.nlo() == A.nlo()+B.nlo());
    TMVAssert(alpha!= T(0));
    TMVAssert(A.rowsize() > 0);
    TMVAssert(C.rowsize() > 0);
    TMVAssert(C.colsize() > 0);
    TMVAssert(C.ct()==NonConj);

    size_t i1=0;
    size_t k=A.nhi();
    size_t i2=A.nlo()+1;

    size_t m1=0;
    size_t n=B.nlo();
    size_t m2=B.nhi()+1;

    C *= beta;
    for(size_t j=0;j<A.rowsize(); ++j) {
      C.SubMatrix(i1,i2,m1,m2) += alpha * A.col(j,i1,i2) ^ B.row(j,m1,m2);
      if (k>0) --k; else ++i1;
      if (i2<A.colsize()) ++i2;
      else if (i1==A.colsize()) break;
      if (n>0) --n; else ++m1;
      if (m2<B.rowsize()) ++m2;
      else if (m1==B.rowsize()) break;
    }
  }

  template <int alpha, class Ta, class Tb, class Tc> 
    inline void DoDiagMajorMultMM(
	const GenBandMatrix<Ta>& A,
	const GenBandMatrix<Tb>& B, const BandMatrixView<Tc>& C)
    // C += alpha * A * B
    {
      //cerr<<"DiagMajor:\n";
      //cerr<<"A = "<<Type(A)<<"  "<<A<<endl;
      //cerr<<"B = "<<Type(B)<<"  "<<B<<endl;
      //cerr<<"C = "<<Type(C)<<"  "<<C<<endl;
      TMVAssert(A.colsize() == C.colsize());
      TMVAssert(A.rowsize() == B.colsize());
      TMVAssert(B.rowsize() == C.rowsize());
      TMVAssert(C.nhi() == A.nhi()+B.nhi());
      TMVAssert(C.nlo() == A.nlo()+B.nlo());
      TMVAssert(alpha == 1 || alpha == -1);
      TMVAssert(A.rowsize() > 0);
      TMVAssert(C.rowsize() > 0);
      TMVAssert(C.colsize() > 0);
      TMVAssert(C.ct()==NonConj);
      TMVAssert(A.isdm());
      TMVAssert(B.isdm());
      TMVAssert(C.isdm());

      for(int kA=-A.nlo();kA<=A.nhi();++kA) {
	int kC = kA-B.nlo();

	size_t len = kC < 0 ?
	  min(min(C.rowsize(),C.colsize()+kC),A.rowsize()-B.nlo()) :
	  min(min(C.rowsize()-kC,C.colsize()),A.rowsize()-kA);
	size_t j2C = kC < 0 ? len : len + kC;

	const Ta* Aptr = kA<0 ?  A.cptr()-kA*A.stepi() : A.cptr()+kA*A.stepj();
	if (kA<0) Aptr += B.nlo();
	else if (kC<0) Aptr -= kC;
	const Tb* Bptr = B.cptr()+B.nlo()*B.stepi();
	if (kC>=0) Bptr += kC;
	Tc* Cptr = kC<0 ?  C.ptr()-kC*C.stepi() : C.ptr()+kC*C.stepj();

	for(int kB=-B.nlo();kB<=B.nhi();++kB,++kC) {

	  // C.diag(kC,m1C,m2C) += alpha * DiagMatrixViewOf(A.diag(kA,m1A,m2A))
	  //      * B.diag(kB,m1B,m2B);
	  if (A.isconj()) 
	    if (B.isconj()) 
	      DoAddElementProd(RealType(Tc)(alpha),
		  CVIt<Ta,Unit,Conj>(Aptr,1), CVIt<Tb,Unit,Conj>(Bptr,1),
		  VIt<Tc,Unit,NonConj>(Cptr,1 FIRSTLAST1(C.first,C.last) ),len);
	    else 
	      DoAddElementProd(RealType(Tc)(alpha),
		  CVIt<Ta,Unit,Conj>(Aptr,1), CVIt<Tb,Unit,NonConj>(Bptr,1),
		  VIt<Tc,Unit,NonConj>(Cptr,1 FIRSTLAST1(C.first,C.last) ),len);
	  else 
	    if (B.isconj()) 
	      DoAddElementProd(RealType(Tc)(alpha),
		  CVIt<Ta,Unit,Conj>(Aptr,1), CVIt<Tb,Unit,Conj>(Bptr,1),
		  VIt<Tc,Unit,NonConj>(Cptr,1 FIRSTLAST1(C.first,C.last) ),len);
	    else 
	      DoAddElementProd(RealType(Tc)(alpha),
		  CVIt<Ta,Unit,Conj>(Aptr,1), CVIt<Tb,Unit,NonConj>(Bptr,1),
		  VIt<Tc,Unit,NonConj>(Cptr,1 FIRSTLAST1(C.first,C.last) ),len);

	  if (kC < 0) {
	    Cptr -= C.stepi();
	    if (kB < 0) {
	      Bptr -= B.stepi();
	      --Aptr;
	      if (j2C != C.rowsize()) { ++j2C; ++len; }
	    }
	    else {
	      Bptr += B.stepj();
	      ++Cptr;
	      if (j2C == C.rowsize()) --len; else ++j2C;
	    }
	  }
	  else {
	    Cptr += C.stepj();
	    if (kB < 0)  {
	      Bptr -= B.stepi();
	      ++Bptr;
	      if (j2C == C.rowsize()) --len; else ++j2C;
	    }
	    else {
	      Bptr += B.stepj();
	      if (j2C == C.rowsize()) --len; else ++j2C;
	    }
	  }
	}
      }
    }
 
  template <int alpha, class Ta, class Tb, class Tc> inline void DoDiagMultMM(
      const GenBandMatrix<Ta>& A,
      const GenBandMatrix<Tb>& B, const BandMatrixView<Tc>& C)
    // C += alpha * A * B
  {
    TMVAssert(A.colsize() == C.colsize());
    TMVAssert(A.rowsize() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(C.nhi() == A.nhi()+B.nhi());
    TMVAssert(C.nlo() == A.nlo()+B.nlo());
    TMVAssert(alpha == 1 || alpha == -1);
    TMVAssert(A.rowsize() > 0);
    TMVAssert(C.rowsize() > 0);
    TMVAssert(C.colsize() > 0);
    TMVAssert(C.ct()==NonConj);

    //cerr<<"Diag MultMM:\n";
    //cerr<<"   A = "<<Type(A)<<endl;
    //cerr<<"   B = "<<Type(B)<<endl;
    //cerr<<"   C = "<<Type(C)<<endl;
    if (A.isdm() && B.isdm() && C.isdm()) DoDiagMajorMultMM<alpha>(A,B,C);
    else {

      // The indices for each diagonal are identified by A,B,C at the end.
      // X.diag(k,m1,m2) extends from (i1,j1) to (i2,j2)
      // i?A = i?C
      // j?A = i?B
      // j?B = j?C

      for(int kA=-A.nlo();kA<=A.nhi();++kA) {
	int kC = kA-B.nlo();
	size_t m1A = kA < 0 ? B.nlo() : kC < 0 ? -kC : 0;
	size_t m1B = kC < 0 ? 0 : kC;
	size_t m1C = 0;
	//size_t i1B = kC < 0 ? B.nlo() : kA; // = j1A
	//size_t i1C = kC < 0 ? -kC : 0; // = i1A
	//size_t j1C = kC < 0 ? 0 : kC;  // = j1B

	size_t m2C = kC < 0 ?
	  min(min(C.rowsize(),C.colsize()+kC),A.rowsize()-B.nlo()) :
	    min(min(C.rowsize()-kC,C.colsize()),A.rowsize()-kA);
	size_t len = m2C;
	size_t m2A = kC < 0 ? m1A + m2C : m2C;
	size_t m2B = kC < 0 ? m2C : m1B + m2C;
	//size_t i2B = i1B + m2C;
	//size_t i2C = i1C + m2C;
	//size_t j2C = j1C + m2C;
	size_t j2C = kC < 0 ? m2C : m2C + kC;
	//TMVAssert(i2B <= B.colsize());
	//TMVAssert(i2C <= C.colsize());
	TMVAssert(j2C <= C.rowsize());


	for(int kB=-B.nlo();kB<=B.nhi();++kB,++kC) {

	  // C.diag(kC,m1C,m2C) += alpha * DiagMatrixViewOf(A.diag(kA,m1A,m2A))
	  //      * B.diag(kB,m1B,m2B);
	  AddElementProd(Tc(alpha),A.diag(kA,m1A,m2A),
	      B.diag(kB,m1B,m2B),C.diag(kC,m1C,m2C));

	  TMVAssert(m2A-m1A == len);
	  TMVAssert(m2B-m1B == len);
	  TMVAssert(m2C-m1C == len);
	  //TMVAssert(i2B-i1B == m2C-m1C);
	  //TMVAssert(i2C-i1C == m2C-m1C);
	  //TMVAssert(j2C-j1C == m2C-m1C);
	  //TMVAssert(i2B <= B.colsize());
	  //TMVAssert(i2C <= C.colsize());
	  TMVAssert(j2C <= C.rowsize());
	  if (kC < 0) {
	    if (kB >= 0) { 
	      TMVAssert(kA<0);
	      //TMVAssert(i1B==0);
	      //TMVAssert(i1C==-kA);
	      //TMVAssert(j1C==kB); ++j1C;
	      TMVAssert(m1A==0);
	      TMVAssert(m1B==0);
	      TMVAssert(int(m1C)==kB); ++m1C;
	      //TMVAssert(i2B == m2C-kB);
	      //TMVAssert(i2C == m2C-kC);
	      TMVAssert(j2C == m2C);
	      TMVAssert(m2A == m2C-kB);
	      TMVAssert(m2B == m2C-kB);
	      if (m2C == C.rowsize()) {
		--m2B; --m2A; --len;
		//--i2C; --i2B;
	      } else {
		++m2C; 
		++j2C;
	      }
	    } else if (kA >= 0) {
	      //TMVAssert(i1B==-kB); --i1B;
	      //TMVAssert(i1C==-kC); --i1C;
	      //TMVAssert(j1C==0);
	      TMVAssert(int(m1A)==-kC); --m1A; 
	      TMVAssert(m1B==0);
	      TMVAssert(m1C==0);
	      //TMVAssert(i2B == m2C-kB);
	      //TMVAssert(i2C == m2C-kC);
	      TMVAssert(j2C == m2C);
	      TMVAssert(m2A == m2C-kC);
	      TMVAssert(m2B == m2C);
	      if (m2C == C.rowsize()) {
		--m2A;
		//--i2C; --i2B;
	      } else {
		++m2B; ++m2C; ++len;
		++j2C;
	      }
	    } else {
	      //TMVAssert(i1B==-kB); --i1B;
	      //TMVAssert(i1C==-kC); --i1C;
	      //TMVAssert(j1C==0);
	      TMVAssert(int(m1A)==-kB); --m1A;
	      TMVAssert(m1B==0);
	      TMVAssert(m1C==0);
	      //TMVAssert(i2B == m2C-kB);
	      //TMVAssert(i2C == m2C-kC);
	      TMVAssert(j2C == m2C);
	      TMVAssert(m2A == m2C-kB);
	      TMVAssert(m2B == m2C);
	      if (m2C == C.rowsize()) {
		--m2A; 
		//--i2C; --i2B;
	      } else {
		++m2B; ++m2C; ++len;
		++j2C;
	      }
	    }
	  } else {
	    if (kB < 0) { 
	      TMVAssert(kA>0);
	      //TMVAssert(i1B==kA);
	      //TMVAssert(i1C==0);
	      //TMVAssert(j1C==kC); ++j1C;
	      TMVAssert(m1A==0);
	      TMVAssert(int(m1B)==kC); ++m1B;
	      TMVAssert(m1C==0);
	      //TMVAssert(i2B == m2C+kA);
	      //TMVAssert(i2C == m2C);
	      TMVAssert(j2C == m2C+kC);
	      TMVAssert(m2A == m2C);
	      TMVAssert(m2B == m2C+kC);
	      if (m2B == C.rowsize()) {
		--m2A; --m2C; --len;
		//--i2C; --i2B;
	      } else {
		++m2B;
		++j2C;
	      }
	    } else if (kA < 0) {
	      //TMVAssert(i1B==0);
	      //TMVAssert(i1C==-kA);
	      //TMVAssert(j1C==kB); ++j1C;
	      TMVAssert(m1A==0);
	      TMVAssert(m1B==0);
	      TMVAssert(int(m1C)==-kA);
	      //TMVAssert(i2B == m2C+kA);
	      //TMVAssert(i2C == m2C);
	      TMVAssert(j2C == m2C+kC);
	      TMVAssert(m2A == m2C+kA);
	      TMVAssert(m2B == m2C+kA);
	      if (j2C == C.rowsize()) {
		--m2A; --m2B; --m2C; --len;
		//--i2C; --i2B;
	      } else {
		++j2C;
	      }
	    } else {
	      //TMVAssert(i1B==kA);
	      //TMVAssert(i1C==0);
	      //TMVAssert(j1C==kC); ++j1C;
	      TMVAssert(m1A==0);
	      TMVAssert(int(m1B)==kA);
	      TMVAssert(m1C==0);
	      //TMVAssert(i2B == m2C+kA);
	      //TMVAssert(i2C == m2C);
	      TMVAssert(j2C == m2C+kC);
	      TMVAssert(m2A == m2C);
	      TMVAssert(m2B == m2C+kA);
	      if (j2C == C.rowsize()) {
		--m2A; --m2B; --m2C; --len;
		//--i2C; --i2B;
	      } else {
		++j2C;
	      }
	    }
	  }
	  /* The distilled version: 
	   * (But I keep the above version since it isn't slower,
	   * and it is a bit more descriptive of how the indices change.)
	   if (kC < 0) 
	   if (kB >= 0) 
	   if (j2C == C.rowsize()) { ++m1C; --m2B; --m2A; --len; }
	   else { ++m1C; ++m2C; ++j2C; }
	   else 
	   if (j2C == C.rowsize()) { --m1A; --m2A; }
	   else { --m1A; ++m2B; ++m2C; ++j2C; ++len; }
	   else 
	   if (kB < 0)  
	   if (j2C == C.rowsize()) { ++m1B; --m2A; --m2C; --len; } 
	   else { ++m1B; ++m2B; ++j2C; }
	   else 
	   if (j2C == C.rowsize()) { --m2A; --m2B; --m2C; --len; }
	   else ++j2C;
	   */
	}
      }
    }
  }
 
  template <class T, class Ta, class Tb> inline void DiagMultMM(const T alpha,
      const GenBandMatrix<Ta>& A, const GenBandMatrix<Tb>& B,
      const T beta, const BandMatrixView<T>& C)
    // C = alpha * A * B + beta * C
  {
    TMVAssert(A.colsize() == C.colsize());
    TMVAssert(A.rowsize() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(C.nhi() == A.nhi()+B.nhi());
    TMVAssert(C.nlo() == A.nlo()+B.nlo());
    TMVAssert(beta == T(0) || alpha == T(1) || alpha == T(-1));
    TMVAssert(A.rowsize() > 0);
    TMVAssert(C.rowsize() > 0);
    TMVAssert(C.colsize() > 0);
    TMVAssert(C.ct()==NonConj);

    if (beta == T(0)) {
      C.Zero();
      DoDiagMultMM<1>(A,B,C);
      C *= alpha;
    } else if (alpha == T(1)) {
      DoDiagMultMM<1>(A,B,C);
      C *= beta;
    } else {
      TMVAssert(alpha == T(-1));
      DoDiagMultMM<-1>(A,B,C);
      C *= beta;
      // If alpha != +- 1, then requires temporary.
    }
  }

  template <class T, class Ta, class Tb> inline void DoMultMM(const T alpha,
      const GenBandMatrix<Ta>& A, const GenBandMatrix<Tb>& B,
      const T beta, const BandMatrixView<T>& C)
  {
    TMVAssert(A.colsize() == C.colsize());
    TMVAssert(A.rowsize() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(C.nhi() == A.nhi()+B.nhi());
    TMVAssert(C.nlo() == A.nlo()+B.nlo());
    TMVAssert(alpha!= T(0));
    TMVAssert(A.rowsize() > 0);
    TMVAssert(C.rowsize() > 0);
    TMVAssert(C.colsize() > 0);
    TMVAssert(C.ct()==NonConj);

    if (A.isrm() && C.isrm()) RowMultMM(alpha,A,B,beta,C);
    else if (A.iscm() && B.isrm()) OPMultMM(alpha,A,B,beta,C);
    else if (B.iscm() && C.iscm()) 
      RowMultMM(alpha,B.Transpose(),A.Transpose(),beta,C.Transpose());
    else if (beta == T(0) || alpha == T(1) || alpha == T(-1)) 
      DiagMultMM(alpha,A,B,beta,C);
    else if (A.isrm() || C.isrm()) RowMultMM(alpha,A,B,beta,C);
    else if (A.iscm() || B.isrm()) OPMultMM(alpha,A,B,beta,C);
    else RowMultMM(alpha,A.Transpose(),B.Transpose(),beta,C.Transpose());
  }

  template <class T, class Ta, class Tb> void TempMultMM(const T alpha,
      const GenBandMatrix<Ta>& A, const GenBandMatrix<Tb>& B,
      const T beta, const BandMatrixView<T>& C)
  {
    if (C.isrm()) {
      BandMatrix<T,RowMajor> tempC = C;
      DoMultMM(alpha,A,B,beta,tempC.View());
      C = tempC;
    } else if (C.iscm()) {
      BandMatrix<T,ColMajor> tempC = C;
      DoMultMM(alpha,A,B,beta,tempC.View());
      C = tempC;
    } else {
      BandMatrix<T,DiagMajor> tempC = C;
      DoMultMM(alpha,A,B,beta,tempC.View());
      C = tempC;
    }
  }

  template <class T, class Ta, class Tb> void MultMM(const T alpha,
      const GenBandMatrix<Ta>& A, const GenBandMatrix<Tb>& B,
      const T beta, const BandMatrixView<T>& C)
    // C = alpha * A * B + beta * C
  {
    //cerr<<"Start Band MultMM\n";
    //cerr<<"A = "<<Type(A)<<"  "<<A.nlo()<<','<<A.nhi()<<"  "<<A<<endl;
    //cerr<<"B = "<<Type(B)<<"  "<<B.nlo()<<','<<B.nhi()<<"  "<<B<<endl;
    //cerr<<"C = "<<Type(C)<<"  "<<C.nlo()<<','<<C.nhi()<<"  "<<C<<endl;
#ifdef XDEBUG
    Matrix<T> C2 = beta*Matrix<T>(C)+alpha*Matrix<T>(A)*Matrix<T>(B);
    Matrix<T> C0 = C;
#endif
    TMVAssert(A.colsize() == C.colsize());
    TMVAssert(A.rowsize() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(C.nhi() == A.nhi()+B.nhi());
    TMVAssert(C.nlo() == A.nlo()+B.nlo());

    if (C.colsize() > 0 && C.rowsize() > 0) {
      if (C.isconj()) MultMM(CONJ(alpha),A.Conjugate(),
	  B.Conjugate(),CONJ(beta),C.Conjugate());
      else if (A.rowsize() == 0 || alpha == T(0)) 
	C *= beta;
      else if (C.SameStorageAs(A) || C.SameStorageAs(B)) 
	TempMultMM(alpha,A,B,beta,C);
      else DoMultMM(alpha, A, B, beta, C);
    }
    //cerr<<"Done: C = "<<Type(C)<<"  "<<C.nlo()<<','<<C.nhi()<<"  "<<C<<endl;
#ifdef XDEBUG
    if (Norm(C2-C) > 0.001*Norm(C)) {
      cerr<<"Band MultMM alpha,beta = "<<alpha<<"  "<<beta<<endl;
      cerr<<"A = "<<Type(A)<<"  "<<A.nlo()<<','<<A.nhi()<<"  "<<A<<endl;
      cerr<<"B = "<<Type(B)<<"  "<<B.nlo()<<','<<B.nhi()<<"  "<<B<<endl;
      cerr<<"C = "<<Type(C)<<"  "<<C.nlo()<<','<<C.nhi()<<"  "<<C0<<endl;
      cerr<<"--> C = "<<C0<<endl;
      cerr<<"C2 = "<<C2<<endl;
      abort();
    }
#endif
  }

#define InstFile "TMV_BandMatrixArith_D.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


