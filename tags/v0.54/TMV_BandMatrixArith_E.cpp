
#include "TMV.h"
#include "TMV_Band.h"

//#define XDEBUG

namespace tmv {

  //
  // MultMM (Band * Band)
  //

  template <class T, class Ta, class Tb> inline void RowMultMM(
      const T alpha, const GenBandMatrix<Ta>& A, const GenBandMatrix<Tb>& B,
      const T beta, const BandMatrixView<T>& C)
  {
    //cerr<<"RowMultMM: alpha,beta = "<<alpha<<','<<beta<<endl;
    //cerr<<"A = "<<Type(A)<<"  "<<A<<endl;
    //cerr<<"B = "<<Type(B)<<"  "<<B<<endl;
    //cerr<<"C = "<<Type(C)<<"  "<<C<<endl;

    TMVAssert(A.colsize() == C.colsize());
    TMVAssert(A.rowsize() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(C.nhi() == min(int(C.rowsize()-1),A.nhi()+B.nhi()));
    TMVAssert(C.nlo() == min(int(C.colsize()-1),A.nlo()+B.nlo()));
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
      //cerr<<"i = "<<i<<", j1,j2 = "<<j1<<','<<j2<<", m1,m2 = "<<m1<<','<<m2<<endl;
      //cerr<<"subnlo,subnhi = "<<subnlo<<','<<subnhi<<endl;

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
    //cerr<<"Done: C = "<<C<<endl;
  }

  template <class T, class Ta, class Tb> inline void OPMultMM(
      const T alpha, const GenBandMatrix<Ta>& A, const GenBandMatrix<Tb>& B,
      const T beta, const BandMatrixView<T>& C)
  {
    //cerr<<"OPMultMM: alpha,beta = "<<alpha<<','<<beta<<endl;
    //cerr<<"A = "<<Type(A)<<"  "<<A<<endl;
    //cerr<<"B = "<<Type(B)<<"  "<<B<<endl;
    //cerr<<"C = "<<Type(C)<<"  "<<C<<endl;

    TMVAssert(A.colsize() == C.colsize());
    TMVAssert(A.rowsize() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(C.nhi() == min(int(C.rowsize()-1),A.nhi()+B.nhi()));
    TMVAssert(C.nlo() == min(int(C.colsize()-1),A.nlo()+B.nlo()));
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
    //cerr<<"Done: C = "<<C<<endl;
  }

  template <int alpha, class Ta, class Tb, class Tc> inline void DoDiagMultMM(
      const GenBandMatrix<Ta>& A,
      const GenBandMatrix<Tb>& B, const BandMatrixView<Tc>& C)
    // C += alpha * A * B
  {
    //cerr<<"DiagMultMM: alpha = "<<alpha<<endl;
    //cerr<<"A = "<<Type(A)<<"  "<<A<<endl;
    //cerr<<"B = "<<Type(B)<<"  "<<B<<endl;
    //cerr<<"C = "<<Type(C)<<"  "<<C<<endl;

    TMVAssert(A.colsize() == C.colsize());
    TMVAssert(A.rowsize() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(C.nhi() == min(int(C.rowsize()-1),A.nhi()+B.nhi()));
    TMVAssert(C.nlo() == min(int(C.colsize()-1),A.nlo()+B.nlo()));
    TMVAssert(alpha == 1 || alpha == -1);
    TMVAssert(A.rowsize() > 0);
    TMVAssert(C.rowsize() > 0);
    TMVAssert(C.colsize() > 0);
    TMVAssert(C.ct()==NonConj);

    // The indices for each diagonal are identified by A,B,C at the end.
    // X.diag(k,m1,m2) extends from (i1,j1) to (i2,j2)
    // i?A = i?C
    // j?A = i?B
    // j?B = j?C

    for(int kA=-A.nlo();kA<=A.nhi();++kA) {
      int kC = kA-B.nlo();
      int kB1 = -B.nlo();
      if (kC < -C.nlo()) { kC = -C.nlo(); kB1 = kC-kA; }

      size_t m1A = kA < 0 ? -kB1 : kC < 0 ? -kC : 0;
      size_t m1B = kC < 0 ? 0 : kC;
      size_t m1C = 0;
      //size_t i1B = kC < 0 ? B.nlo() : kA; // = j1A
      //size_t i1C = kC < 0 ? -kC : 0; // = i1A
      //size_t j1C = kC < 0 ? 0 : kC;  // = j1B

      size_t len = kC < 0 ?
	min(min(C.rowsize(),C.colsize()+kC),A.rowsize()+kB1) :
	min(min(C.rowsize()-kC,C.colsize()),A.rowsize()-kA);
      
      size_t m2C = len;
      size_t m2A = kC < 0 ? m1A + m2C : m2C;
      size_t m2B = kC < 0 ? m2C : m1B + m2C;
      //size_t i2B = i1B + m2C;
      //size_t i2C = i1C + m2C;
      //size_t j2C = j1C + m2C;
      size_t j2C = kC < 0 ? m2C : m2C + kC;
      //TMVAssert(i2B <= B.colsize());
      //TMVAssert(i2C <= C.colsize());
      TMVAssert(j2C <= C.rowsize());

      for(int kB=kB1;kB<=B.nhi();++kB,++kC) {
	if (kC > C.nhi()) break;
	//cerr<<"kA, kB, kC = "<<kA<<','<<kB<<','<<kC<<endl;
	//cerr<<"m1A, m2A = "<<m1A<<','<<m2A<<endl;
	//cerr<<"m1B, m2B = "<<m1B<<','<<m2B<<endl;
	//cerr<<"m1C, m2C = "<<m1C<<','<<m2C<<endl;

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
    //cerr<<"Done: C = "<<C<<endl;
  }
 
  template <class T, class Ta, class Tb> inline void DiagMultMM(
      const T alpha, const GenBandMatrix<Ta>& A, const GenBandMatrix<Tb>& B,
      const T beta, const BandMatrixView<T>& C)
    // C = alpha * A * B + beta * C
  {
    TMVAssert(A.colsize() == C.colsize());
    TMVAssert(A.rowsize() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(C.nhi() == min(int(C.rowsize()-1),A.nhi()+B.nhi()));
    TMVAssert(C.nlo() == min(int(C.colsize()-1),A.nlo()+B.nlo()));
    // If alpha != +- 1 and beta != 0, then requires temporary.
    TMVAssert(beta == T(0) || alpha == T(1) || alpha == T(-1));
    TMVAssert(A.rowsize() > 0);
    TMVAssert(C.rowsize() > 0);
    TMVAssert(C.colsize() > 0);
    TMVAssert(C.ct()==NonConj);

    if (beta == T(0)) {
      C.Zero();
      DoDiagMultMM<1>(A,B,C);
      if (alpha != T(1)) C *= alpha;
    } else if (alpha == T(1)) {
      DoDiagMultMM<1>(A,B,C);
      if (beta != T(1)) C *= beta;
    } else {
      TMVAssert(alpha == T(-1));
      DoDiagMultMM<-1>(A,B,C);
      if (beta != T(1)) C *= beta;
    }
  }

  template <class T, class Ta, class Tb> inline void DoMultMM(
      const T alpha, const GenBandMatrix<Ta>& A, const GenBandMatrix<Tb>& B,
      const T beta, const BandMatrixView<T>& C)
  {
    TMVAssert(A.colsize() == C.colsize());
    TMVAssert(A.rowsize() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(C.nhi() == min(int(C.rowsize()-1),A.nhi()+B.nhi()));
    TMVAssert(C.nlo() == min(int(C.colsize()-1),A.nlo()+B.nlo()));
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
    else if (A.isdm() && B.isdm() && C.isdm()) {
      BandMatrix<T,DiagMajor> CC(C.colsize(),C.rowsize(),C.nlo(),C.nhi());
      DiagMultMM(T(1),A,B,T(0),CC.View());
      AddMM(alpha,CC,beta,C);
    }
    else if (A.isrm() || C.isrm()) RowMultMM(alpha,A,B,beta,C);
    else if (A.iscm() || B.isrm()) OPMultMM(alpha,A,B,beta,C);
    else RowMultMM(alpha,A.Transpose(),B.Transpose(),beta,C.Transpose());
  }

  template <class T, class Ta, class Tb> inline void TempMultMM(
      const T alpha, const GenBandMatrix<Ta>& A, const GenBandMatrix<Tb>& B,
      const T beta, const BandMatrixView<T>& C)
  {
    if (C.isrm()) {
      BandMatrix<T,RowMajor> C2(C.colsize(),C.rowsize(),C.nlo(),C.nhi());
      DoMultMM(T(1),A,B,T(0),C2.View());
      C = alpha*C2+beta*C;
    } else if (C.iscm()) {
      BandMatrix<T,ColMajor> C2(C.colsize(),C.rowsize(),C.nlo(),C.nhi());
      DoMultMM(T(1),A,B,T(0),C2.View());
      C = alpha*C2+beta*C;
    } else {
      BandMatrix<T,DiagMajor> C2(C.colsize(),C.rowsize(),C.nlo(),C.nhi());
      DoMultMM(T(1),A,B,T(0),C2.View());
      C = alpha*C2+beta*C;
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
    Matrix<T> A0 = A;
    Matrix<T> B0 = B;
    Matrix<T> C0 = C;
    Matrix<T> C2 = beta*C0+alpha*A0*B0;
#endif
    TMVAssert(A.colsize() == C.colsize());
    TMVAssert(A.rowsize() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(C.nhi() >= min(int(C.rowsize()-1),A.nhi()+B.nhi()));
    TMVAssert(C.nlo() >= min(int(C.colsize()-1),A.nlo()+B.nlo()));

    if (C.colsize() > 0 && C.rowsize() > 0) {
      if (A.rowsize() == 0 || alpha == T(0)) {
	C *= beta;
      } else if (A.rowsize() > A.colsize()+A.nhi()) {
	MultMM(alpha,A.Cols(0,A.colsize()+A.nhi()),
	    B.Rows(0,A.colsize()+A.nhi()),beta,C);
      } else if (A.colsize() > A.rowsize()+A.nlo()) {
	MultMM(alpha,A.Rows(0,A.rowsize()+A.nlo()),
	    B,beta,C.Rows(0,A.rowsize()+A.nlo()));
	C.Rows(A.rowsize()+A.nlo(),A.colsize()) *= beta;
      } else if (B.colsize() > B.rowsize()+B.nlo()) {
	MultMM(alpha,A.Cols(0,B.rowsize()+B.nlo()),
	    B.Rows(0,B.rowsize()+B.nlo()),beta,C);
      } else if (B.rowsize() > B.colsize()+B.nhi()) {
	MultMM(alpha,A,B.Cols(0,B.colsize()+B.nhi()),
	    beta,C.Cols(0,B.colsize()+B.nhi()));
	C.Cols(B.colsize()+B.nhi(),B.rowsize()) *= beta;
      } else {
	int nhi = min(int(C.rowsize()-1),A.nhi()+B.nhi());
	int nlo = min(int(C.colsize()-1),A.nlo()+B.nlo());
	if (C.nhi() > nhi || C.nlo() > nlo) {
	  MultMM(alpha,A,B,beta,C.Diags(-nlo,nhi+1));
	  if (beta != T(1)) {
	    C.Diags(-C.nlo(),-nlo) *= beta;
	    C.Diags(nhi+1,C.nhi()+1) *= beta;
	  }
	}
	else if (C.isconj()) MultMM(CONJ(alpha),A.Conjugate(),B.Conjugate(),
	    CONJ(beta),C.Conjugate());
	else if (C.SameStorageAs(A) || C.SameStorageAs(B)) 
	  TempMultMM(alpha,A,B,beta,C);
	else DoMultMM(alpha, A, B, beta, C);
      }
    }
    //cerr<<"Done: C = "<<Type(C)<<"  "<<C.nlo()<<','<<C.nhi()<<"  "<<C<<endl;
#ifdef XDEBUG
    if (Norm(C2-C) > 0.001*(abs(alpha)*Norm(A0)*Norm(B0)+
	  (beta==T(0)?RealType(T)(0):abs(beta)*Norm(C0)))) {
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

#define InstFile "TMV_BandMatrixArith_E.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


