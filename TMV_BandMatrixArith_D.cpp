
#include "TMV_VectorArith_Inline.h"
#include "TMV.h"
#include "TMV_Band.h"

namespace tmv {

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

    if (beta == T(0)) {
      for(size_t i=0;i<A.colsize(); ++i) {

	// C.row(i) = alpha * A.row(i) * B + beta * C.row(i)
	// C.row(i) = alpha * B.Transpose() * A.row(i) + beta * C.row(i)
	MultMV(alpha,B.Rows(j1,j2).QuickTranspose(),A.row(i,j1,j2),
	    beta,C.row(i));

	if (k>0) --k; else ++j1;
	if (j2<A.rowsize()) ++j2;
	else if (j1==A.rowsize()) {
	  C.Rows(i+1,C.colsize()).Zero();
	  break;
	}
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

    if (beta == T(0)) {
      C.Zero();
      if (alpha != T(0)) {
	for(size_t j=0;j<A.rowsize(); ++j) {
	  C.Rows(i1,i2) += A.col(j,i1,i2) ^ B.row(j);
	  if (k>0) --k; else ++i1;
	  if (i2<A.colsize()) ++i2;
	  else if (i1==A.colsize()) break;
	}
	if (alpha != T(1)) C *= alpha;
      }
    } else {
      if (alpha == T(1)) {
	if (beta != T(1)) C *= beta;
	for(size_t j=0;j<A.rowsize(); ++j) {
	  C.Rows(i1,i2) += A.col(j,i1,i2) ^ B.row(j);
	  if (k>0) --k; else ++i1;
	  if (i2<A.colsize()) ++i2;
	  else if (i1==A.colsize()) break;
	}
      } else if (alpha == T(-1)) {
	if (beta != T(1)) C *= beta;
	for(size_t j=0;j<A.rowsize(); ++j) {
	  C.Rows(i1,i2) -= A.col(j,i1,i2) ^ B.row(j);
	  if (k>0) --k; else ++i1;
	  if (i2<A.colsize()) ++i2;
	  else if (i1==A.colsize()) break;
	}
      } else {
	// Requires temporary
	// (or extra multiplications by alpha - 
	//  extra storage is usually preferred.)
	if (C.isrm()) {
	  Matrix<T,RowMajor> betaC = beta*C;
	  OPMultMM(alpha,A,B,T(0),C);
	  C += betaC;
	} else {
	  Matrix<T,ColMajor> betaC = beta*C;
	  OPMultMM(alpha,A,B,T(0),C);
	  C += betaC;
	}
      }
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
    TMVAssert(C.ct()==NonConj);

    for(size_t j=0;j<B.rowsize();j++)
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

    if (A.isrm() && C.isrm())
      RowMultMM(alpha,A,B,beta,C);
    else if (A.iscm() && B.isrm())
      OPMultMM(alpha,A,B,beta,C);
    else if (B.iscm() && C.iscm())
      ColMultMM(alpha,A,B,beta,C);
    else {
      if (C.colsize() < C.rowsize()) RowMultMM(alpha,A,B,beta,C);
      else ColMultMM(alpha,A,B,beta,C);
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
      if (A.rowsize() == 0) C.Zero(); 
      else if (C.isconj()) MultMM(CONJ(alpha),A.QuickConjugate(),
	  B.QuickConjugate(),CONJ(beta),C.QuickConjugate());
      else if (alpha == T(0)) {
	if (beta != T(1)) C *= beta;
      } else if (C.SameStorageAs(A)) {
	if (C.isrm()) {
	  Matrix<T,RowMajor> tempC = C;
	  DoMultMM(alpha,A,B,beta,tempC.QuickView());
	  C = tempC;
	} else {
	  Matrix<T,ColMajor> tempC = C;
	  DoMultMM(alpha,A,B,beta,tempC.QuickView());
	  C = tempC;
	}
      } else if (C.SameStorageAs(B)) {
	if (C.SameAs(B) && C.iscm()) {
	  ColMultMM(alpha,A,B,beta,C);
	} else {
	  if (C.isrm()) {
	    Matrix<T,RowMajor> tempC = C;
	    DoMultMM(alpha,A,B,beta,tempC.QuickView());
	    C = tempC;
	  } else {
	    Matrix<T,ColMajor> tempC = C;
	    DoMultMM(alpha,A,B,beta,tempC.QuickView());
	    C = tempC;
	  }
	}
      } 
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

      // C.row(i) = alpha * A.row(i) * B + beta * C.row(i)
      // C.row(i) = alpha * B.Transpose() * A.row(i) + beta * C.row(i)
      MultMV(alpha,B.SubBandMatrix(m1,m2,j1,j2,subnlo,subnhi).QuickTranspose(),
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
	  C.SubBandMatrix(i,C.colsize(),j1,C.rowsize(),0,j2-j1-1).Zero();
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

    if (beta == T(0)) {
      C.Zero();
      if (alpha != T(0)) {
	for(size_t j=0;j<A.rowsize(); ++j) {
	  C.SubMatrix(i1,i2,m1,m2) += A.col(j,i1,i2) ^ B.row(j,m1,m2);
	  if (k>0) --k; else ++i1;
	  if (i2<A.colsize()) ++i2;
	  else if (i1==A.colsize()) break;
	  if (n>0) --n; else ++m1;
	  if (m2<B.rowsize()) ++m2;
	  else if (m1==B.rowsize()) break;
	}
	if (alpha != T(1)) C *= alpha;
      }
    } else {
      if (alpha == T(1)) {
	if (beta != T(1)) C *= beta;
	for(size_t j=0;j<A.rowsize(); ++j) {
	  C.SubMatrix(i1,i2,m1,m2) += A.col(j,i1,i2) ^ B.row(j,m1,m2);
	  if (k>0) --k; else ++i1;
	  if (i2<A.colsize()) ++i2;
	  else if (i1==A.colsize()) break;
	  if (n>0) --n; else ++m1;
	  if (m2<B.rowsize()) ++m2;
	  else if (m1==B.rowsize()) break;
	}
      } else if (alpha == T(-1)) {
	if (beta != T(1)) C *= beta;
	for(size_t j=0;j<A.rowsize(); ++j) {
	  C.SubMatrix(i1,i2,m1,m2) -= A.col(j,i1,i2) ^ B.row(j,m1,m2);
	  if (k>0) --k; else ++i1;
	  if (i2<A.colsize()) ++i2;
	  else if (i1==A.colsize()) break;
	  if (n>0) --n; else ++m1;
	  if (m2<B.rowsize()) ++m2;
	  else if (m1==B.rowsize()) break;
	}
      } else {
	// Requires temporary
	// (or extra multiplications by alpha - 
	//  extra storage is usually preferred.)
	if (C.isrm()) {
	  BandMatrix<T,RowMajor> betaC = beta*C;
	  OPMultMM(alpha,A,B,T(0),C);
	  C += betaC;
	} else { 
	  BandMatrix<T,ColMajor> betaC = beta*C;
	  OPMultMM(alpha,A,B,T(0),C);
	  C += betaC;
	}
      }
    }
  }

  template <class T, class Ta, class Tb> inline void ColMultMM(const T alpha,
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
    size_t i2=C.nlo()+1;
    size_t k=C.nhi();

    size_t m1=0;
    size_t n=B.nhi();
    size_t m2=B.nlo()+1;

    size_t k2=A.colsize()-A.nlo();

    int subnlo = A.nlo();
    int subnhi = min(A.nhi(),B.nlo());
    for(size_t j=0;j<B.rowsize();j++) {
      MultMV(alpha, A.SubBandMatrix(i1,i2,m1,m2,subnlo,subnhi),
	  B.col(j,m1,m2), beta, C.col(j,i1,i2));

      if (k==0) { ++m1; ++i1; }
      else if (n==0) { ++m1; --k; ++subnlo; if(int(m2)>A.nhi()) --subnhi; }
      else { --n; --k; if(subnhi<A.nhi()) ++subnhi; }

      if (i2<C.colsize()) ++i2;
      else if (i1==C.colsize()) break;
      else if (m1>=k2) --subnlo;

      if (m2<B.colsize()) ++m2;
      else if (m1==B.colsize()) {
	if (++j < C.rowsize()) 
	  C.SubBandMatrix(i1,C.colsize(),j,C.rowsize(),i2-i1-1,0).Zero();
	break;
      }
    }
  }

  template <class T, class Ta, class Tb> inline void DoDiagMajorMultMM(
      const RealType(T) alpha, const GenBandMatrix<Ta>& A,
      const GenBandMatrix<Tb>& B, const BandMatrixView<T>& C)
    // C += alpha * A * B
  {
    TMVAssert(A.colsize() == C.colsize());
    TMVAssert(A.rowsize() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(C.nhi() == A.nhi()+B.nhi());
    TMVAssert(C.nlo() == A.nlo()+B.nlo());
    TMVAssert(alpha == RealType(T)(1) || alpha == RealType(T)(-1));
    TMVAssert(A.rowsize() > 0);
    TMVAssert(C.rowsize() > 0);
    TMVAssert(C.colsize() > 0);
    TMVAssert(C.ct()==NonConj);
    TMVAssert(A.isdm());
    TMVAssert(B.isdm());
    TMVAssert(C.isdm());

    for(int kA=-A.nlo();kA<=A.nhi();++kA) {
      int kC = kA-B.nlo();
      size_t m1A = kA < 0 ? B.nlo() : kC < 0 ? -kC : 0;
      size_t m1B = kC < 0 ? 0 : kC;
      size_t m1C = 0;

      size_t m2C = kC < 0 ?
	min(min(C.rowsize(),C.colsize()+kC),A.rowsize()-B.nlo()) :
	  min(min(C.rowsize()-kC,C.colsize()),A.rowsize()-kA);
      size_t len = m2C;
      size_t m2A = kC < 0 ? m1A + m2C : m2C;
      size_t m2B = kC < 0 ? m2C : m1B + m2C;
      size_t j2C = kC < 0 ? m2C : m2C + kC;
      
      // MJ: The next improvement is to keep track of where 
      // the diag.begin's are and call CVIt<Ta...>(Ait,1)
      for(int kB=-B.nlo();kB<=B.nhi();++kB,++kC) {

	if (A.isconj()) 
	  if (B.isconj()) 
	    DoAddElementProd(alpha,
		CVIt<Ta,Unit,Conj>(A.diag(kA,m1A,m2A).begin()),
		CVIt<Tb,Unit,Conj>(B.diag(kB,m1B,m2B).begin()),
		RealType(T)(1),
		VIt<T,Unit,NonConj>(C.diag(kC,m1C,m2C).begin()),len);
	  else 
	    DoAddElementProd(alpha,
		CVIt<Ta,Unit,Conj>(A.diag(kA,m1A,m2A).begin()),
		CVIt<Tb,Unit,NonConj>(B.diag(kB,m1B,m2B).begin()),
		RealType(T)(1),
		VIt<T,Unit,NonConj>(C.diag(kC,m1C,m2C).begin()),len);
	else 
	  if (B.isconj()) 
	    DoAddElementProd(alpha,
		CVIt<Ta,Unit,Conj>(A.diag(kA,m1A,m2A).begin()),
		CVIt<Tb,Unit,Conj>(B.diag(kB,m1B,m2B).begin()),
		RealType(T)(1),
		VIt<T,Unit,NonConj>(C.diag(kC,m1C,m2C).begin()),len);
	  else 
	    DoAddElementProd(alpha,
		CVIt<Ta,Unit,Conj>(A.diag(kA,m1A,m2A).begin()),
		CVIt<Tb,Unit,NonConj>(B.diag(kB,m1B,m2B).begin()),
		RealType(T)(1),
		VIt<T,Unit,NonConj>(C.diag(kC,m1C,m2C).begin()),len);

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
      }
    }
  }
 
  template <class T, class Ta, class Tb> inline void DoDiagMultMM(
      const RealType(T) alpha, const GenBandMatrix<Ta>& A,
      const GenBandMatrix<Tb>& B, const BandMatrixView<T>& C)
    // C += alpha * A * B
  {
    TMVAssert(A.colsize() == C.colsize());
    TMVAssert(A.rowsize() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(C.nhi() == A.nhi()+B.nhi());
    TMVAssert(C.nlo() == A.nlo()+B.nlo());
    TMVAssert(alpha == RealType(T)(1) || alpha == RealType(T)(-1));
    TMVAssert(A.rowsize() > 0);
    TMVAssert(C.rowsize() > 0);
    TMVAssert(C.colsize() > 0);
    TMVAssert(C.ct()==NonConj);

    if (A.isdm() && B.isdm() && C.isdm()) DoDiagMajorMultMM(alpha,A,B,C);
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

	  AddElementProd(T(alpha),A.diag(kA,m1A,m2A),
	      B.diag(kB,m1B,m2B),T(1),C.diag(kC,m1C,m2C));

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
    TMVAssert(alpha!= T(0));
    TMVAssert(A.rowsize() > 0);
    TMVAssert(C.rowsize() > 0);
    TMVAssert(C.colsize() > 0);
    TMVAssert(C.ct()==NonConj);

    if (beta == T(0)) {
      C.Zero();
      DoDiagMultMM(RealType(T)(1),A,B,C);
      if (alpha != T(1)) C *= alpha;
    } else {
      if (alpha == T(1)) {
	DoDiagMultMM(RealType(T)(1),A,B,C);
	if (beta != T(1)) C *= beta;
      } else if (alpha == T(-1)) {
	DoDiagMultMM(RealType(T)(-1),A,B,C);
	if (beta != T(1)) C *= beta;
      } else {
	// Requires temporary
	// (or extra multiplications by alpha - 
	//  extra storage is usually preferred.)
	BandMatrix<T,DiagMajor> betaC = beta*C;
	C.Zero();
	DoDiagMultMM(RealType(T)(1),A,B,C);
	if (alpha != T(1)) C *= alpha;
	C += betaC;
      }
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
    else if (B.iscm() && C.iscm()) ColMultMM(alpha,A,B,beta,C);
    else DiagMultMM(alpha,A,B,beta,C);
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
    //Matrix<T,RowMajor> C2 = beta*C;
    //C2 += alpha*Matrix<T,RowMajor>(A)*Matrix<T,RowMajor>(B);
    TMVAssert(A.colsize() == C.colsize());
    TMVAssert(A.rowsize() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(C.nhi() == A.nhi()+B.nhi());
    TMVAssert(C.nlo() == A.nlo()+B.nlo());

    if (C.colsize() > 0 && C.rowsize() > 0) {
      if (A.rowsize() == 0) C.Zero(); 
      else if (C.isconj()) MultMM(CONJ(alpha),A.QuickConjugate(),
	  B.QuickConjugate(),CONJ(beta),C.QuickConjugate());
      else if (alpha == T(0)) {
	if (beta != T(1)) C *= beta;
      } else if (C.SameStorageAs(A)) {
	if (C.SameAs(A) && !C.SameStorageAs(B) && C.isrm())
	  ColMultMM(alpha,B.QuickTranspose(),A.QuickTranspose(),beta,
	      C.QuickTranspose());
	else {
	  if (C.isrm()) {
	    BandMatrix<T,RowMajor> tempC = C;
	    DoMultMM(alpha,A,B,beta,tempC.QuickView());
	    C = tempC;
	  } else if (C.iscm()) {
	    BandMatrix<T,ColMajor> tempC = C;
	    DoMultMM(alpha,A,B,beta,tempC.QuickView());
	    C = tempC;
	  } else {
	    BandMatrix<T,DiagMajor> tempC = C;
	    DoMultMM(alpha,A,B,beta,tempC.QuickView());
	    C = tempC;
	  }
	}
      } else if (C.SameStorageAs(B)) {
	if (C.SameAs(B) && C.iscm()) 
	  ColMultMM(alpha,A,B,beta,C);
	else {
	  if (C.isrm()) {
	    BandMatrix<T,RowMajor> tempC = C;
	    DoMultMM(alpha,A,B,beta,tempC.QuickView());
	    C = tempC;
	  } else if (C.iscm()) {
	    BandMatrix<T,ColMajor> tempC = C;
	    DoMultMM(alpha,A,B,beta,tempC.QuickView());
	    C = tempC;
	  } else {
	    BandMatrix<T,DiagMajor> tempC = C;
	    DoMultMM(alpha,A,B,beta,tempC.QuickView());
	    C = tempC;
	  }
	}
      } 
      else DoMultMM(alpha, A, B, beta, C);
    }
    //cerr<<"Done: C = "<<Type(C)<<"  "<<C.nlo()<<','<<C.nhi()<<"  "<<C<<endl;
    //cerr<<"C2 = "<<C2<<endl;
  }

#define InstFile "TMV_BandMatrixArith_D.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


