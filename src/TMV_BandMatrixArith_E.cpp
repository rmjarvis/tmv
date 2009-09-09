///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// The Template Matrix/Vector Library for C++ was created by Mike Jarvis     //
// Copyright (C) 2007                                                        //
//                                                                           //
// The project is hosted at http://sourceforge.net/projects/tmv-cpp/         //
// where you can find the current version and current documention.           //
//                                                                           //
// For concerns or problems with the software, Mike may be contacted at      //
// mike_jarvis@users.sourceforge.net                                         //
//                                                                           //
// This program is free software; you can redistribute it and/or             //
// modify it under the terms of the GNU General Public License               //
// as published by the Free Software Foundation; either version 2            //
// of the License, or (at your option) any later version.                    //
//                                                                           //
// This program is distributed in the hope that it will be useful,           //
// but WITHOUT ANY WARRANTY; without even the implied warranty of            //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             //
// GNU General Public License for more details.                              //
//                                                                           //
// You should have received a copy of the GNU General Public License         //
// along with this program in the file gpl.txt.                              //
//                                                                           //
// If not, write to:                                                         //
// The Free Software Foundation, Inc.                                        //
// 51 Franklin Street, Fifth Floor,                                          //
// Boston, MA  02110-1301, USA.                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////



#include "TMV_BandMatrix.h"
#include "TMV_VectorArith.h"
#include "TMV_MatrixArith.h"
#include "TMV_BandMatrixArith.h"

//#define XDEBUG

#ifdef XDEBUG
#include <iostream>
using std::cerr;
using std::endl;
#endif

namespace tmv {

  //
  // MultMM (Band * Band)
  //

  template <bool add, class T, class Ta, class Tb> inline void RowMultMM(
      const T alpha, const GenBandMatrix<Ta>& A, const GenBandMatrix<Tb>& B,
      const BandMatrixView<T>& C)
  {
    //cerr<<"RowMultMM: alpha = "<<alpha<<endl;
    //cerr<<"add = "<<add<<endl;
    //cerr<<"A = "<<Type(A)<<"  "<<A<<endl;
    //cerr<<"B = "<<Type(B)<<"  "<<B<<endl;
    //cerr<<"C = "<<Type(C)<<"  "<<C<<endl;

    TMVAssert(A.colsize() == C.colsize());
    TMVAssert(A.rowsize() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(C.nhi() == std::min(int(C.rowsize()-1),A.nhi()+B.nhi()));
    TMVAssert(C.nlo() == std::min(int(C.colsize()-1),A.nlo()+B.nlo()));
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

    int subnlo = std::min(A.nhi(),B.nlo());
    int subnhi = B.nhi();
    for(size_t i=0;i<C.colsize(); ++i) {
      //cerr<<"i = "<<i<<", j1,j2 = "<<j1<<','<<j2<<", m1,m2 = "<<m1<<','<<m2<<endl;
      //cerr<<"subnlo,subnhi = "<<subnlo<<','<<subnhi<<endl;

      // C.row(i,j1,j2) (+)= alpha * A.row(i,m1,m2) * 
      //     B.SubBandMatrix(m1,m2,j1,j2,subnlo,subnhi);
      MultMV<add>(alpha,B.SubBandMatrix(m1,m2,j1,j2,subnlo,subnhi).Transpose(),
	  A.row(i,m1,m2),C.row(i,j1,j2));
      //cerr<<"C.row("<<i<<','<<j1<<','<<j2<<") = "<<C.row(i,j1,j2)<<endl;

      if (k==0) { ++m1; ++j1; }
      else if (n==0) { --k; ++m1; ++subnhi; if(int(m2)>B.nlo()) --subnlo; }
      else { --n; --k; if(subnlo<B.nlo()) ++subnlo; }

      if (j2<C.rowsize()) ++j2;
      else if (j1==C.rowsize()) break;
      else if (m1>=k2) --subnhi; 

      if (m2<A.rowsize()) ++m2;
      else if (m1==A.rowsize()) {
	if (!add && ++i < C.colsize()) 
	  C.SubBandMatrix(i,C.colsize(),j1,C.rowsize(),0,j2-j1-1).Zero();
	break;
      }
    }
    //cerr<<"Done: C = "<<C<<endl;
  }

  template <bool add, class T, class Ta, class Tb> inline void OPMultMM(
      const T alpha, const GenBandMatrix<Ta>& A, const GenBandMatrix<Tb>& B,
      const BandMatrixView<T>& C)
  {
    //cerr<<"OPMultMM: alpha = "<<alpha<<endl;
    //cerr<<"add = "<<add<<endl;
    //cerr<<"A = "<<Type(A)<<"  "<<A<<endl;
    //cerr<<"B = "<<Type(B)<<"  "<<B<<endl;
    //cerr<<"C = "<<Type(C)<<"  "<<C<<endl;

    TMVAssert(A.colsize() == C.colsize());
    TMVAssert(A.rowsize() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(C.nhi() == std::min(int(C.rowsize()-1),A.nhi()+B.nhi()));
    TMVAssert(C.nlo() == std::min(int(C.colsize()-1),A.nlo()+B.nlo()));
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

    if (!add) C.Zero();
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
    TMVAssert(C.nhi() == std::min(int(C.rowsize()-1),A.nhi()+B.nhi()));
    TMVAssert(C.nlo() == std::min(int(C.colsize()-1),A.nlo()+B.nlo()));
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
	std::min(std::min(C.rowsize(),C.colsize()+kC),A.rowsize()+kB1) :
	std::min(std::min(C.rowsize()-kC,C.colsize()),A.rowsize()-kA);
      
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
 
  template <bool add, class T, class Ta, class Tb> inline void DiagMultMM(
      const T alpha, const GenBandMatrix<Ta>& A, const GenBandMatrix<Tb>& B,
      const BandMatrixView<T>& C)
    // C (+)= alpha * A * B
  {
    TMVAssert(A.colsize() == C.colsize());
    TMVAssert(A.rowsize() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(C.nhi() == std::min(int(C.rowsize()-1),A.nhi()+B.nhi()));
    TMVAssert(C.nlo() == std::min(int(C.colsize()-1),A.nlo()+B.nlo()));
    // If alpha != +- 1 and add, then requires temporary.
    TMVAssert(!add || alpha == T(1) || alpha == T(-1));
    TMVAssert(A.rowsize() > 0);
    TMVAssert(C.rowsize() > 0);
    TMVAssert(C.colsize() > 0);
    TMVAssert(C.ct()==NonConj);

    if (!add) {
      C.Zero();
      DoDiagMultMM<1>(A,B,C);
      if (alpha != T(1)) C *= alpha;
    } else if (alpha == T(1)) {
      DoDiagMultMM<1>(A,B,C);
    } else {
      TMVAssert(alpha == T(-1));
      DoDiagMultMM<-1>(A,B,C);
    }
  }

  template <bool add, class T, class Ta, class Tb> inline void DoMultMM(
      const T alpha, const GenBandMatrix<Ta>& A, const GenBandMatrix<Tb>& B,
      const BandMatrixView<T>& C)
  {
    TMVAssert(A.colsize() == C.colsize());
    TMVAssert(A.rowsize() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(C.nhi() == std::min(int(C.rowsize()-1),A.nhi()+B.nhi()));
    TMVAssert(C.nlo() == std::min(int(C.colsize()-1),A.nlo()+B.nlo()));
    TMVAssert(alpha!= T(0));
    TMVAssert(A.rowsize() > 0);
    TMVAssert(C.rowsize() > 0);
    TMVAssert(C.colsize() > 0);
    TMVAssert(C.ct()==NonConj);

#ifdef ACML
    // This is really strange.  The ACML version gets this wrong sometimes
    // when A and B have the same storage.
    // It seems that the error happens in the blas gemv routine.
    // Specifically, I found that:
    //
    // [ a b c d e f g h ]   [ i ]
    // [ i j k l m n o p ] * [ j ]
    // [ q r s t u v w x ]   [ k ]
    //                       [ l ]
    //                       [ m ]
    //                       [ n ]
    //                       [ o ]
    //                       [ p ]
    // 
    // with the matrix being Conj and the vector being NonConj
    //
    // gets the 2nd element in the resulting vector wrong.
    // This kind of thing can happen sometimes (depending on nlo, nhi, etc.)
    // in this function.  So to be safe, we just copy B for ACML compiles.
    //
    if (SameStorage(A,B)) {
      if (B.isrm()) {
	BandMatrix<Tb,RowMajor> B2 = B;
	return DoMultMM<add>(alpha,A,B2,C);
      } else if (B.iscm()) {
	BandMatrix<Tb,ColMajor> B2 = B;
	return DoMultMM<add>(alpha,A,B2,C);
      } else {
	BandMatrix<Tb,DiagMajor> B2 = B;
	return DoMultMM<add>(alpha,A,B2,C);
      }
    } 
#endif

    if (A.isrm() && C.isrm()) RowMultMM<add>(alpha,A,B,C);
    else if (A.iscm() && B.isrm()) OPMultMM<add>(alpha,A,B,C);
    else if (B.iscm() && C.iscm()) 
      RowMultMM<add>(alpha,B.Transpose(),A.Transpose(),C.Transpose());
    else if (!add || alpha == T(1) || alpha == T(-1)) 
      DiagMultMM<add>(alpha,A,B,C);
    else if (A.isdm() && B.isdm() && C.isdm()) {
      BandMatrix<T,DiagMajor> CC(C.colsize(),C.rowsize(),C.nlo(),C.nhi());
      DiagMultMM<false>(T(1),A,B,CC.View());
      if (add) C += alpha*CC;
      else C = alpha*CC;
    }
    else if (A.isrm() || C.isrm()) RowMultMM<add>(alpha,A,B,C);
    else if (A.iscm() || B.isrm()) OPMultMM<add>(alpha,A,B,C);
    else RowMultMM<add>(alpha,B.Transpose(),A.Transpose(),C.Transpose());
  }

  template <bool add, class T, class Ta, class Tb> inline void TempMultMM(
      const T alpha, const GenBandMatrix<Ta>& A, const GenBandMatrix<Tb>& B,
      const BandMatrixView<T>& C)
  {
    if (C.isrm()) {
      BandMatrix<T,RowMajor> C2(C.colsize(),C.rowsize(),C.nlo(),C.nhi());
      DoMultMM<false>(T(1),A,B,C2.View());
      if (add) C += alpha*C2;
      else C = alpha*C2;
    } else if (C.iscm()) {
      BandMatrix<T,ColMajor> C2(C.colsize(),C.rowsize(),C.nlo(),C.nhi());
      DoMultMM<false>(T(1),A,B,C2.View());
      if (add) C += alpha*C2;
      else C = alpha*C2;
    } else {
      BandMatrix<T,DiagMajor> C2(C.colsize(),C.rowsize(),C.nlo(),C.nhi());
      DoMultMM<false>(T(1),A,B,C2.View());
      if (add) C += alpha*C2;
      else C = alpha*C2;
    }
  }

  template <bool add, class T, class Ta, class Tb> void MultMM(const T alpha,
      const GenBandMatrix<Ta>& A, const GenBandMatrix<Tb>& B,
      const BandMatrixView<T>& C)
    // C (+)= alpha * A * B
  {
#ifdef XDEBUG
    //cerr<<"Start Band MultMM\n";
    //cerr<<"A = "<<A.cptr()<<"  "<<Type(A)<<"  "<<A.cptr()<<"  "<<A.nlo()<<','<<A.nhi()<<"  "<<A<<endl;
    //cerr<<"B = "<<B.cptr()<<"  "<<Type(B)<<"  "<<B.cptr()<<"  "<<B.nlo()<<','<<B.nhi()<<"  "<<B<<endl;
    //cerr<<"C = "<<C.cptr()<<"  "<<Type(C)<<"  "<<C.cptr()<<"  "<<C.nlo()<<','<<C.nhi()<<"  "<<C<<endl;
    Matrix<Ta> A0 = A;
    Matrix<Tb> B0 = B;
    Matrix<T> C0 = C;
    Matrix<T> C2 = alpha*A0*B0;
    if (add) C2 += C0;
#endif
    TMVAssert(A.colsize() == C.colsize());
    TMVAssert(A.rowsize() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(C.nhi() >= std::min(int(C.rowsize()-1),A.nhi()+B.nhi()));
    TMVAssert(C.nlo() >= std::min(int(C.colsize()-1),A.nlo()+B.nlo()));

    if (C.colsize() > 0 && C.rowsize() > 0) {
      if (A.rowsize() == 0 || alpha == T(0)) {
	if (!add) C.Zero();
      } else if (A.rowsize() > A.colsize()+A.nhi()) {
	ConstBandMatrixView<Ta> AA = A.Cols(0,A.colsize()+A.nhi());
	ConstBandMatrixView<Tb> BB = B.SubBandMatrix(0,AA.rowsize(),0,B.rowsize(),
	    std::min(B.nlo(),int(AA.rowsize())-1),B.nhi());
	MultMM<add>(alpha,AA,BB,C);
      } else if (A.colsize() > A.rowsize()+A.nlo()) {
	ConstBandMatrixView<Ta> AA = A.Rows(0,A.rowsize()+A.nlo());
	BandMatrixView<T> CC = C.SubBandMatrix(0,AA.colsize(),0,C.rowsize(),
	    std::min(C.nlo(),int(AA.colsize())-1),C.nhi());
	MultMM<add>(alpha,AA,B,CC);
	if (!add) C.Rows(A.rowsize()+A.nlo(),A.colsize()).Zero();
      } else if (B.colsize() > B.rowsize()+B.nlo()) {
	ConstBandMatrixView<Tb> BB = B.Rows(0,B.rowsize()+B.nlo());
	ConstBandMatrixView<Ta> AA = A.SubBandMatrix(0,A.colsize(),0,BB.rowsize(),
	    A.nlo(),std::min(A.nhi(),int(BB.rowsize())-1));
	MultMM<add>(alpha,AA,BB,C);
      } else if (B.rowsize() > B.colsize()+B.nhi()) {
	ConstBandMatrixView<Tb> BB = B.Cols(0,B.colsize()+B.nhi());
	BandMatrixView<T> CC = C.SubBandMatrix(0,C.colsize(),0,BB.rowsize(),
	    C.nlo(),std::min(C.nhi(),int(BB.rowsize())-1));
	MultMM<add>(alpha,A,BB,CC);
	if (!add) C.Cols(B.colsize()+B.nhi(),B.rowsize()).Zero();
      } else {
	int nhi = std::min(int(C.rowsize()-1),A.nhi()+B.nhi());
	int nlo = std::min(int(C.colsize()-1),A.nlo()+B.nlo());
	if (C.nhi() > nhi || C.nlo() > nlo) {
	  MultMM<add>(alpha,A,B,C.Diags(-nlo,nhi+1));
	  if (!add) {
	    if (C.nlo() > nlo)
	      C.Diags(-C.nlo(),-nlo).Zero();
	    if (C.nhi() > nhi)
	      C.Diags(nhi+1,C.nhi()+1).Zero();
	  }
	}
	else if (C.isconj()) 
	  MultMM<add>(CONJ(alpha),A.Conjugate(),B.Conjugate(),C.Conjugate());
	else if (SameStorage(A,C) || SameStorage(B,C)) 
	  TempMultMM<add>(alpha,A,B,C);
	else DoMultMM<add>(alpha, A, B, C);
      }
    }
#ifdef XDEBUG
    //cerr<<"Done: C = "<<Type(C)<<"  "<<C.nlo()<<','<<C.nhi()<<"  "<<C<<endl;
    //cerr<<"C2-C = "<<(C2-C)<<endl;
    //cerr<<"Norm(C2-C) = "<<Norm(C2-C)<<endl;
    if (Norm(C2-C) > 0.001*(ABS(alpha)*Norm(A0)*Norm(B0)+
	  (add?Norm(C0):RealType(T)(0)))) {
      cerr<<"Band MultMM alpha = "<<alpha<<endl;
      cerr<<"add = "<<add<<endl;
      cerr<<"A = "<<Type(A)<<"  "<<A.nlo()<<','<<A.nhi()<<"  "<<A<<endl;
      cerr<<"B = "<<Type(B)<<"  "<<B.nlo()<<','<<B.nhi()<<"  "<<B<<endl;
      cerr<<"C = "<<Type(C)<<"  "<<C.nlo()<<','<<C.nhi()<<"  "<<C0<<endl;
      cerr<<"--> C = "<<C0<<endl;
      cerr<<"C2 = "<<C2<<endl;
      abort();
    }
    //cerr<<"End MultMM\n";
#endif
  }

#define InstFile "TMV_BandMatrixArith_E.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


