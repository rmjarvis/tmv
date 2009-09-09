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



#include "TMV_VIt.h"
#include "TMV_Matrix.h"
#include "TMV_TriMatrix.h"
#include "TMV_DiagMatrix.h"
#include "TMV_SymMatrix.h"
#include "TMV_SymLUDiv.h"
#include "TMV_SymLUDiv_A.h"
#include "TMV_SymCHDiv_A.h"
#include "TMV_TriDiv.h"
#include "TMV_TriMatrixArith.h"
#include "TMV_SymMatrixArith.h"
#include "TMV_MatrixArith.h"

//#define XDEBUG

#ifdef XDEBUG
#include <iostream>
using std::cerr;
using std::endl;
#endif

namespace tmv {

#ifdef TMV_BLOCKSIZE
#define SYM_LU_BLOCKSIZE2 (TMV_BLOCKSIZE/2)
#else
#define SYM_LU_BLOCKSIZE2 32
#endif

  template <bool herm, class T> inline void LDLt_InvertComponents(
      const SymMatrixView<T>& sinv, const VectorView<T>& xD)
  {
    TMVAssert(IsReal(T()) || herm == sinv.isherm());
    TMVAssert(sinv.uplo() == Lower);
    TMVAssert(sinv.size() == xD.size()+1);
    TMVAssert(xD.step() == 1);

#ifdef XDEBUG
    Matrix<T> sinv0(sinv);
    Matrix<T> D0(DiagMatrixViewOf(sinv.diag()));
    D0.diag(-1) = xD;
    D0.diag(1) = sinv.isherm() ? xD.Conjugate() : xD.View();
    Matrix<T> L0(sinv.LowerTri(UnitDiag));
#endif

    LowerTriMatrixView<T> L = sinv.LowerTri(UnitDiag);
    L.InvertSelf();

    T* Di = sinv.diag().ptr();
    const int Dstep = sinv.diag().step();
    T* xDi = xD.ptr();
    const size_t N = sinv.size();
    for(size_t i=0;i<N;) {
      if (i==N-1 || *xDi == T(0)) { // Then 1x1
	if (herm)
	  *Di = RealType(T)(1)/REAL(*Di);
	else
	  *Di = RealType(T)(1)/(*Di);
	Di+=Dstep; ++xDi; ++i;
      } else { // 2x2
	SymInvert_2x2<herm>(*Di,*(Di+Dstep),*xDi);
	Di+=2*Dstep; xDi+=2; i+=2;
      }
    }
#ifdef XDEBUG
    RealType(T) normdiff1 = Norm(L0 * L - T(1));
    Matrix<T> D2(DiagMatrixViewOf(sinv.diag()));
    D2.diag(-1) = xD;
    D2.diag(1) = sinv.isherm() ? xD.Conjugate() : xD.View();
    RealType(T) normdiff2 = Norm(D0 * D2 - T(1));
    RealType(T) kappa1 = Norm(L0) * Norm(L);
    RealType(T) kappa2 = Norm(D0) * Norm(D2);
    if (normdiff1 > 0.0001*kappa1*sinv.size() || 
	normdiff2 > 0.0001*kappa2*sinv.size()) {
      cerr<<"SymLUInverse - Invert Components\n";
      cerr<<"sinv input = "<<sinv0<<endl;
      cerr<<"sinv output = "<<sinv<<endl;
      cerr<<"A = "<<L0 * D0 * (sinv.issym() ? L0.Transpose() : L0.Adjoint())<<endl;
      cerr<<"xD = "<<xD<<endl;
      cerr<<"sinv = "<<sinv<<endl;
      cerr<<"D = "<<D0<<endl;
      cerr<<"D.inverse = "<<D2<<endl;
      cerr<<"L = "<<L0<<endl;
      cerr<<"L.inverse = "<<L<<endl;
      cerr<<"Norm(D*Dinv-1) = "<<normdiff1<<endl;
      cerr<<"Norm(L*Linv-1) = "<<normdiff2<<endl;
      cerr<<"kappa = "<<kappa1<<"  "<<kappa2<<endl;
      abort();
    }
#endif
  }

  template <bool herm, class T> inline void LDLt_AddXtDX(
      const SymMatrixView<T>& sinv, const GenMatrix<T>& X,
      const GenVector<T>& D, const GenVector<T>& xD)
    // sinv += Xt D X
  {
    TMVAssert(IsReal(T()) || herm == sinv.isherm());
    TMVAssert(sinv.uplo() == Lower);
    TMVAssert(xD.step() == 1);
    TMVAssert(sinv.size() == X.rowsize());
    TMVAssert(X.colsize() == xD.size()+1);
    TMVAssert(X.colsize() == D.size());

#ifdef XDEBUG
    Matrix<T> sinv0(sinv);
    Matrix<T> X0(X);
    Vector<T> D0(D);
    Vector<T> xD0(xD);
    Matrix<T> DD(DiagMatrixViewOf(D));
    if (D.size() > 1) {
      DD.diag(-1) = xD;
      DD.diag(1) = herm ? xD.Conjugate() : xD.View();
    }
    Matrix<T> m1 = sinv + (herm ? X.Adjoint() : X.Transpose())*DD*X;
#endif

    const size_t N = D.size();
    if (N == 1) {
      sinv += D(0) * ((herm ? X.row(0).Conjugate() : X.row(0)) ^ X.row(0));
    } else if ((N == 2 && xD(0) != T(0)) || N <= SYM_LU_BLOCKSIZE2) {
      //cerr<<"Non recurse, N>1\n";
      Matrix<T> DX = X;
      //cerr<<"D = "<<D<<endl;
      //cerr<<"xD = "<<xD<<endl;
      //cerr<<"DD = "<<DD<<endl;
      //cerr<<"X = "<<X<<endl;
      //cerr<<"DX = "<<DX<<endl;
      PseudoDiag_LMultEq<herm>(D,xD,DX.View());
      //cerr<<"DX = "<<DX<<endl;
      //cerr<<"DD*X = "<<DD*X<<endl;
      //SymMatrix<T,Lower,ColMajor> sinv2 = sinv;
      //SymMatrix<T,Lower,RowMajor> sinv3 = sinv;
      //sinv2.Zero();
      //sinv3.Zero();
      //SymMultMM<true>(T(1),herm?X.Adjoint():X.Transpose(),DX,sinv2.View());
      //SymMultMM<true>(T(1),herm?X.Adjoint():X.Transpose(),DX,sinv3.View());
      //cerr<<"sinv2 = "<<sinv2<<endl;
      //cerr<<"sinv3 = "<<sinv3<<endl;
      //cerr<<"Xt*DD*X = "<<(herm?X.Adjoint():X.Transpose())*DD*X<<endl;
      SymMultMM<true>(T(1),herm?X.Adjoint():X.Transpose(),DX,sinv);
      //cerr<<"sinv => "<<sinv<<endl;
      //cerr<<"sinv0 + sinv2 = "<<sinv0+sinv2<<endl;
      //cerr<<"sinv0 + sinv3 = "<<sinv0+sinv3<<endl;
      //cerr<<"m1 = "<<m1<<endl;
    } else {
      size_t No2 = N/2;
      if (xD(No2-1) != T(0)) ++No2;
      TMVAssert(xD(No2-1) == T(0));
      LDLt_AddXtDX<herm>(sinv,X.Rows(0,No2),
	  D.SubVector(0,No2),xD.SubVector(0,No2-1));
      LDLt_AddXtDX<herm>(sinv,X.Rows(No2,N),
	  D.SubVector(No2,N),xD.SubVector(No2,N-1));
    }
#ifdef XDEBUG
    if (Norm(m1-sinv) > 0.0001*Norm(m1)) {
      cerr<<"AddXtDX: \n";
      cerr<<"init sinv = "<<sinv0<<endl;
      cerr<<"xD = "<<xD0<<endl;
      cerr<<"D = "<<D0<<endl;
      cerr<<"DD = "<<DD<<endl;
      cerr<<"X = "<<X0<<endl;
      cerr<<"right answer = "<<m1<<endl;
      cerr<<"sinv => "<<sinv<<endl;
      abort();
    }
#endif
  }

  template <bool herm, class T> inline void LDLt_CombineInverse(
      const SymMatrixView<T>& sinv, const VectorView<T>& xDinv)
    // Do LDLt inverse recursively:
    // A = [ L00  0  ] [ D00  0  ] [ L00t L10t ]
    //     [ L10 L11 ] [  0  D11 ] [  0   L11t ]
    //
    // A^-1 = [ (L^-1)00t (L^-1)10t ] [ D00^-1  0   ] [ (L^-1)00    0     ]
    //        [     0     (L^-1)11t ] [  0   D11^-1 ] [ (L^-1)10 (L^-1)11 ]
    // Let L^-1 = X
    //      = [ X00t D00^-1 X00 + X10t D11^-1 X10    X10t D11^-1 X11  ]
    //        [       X11t D11^-1 X10                X11t D11^-1 X11  ]
  {
    TMVAssert(IsReal(T()) || herm == sinv.isherm());
    TMVAssert(sinv.uplo() == Lower);
    TMVAssert(sinv.size() == xDinv.size()+1);
    TMVAssert(xDinv.step() == 1);

    const size_t N = sinv.size();

#ifdef XDEBUG
    Matrix<T> sinv0(sinv);
    Matrix<T> Dinv(DiagMatrixViewOf(sinv.diag()));
    Dinv.diag(-1) = xDinv;
    Dinv.diag(1) = herm ? xDinv.Conjugate() : xDinv;
    Matrix<T> minv1 = 
      (herm ? sinv.LowerTri(UnitDiag).Adjoint() :
       sinv.LowerTri(UnitDiag).Transpose())
      * Dinv * sinv.LowerTri(UnitDiag);
#endif

    if (N > 1) {
      if (N == 2) {
	if (xDinv(0) != T(0)) {
	  sinv(1,0) = xDinv(0);
	} else {
	  if (herm) {
	    RealType(T) D11inv = REAL(sinv(1,1));
	    T X10 = sinv(1,0);
	    T temp = D11inv * X10;
	    sinv(1,0) = temp;
	    sinv(0,0) += REAL(CONJ(X10) * temp);
	  } else {
	    T D11inv = sinv(1,1);
	    T X10 = sinv(1,0);
	    T temp = D11inv * X10;
	    sinv(1,0) = temp;
	    sinv(0,0) += X10 * temp;
	  }
	}
      } else {
	size_t No2 = N/2;
	if (xDinv(No2-1) != T(0)) ++No2;
	TMVAssert(xDinv(No2-1) == T(0));
	TMVAssert(No2 < N);
	MatrixView<T> X10 = sinv.SubMatrix(No2,N,0,No2);
	LowerTriMatrixView<T> X11 = 
	  sinv.LowerTri(UnitDiag).SubTriMatrix(No2,N);
	VectorView<T> D = sinv.diag();

	LDLt_CombineInverse<herm>(sinv.SubSymMatrix(0,No2),
	    xDinv.SubVector(0,No2-1));
	LDLt_AddXtDX<herm>(sinv.SubSymMatrix(0,No2),X10,
	    D.SubVector(No2,N),xDinv.SubVector(No2,N-1));
	PseudoDiag_LMultEq<herm>(D.SubVector(No2,N),
	    xDinv.SubVector(No2,N-1),X10);
	if (herm) 
	  X10 = X11.Adjoint() * X10;
	else
	  X10 = X11.Transpose() * X10;
	LDLt_CombineInverse<herm>(sinv.SubSymMatrix(No2,N),
	    xDinv.SubVector(No2,N-1));
      }
    }
#ifdef XDEBUG
    if (Norm(minv1-sinv) > 0.0001*Norm(sinv)) {
      cerr<<"Combine inverse\n";
      cerr<<"Init sinv = "<<sinv0<<endl;
      cerr<<"xDinv = "<<xDinv<<endl;
      cerr<<"Right answer is "<<minv1<<endl;
      cerr<<"sinv = "<<sinv<<endl;
      cerr<<"Norm(diff) = "<<Norm(minv1-sinv)<<endl;
      abort();
    }
#endif
  }

  template <class T, class T1> void SymLU_Inverse(
      const GenSymMatrix<T>& LLx, const GenVector<T>& xD, const size_t* P,
      const SymMatrixView<T1>& sinv)
  {
#ifdef XTEST
    TMVAssert(LLx.HermOK());
#endif
    TMVAssert(sinv.size() == LLx.size());
    TMVAssert(IsReal(T()) || LLx.isherm() == sinv.isherm());
    TMVAssert(IsReal(T()) || LLx.issym() == sinv.issym());

#ifdef XDEBUG
    Matrix<T> D(DiagMatrixViewOf(LLx.diag()));
    D.diag(-1) = xD;
    D.diag(1) = LLx.isherm() ? xD.Conjugate() : xD.View();
    Matrix<T> A = LLx.LowerTri(UnitDiag) * D * LLx.UpperTri(UnitDiag);
    A.ReversePermuteRows(P);
    A.ReversePermuteCols(P);
    //cerr<<"Start SymLU_Inverse: \n";
    //cerr<<"sinv = "<<tmv::Type(sinv)<<endl;
    //cerr<<"A = "<<A<<endl;
    //cerr<<"L = "<<LLx.LowerTri(UnitDiag)<<endl;
    //cerr<<"D = "<<D<<endl;
#endif

    if (sinv.uplo() == Upper) 
      SymLU_Inverse(LLx,xD,P,
	  sinv.issym() ? sinv.Transpose() : sinv.Adjoint());
    else {
      if (sinv.size() > 0) {
	if (!(sinv.iscm() || sinv.isrm())) {
	  Matrix<T,ColMajor> minv1(LLx.size(),LLx.size());
	  SymMatrixView<T> sinv1 = LLx.isherm() ?
	    HermMatrixViewOf(minv1,Lower) :
	      SymMatrixViewOf(minv1,Lower);
	  SymLU_Inverse(LLx,xD,P,sinv1);
	  sinv = sinv1;
	}
	else {
	  sinv = LLx;
	  Vector<T1> xDinv = xD;
	  if (sinv.isherm()) {
	    LDLt_InvertComponents<true>(sinv,xDinv.View());
	    LDLt_CombineInverse<true>(sinv,xDinv.View());
	  } else {
	    LDLt_InvertComponents<false>(sinv,xDinv.View());
	    LDLt_CombineInverse<false>(sinv,xDinv.View());
	  }
	  sinv.ReversePermuteRowsCols(P);
	}
      }
    }

#ifdef XDEBUG
    RealType(T) normdiff = Norm(A * sinv - T(1));
    RealType(T) kappa = Norm(A) * Norm(sinv);
    if (normdiff > 0.0001*kappa*sinv.size()) {
      cerr<<"SymLUInverse\n";
      cerr<<"A = "<<A<<endl;
      cerr<<"LLx = "<<LLx<<endl;
      cerr<<"xD = "<<xD<<endl;
      cerr<<"sinv = "<<sinv<<endl;
      cerr<<"A.inverse = "<<A.Inverse()<<endl;
      cerr<<"A*sinv = "<<A*sinv<<endl;
      cerr<<"sinv*A = "<<sinv*A<<endl;
      cerr<<"Norm(A*sinv-1) = "<<Norm(A*sinv-T(1))<<endl;
      cerr<<"Norm(sinv*A-1) = "<<Norm(sinv*A-T(1))<<endl;
      cerr<<"kappa = "<<kappa<<endl;
      abort();
    }
#endif
#ifdef XTEST
    TMVAssert(sinv.HermOK());
#endif
  }

#define InstFile "TMV_SymLUDiv_C.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


