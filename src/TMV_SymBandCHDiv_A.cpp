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



#include "TMV_Blas.h"
#include "TMV_SymBandMatrix.h"
#include "TMV_SymBandCHDiv.h"
#include "TMV_SymBandCHDiv_A.h"
#include "TMV_VectorArith.h"
#include "TMV_MatrixArith.h"
#include "TMV_BandMatrixArith.h"
#include "TMV_SymMatrixArith.h"
#include "TMV_SymBandMatrixArith.h"

//#define XDEBUG

#ifdef XDEBUG
#include <iostream>
using std::cerr;
using std::endl;
#endif

namespace tmv {

  //
  // Decompose
  //

  template <bool cm, class T> inline void DoNonLapHermBandCH_Decompose(
      const SymBandMatrixView<T>& A)
  {
    // Cholesky decompostion for a banded Hermitian matrix follows the 
    // same structure as a regular Cholesky decomposition, but we 
    // take advantage of the fact that the column or row lengths are 
    // shorter.
    //
    TMVAssert(A.uplo() == Lower);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(IsReal(T()) || A.isherm());
    TMVAssert(A.iscm() || A.isrm());
    TMVAssert(cm == A.iscm());
    const size_t N = A.size();
#ifdef XDEBUG
    Matrix<T> A0(A);
#endif

    const VectorView<RealType(T)> Adiag = A.Real().diag();
    const size_t nlo = A.nlo();
    
    if (cm) {
      RealType(T)* Ajj= Adiag.ptr();
      const int ds = Adiag.step();
      size_t endcol = nlo+1;
      for(size_t j=0;j<N-1;++j,Ajj+=ds) {
	if (*Ajj <= RealType(T)(0)) 
	  throw NonPosDefHermBandMatrix<T>(A);
	*Ajj = SQRT(*Ajj);
	A.col(j,j+1,endcol) /= *Ajj;
	A.SubSymMatrix(j+1,endcol) -= 
	  A.col(j,j+1,endcol) ^ A.col(j,j+1,endcol).Conjugate();
	if (endcol < N) ++endcol;
      }
      if (*Ajj <= RealType(T)(0)) 
	throw NonPosDefHermBandMatrix<T>(A);
      *Ajj = SQRT(*Ajj);
    } else {
      RealType(T)* Aii = Adiag.ptr();
      const int ds = Adiag.step();
      size_t startrow = 0;
      if (*Aii <= RealType(T)(0))
	throw NonPosDefHermBandMatrix<T>(A);
      *Aii = SQRT(*Aii);
      for(size_t i=1;i<N;++i) {
	if (i > nlo) ++startrow;
	Aii+=ds;
	A.row(i,startrow,i) %= A.SubSymBandMatrix(0,i).LowerBand().Adjoint();
	*Aii -= NormSq(A.row(i,startrow,i));
	if (*Aii <= RealType(T)(0)) 
	  throw NonPosDefHermBandMatrix<T>(A);
	*Aii = SQRT(*Aii);
      }
    }

#ifdef XDEBUG
    BandMatrix<T> L = A.LowerBand();
    BandMatrix<T> A2 = L * L.Adjoint();
    RealType(T) normll = SQR(Norm(L));
    if (Norm(A2-A0) > 0.001*normll) {
      cerr<<"HermBandCHDecomp: A = "<<Type(A)<<"  "<<A0<<endl;
      cerr<<"Done: A = "<<A<<endl;
      cerr<<"L = "<<L<<endl;
      cerr<<"Lt = "<<L.Adjoint()<<endl;
      cerr<<"L*Lt = "<<A2<<endl;
      cerr<<"norm(diff) = "<<Norm(A2-A0);
      cerr<<"  Norm(L)^2 = "<<normll<<endl;
      abort();
    }
#endif
  }

  template <class T> inline void NonLapHermBandCH_Decompose(
      const SymBandMatrixView<T>& A)
  {
    TMVAssert(A.iscm() || A.isrm());
    if (A.iscm()) DoNonLapHermBandCH_Decompose<true>(A);
    else DoNonLapHermBandCH_Decompose<false>(A);
  }

  template <class T> inline void NonLapHermTriDiagCH_Decompose(
      const SymBandMatrixView<T>& A)
  {
    // For tridiagonal Hermitian band matrices, the sqrt can 
    // become a significant fraction of the calculation.
    // This is unnecessary if we instead decompose A into L D Lt
    // where L is unit diagonal.
    //
    // In this case, the equations become:
    //
    // A = L0 D L0t
    //
    // ( A00 Ax0t ) = (  1  0 ) ( D0  0 ) ( 1 Lx0t )
    // ( Ax0 Axx  ) = ( Lx0 1 ) (  0 Dx ) ( 0   1  )
    //              = (  1  0 ) ( D0  D0 Lx0t )
    //              = ( Lx0 1 ) ( 0   Dx      )
    //
    // The equations from this are:
    //
    // A00 = D0
    // Ax0 = Lx0 D0
    // Axx = Lx0 D0 Lx0t + Axx'
    //
    // where Axx' = Dx is the sub-matrix for the next step.
    //
    TMVAssert(A.uplo() == Lower);
    TMVAssert(A.isdm());
    TMVAssert(A.ct() == NonConj);
    TMVAssert(IsReal(T()) || A.isherm());
    TMVAssert(A.nlo() == 1);
    const size_t N = A.size();
#ifdef XDEBUG
    Matrix<T> A0(A);
#endif

    RealType(T)* Dj = A.Real().ptr();
    T* Lj = A.diag(-1).ptr();
    
    for(size_t j=0;j<N-1;++j) {
      if (*Dj <= RealType(T)(0)) throw NonPosDefHermBandMatrix<T>(A);
      T Ax0 = *Lj;
      *Lj /= *Dj;
      if (IsReal(T())) ++Dj; else Dj+=2; 
      *Dj -= REAL(CONJ(*Lj) * Ax0);
      ++Lj;
    }
    if (*Dj <= RealType(T)(0)) throw NonPosDefHermBandMatrix<T>(A);

#ifdef XDEBUG
    BandMatrix<T> L = A.LowerBand();
    L.diag().SetAllTo(T(1));
    DiagMatrix<T> D = DiagMatrixViewOf(A.diag());
    BandMatrix<T> A2 = L * BandMatrixViewOf(D) * L.Adjoint();
    RealType(T) normldl = SQR(Norm(L))*Norm(D);
    if (Norm(A2-A0) > 0.001*normldl) {
      cerr<<"HermBandCHDecomp: A = "<<Type(A)<<"  "<<A0<<endl;
      cerr<<"Done: A = "<<A<<endl;
      cerr<<"L = "<<L<<endl;
      cerr<<"D = "<<D<<endl;
      cerr<<"Lt = "<<L.Adjoint()<<endl;
      cerr<<"L*D*Lt = "<<A2<<endl;
      cerr<<"norm(diff) = "<<Norm(A2-A0);
      cerr<<"  Norm(L)^2*Norm(D) = "<<normldl<<endl;
      abort();
    }
#endif
  }

#ifdef LAP
  template <class T> inline void LapHermBandCH_Decompose(
      const SymBandMatrixView<T>& A)
  { NonLapHermBandCH_Decompose(A); }
#ifdef INST_DOUBLE
  template <> inline void LapHermBandCH_Decompose(
      const SymBandMatrixView<double>& A)
  {
    TMVAssert(A.uplo() == Lower);
    TMVAssert(A.iscm() || A.isrm());
    TMVAssert(A.ct()==NonConj);

    int n = A.size();
    int kl = A.nlo();
    int lda = (A.iscm() ? A.stepj() : A.stepi()) + 1;
    double* Aptr = A.ptr();
    if (A.isrm()) Aptr -= A.nlo();

    LAPNAME(dpbtrf) (LAPCM A.iscm()?LAPCH_LO:LAPCH_UP, LAPV(n), LAPV(kl),
	LAPP(Aptr),LAPV(lda) LAPINFO LAP1);
    if (Lap_info > 0) throw NonPosDefHermBandMatrix<double>(A);
    LAP_Results("dpbtrf");
  }
  template <> inline void LapHermBandCH_Decompose(
      const SymBandMatrixView<std::complex<double> >& A)
  {
    TMVAssert(A.uplo() == Lower);
    TMVAssert(A.iscm() || A.isrm());
    TMVAssert(A.ct()==NonConj);
    TMVAssert(A.isherm());

    int n = A.size();
    int kl = A.nlo();
    int lda = (A.iscm() ? A.stepj() : A.stepi()) + 1;
    std::complex<double>* Aptr = A.ptr();
    if (A.isrm()) Aptr -= A.nlo();

    LAPNAME(zpbtrf) (LAPCM A.iscm()?LAPCH_LO:LAPCH_UP, LAPV(n), LAPV(kl),
	LAPP(Aptr),LAPV(lda) LAPINFO LAP1);
    if (Lap_info > 0) throw NonPosDefHermBandMatrix<std::complex<double> >(A);
    LAP_Results("zpbtrf");
  }
#endif
#ifdef INST_FLOAT
  template <> inline void LapHermBandCH_Decompose(
      const SymBandMatrixView<float>& A)
  {
    TMVAssert(A.uplo() == Lower);
    TMVAssert(A.iscm() || A.isrm());
    TMVAssert(A.ct()==NonConj);

    int n = A.size();
    int kl = A.nlo();
    int lda = (A.iscm() ? A.stepj() : A.stepi()) + 1;
    float* Aptr = A.ptr();
    if (A.isrm()) Aptr -= A.nlo();

    LAPNAME(spbtrf) (LAPCM A.iscm()?LAPCH_LO:LAPCH_UP, LAPV(n), LAPV(kl),
	LAPP(Aptr),LAPV(lda) LAPINFO LAP1);
    if (Lap_info > 0) throw NonPosDefHermBandMatrix<float>(A);
    LAP_Results("spbtrf");
  }
  template <> inline void LapHermBandCH_Decompose(
      const SymBandMatrixView<std::complex<float> >& A)
  {
    TMVAssert(A.uplo() == Lower);
    TMVAssert(A.iscm() || A.isrm());
    TMVAssert(A.ct()==NonConj);
    TMVAssert(A.isherm());

    int n = A.size();
    int kl = A.nlo();
    int lda = (A.iscm() ? A.stepj() : A.stepi()) + 1;
    std::complex<float>* Aptr = A.ptr();
    if (A.isrm()) Aptr -= A.nlo();

    LAPNAME(cpbtrf) (LAPCM A.iscm()?LAPCH_LO:LAPCH_UP, LAPV(n), LAPV(kl),
	LAPP(Aptr),LAPV(lda) LAPINFO LAP1);
    if (Lap_info > 0) throw NonPosDefHermBandMatrix<std::complex<float> >(A);
    LAP_Results("cpbtrf");
  }
#endif 

  template <class T> inline void LapHermTriDiagCH_Decompose(
      const SymBandMatrixView<T>& A)
  { NonLapHermTriDiagCH_Decompose(A); }
#ifdef INST_DOUBLE
  template <> inline void LapHermTriDiagCH_Decompose(
      const SymBandMatrixView<double>& A)
  {
    TMVAssert(A.uplo() == Lower);
    TMVAssert(A.isdm() && A.nlo()==1);
    TMVAssert(A.ct()==NonConj);

    int n = A.size();
    LAPNAME(dpttrf) (LAPCM LAPV(n), LAPP(A.ptr()),
	LAPP(A.ptr()+A.stepi()) LAPINFO);
    if (Lap_info > 0) throw NonPosDefHermBandMatrix<double>(A);
    LAP_Results("dpttrf");
  }
  template <> inline void LapHermTriDiagCH_Decompose(
      const SymBandMatrixView<std::complex<double> >& A)
  {
    TMVAssert(A.uplo() == Lower);
    TMVAssert(A.isdm() && A.nlo()==1);
    TMVAssert(A.ct()==NonConj);
    TMVAssert(A.isherm());

    int n = A.size();
    // Need A.diag.Real to have unit step...
    Vector<double> Ad = A.diag().Real();
    LAPNAME(zpttrf) (LAPCM LAPV(n), LAPP(Ad.ptr()),
	LAPP(A.ptr()+A.stepi()) LAPINFO);
    A.diag().Real() = Ad;
    if (Lap_info > 0) throw NonPosDefHermBandMatrix<std::complex<double> >(A);
    LAP_Results("zpttrf");
  }
#endif
#ifdef INST_FLOAT
  template <> inline void LapHermTriDiagCH_Decompose(
      const SymBandMatrixView<float>& A)
  {
    TMVAssert(A.uplo() == Lower);
    TMVAssert(A.isdm() && A.nlo()==1);
    TMVAssert(A.ct()==NonConj);

    int n = A.size();
    LAPNAME(spttrf) (LAPCM LAPV(n), LAPP(A.ptr()),
	LAPP(A.ptr()+A.stepi()) LAPINFO LAP1);
    if (Lap_info > 0) throw NonPosDefHermBandMatrix<float>(A);
    LAP_Results("spttrf");
  }
  template <> inline void LapHermTriDiagCH_Decompose(
      const SymBandMatrixView<std::complex<float> >& A)
  {
    TMVAssert(A.uplo() == Lower);
    TMVAssert(A.isdm() && A.nlo()==1);
    TMVAssert(A.ct()==NonConj);
    TMVAssert(A.isherm());

    int n = A.size();
    Vector<float> Ad = A.diag().Real();
    LAPNAME(cpttrf) (LAPCM LAPV(n), LAPP(Ad.ptr()),
	LAPP(A.ptr()+A.stepi()) LAPINFO);
    A.diag().Real() = Ad;
    if (Lap_info > 0) throw NonPosDefHermBandMatrix<std::complex<float> >(A);
    LAP_Results("cpttrf");
  }
#endif 
#endif // LAP

  template <class T> void HermBandCH_Decompose(
      const SymBandMatrixView<T>& A)
  {
#ifdef XTEST
    TMVAssert(A.HermOK());
#endif
    TMVAssert(IsReal(T()) || A.isherm());
    TMVAssert((A.nlo() == 1 && A.isdm()) || A.iscm() || A.isrm());

    if (A.uplo() == Upper) HermBandCH_Decompose(A.Adjoint());
    else if (A.isconj()) HermBandCH_Decompose(A.Conjugate());
    else if (A.size() > 0) {
#ifdef LAP 
      if (A.isdm()) 
	LapHermTriDiagCH_Decompose(A);
      else
	LapHermBandCH_Decompose(A);
#else
      if (A.isdm())
	NonLapHermTriDiagCH_Decompose(A);
      else
	NonLapHermBandCH_Decompose(A);
#endif
    }
#ifdef XTEST
    TMVAssert(A.HermOK());
#endif
  }

#define InstFile "TMV_SymBandCHDiv_A.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


