///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// The Template Matrix/Vector Library for C++ was created by Mike Jarvis     //
// Copyright (C) 1998 - 2009                                                 //
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
// along with this program in the file LICENSE.                              //
//                                                                           //
// If not, write to:                                                         //
// The Free Software Foundation, Inc.                                        //
// 51 Franklin Street, Fifth Floor,                                          //
// Boston, MA  02110-1301, USA.                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


//#define XDEBUG


#include "TMV_Blas.h"
#include "TMV_SymBandCHDiv.h"
#include "tmv/TMV_SymBandCHD.h"
#include "tmv/TMV_SymBandMatrix.h"
#include "tmv/TMV_BandMatrixArith.h"
#include "tmv/TMV_SymMatrixArith.h"

#ifdef NOTHROW
#include <iostream>
#endif

#ifdef XDEBUG
#include "tmv/TMV_BandMatrixArith.h"
#include "tmv/TMV_MatrixArith.h"
#include "tmv/TMV_DiagMatrix.h"
#include "tmv/TMV_DiagMatrixArith.h"
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#endif

namespace tmv {

    //
    // Decompose
    //

    template <bool cm, class T> 
    static void DoNonLapCH_Decompose(const SymBandMatrixView<T>& A)
    {
        // Cholesky decompostion for a banded Hermitian matrix follows the 
        // same structure as a regular Cholesky decomposition, but we 
        // take advantage of the fact that the column or row lengths are 
        // shorter.
        //
        TMVAssert(A.uplo() == Lower);
        TMVAssert(A.ct() == NonConj);
        TMVAssert(isReal(T()) || A.isherm());
        TMVAssert(A.iscm() || A.isrm());
        TMVAssert(cm == A.iscm());
        const int N = A.size();
#ifdef XDEBUG
        Matrix<T> A0(A);
        //cout<<"CHDecompose:\n";
        //cout<<"A = "<<TMV_Text(A)<<"  "<<A<<endl;
#endif

        const VectorView<TMV_RealType(T)> Adiag = A.diag().realPart();
        const int nlo = A.nlo();

        if (nlo == 0) {
            TMV_RealType(T)* Ajj= Adiag.ptr();
            const int ds = Adiag.step();
            for(int j=0;j<N;++j,Ajj+=ds) {
#ifdef TMVFLDEBUG
                TMVAssert(Ajj >= A.realPart().first);
                TMVAssert(Ajj < A.realPart().last);
#endif
                if (*Ajj <= TMV_RealType(T)(0))  {
#ifdef NOTHROW
                    std::cerr<<"Non Posdef HermBandMatrix found\n"; 
                    exit(1); 
#else
                    throw NonPosDefHermBandMatrix<T>(A);
#endif
                }
                *Ajj = TMV_SQRT(*Ajj);
            }
        } else if (cm) {
            TMV_RealType(T)* Ajj= Adiag.ptr();
            const int ds = Adiag.step();
            int endcol = nlo+1;
            for(int j=0;j<N-1;++j,Ajj+=ds) {
#ifdef TMVFLDEBUG
                TMVAssert(Ajj >= A.realPart().first);
                TMVAssert(Ajj < A.realPart().last);
#endif
                if (*Ajj <= TMV_RealType(T)(0))  {
#ifdef NOTHROW
                    std::cerr<<"Non Posdef HermBandMatrix found\n"; 
                    exit(1); 
#else
                    throw NonPosDefHermBandMatrix<T>(A);
#endif
                }
                *Ajj = TMV_SQRT(*Ajj);
                A.col(j,j+1,endcol) /= *Ajj;
                A.subSymMatrix(j+1,endcol) -= 
                    A.col(j,j+1,endcol) ^ A.col(j,j+1,endcol).conjugate();
                if (endcol < N) ++endcol;
            }
#ifdef TMVFLDEBUG
            TMVAssert(Ajj >= A.realPart().first);
            TMVAssert(Ajj < A.realPart().last);
#endif
            if (*Ajj <= TMV_RealType(T)(0))  {
#ifdef NOTHROW
                std::cerr<<"Non Posdef HermBandMatrix found\n"; 
                exit(1); 
#else
                throw NonPosDefHermBandMatrix<T>(A);
#endif
            }
            *Ajj = TMV_SQRT(*Ajj);
        } else {
            TMV_RealType(T)* Aii = Adiag.ptr();
            const int ds = Adiag.step();
            int startrow = 0;
#ifdef TMVFLDEBUG
            TMVAssert(Aii >= A.realPart().first);
            TMVAssert(Aii < A.realPart().last);
#endif
            if (*Aii <= TMV_RealType(T)(0)) {
#ifdef NOTHROW
                std::cerr<<"Non Posdef HermBandMatrix found\n"; 
                exit(1); 
#else
                throw NonPosDefHermBandMatrix<T>(A);
#endif
            }
            *Aii = TMV_SQRT(*Aii);
            for(int i=1;i<N;++i) {
                if (i > nlo) ++startrow;
                Aii+=ds;
                A.row(i,startrow,i) %= 
                    A.subSymBandMatrix(startrow,i).lowerBand().adjoint();
#ifdef TMVFLDEBUG
                TMVAssert(Aii >= A.realPart().first);
                TMVAssert(Aii < A.realPart().last);
#endif
                *Aii -= NormSq(A.row(i,startrow,i));
                if (*Aii <= TMV_RealType(T)(0))  {
#ifdef NOTHROW
                    std::cerr<<"Non Posdef HermBandMatrix found\n"; 
                    exit(1); 
#else
                    throw NonPosDefHermBandMatrix<T>(A);
#endif
                }
                *Aii = TMV_SQRT(*Aii);
            }
        }

#ifdef XDEBUG
        //cout<<"Done CHDecompose\n";
        BandMatrix<T> L = A.lowerBand();
        BandMatrix<T> A2 = L * L.adjoint();
        TMV_RealType(T) normll = TMV_SQR(Norm(L));
        if (Norm(A2-A0) > 0.001*normll) {
            cerr<<"CH_Decompose: A = "<<TMV_Text(A)<<"  "<<A0<<endl;
            cerr<<"Done: A = "<<A<<endl;
            cerr<<"L = "<<L<<endl;
            cerr<<"Lt = "<<L.adjoint()<<endl;
            cerr<<"L*Lt = "<<A2<<endl;
            cerr<<"Norm(diff) = "<<Norm(A2-A0);
            cerr<<"  Norm(L)^2 = "<<normll<<endl;
            abort();
        }
#endif
    }

    template <class T> 
    static inline void NonLapCH_Decompose(const SymBandMatrixView<T>& A)
    {
        TMVAssert(A.iscm() || A.isrm());
        if (A.iscm()) DoNonLapCH_Decompose<true>(A);
        else DoNonLapCH_Decompose<false>(A);
    }

    template <bool dm, class T> 
    static void NonLapLDL_Decompose(const SymBandMatrixView<T>& A)
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
        TMVAssert(dm == A.isdm());
        TMVAssert(A.ct() == NonConj);
        TMVAssert(A.isherm());
        TMVAssert(A.nlo() == 1);
        const int N = A.size();
#ifdef XDEBUG
        Matrix<T> A0(A);
        //cout<<"LDLDecompose:\n";
        //cout<<"A = "<<TMV_Text(A)<<"  "<<A<<endl;
#endif

        const int Dstep = dm ? 0 : A.diag().realPart().step();
        const int Lstep = dm ? 0 : A.diag().step();
        TMV_RealType(T)* Dj = A.realPart().ptr();
        T* Lj = A.diag(-1).ptr();

        for(int j=0;j<N-1;++j) {
#ifdef TMVFLDEBUG
            TMVAssert(Lj >= A.first);
            TMVAssert(Lj < A.last);
            TMVAssert(Dj >= A.realPart().first);
            TMVAssert(Dj < A.realPart().last);
#endif
            T Ax0 = *Lj;
            if (*Dj == TMV_RealType(T)(0))  {
#ifdef NOTHROW
                std::cerr<<"Non Posdef HermBandMatrix found\n"; 
                exit(1); 
#else
                throw NonPosDefHermBandLDL<T>(A);
#endif
            }
            *Lj /= *Dj;
            if (dm) { if (isReal(T())) ++Dj; else Dj+=2; }
            else Dj += Dstep;
#ifdef TMVFLDEBUG
            TMVAssert(Dj >= A.realPart().first);
            TMVAssert(Dj < A.realPart().last);
#endif
            *Dj -= TMV_REAL(TMV_CONJ(*Lj) * Ax0);
            if (dm) ++Lj;
            else Lj += Lstep;
        }
        if (*Dj == TMV_RealType(T)(0))  {
#ifdef NOTHROW
            std::cerr<<"Non Posdef HermBandMatrix found\n"; 
            exit(1); 
#else
            throw NonPosDefHermBandLDL<T>(A);
#endif
        }

#ifdef XDEBUG
        //cout<<"Done LDLDecompose\n";
        BandMatrix<T> L = A.lowerBand();
        L.diag().SetAllTo(T(1));
        DiagMatrix<T> D = DiagMatrixViewOf(A.diag());
        BandMatrix<T> A2 = L * D * L.adjoint();
        TMV_RealType(T) normldl = TMV_SQR(Norm(L))*Norm(D);
        //cout<<"Done: A = "<<A<<endl;
        //cout<<"L = "<<L<<endl;
        //cout<<"D = "<<D<<endl;
        //cout<<"A2 = "<<A2<<endl;
        //cout<<"cf A0 = "<<A0<<endl;
        if (Norm(A2-A0) > 0.001*normldl) {
            cerr<<"LDL_Decompose: A = "<<TMV_Text(A)<<"  "<<A0<<endl;
            cerr<<"Done: A = "<<A<<endl;
            cerr<<"L = "<<L<<endl;
            cerr<<"D = "<<D<<endl;
            cerr<<"Lt = "<<L.adjoint()<<endl;
            cerr<<"L*D*Lt = "<<A2<<endl;
            cerr<<"Norm(diff) = "<<Norm(A2-A0);
            cerr<<"  Norm(L)^2*Norm(D) = "<<normldl<<endl;
            abort();
        }
#endif
    }

    // Same thing, but A is symmetric, rather than hermitian
    template <bool dm, class T> 
    static void SymLDL_Decompose(const SymBandMatrixView<T>& A)
    {
        TMVAssert(A.uplo() == Lower);
        TMVAssert(dm == A.isdm());
        TMVAssert(A.ct() == NonConj);
        TMVAssert(A.issym());
        TMVAssert(A.nlo() == 1);
        const int N = A.size();
#ifdef XDEBUG
        Matrix<T> A0(A);
        //cout<<"SymLDLDecompose:\n";
        //cout<<"A = "<<TMV_Text(A)<<"  "<<A<<endl;
#endif

        int step = dm ? 0 : A.diagstep();
        T* Dj = A.ptr();
        T* Lj = A.diag(-1).ptr();

        for(int j=0;j<N-1;++j) {
            T Ax0 = *Lj;
#ifdef TMVFLDEBUG
            TMVAssert(Lj >= A.first);
            TMVAssert(Lj < A.last);
#endif
            if (*Dj == T(0))  {
#ifdef NOTHROW
                std::cerr<<"Non Posdef HermBandMatrix found\n"; 
                exit(1); 
#else
                throw NonPosDefSymBandLDL<T>(A);
#endif
            }
            *Lj /= *Dj;
            if (dm) ++Dj;
            else Dj += step;
#ifdef TMVFLDEBUG
            TMVAssert(Dj >= A.first);
            TMVAssert(Dj < A.last);
#endif
            *Dj -= *Lj * Ax0;
            if (dm) ++Lj;
            else Lj += step;
        }
        if (*Dj == T(0))  {
#ifdef NOTHROW
            std::cerr<<"Non Posdef HermBandMatrix found\n"; 
            exit(1); 
#else
            throw NonPosDefSymBandLDL<T>(A);
#endif
        }

#ifdef XDEBUG
        //cout<<"Done SymLDLDecom\n";
        BandMatrix<T> L = A.lowerBand();
        L.diag().SetAllTo(T(1));
        DiagMatrix<T> D = DiagMatrixViewOf(A.diag());
        BandMatrix<T> A2 = L * D * L.transpose();
        TMV_RealType(T) normldl = TMV_SQR(Norm(L))*Norm(D);
        //cout<<"L = "<<L<<endl;
        //cout<<"D = "<<D<<endl;
        //cout<<"A2 = "<<A2<<endl;
        //cout<<"cf. A0 = "<<A0<<endl;
        if (Norm(A2-A0) > 0.001*normldl) {
            cerr<<"SymLDLDecomp: A = "<<TMV_Text(A)<<"  "<<A0<<endl;
            cerr<<"Done: A = "<<A<<endl;
            cerr<<"L = "<<L<<endl;
            cerr<<"D = "<<D<<endl;
            cerr<<"LT = "<<L.transpose()<<endl;
            cerr<<"L*D*LT = "<<A2<<endl;
            cerr<<"Norm(diff) = "<<Norm(A2-A0);
            cerr<<"  Norm(L)^2*Norm(D) = "<<normldl<<endl;
            abort();
        }
#endif
    }

#ifdef LAP
    template <class T> 
    static inline void LapCH_Decompose(const SymBandMatrixView<T>& A)
    { NonLapCH_Decompose(A); }
    template <class T> 
    static inline void LapLDL_Decompose(const SymBandMatrixView<T>& A)
    { NonLapLDL_Decompose<true>(A); }
#ifdef INST_DOUBLE
    template <> 
    void LapCH_Decompose(const SymBandMatrixView<double>& A)
    {
        TMVAssert(A.uplo() == Lower);
        TMVAssert(A.iscm() || A.isrm());
        TMVAssert(A.ct()==NonConj);

        int n = A.size();
        int kl = A.nlo();
        int lda = (A.iscm() ? A.stepj() : A.stepi()) + 1;
        double* Aptr = A.ptr();
        if (A.isrm()) Aptr -= A.nlo();

        LAPNAME(dpbtrf) (
            LAPCM A.iscm()?LAPCH_LO:LAPCH_UP, LAPV(n), LAPV(kl),
            LAPP(Aptr),LAPV(lda) LAPINFO LAP1);
        if (Lap_info > 0)  {
#ifdef NOTHROW
            std::cerr<<"Non Posdef HermBandMatrix found\n"; 
            exit(1); 
#else
            throw NonPosDefHermBandMatrix<double>(A);
#endif
        }
        LAP_Results("dpbtrf");
    }
    template <> 
    void LapCH_Decompose(const SymBandMatrixView<std::complex<double> >& A)
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

        LAPNAME(zpbtrf) (
            LAPCM A.iscm()?LAPCH_LO:LAPCH_UP, LAPV(n), LAPV(kl),
            LAPP(Aptr),LAPV(lda) LAPINFO LAP1);
        if (Lap_info > 0)  {
#ifdef NOTHROW
            std::cerr<<"Non Posdef HermBandMatrix found\n"; 
            exit(1); 
#else
            throw NonPosDefHermBandMatrix<std::complex<double> >(A);
#endif
        }
        LAP_Results("zpbtrf");
    }
    template <> 
    void LapLDL_Decompose(const SymBandMatrixView<double>& A)
    {
        TMVAssert(A.uplo() == Lower);
        TMVAssert(A.isdm() && A.nlo()==1);
        TMVAssert(A.ct()==NonConj);

        int n = A.size();
        LAPNAME(dpttrf) (
            LAPCM LAPV(n), LAPP(A.ptr()),
            LAPP(A.ptr()+A.stepi()) LAPINFO);
        if (Lap_info > 0)  {
#ifdef NOTHROW
            std::cerr<<"Non Posdef HermBandMatrix found\n"; 
            exit(1); 
#else
            throw NonPosDefHermBandMatrix<double>(A);
#endif
        }
        LAP_Results("dpttrf");
    }
    template <> 
    void LapLDL_Decompose(const SymBandMatrixView<std::complex<double> >& A)
    {
        TMVAssert(A.uplo() == Lower);
        TMVAssert(A.isdm() && A.nlo()==1);
        TMVAssert(A.ct()==NonConj);
        TMVAssert(A.isherm());

        int n = A.size();
        // Need A.diag.Real to have unit step...
        Vector<double> Ad = A.diag().realPart();
        LAPNAME(zpttrf) (
            LAPCM LAPV(n), LAPP(Ad.ptr()),
            LAPP(A.ptr()+A.stepi()) LAPINFO);
        A.diag().realPart() = Ad;
        if (Lap_info > 0)  {
#ifdef NOTHROW
            std::cerr<<"Non Posdef HermBandMatrix found\n"; 
            exit(1); 
#else
            throw NonPosDefHermBandMatrix<std::complex<double> >(A);
#endif
        }
        LAP_Results("zpttrf");
    }
#endif
#ifdef INST_FLOAT
    template <> 
    void LapLDL_Decompose(const SymBandMatrixView<float>& A)
    {
        TMVAssert(A.uplo() == Lower);
        TMVAssert(A.isdm() && A.nlo()==1);
        TMVAssert(A.ct()==NonConj);

        int n = A.size();
        LAPNAME(spttrf) (
            LAPCM LAPV(n), LAPP(A.ptr()),
            LAPP(A.ptr()+A.stepi()) LAPINFO LAP1);
        if (Lap_info > 0)  {
#ifdef NOTHROW
            std::cerr<<"Non Posdef HermBandMatrix found\n"; 
            exit(1); 
#else
            throw NonPosDefHermBandMatrix<float>(A);
#endif
        }
        LAP_Results("spttrf");
    }
    template <> 
    void LapLDL_Decompose(const SymBandMatrixView<std::complex<float> >& A)
    {
        TMVAssert(A.uplo() == Lower);
        TMVAssert(A.isdm() && A.nlo()==1);
        TMVAssert(A.ct()==NonConj);
        TMVAssert(A.isherm());

        int n = A.size();
        Vector<float> Ad = A.diag().realPart();
        LAPNAME(cpttrf) (
            LAPCM LAPV(n), LAPP(Ad.ptr()),
            LAPP(A.ptr()+A.stepi()) LAPINFO);
        A.diag().realPart() = Ad;
        if (Lap_info > 0)  {
#ifdef NOTHROW
            std::cerr<<"Non Posdef HermBandMatrix found\n"; 
            exit(1); 
#else
            throw NonPosDefHermBandMatrix<std::complex<float> >(A);
#endif
        }
        LAP_Results("cpttrf");
    }
    template <> 
    void LapCH_Decompose(const SymBandMatrixView<float>& A)
    {
        TMVAssert(A.uplo() == Lower);
        TMVAssert(A.iscm() || A.isrm());
        TMVAssert(A.ct()==NonConj);

        int n = A.size();
        int kl = A.nlo();
        int lda = (A.iscm() ? A.stepj() : A.stepi()) + 1;
        float* Aptr = A.ptr();
        if (A.isrm()) Aptr -= A.nlo();

        LAPNAME(spbtrf) (
            LAPCM A.iscm()?LAPCH_LO:LAPCH_UP, LAPV(n), LAPV(kl),
            LAPP(Aptr),LAPV(lda) LAPINFO LAP1);
        if (Lap_info > 0)  {
#ifdef NOTHROW
            std::cerr<<"Non Posdef HermBandMatrix found\n"; 
            exit(1); 
#else
            throw NonPosDefHermBandMatrix<float>(A);
#endif
        }
        LAP_Results("spbtrf");
    }
    template <> 
    void LapCH_Decompose(const SymBandMatrixView<std::complex<float> >& A)
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

        LAPNAME(cpbtrf) (
            LAPCM A.iscm()?LAPCH_LO:LAPCH_UP, LAPV(n), LAPV(kl),
            LAPP(Aptr),LAPV(lda) LAPINFO LAP1);
        if (Lap_info > 0)  {
#ifdef NOTHROW
            std::cerr<<"Non Posdef HermBandMatrix found\n"; 
            exit(1); 
#else
            throw NonPosDefHermBandMatrix<std::complex<float> >(A);
#endif
        }
        LAP_Results("cpbtrf");
    }
#endif 
#endif // LAP

    template <class T> 
    void LDL_Decompose(const SymBandMatrixView<T>& A)
    {
        TMVAssert(A.nlo() == 1);

        if (A.uplo() == Upper) LDL_Decompose(A.adjoint());
        else if (A.isconj()) LDL_Decompose(A.conjugate());
        else if (A.size() > 0) {
            if (A.isherm()) {
#ifdef LAP 
                if (A.isdm()) LapLDL_Decompose(A);
                else {
                    HermBandMatrix<T,Lower,DiagMajor> A2 = A;
                    LapLDL_Decompose(A2.view());
                    A = A2;
                }
#else
                if (A.isdm()) NonLapLDL_Decompose<true>(A);
                else NonLapLDL_Decompose<false>(A);
#endif
            } else {
                if (A.isdm()) SymLDL_Decompose<true>(A);
                else SymLDL_Decompose<false>(A);
            }
        }
    }

    template <class T> 
    void CH_Decompose(const SymBandMatrixView<T>& A)
    {
        TMVAssert(isReal(T()) || A.isherm());
        TMVAssert(A.iscm() || A.isrm());

        if (A.uplo() == Upper) CH_Decompose(A.adjoint());
        else if (A.isconj()) CH_Decompose(A.conjugate());
        else if (A.size() > 0) {
#ifdef LAP 
            LapCH_Decompose(A);
#else
            NonLapCH_Decompose(A);
#endif
        }
    }

#define InstFile "TMV_SymBandCHDecompose.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv

