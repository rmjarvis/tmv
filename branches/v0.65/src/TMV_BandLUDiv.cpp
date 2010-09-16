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
#include "TMV_BandLUDiv.h"
#include "tmv/TMV_BandMatrix.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_MatrixArith.h"

#ifdef XDEBUG
#include "tmv/TMV_TriDiv.h"
#include "tmv/TMV_TriMatrixArith.h"
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#endif

namespace tmv {

    //
    // Unpack packed PL
    //

    template <class T, class T1> 
    void LU_PackedPL_Unpack(
        const GenBandMatrix<T1>& LUx, const int* p,
        const LowerTriMatrixView<T>& L)
    {
        TMVAssert(L.isunit());
        TMVAssert(LUx.colsize() == L.size());
        TMVAssert(LUx.rowsize() == L.size());

        // Unpack a packed PL BandMatrix
        //
        // L = (  1   0 ) ( 1  0   0 ) ( I  0   0 )
        //     ( P0L0 I ) ( 0  1   0 ) ( 0  1   0 ) ...
        //                ( 0 P1L1 I ) ( 0 P2L2 I )

        int N = L.size();
        int nlo = LUx.nlo();
        if (nlo == 0) L.setToIdentity();
        else {
            L.setZero();
            for(int i=0;i<N;++i) {
                Swap(L.row(i,0,i),L.row(p[i],0,i));
                int end = TMV_MIN(i+nlo+1,N);
                L.col(i,i+1,end) = LUx.col(i,i+1,end);
            }
        }
    }

    //
    // LDivEq
    //

    template <class T, class T1> 
    void LU_PackedPL_LDivEq(
        const GenBandMatrix<T1>& LUx,
        const int* p, const MatrixView<T>& m) 
    {
#ifdef XDEBUG
        Matrix<T> m0 = m;
        Matrix<T> PL0(LUx.colsize(),LUx.colsize());
        PL0 = T(1);
        LU_PackedPL_Unpack(LUx,p,PL0.lowerTri(UnitDiag));
        if (LUx.nlo() > 0) PL0.reversePermuteRows(p);
        Matrix<T> m2 = m;
        m2 /= PL0;
#endif

        // Solve L y = m by forward substitution
        // Remember L is really:
        //
        // L = (  1   0 ) ( 1  0   0 ) ( I  0   0 )
        //     ( P0L0 I ) ( 0  1   0 ) ( 0  1   0 ) ...
        //                ( 0 P1L1 I ) ( 0 P2L2 I )
        //
        // where the Li are columns of nlo length which are
        // stored in the lower band of LUx,
        // and each Pi is a row swap of i with p[i]
        //
        const int N = LUx.colsize();
        const int nlo = LUx.nlo();
        if (nlo > 0) {
            int jn=nlo+1;  // jn = j+nlo+1
            const int* pj = p;
            for(int j=0; j+1<N; ++j,++pj) {
                TMVAssert(*pj<int(m.colsize()));
                m.swapRows(j,*pj);
                m.rowRange(j+1,jn) -= LUx.col(j,j+1,jn) ^ m.row(j);
                if (jn<N) ++jn;
            }
        }
#ifdef XDEBUG
        TMV_RealType(T) kappa = Norm(PL0)*Norm(PL0.inverse());
        TMV_RealType(T) normdiff = Norm(m2-m);
        if (normdiff > 0.0001*kappa*Norm(m2)) {
            cerr<<"LU_PackedPL_LDivEq\n";
            cerr<<"LUx = "<<LUx<<endl;
            cerr<<"PL = "<<PL0<<endl;
            cerr<<"p = ";
            for(size_t i=0;i<LUx.colsize();i++) cerr<<p[i]<<" ";
            cerr<<endl;
            cerr<<"m = "<<m0<<endl;
            cerr<<"m => "<<m<<endl;
            cerr<<"m2 = "<<m2<<endl;
            cerr<<"Norm(m-m2) = "<<normdiff<<endl;
            cerr<<"kappa = "<<kappa<<endl;
            abort();
        }
#endif
    } 

    template <class T, class T1> 
    static void NonLapLU_LDivEq(
        const GenBandMatrix<T1>& LUx,
        const int* p, const MatrixView<T>& m) 
    { 
#ifdef XDEBUG
        Matrix<T> m0 = m;
        LowerTriMatrix<T,UnitDiag> L0(LUx.colsize());
        LU_PackedPL_Unpack(LUx,p,L0.view());
        UpperTriMatrix<T> U0 = BandMatrixViewOf(LUx,0,LUx.nhi());
        Matrix<T> PLU = L0*U0;
        if (LUx.nlo() > 0) PLU.reversePermuteRows(p);
        Matrix<T> m2 = m;
        m2 /= PLU;
#endif
        // Solve A x = m given that A = L U
        // L U x = m
        TMVAssert(m.colsize() == LUx.colsize());
        TMVAssert(LUx.isSquare());

        // Solve L y = m by forward substitution
        LU_PackedPL_LDivEq(LUx,p,m);

        // Next solve U x = y by back substitution
        TriLDivEq(BandMatrixViewOf(LUx,0,LUx.nhi()),m,NonUnitDiag);

#ifdef XDEBUG
        TMV_RealType(T) kappa = Norm(PLU)*Norm(PLU.inverse());
        TMV_RealType(T) normdiff = Norm(m2-m);
        if (normdiff > 0.0001*kappa*Norm(m2)) {
            cerr<<"LU_PackedPL_LDivEq\n";
            cerr<<"LUx = "<<LUx<<endl;
            cerr<<"L = "<<L0<<endl;
            cerr<<"U = "<<U0<<endl;
            cerr<<"p = ";
            for(size_t i=0;i<LUx.colsize();i++) cerr<<p[i]<<" ";
            cerr<<endl;
            cerr<<"PLU = "<<PLU<<endl;
            cerr<<"m = "<<m0<<endl;
            cerr<<"m => "<<m<<endl;
            cerr<<"m2 = "<<m2<<endl;
            cerr<<"Norm(m-m2) = "<<normdiff<<endl;
            cerr<<"kappa = "<<kappa<<endl;
            abort();
        }
#endif
    }

    //
    // LDivEq
    //

    template <class T, class T1> 
    void LU_PackedPL_RDivEq(
        const GenBandMatrix<T1>& LUx,
        const int* p, const MatrixView<T>& m) 
    {
#ifdef XDEBUG
        Matrix<T> m0 = m;
        Matrix<T> PL0(LUx.colsize(),LUx.colsize());
        PL0 = T(1);
        LU_PackedPL_Unpack(LUx,p,PL0.lowerTri(UnitDiag));
        if (LUx.nlo() > 0) PL0.reversePermuteRows(p);
        Matrix<T> m2 = m;
        m2 %= PL0;
#endif
        // Solve z L = y by back substitution with L = :
        //
        // L = (  1   0 ) ( 1  0   0 ) ( I  0   0 )     ( I    0     0 ) ( I 0 )
        //     ( P0L0 I ) ( 0  1   0 ) ( 0  1   0 ) ... ( 0    1     0 ) ( 0 1 )
        //                ( 0 P1L1 I ) ( 0 P2L2 I )     ( 0 Pn-1Ln-1 1 )
        //
        const int N = LUx.colsize();
        const int nlo = LUx.nlo();
        if (nlo > 0) {
            int jn=N;
            int k=nlo-1;
            const int* pj = p+N-1;
            for(int j=N-1;j>0;) {
                --j; --pj;
                m.col(j) -= m.colRange(j+1,jn) * LUx.col(j,j+1,jn);
                TMVAssert(*pj<int(m.rowsize()));
                m.swapCols(j,*pj);
                if (k>0) --k; else --jn;
            }
        }
#ifdef XDEBUG
        TMV_RealType(T) kappa = Norm(PL0)*Norm(PL0.inverse());
        TMV_RealType(T) normdiff = Norm(m2-m);
        if (normdiff > 0.0001*kappa*Norm(m2)) {
            cerr<<"LU_PackedPL_RDivEq\n";
            cerr<<"LUx = "<<LUx<<endl;
            cerr<<"PL = "<<PL0<<endl;
            cerr<<"p = ";
            for(size_t i=0;i<LUx.colsize();i++) cerr<<p[i]<<" ";
            cerr<<endl;
            cerr<<"m = "<<m0<<endl;
            cerr<<"m => "<<m<<endl;
            cerr<<"m2 = "<<m2<<endl;
            cerr<<"Norm(m-m2) = "<<normdiff<<endl;
            cerr<<"kappa = "<<kappa<<endl;
            abort();
        }
#endif
    }

    template <class T, class T1> 
    static void NonLapLU_RDivEq(
        const GenBandMatrix<T1>& LUx,
        const int* p, const MatrixView<T>& m) 
    { 
#ifdef XDEBUG
        Matrix<T> m0 = m;
        LowerTriMatrix<T,UnitDiag> L0(LUx.colsize());
        LU_PackedPL_Unpack(LUx,p,L0.view());
        UpperTriMatrix<T> U0 = BandMatrixViewOf(LUx,0,LUx.nhi());
        Matrix<T> PLU = L0*U0;
        if (LUx.nlo() > 0) PLU.reversePermuteRows(p);
        Matrix<T> m2 = m;
        m2 %= PLU;
#endif
        // Solve x A = m given that A = L U
        // x L U = m
        TMVAssert(m.rowsize() == LUx.colsize());
        TMVAssert(LUx.isSquare());

        // First solve y U = m by forward substitution
        // Or: UT yT = mT
        ConstBandMatrixView<T1> U = BandMatrixViewOf(LUx,0,LUx.nhi());
        TriLDivEq(Transpose(U),Transpose(m),NonUnitDiag);

        // Next solve z L = y by back substitution with L 
        LU_PackedPL_RDivEq(LUx,p,m);
#ifdef XDEBUG
        TMV_RealType(T) kappa = Norm(PLU)*Norm(PLU.inverse());
        TMV_RealType(T) normdiff = Norm(m2-m);
        if (normdiff > 0.0001*kappa*Norm(m2)) {
            cerr<<"LU_PackedPL_RDivEq\n";
            cerr<<"LUx = "<<LUx<<endl;
            cerr<<"L = "<<L0<<endl;
            cerr<<"U = "<<U0<<endl;
            cerr<<"p = ";
            for(size_t i=0;i<LUx.colsize();i++) cerr<<p[i]<<" ";
            cerr<<endl;
            cerr<<"PLU = "<<PLU<<endl;
            cerr<<"m = "<<m0<<endl;
            cerr<<"m => "<<m<<endl;
            cerr<<"m2 = "<<m2<<endl;
            cerr<<"Norm(m-m2) = "<<normdiff<<endl;
            cerr<<"kappa = "<<kappa<<endl;
            abort();
        }
#endif
    }

#ifdef LAP
    template <class T, class T1> 
    static inline void LapLU_LDivEq(
        const GenBandMatrix<T1>& LUx, const int* P, const MatrixView<T>& m)
    { NonLapLU_LDivEq(LUx,P,m); }
    template <class T, class T1> 
    static inline void LapTriDiagLU_LDivEq(
        const GenBandMatrix<T1>& LUx, const int* P, const MatrixView<T>& m)
    { NonLapLU_LDivEq(LUx,P,m); }
    template <class T, class T1> 
    static inline void LapLU_RDivEq(
        const GenBandMatrix<T1>& LUx, const int* P, const MatrixView<T>& m)
    { NonLapLU_RDivEq(LUx,P,m); }
    template <class T, class T1> 
    static inline void LapTriDiagLU_RDivEq(
        const GenBandMatrix<T1>& LUx, const int* P, const MatrixView<T>& m)
    { NonLapLU_RDivEq(LUx,P,m); }
#ifdef INST_DOUBLE
    template <> 
    void LapLU_LDivEq(
        const GenBandMatrix<double>& LUx,
        const int* P, const MatrixView<double>& m) 
    {
        TMVAssert(m.colsize() == LUx.colsize());
        TMVAssert(LUx.isSquare());
        TMVAssert(m.rowsize() > 0);
        TMVAssert(m.colsize() > 0);
        TMVAssert(LUx.nhi() > LUx.nlo());
        TMVAssert(LUx.nlo() > 0);
        TMVAssert(LUx.nhi() < int(LUx.colsize()));
        TMVAssert(LUx.nlo() < int(LUx.colsize()));
        TMVAssert(LUx.iscm());
        TMVAssert(m.iscm());

        int n = LUx.colsize();
        int kl = LUx.nlo();
        int ku = LUx.nhi()-LUx.nlo();
        int nrhs = m.rowsize();
        int lda = LUx.stepj()+1;
        int ldm = m.stepj();
        AlignedArray<int> lap_p(n);
        for(int i=0;i<n;i++) (lap_p.get())[i] = P[i] LAPPLUS1;
        LAPNAME(dgbtrs) (
            LAPCM LAPCH_NT,LAPV(n),LAPV(kl),LAPV(ku),LAPV(nrhs),
            LAPP(LUx.cptr()-LUx.nhi()),LAPV(lda),
            LAPP(lap_p.get()),LAPP(m.ptr()),LAPV(ldm) LAPINFO LAP1);
        LAP_Results("dgbtrs");
    }
    template <> 
    void LapLU_LDivEq(
        const GenBandMatrix<std::complex<double> >& LUx,
        const int* P, const MatrixView<std::complex<double> >& m)
    {
        TMVAssert(m.colsize() == LUx.colsize());
        TMVAssert(LUx.isSquare());
        TMVAssert(m.rowsize() > 0);
        TMVAssert(m.colsize() > 0);
        TMVAssert(LUx.nhi() > LUx.nlo());
        TMVAssert(LUx.nlo() > 0);
        TMVAssert(LUx.nhi() < int(LUx.colsize()));
        TMVAssert(LUx.nlo() < int(LUx.colsize()));
        TMVAssert(LUx.ct() == NonConj);
        TMVAssert(m.ct() == NonConj);
        TMVAssert(LUx.iscm());
        TMVAssert(m.iscm());

        int n = LUx.colsize();
        int kl = LUx.nlo();
        int ku = LUx.nhi()-LUx.nlo();
        int nrhs = m.rowsize();
        int lda = LUx.diagstep();
        int ldm = m.stepj();
        AlignedArray<int> lap_p(n);
        for(int i=0;i<n;i++) (lap_p.get())[i] = P[i] LAPPLUS1;
        LAPNAME(zgbtrs) (
            LAPCM LAPCH_NT,LAPV(n),LAPV(kl),LAPV(ku),LAPV(nrhs),
            LAPP(LUx.cptr()-LUx.nhi()),LAPV(lda),
            LAPP(lap_p.get()),LAPP(m.ptr()),LAPV(ldm) LAPINFO LAP1);
        LAP_Results("zgbtrs");
    }
    template <> 
    void LapTriDiagLU_LDivEq(
        const GenBandMatrix<double>& LUx,
        const int* P, const MatrixView<double>& m) 
    {
        TMVAssert(m.colsize() == LUx.colsize());
        TMVAssert(LUx.isSquare());
        TMVAssert(m.rowsize() > 0);
        TMVAssert(m.colsize() > 0);
        TMVAssert(LUx.nhi() == 2);
        TMVAssert(LUx.nlo() == 1);
        TMVAssert(LUx.nhi() < int(LUx.colsize()));
        TMVAssert(LUx.nlo() < int(LUx.colsize()));
        TMVAssert(LUx.isdm());
        TMVAssert(m.iscm());

        int n = LUx.colsize();
        int nrhs = m.rowsize();
        int ldm = m.stepj();
        AlignedArray<int> lap_p(n);
        for(int i=0;i<n;i++) (lap_p.get())[i] = P[i] LAPPLUS1;
        LAPNAME(dgttrs) (
            LAPCM LAPCH_NT,LAPV(n),LAPV(nrhs),
            LAPP(LUx.diag(-1).cptr()),LAPP(LUx.diag().cptr()),
            LAPP(LUx.diag(1).cptr()),LAPP(LUx.diag(2).cptr()),
            LAPP(lap_p.get()),LAPP(m.ptr()),LAPV(ldm) LAPINFO LAP1);
        LAP_Results("dgttrs");
    }
    template <> 
    void LapTriDiagLU_LDivEq(
        const GenBandMatrix<std::complex<double> >& LUx,
        const int* P, const MatrixView<std::complex<double> >& m)
    {
        TMVAssert(m.colsize() == LUx.colsize());
        TMVAssert(LUx.isSquare());
        TMVAssert(m.rowsize() > 0);
        TMVAssert(m.colsize() > 0);
        TMVAssert(LUx.nhi() == 2);
        TMVAssert(LUx.nlo() == 1);
        TMVAssert(LUx.nhi() < int(LUx.colsize()));
        TMVAssert(LUx.nlo() < int(LUx.colsize()));
        TMVAssert(LUx.ct() == NonConj);
        TMVAssert(m.ct() == NonConj);
        TMVAssert(LUx.isdm());
        TMVAssert(m.iscm());

        int n = LUx.colsize();
        int nrhs = m.rowsize();
        int ldm = m.stepj();
        AlignedArray<int> lap_p(n);
        for(int i=0;i<n;i++) (lap_p.get())[i] = P[i] LAPPLUS1;
        LAPNAME(zgttrs) (
            LAPCM LAPCH_NT,LAPV(n),LAPV(nrhs),
            LAPP(LUx.diag(-1).cptr()), LAPP(LUx.diag().cptr()),
            LAPP(LUx.diag(1).cptr()), LAPP(LUx.diag(2).cptr()), 
            LAPP(lap_p.get()), LAPP(m.ptr()),LAPV(ldm) LAPINFO LAP1);
        LAP_Results("zgttrs");
    }
    template <> 
    void LapLU_RDivEq(
        const GenBandMatrix<double>& LUx,
        const int* P, const MatrixView<double>& m) 
    {
        TMVAssert(m.rowsize() == LUx.colsize());
        TMVAssert(LUx.isSquare());
        TMVAssert(m.rowsize() > 0);
        TMVAssert(m.colsize() > 0);
        TMVAssert(LUx.nhi() > LUx.nlo());
        TMVAssert(LUx.nlo() > 0);
        TMVAssert(LUx.nhi() < int(LUx.colsize()));
        TMVAssert(LUx.nlo() < int(LUx.colsize()));
        TMVAssert(LUx.iscm());
        TMVAssert(m.isrm());

        int n = LUx.colsize();
        int kl = LUx.nlo();
        int ku = LUx.nhi()-LUx.nlo();
        int nrhs = m.colsize();
        int lda = LUx.stepj()+1;
        int ldm = m.stepi();
        AlignedArray<int> lap_p(n);
        for(int i=0;i<n;i++) (lap_p.get())[i] = P[i]+1;
        LAPNAME(dgbtrs) (
            LAPCM LAPCH_T,LAPV(n),LAPV(kl),LAPV(ku),LAPV(nrhs),
            LAPP(LUx.cptr()-LUx.nhi()),LAPV(lda),
            LAPP(lap_p.get()),LAPP(m.ptr()),LAPV(ldm) LAPINFO LAP1);
        LAP_Results("dgbtrs");
    }
    template <> 
    void LapLU_RDivEq(
        const GenBandMatrix<std::complex<double> >& LUx,
        const int* P, const MatrixView<std::complex<double> >& m)
    {
        TMVAssert(m.rowsize() == LUx.colsize());
        TMVAssert(LUx.isSquare());
        TMVAssert(m.rowsize() > 0);
        TMVAssert(m.colsize() > 0);
        TMVAssert(LUx.nhi() > LUx.nlo());
        TMVAssert(LUx.nlo() > 0);
        TMVAssert(LUx.nhi() < int(LUx.colsize()));
        TMVAssert(LUx.nlo() < int(LUx.colsize()));
        TMVAssert(m.ct() == NonConj);
        TMVAssert(LUx.iscm());
        TMVAssert(m.isrm());

        int n = LUx.colsize();
        int kl = LUx.nlo();
        int ku = LUx.nhi()-LUx.nlo();
        int nrhs = m.colsize();
        int lda = LUx.stepj()+1;
        int ldm = m.stepi();
        AlignedArray<int> lap_p(n);
        for(int i=0;i<n;i++) (lap_p.get())[i] = P[i]+1;
        LAPNAME(zgbtrs) (
            LAPCM LUx.isconj()?LAPCH_CT:LAPCH_T,
            LAPV(n),LAPV(kl),LAPV(ku),LAPV(nrhs),
            LAPP(LUx.cptr()-LUx.nhi()),LAPV(lda),
            LAPP(lap_p.get()),LAPP(m.ptr()),LAPV(ldm) LAPINFO LAP1);
        LAP_Results("zgbtrs");
    }
    template <> 
    void LapTriDiagLU_RDivEq(
        const GenBandMatrix<double>& LUx,
        const int* P, const MatrixView<double>& m) 
    {
        TMVAssert(m.rowsize() == LUx.colsize());
        TMVAssert(LUx.isSquare());
        TMVAssert(m.rowsize() > 0);
        TMVAssert(m.colsize() > 0);
        TMVAssert(LUx.nhi() == 2);
        TMVAssert(LUx.nlo() == 1);
        TMVAssert(LUx.nhi() < int(LUx.colsize()));
        TMVAssert(LUx.nlo() < int(LUx.colsize()));
        TMVAssert(LUx.isdm());
        TMVAssert(m.isrm());

        int n = LUx.colsize();
        int nrhs = m.colsize();
        int ldm = m.stepi();
        AlignedArray<int> lap_p(n);
        for(int i=0;i<n;i++) (lap_p.get())[i] = P[i]+1;
        LAPNAME(dgttrs) (
            LAPCM LAPCH_T,LAPV(n),LAPV(nrhs),
            LAPP(LUx.diag(-1).cptr()), LAPP(LUx.diag().cptr()),
            LAPP(LUx.diag(1).cptr()), LAPP(LUx.diag(2).cptr()),
            LAPP(lap_p.get()),LAPP(m.ptr()),LAPV(ldm) LAPINFO LAP1);
        LAP_Results("dgttrs");
    }
    template <> 
    void LapTriDiagLU_RDivEq(
        const GenBandMatrix<std::complex<double> >& LUx,
        const int* P, const MatrixView<std::complex<double> >& m)
    {
        TMVAssert(m.rowsize() == LUx.colsize());
        TMVAssert(LUx.isSquare());
        TMVAssert(m.rowsize() > 0);
        TMVAssert(m.colsize() > 0);
        TMVAssert(LUx.nhi() == 2);
        TMVAssert(LUx.nlo() == 1);
        TMVAssert(LUx.nhi() < int(LUx.colsize()));
        TMVAssert(LUx.nlo() < int(LUx.colsize()));
        TMVAssert(m.ct() == NonConj);
        TMVAssert(LUx.isdm());
        TMVAssert(m.isrm());

        int n = LUx.colsize();
        int nrhs = m.colsize();
        int ldm = m.stepi();
        AlignedArray<int> lap_p(n);
        for(int i=0;i<n;i++) (lap_p.get())[i] = P[i]+1;
        LAPNAME(zgttrs) (
            LAPCM LUx.isconj()?LAPCH_CT:LAPCH_T, LAPV(n),LAPV(nrhs),
            LAPP(LUx.diag(-1).cptr()), LAPP(LUx.diag().cptr()),
            LAPP(LUx.diag(1).cptr()), LAPP(LUx.diag(2).cptr()),
            LAPP(lap_p.get()), LAPP(m.ptr()),LAPV(ldm) LAPINFO LAP1);
        LAP_Results("zgttrs");
    }
#endif // DOUBLE
#ifdef INST_FLOAT
    template <> 
    void LapLU_LDivEq(
        const GenBandMatrix<float>& LUx,
        const int* P, const MatrixView<float>& m) 
    {
        TMVAssert(m.colsize() == LUx.colsize());
        TMVAssert(LUx.isSquare());
        TMVAssert(m.rowsize() > 0);
        TMVAssert(m.colsize() > 0);
        TMVAssert(LUx.nhi() > LUx.nlo());
        TMVAssert(LUx.nlo() > 0);
        TMVAssert(LUx.nhi() < int(LUx.colsize()));
        TMVAssert(LUx.nlo() < int(LUx.colsize()));
        TMVAssert(LUx.iscm());
        TMVAssert(m.iscm());

        int n = LUx.colsize();
        int kl = LUx.nlo();
        int ku = LUx.nhi()-LUx.nlo();
        int nrhs = m.rowsize();
        int lda = LUx.stepj()+1;
        int ldm = m.stepj();
        AlignedArray<int> lap_p(n);
        for(int i=0;i<n;i++) (lap_p.get())[i] = P[i] LAPPLUS1;
        LAPNAME(sgbtrs) (
            LAPCM LAPCH_NT,LAPV(n),LAPV(kl),LAPV(ku),LAPV(nrhs),
            LAPP(LUx.cptr()-LUx.nhi()),LAPV(lda),
            LAPP(lap_p.get()),LAPP(m.ptr()),LAPV(ldm) LAPINFO LAP1);
        LAP_Results("sgbtrs");
    }
    template <> 
    void LapLU_LDivEq(
        const GenBandMatrix<std::complex<float> >& LUx,
        const int* P, const MatrixView<std::complex<float> >& m)
    {
        TMVAssert(m.colsize() == LUx.colsize());
        TMVAssert(LUx.isSquare());
        TMVAssert(m.rowsize() > 0);
        TMVAssert(m.colsize() > 0);
        TMVAssert(LUx.nhi() > LUx.nlo());
        TMVAssert(LUx.nlo() > 0);
        TMVAssert(LUx.nhi() < int(LUx.colsize()));
        TMVAssert(LUx.nlo() < int(LUx.colsize()));
        TMVAssert(LUx.ct() == NonConj);
        TMVAssert(m.ct() == NonConj);
        TMVAssert(LUx.iscm());
        TMVAssert(m.iscm());

        int n = LUx.colsize();
        int kl = LUx.nlo();
        int ku = LUx.nhi()-LUx.nlo();
        int nrhs = m.rowsize();
        int lda = LUx.stepj()+1;
        int ldm = m.stepj();
        AlignedArray<int> lap_p(n);
        for(int i=0;i<n;i++) (lap_p.get())[i] = P[i] LAPPLUS1;
        LAPNAME(cgbtrs) (
            LAPCM LAPCH_NT,LAPV(n),LAPV(kl),LAPV(ku),LAPV(nrhs),
            LAPP(LUx.cptr()-LUx.nhi()),LAPV(lda),
            LAPP(lap_p.get()),LAPP(m.ptr()),LAPV(ldm) LAPINFO LAP1);
        LAP_Results("cgbtrs");
    }
    template <> 
    void LapTriDiagLU_LDivEq(
        const GenBandMatrix<float>& LUx,
        const int* P, const MatrixView<float>& m) 
    {
        TMVAssert(m.colsize() == LUx.colsize());
        TMVAssert(LUx.isSquare());
        TMVAssert(m.rowsize() > 0);
        TMVAssert(m.colsize() > 0);
        TMVAssert(LUx.nhi() == 2);
        TMVAssert(LUx.nlo() == 1);
        TMVAssert(LUx.nhi() < int(LUx.colsize()));
        TMVAssert(LUx.nlo() < int(LUx.colsize()));
        TMVAssert(LUx.isdm());
        TMVAssert(m.iscm());

        int n = LUx.colsize();
        int nrhs = m.rowsize();
        int ldm = m.stepj();
        AlignedArray<int> lap_p(n);
        for(int i=0;i<n;i++) (lap_p.get())[i] = P[i] LAPPLUS1;
        LAPNAME(sgttrs) (
            LAPCM LAPCH_NT,LAPV(n),LAPV(nrhs),
            LAPP(LUx.diag(-1).cptr()), LAPP(LUx.diag().cptr()),
            LAPP(LUx.diag(1).cptr()), LAPP(LUx.diag(2).cptr()),
            LAPP(lap_p.get()),LAPP(m.ptr()),LAPV(ldm) LAPINFO LAP1);
        LAP_Results("sgttrs");
    }
    template <> 
    void LapTriDiagLU_LDivEq(
        const GenBandMatrix<std::complex<float> >& LUx,
        const int* P, const MatrixView<std::complex<float> >& m)
    {
        TMVAssert(m.colsize() == LUx.colsize());
        TMVAssert(LUx.isSquare());
        TMVAssert(m.rowsize() > 0);
        TMVAssert(m.colsize() > 0);
        TMVAssert(LUx.nhi() == 2);
        TMVAssert(LUx.nlo() == 1);
        TMVAssert(LUx.nhi() < int(LUx.colsize()));
        TMVAssert(LUx.nlo() < int(LUx.colsize()));
        TMVAssert(LUx.ct() == NonConj);
        TMVAssert(m.ct() == NonConj);
        TMVAssert(LUx.isdm());
        TMVAssert(m.iscm());

        int n = LUx.colsize();
        int nrhs = m.rowsize();
        int ldm = m.stepj();
        AlignedArray<int> lap_p(n);
        for(int i=0;i<n;i++) (lap_p.get())[i] = P[i] LAPPLUS1;
        LAPNAME(cgttrs) (
            LAPCM LAPCH_NT,LAPV(n),LAPV(nrhs),
            LAPP(LUx.diag(-1).cptr()), LAPP(LUx.diag().cptr()),
            LAPP(LUx.diag(1).cptr()), LAPP(LUx.diag(2).cptr()),
            LAPP(lap_p.get()),LAPP(m.ptr()),LAPV(ldm) LAPINFO LAP1);
        LAP_Results("cgttrs");
    }
    template <> 
    void LapLU_RDivEq(
        const GenBandMatrix<float>& LUx,
        const int* P, const MatrixView<float>& m) 
    {
        TMVAssert(m.rowsize() == LUx.colsize());
        TMVAssert(LUx.isSquare());
        TMVAssert(m.rowsize() > 0);
        TMVAssert(m.colsize() > 0);
        TMVAssert(LUx.nhi() > LUx.nlo());
        TMVAssert(LUx.nlo() > 0);
        TMVAssert(LUx.nhi() < int(LUx.colsize()));
        TMVAssert(LUx.nlo() < int(LUx.colsize()));
        TMVAssert(LUx.iscm());
        TMVAssert(m.isrm());

        int n = LUx.colsize();
        int kl = LUx.nlo();
        int ku = LUx.nhi()-LUx.nlo();
        int nrhs = m.colsize();
        int lda = LUx.stepj()+1;
        int ldm = m.stepi();
        AlignedArray<int> lap_p(n);
        for(int i=0;i<n;i++) (lap_p.get())[i] = P[i]+1;
        LAPNAME(sgbtrs) (
            LAPCM LAPCH_T,LAPV(n),LAPV(kl),LAPV(ku),LAPV(nrhs),
            LAPP(LUx.cptr()-LUx.nhi()),LAPV(lda),
            LAPP(lap_p.get()),LAPP(m.ptr()),LAPV(ldm) LAPINFO LAP1);
        LAP_Results("sgbtrs");
    }
    template <> 
    void LapLU_RDivEq(
        const GenBandMatrix<std::complex<float> >& LUx,
        const int* P, const MatrixView<std::complex<float> >& m)
    {
        TMVAssert(m.rowsize() == LUx.colsize());
        TMVAssert(LUx.isSquare());
        TMVAssert(m.rowsize() > 0);
        TMVAssert(m.colsize() > 0);
        TMVAssert(LUx.nhi() > LUx.nlo());
        TMVAssert(LUx.nlo() > 0);
        TMVAssert(LUx.nhi() < int(LUx.colsize()));
        TMVAssert(LUx.nlo() < int(LUx.colsize()));
        TMVAssert(m.ct() == NonConj);
        TMVAssert(LUx.iscm());
        TMVAssert(m.isrm());

        int n = LUx.colsize();
        int kl = LUx.nlo();
        int ku = LUx.nhi()-LUx.nlo();
        int nrhs = m.colsize();
        int lda = LUx.stepj()+1;
        int ldm = m.stepi();
        AlignedArray<int> lap_p(n);
        for(int i=0;i<n;i++) (lap_p.get())[i] = P[i]+1;
        LAPNAME(cgbtrs) (
            LAPCM LUx.isconj()?LAPCH_CT:LAPCH_T,
            LAPV(n),LAPV(kl),LAPV(ku),LAPV(nrhs),
            LAPP(LUx.cptr()-LUx.nhi()),LAPV(lda),
            LAPP(lap_p.get()),LAPP(m.ptr()),LAPV(ldm) LAPINFO LAP1);
        LAP_Results("cgbtrs");
    }
    template <> 
    void LapTriDiagLU_RDivEq(
        const GenBandMatrix<float>& LUx,
        const int* P, const MatrixView<float>& m) 
    {
        TMVAssert(m.rowsize() == LUx.colsize());
        TMVAssert(LUx.isSquare());
        TMVAssert(m.rowsize() > 0);
        TMVAssert(m.colsize() > 0);
        TMVAssert(LUx.nhi() == 2);
        TMVAssert(LUx.nlo() == 1);
        TMVAssert(LUx.nhi() < int(LUx.colsize()));
        TMVAssert(LUx.nlo() < int(LUx.colsize()));
        TMVAssert(LUx.isdm());
        TMVAssert(m.isrm());

        int n = LUx.colsize();
        int nrhs = m.colsize();
        int ldm = m.stepi();
        AlignedArray<int> lap_p(n);
        for(int i=0;i<n;i++) (lap_p.get())[i] = P[i]+1;
        LAPNAME(sgttrs) (
            LAPCM LAPCH_T,LAPV(n),LAPV(nrhs),
            LAPP(LUx.diag(-1).cptr()), LAPP(LUx.diag().cptr()),
            LAPP(LUx.diag(1).cptr()), LAPP(LUx.diag(2).cptr()),
            LAPP(lap_p.get()),LAPP(m.ptr()),LAPV(ldm) LAPINFO LAP1);
        LAP_Results("sgttrs");
    }
    template <> 
    void LapTriDiagLU_RDivEq(
        const GenBandMatrix<std::complex<float> >& LUx,
        const int* P, const MatrixView<std::complex<float> >& m)
    {
        TMVAssert(m.rowsize() == LUx.colsize());
        TMVAssert(LUx.isSquare());
        TMVAssert(m.rowsize() > 0);
        TMVAssert(m.colsize() > 0);
        TMVAssert(LUx.nhi() == 2);
        TMVAssert(LUx.nlo() == 1);
        TMVAssert(LUx.nhi() < int(LUx.colsize()));
        TMVAssert(LUx.nlo() < int(LUx.colsize()));
        TMVAssert(m.ct() == NonConj);
        TMVAssert(LUx.isdm());
        TMVAssert(m.isrm());

        int n = LUx.colsize();
        int nrhs = m.colsize();
        int ldm = m.stepi();
        AlignedArray<int> lap_p(n);
        for(int i=0;i<n;i++) (lap_p.get())[i] = P[i]+1;
        LAPNAME(cgttrs) (
            LAPCM LUx.isconj()?LAPCH_CT:LAPCH_T,LAPV(n),LAPV(nrhs),
            LAPP(LUx.diag(-1).cptr()), LAPP(LUx.diag().cptr()),
            LAPP(LUx.diag(1).cptr()), LAPP(LUx.diag(2).cptr()),
            LAPP(lap_p.get()), LAPP(m.ptr()),LAPV(ldm) LAPINFO LAP1);
        LAP_Results("cgttrs");
    }
#endif // FLOAT
#endif // LAP

    template <class T, class T1> 
    void LU_LDivEq(
        const GenBandMatrix<T1>& LUx,
        const int* P, const MatrixView<T>& m) 
    {
        TMVAssert(m.colsize() == LUx.colsize());
        TMVAssert(LUx.isSquare());

        if (m.colsize() > 0 && m.rowsize() > 0) {
#ifdef LAP
            if (m.isconj())
                LU_LDivEq(LUx.conjugate(),P,m.conjugate());
            else if (m.iscm() && LUx.iscm() && !LUx.isconj() && LUx.nlo() > 0)
                LapLU_LDivEq(LUx,P,m);
            else if (m.iscm() && LUx.isdm() && 
                     LUx.nlo() == 1 && LUx.nhi() == 2 && !LUx.isconj())
                LapTriDiagLU_LDivEq(LUx,P,m);
            else
#endif
                NonLapLU_LDivEq(LUx,P,m);
        }
    }

    template <class T, class T1> 
    void LU_RDivEq(
        const GenBandMatrix<T1>& LUx,
        const int* P, const MatrixView<T>& m) 
    {
        TMVAssert(m.rowsize() == LUx.colsize());
        TMVAssert(LUx.isSquare());

        if (m.colsize() > 0 && m.rowsize() > 0) {
#ifdef LAP
            if (m.isconj())
                LU_RDivEq(LUx.conjugate(),P,m.conjugate());
            else if (m.isrm() && LUx.iscm() && LUx.nlo()>0)
                LapLU_RDivEq(LUx,P,m);
            else if (m.isrm() && LUx.isdm() && 
                     LUx.nlo() == 1 && LUx.nhi() == 2)
                LapTriDiagLU_RDivEq(LUx,P,m);
            else
#endif
                NonLapLU_RDivEq(LUx,P,m);
        }
    }

#ifdef INST_INT
#undef INST_INT
#endif

#define InstFile "TMV_BandLUDiv.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


