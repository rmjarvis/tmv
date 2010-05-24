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
#include "tmv/TMV_BandLUD.h"
#include "tmv/TMV_BandMatrix.h"
#include "tmv/TMV_TriMatrix.h"
#include "tmv/TMV_VectorArith.h"
#include "tmv/TMV_MatrixArith.h"

#ifdef XDEBUG
#include <iostream>
#include "tmv/TMV_BandMatrixArith.h"
using std::cout;
using std::cerr;
using std::endl;
#endif

namespace tmv {

    //
    // Decompose
    //

    template <class T> 
    static void NonLapBandLU_Decompose(
        const BandMatrixView<T>& A, int* P, T& det)
    {
        // LU Decompostion with partial pivoting.
        //
        // For band matrices, we use a somewhat different algorithm than for
        // regular matrices.  With regular matrices, the main operations were 
        // Matrix * Vector.  With the band matrix, these become difficult
        // to implement, since the submatrix that is doing the multiplying
        // goes off the edge of the bands.  So we choose to implement an
        // algorithm based on outerproduct updates which works better, since
        // the matrix being updated is completely within the band.
        //
        // On input, A contains the original band matrix with nlo subdiagonals
        // and nhi superdiagonals.  A must be created with nlo subdiagonals
        // and nlo+nhi superdiagonals.
        //
        // On each step, we calculate the j column of L, the j row of U
        // and the diagonal Ujj.  here is the first step:
        // 
        // A = L0 U0
        // ( A00 A0x ) = (  1  0 ) ( U00 U0x )
        // ( Ax0 Axx )   ( Lx0 I ) (  0  A'  )
        //
        // In Ax0, only A_10..A_nlo,0 are nonzero.
        // In A0x, only A_01..A_0,nhi are nonzero.
        // Axx is also band diagonal with the same nlo,nhi
        //
        // The formulae for L,U components are:
        //
        // U00 = A00
        // Lx0 = Ax0/U00
        // U0x = A0x
        // Uxx = Axx - Lx0 U0x
        //
        // It is apparent that Lx0 and U0x will have the same nonzero structure 
        // as Ax0 and A0x respectively.  This continues down the recursion,
        // so when we are done, L is lower banded with nlo subdiagonals, and
        // U is upper banded with nhi superdiagonals.
        //  
        // Unfortunately, this gets messed up a bit with pivoting.  
        // If the pivot element is as low as it can be, A_nlo,0, then 
        // swapping rows nlo and 0 will put A_nlo,nlo+nhi into A_0,nlo+nhi.
        // So we need to expand the upper band storage to nlo+nhi.
        //
        // The other problem is a bit more subtle.  Swapping rows j+nlo and j
        // also moves data from the j row down to the j+nlo in some of the 
        // previous columns which store data for L.  This would also screw up
        // the band structure, but we don't actually need to do this swap.
        // If we just keep track of what the swaps would be, we can just swap
        // rows in the remaining parts of A without swapping rows for L.
        // This makes the LDivEq, RDivEq functions a bit more complicated.
        //
        TMVAssert(A.nhi() > A.nlo());
        TMVAssert(A.iscm());
        TMVAssert(A.ct() == NonConj);
        TMVAssert(A.isSquare());
        const int N = A.rowsize();

        typedef TMV_RealType(T) RT;
        const RT halfeps = TMV_Epsilon<RT>()/RT(2);

        const int ds = A.diagstep();
        T* Ujj = A.ptr();
        int* Pj=P;

        int endcol = A.nlo()+1;
        int endrow = A.nhi()+1;
        for (int j=0; j<N-1; ++j,Ujj+=ds,++Pj) {
            // Find the pivot element
            int ip;
            T piv = A.col(j,j,endcol).maxAbsElement(&ip);

            if (piv * halfeps == RT(0)) {
                ip = 0;
                A.col(j,j,endcol).setZero();
            }

            if (ip != 0) {
                ip += j;
                Swap(A.row(ip,j,endrow),A.row(j,j,endrow));
                *Pj = ip;
                det = -det;
            } else *Pj = j;

            // If Ujj is 0, then all of the L's are 0.
            // ie. Ujj Lij = 0 for all i>j
            // Any value for Lij is valid, so leave them 0.
            if (*Ujj != T(0)) A.col(j,j+1,endcol) /= *Ujj;

            A.subMatrix(j+1,endcol,j+1,endrow) -= 
                (A.col(j,j+1,endcol) ^ A.row(j,j+1,endrow));
            if (endcol < N) ++endcol;
            if (endrow < N) ++endrow;
        }
        // j == N-1
        *Pj = N-1;
    }

    template <class T> 
    static void NonLapTriDiagLU_Decompose(
        const BandMatrixView<T>& A, int* P, T& det)
    {
        // LU Decompostion for TriDiagonal BandMatrix
        //
        // For TriDiagonal BandMatrices, there are only two choices for 
        // each pivot, so we can specialize some of the calculations.
        // Otherwise, this is the same algorithm as above.
        //
        TMVAssert(A.isdm());
        TMVAssert(A.nlo() == 1);
        TMVAssert(A.nhi() == 2);
        TMVAssert(A.colsize() >= A.rowsize());
        TMVAssert(A.ct() == NonConj);
        TMVAssert(A.isSquare());
        const int N = A.rowsize();

        typedef TMV_RealType(T) RT;
        const RT halfeps = TMV_Epsilon<RT>()/RT(2);

        T* Ujj = A.ptr(); // = U(j,j)
        T* Lj = Ujj+A.stepi();  // = L(j+1,j)
        T* Uj1 = Ujj+A.stepj(); // = U(j,j+1)
        T* Uj2 = Uj1+A.stepj(); // = U(j,j+2)
        int* Pj=P;

        for (int j=0; j<N-1; ++j,++Ujj,++Lj,++Uj1,++Uj2,++Pj) {
            bool pivot = TMV_ABS(*Lj) > TMV_ABS(*Ujj);

            if (pivot) {
                //swap(A.row(j+1,j,endrow),A.row(j,j,endrow));
                TMV_SWAP(*Lj,*Ujj);
                if (j+1<N) TMV_SWAP(*(Ujj+1),*Uj1);
                if (j+2<N) TMV_SWAP(*(Uj1+1),*Uj2);
                *Pj = j+1;
                det = -det;
            } else *Pj = j;

#ifdef TMVFLDEBUG
            TMVAssert(Lj >= A.first);
            TMVAssert(Lj < A.last);
#endif
            if (*Ujj * halfeps == T(0)) {
                *Ujj = *Lj = T(0);
            } else {
                *Lj /= *Ujj;
            }

            //A.row(j+1,j+1,endrow) -= *Lj * A.row(j,j+1,endrow);
#ifdef TMVFLDEBUG
            TMVAssert(Ujj+1 >= A.first);
            TMVAssert(Ujj+1 < A.last);
#endif
            *(Ujj+1) -= *Lj * (*Uj1);
            if (pivot && j+2<N) *(Uj1+1) -= *Lj * (*Uj2); // (if !pivot, *Uj2 == 0)
        }
        // j == N-1
        *Pj = N-1;
    }

#ifdef LAP
    template <class T> 
    static inline void LapBandLU_Decompose(
        const BandMatrixView<T>& A, int* P, T& det)
    { NonLapBandLU_Decompose(A,P,det); }
#ifdef INST_DOUBLE
    template <> 
    void LapBandLU_Decompose(
        const BandMatrixView<double>& A, int* P, double& det)
    {
        TMVAssert(A.isSquare());
        TMVAssert(A.iscm());
        TMVAssert(A.ct()==NonConj);
        int n = A.rowsize();
        int kl = A.nlo();
        int ku = A.nhi()-kl;
        int lda = A.stepj()+1;
        auto_array<int> lap_p(new int[n]);
        LAPNAME(dgbtrf) (
            LAPCM LAPV(n),LAPV(n),LAPV(kl),LAPV(ku),
            LAPP(A.ptr()-A.nhi()),LAPV(lda),LAPP(lap_p.get()) LAPINFO );
        LAP_Results("dgbtrf");
        const int M = A.colsize();
        for(int i=0;i<M;i++) {
            P[i] = (lap_p.get())[i] LAPMINUS1;
            if (P[i]!=i) det = -det;
        }
    }
    template <> 
    void LapBandLU_Decompose(
        const BandMatrixView<std::complex<double> >& A, int* P,
        std::complex<double>& det)
    {
        TMVAssert(A.isSquare());
        TMVAssert(A.iscm());
        TMVAssert(A.ct()==NonConj);
        int n = A.rowsize();
        int kl = A.nlo();
        int ku = A.nhi()-kl;
        int lda = A.stepj()+1;
        auto_array<int> lap_p(new int[n]);
        LAPNAME(zgbtrf) (
            LAPCM LAPV(n),LAPV(n),LAPV(kl),LAPV(ku),
            LAPP(A.ptr()-A.nhi()),LAPV(lda),LAPP(lap_p.get()) LAPINFO );
        LAP_Results("zgbtrf");
        const int M = A.colsize();
        for(int i=0;i<M;i++) {
            P[i] = (lap_p.get())[i] LAPMINUS1;
            if (P[i]!=i) det = -det;
        }
    }
#endif
#ifdef INST_FLOAT
    template <> 
    void LapBandLU_Decompose(
        const BandMatrixView<float>& A, int* P, float& det)
    {
        TMVAssert(A.isSquare());
        TMVAssert(A.iscm());
        TMVAssert(A.ct()==NonConj);
        int n = A.rowsize();
        int kl = A.nlo();
        int ku = A.nhi()-kl;
        int lda = A.stepj()+1;
        auto_array<int> lap_p(new int[n]);
        LAPNAME(sgbtrf) (
            LAPCM LAPV(n),LAPV(n),LAPV(kl),LAPV(ku),
            LAPP(A.ptr()-A.nhi()),LAPV(lda),LAPP(lap_p.get()) LAPINFO );
        LAP_Results("sgbtrf");
        const int M = A.colsize();
        for(int i=0;i<M;i++) {
            P[i] = (lap_p.get())[i] LAPMINUS1;
            if (P[i]!=i) det = -det;
        }
    }
    template <> 
    void LapBandLU_Decompose(
        const BandMatrixView<std::complex<float> >& A, int* P,
        std::complex<float>& det)
    {
        TMVAssert(A.isSquare());
        TMVAssert(A.iscm());
        TMVAssert(A.ct()==NonConj);
        int n = A.rowsize();
        int kl = A.nlo();
        int ku = A.nhi()-kl;
        int lda = A.stepj()+1;
        auto_array<int> lap_p(new int[n]);
        LAPNAME(cgbtrf) (
            LAPCM LAPV(n),LAPV(n),LAPV(kl),LAPV(ku),
            LAPP(A.ptr()-A.nhi()),LAPV(lda),LAPP(lap_p.get()) LAPINFO );
        LAP_Results("cgbtrf");
        const int M = A.colsize();
        for(int i=0;i<M;i++) {
            P[i] = (lap_p.get())[i] LAPMINUS1;
            if (P[i]!=i) det = -det;
        }
    }
#endif 
    // Now Lap version of DiagMajor Tridiagonal:
    template <class T> 
    static inline void LapTriDiagLU_Decompose(
        const BandMatrixView<T>& A, int* P, T& det)
    { NonLapTriDiagLU_Decompose(A,P,det); }
#ifdef INST_DOUBLE
    template <> 
    void LapTriDiagLU_Decompose(
        const BandMatrixView<double>& A, int* P, double& det)
    {
        TMVAssert(A.isSquare());
        TMVAssert(A.isdm());
        TMVAssert(A.nlo()==1);
        TMVAssert(A.nhi()==2);
        TMVAssert(A.ct()==NonConj);
        int n = A.colsize();
        auto_array<int> lap_p(new int[n]);
        LAPNAME(dgttrf) (
            LAPCM LAPV(n),LAPP(A.diag(-1).ptr()),LAPP(A.diag().ptr()),
            LAPP(A.diag(1).ptr()),LAPP(A.diag(2).ptr()),LAPP(lap_p.get()) LAPINFO );
        LAP_Results("dgttrf");
        const int M = A.colsize();
        for(int i=0;i<M;i++) {
            P[i] = (lap_p.get())[i] LAPMINUS1;
            if (P[i]!=i) det = -det;
        }
    }
    template <> 
    void LapTriDiagLU_Decompose(
        const BandMatrixView<std::complex<double> >& A, int* P,
        std::complex<double>& det)
    {
        TMVAssert(A.isSquare());
        TMVAssert(A.isdm());
        TMVAssert(A.nlo()==1);
        TMVAssert(A.nhi()==2);
        TMVAssert(A.ct()==NonConj);
        int n = A.colsize();
        auto_array<int> lap_p(new int[n]);
        LAPNAME(zgttrf) (
            LAPCM LAPV(n),LAPP(A.diag(-1).ptr()),LAPP(A.diag().ptr()), 
            LAPP(A.diag(1).ptr()),LAPP(A.diag(2).ptr()),LAPP(lap_p.get()) LAPINFO);
        LAP_Results("zgttrf");
        const int M = A.colsize();
        for(int i=0;i<M;i++) {
            P[i] = (lap_p.get())[i] LAPMINUS1;
            if (P[i]!=i) det = -det;
        }
    }
#endif
#ifdef INST_FLOAT
    template <> 
    void LapTriDiagLU_Decompose(
        const BandMatrixView<float>& A, int* P, float& det)
    {
        TMVAssert(A.isSquare());
        TMVAssert(A.isdm());
        TMVAssert(A.nlo()==1);
        TMVAssert(A.nhi()==2);
        TMVAssert(A.ct()==NonConj);
        int n = A.colsize();
        auto_array<int> lap_p(new int[n]);
        LAPNAME(sgttrf) (
            LAPCM LAPV(n),LAPP(A.diag(-1).ptr()),LAPP(A.diag().ptr()),
            LAPP(A.diag(1).ptr()),LAPP(A.diag(2).ptr()),LAPP(lap_p.get()) LAPINFO );
        LAP_Results("sgttrf");
        const int M = A.colsize();
        for(int i=0;i<M;i++) {
            P[i] = (lap_p.get())[i] LAPMINUS1;
            if (P[i]!=i) det = -det;
        }
    }
    template <> 
    void LapTriDiagLU_Decompose(
        const BandMatrixView<std::complex<float> >& A, int* P,
        std::complex<float>& det)
    {
        TMVAssert(A.isSquare());
        TMVAssert(A.isdm());
        TMVAssert(A.nlo()==1);
        TMVAssert(A.nhi()==2);
        TMVAssert(A.ct()==NonConj);
        int n = A.colsize();
        auto_array<int> lap_p(new int[n]);
        LAPNAME(cgttrf) (
            LAPCM LAPV(n),LAPP(A.diag(-1).ptr()),LAPP(A.diag().ptr()), 
            LAPP(A.diag(1).ptr()),LAPP(A.diag(2).ptr()),LAPP(lap_p.get()) LAPINFO);
        LAP_Results("cgttrf");
        const int M = A.colsize();
        for(int i=0;i<M;i++) {
            P[i] = (lap_p.get())[i] LAPMINUS1;
            if (P[i]!=i) det = -det;
        }
    }
#endif 
#endif // LAP

    template <class T> 
    void LU_Decompose(
        const BandMatrixView<T>& A, int* P, T& det, int 
#ifdef LAP
        Anhi
#endif
    )
    {
#ifdef XDEBUG
        BandMatrix<T> A0 = A;
#endif

        TMVAssert(A.isSquare());
        TMVAssert(A.ct()==NonConj);
        TMVAssert(A.nlo()>0);
        TMVAssert(A.iscm() || (A.isdm() && A.nlo()==1 && A.nhi()==2));
        if (A.colsize() > 0 && A.rowsize() > 0) {
            if (A.iscm()) {
#ifdef LAP
                if (A.nlo()+Anhi+1 > int(A.rowsize())) {
                    TMVAssert(A.nhi()+1 == int(A.rowsize()));
                    NonLapBandLU_Decompose(A,P,det);
                } else {
                    LapBandLU_Decompose(A,P,det);
                }
#else
                NonLapBandLU_Decompose(A,P,det);
#endif
            } else {
#ifdef LAP
                LapTriDiagLU_Decompose(A,P,det);
#else
                NonLapTriDiagLU_Decompose(A,P,det);
#endif
            }
        }

#ifdef XDEBUG
        int N = A.colsize();
        int nlo = A.nlo();
        Matrix<T> L(N,N,T(0));
        L.diag().setAllTo(T(1));
        for(int i=0;i<N;++i) {
            Swap(L.row(i,0,i),L.row(P[i],0,i));
            int end = TMV_MIN(i+nlo+1,N);
            L.col(i,i+1,end) = A.col(i,i+1,end);
        }
        Matrix<T> U(BandMatrixViewOf(A,0,A.nhi()));
        Matrix<T> AA = L*U;
        AA.reversePermuteRows(P);
        if (!(Norm(AA-A0) < 0.001 * Norm(A0))) {
            cerr<<"LU_Decompose: A = "<<TMV_Text(A)<<A0<<endl;
            cerr<<"AA = "<<AA<<endl;
            cerr<<"Norm(diff) = "<<Norm(A0-AA)<<endl;
            cerr<<"BandLU = "<<A<<endl;
            cerr<<"L = "<<L<<endl;
            cerr<<"U = "<<U<<endl;
            cerr<<"P = ";
            for(int i=0;i<N;i++) cerr<<P[i]<<" ";
            cerr<<endl;
#ifdef LAP
            BandMatrix<T,ColMajor> A2 = A0;
            auto_array<int> P2(new int[A.colsize()]);
            T det2(0);
            NonLapBandLU_Decompose(A2.view(),P2.get(),det2);
            cerr<<"NonLap version = "<<A2<<endl;
            cerr<<"P2 = ";
            for(int i=0;i<N;i++) cerr<<(P2.get())[i]<<" ";
            cerr<<endl;
#endif
            abort();
        }
#endif
    }

    template <class T> 
    void LU_Decompose(
        const GenBandMatrix<T>& A, const LowerTriMatrixView<T>& L,
        const BandMatrixView<T>& U, int* P)
    {
        TMVAssert(L.size() == A.colsize());
        TMVAssert(L.size() == A.rowsize());
        TMVAssert(L.size() == U.size());
        TMVAssert(U.nhi() >= A.nlo() + A.nhi());
        T d(0);
        BandMatrix<T> LU(A.colsize(),A.rowsize(),A.nlo(),A.nlo()+A.nhi());
        LU = A;
        LU_Decompose(LU.view(),P,d,A.nhi());
        U = LU.upperBand();
        if (L.isunit())
            LU_PackedPL_Unpack(LU,P,L);
        else {
            L.diag().setAllTo(T(1));
            LU_PackedPL_Unpack(LU,P,L.viewAsUnitDiag());
        }
    }

#define InstFile "TMV_BandLUDecompose.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv

