

#ifndef TMV_BandLUDecompose_H
#define TMV_BandLUDecompose_H

#include "TMV_BaseMatrix_Rec.h"
#include "TMV_BaseMatrix_Tri.h"
#include "TMV_MultMV.h"
#include "TMV_MultMM.h"
#include "TMV_DivVU.h"
#include "TMV_DivMU.h"
#include "TMV_Permutation.h"

#ifdef PRINTALGO_BandLU
#include <iostream>
#include "TMV_MatrixIO.h"
#include "TMV_VectorIO.h"
#include "TMV_TriMatrixIO.h"
#include "TMV_ProdMM.h"
#include "TMV_MultUL.h"
#include "TMV_PermuteM.h"
#endif

// INLINE_MV = Inline the MV (MultMV, LDivEqVU) calls.
#if TMV_OPT >= 1
#define TMV_BandLU_INLINE_VVM
#endif

namespace tmv {

    // Defined in TMV_BandLUDecompose.cpp
    template <class T>
    void InstBandLU_Decompose(BandMatrixView<T> m, ptrdiff_t* P);

    template <int algo, ptrdiff_t cs, ptrdiff_t rs, class M>
    struct BandLUDecompose_Helper;

    // algo 0: Trivial, nothing to do (M == 0 or 1, or N == 0)
    template <ptrdiff_t cs, ptrdiff_t rs, class M>
    struct BandLUDecompose_Helper<0,cs,rs,M>
    { static TMV_INLINE void call(M& A, ptrdiff_t* P) {} };

    // algo 11: Regular partial pivoting algorithm
    template <ptrdiff_t cs, ptrdiff_t rs, class M1>
    struct BandLUDecompose_Helper<11,cs,rs,M1>
    {
        static void call(M1& A, ptrdiff_t* P)
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
            typedef typename M1::value_type T;
            typedef typename M1::real_type RT;

            const ptrdiff_t N = rs==Unknown ? A.rowsize() : rs;
            const ptrdiff_t xx = Unknown;
#ifdef PRINTALGO_LU
            std::cout<<"BandLUDecompose algo 11: M,N,cs,rs = "<<A.colsize()<<','<<N<<
                ','<<cs<<','<<rs<<std::endl;
#endif

            typedef typename M1::col_sub_type M1c;
            typedef typename M1::const_col_sub_type M1cc;
            typedef typename M1::const_row_sub_type M1r;
            typedef typename M1::submatrix_type M1s;
            typedef typename M1::diag_type::iterator IT;

            const Scaling<-1,RT> mone;

#ifdef TMV_BandLU_INLINE_VVM
            const int algo2 = -4;
#else
            const int algo2 = -2;
#endif

            ptrdiff_t endcol = A.nlo()+1;
            ptrdiff_t endrow = A.nhi()+1;
            IT Ajj = A.diag().begin().nonConj();
            for (ptrdiff_t j=0; j<N-1; ++j, ++Ajj) {
                M1c Ajb = A.get_col(j,j,endcol);
                M1c Ajc = A.get_col(j,j+1,endcol);
                M1r Ajr = A.get_row(j,j+1,endrow);
                M1s Ajs = A.cSubMatrix(j+1,endcol,j+1,endrow);

                // Find the pivot element
                ptrdiff_t ip;
                // ip is relative to j index, not absolute.
                RT piv = Ajb.maxAbsElement(&ip);

                // Check for underflow:
                if (TMV_Underflow(piv)) {
                    P[j] = j;
                    Ajb.setZero();
                    continue;
                }

                // Swap the pivot row with j if necessary
                if (ip != 0) {
                    ip += j;
                    TMVAssert(ip < A.colsize());
                    TMVAssert(j < A.colsize());
                    // This does both Lkb and A'
                    Swap(A.row(ip,j,endrow),A.row(j,j,endrow));
                    P[j] = ip;
                } else P[j] = j;

                // If A(j,j) is 0, then all of the L's are 0.
                // ie. Ujj Lij = 0 for all i>j
                // Any value for Lij is valid, so leave them 0.
                if (*Ajj != T(0)) Ajc /= *Ajj;

                //A.subMatrix(j+1,endcol,j+1,endrow) -= (A.col(j,j+1,endcol) ^ A.row(j,j+1,endrow));
                Rank1VVM_Helper<algo2,xx,xx,true,-1,RT,M1cc,M1r,M1s>::call(mone,Ajc,Ajr,Ajs);

                if (endcol < N) ++endcol;
                if (endrow < N) ++endrow;
            }
            // j == N-1
            P[N-1] = N-1;
        }
    };

    // algo 31: Tridiagonal bandmatrix: nlo = nhi = 1
    template <ptrdiff_t cs, ptrdiff_t rs, class M1>
    struct BandLUDecompose_Helper<31,cs,rs,M1>
    {
        static void call(M1& A, ptrdiff_t* P)
        {
            // LU Decompostion for TriDiagonal BandMatrix
            //
            // For TriDiagonal BandMatrices, there are only two choices for 
            // each pivot, so we can specialize some of the calculations.
            // Otherwise, this is the same algorithm as above.
            //
            const ptrdiff_t N = rs==Unknown ? A.rowsize() : rs;
            const ptrdiff_t xx = Unknown;
#ifdef PRINTALGO_LU
            std::cout<<"BandLUDecompose algo 11: M,N,cs,rs = "<<A.colsize()<<','<<N<<
                ','<<cs<<','<<rs<<std::endl;
#endif
            typedef typename M1::value_type T;
            typedef typename M1::real_type RT;

            typedef typename M1::col_sub_type M1c;
            typedef typename M1::const_col_sub_type M1cc;
            typedef typename M1::const_row_sub_type M1r;
            typedef typename M1::submatrix_type M1s;
            typedef typename M1::diag_type::iterator IT;

            const Scaling<-1,RT> mone;

            IT Ujj = A.diag().begin().nonConj();  // U(j,j)
            IT Lj = Ujj; Lj.shiftP(A.stepi());    // L(j+1,j)
            IT Uj1 = Ujj; Uj1.shiftP(A.stepj());  // U(j,j+1)
            IT Uj2 = Uj1; Uj2.shiftP(A.stepj());  // U(j,j+2)

            for (ptrdiff_t j=0; j<N-1; ++j,++Ujj,++Lj,++Uj1,++Uj2) {
                bool pivot = TMV_ABS(*Lj) > TMV_ABS(*Ujj);

                if (pivot) {
                    //swap(A.row(j+1,j,endrow),A.row(j,j,endrow));
                    TMV_SWAP(*Lj,*Ujj);
                    if (j+1<N) TMV_SWAP(*(Ujj+1),*Uj1);
                    if (j+2<N) TMV_SWAP(*(Uj1+1),*Uj2);
                    P[j] = j+1;
                } else P[j] = j;

                if (TMV_Underflow(*Ujj)) {
                    *Ujj = *Lj = T(0);
                } else {
                    *Lj /= *Ujj;
                }

                //A.row(j+1,j+1,endrow) -= *Lj * A.row(j,j+1,endrow);
                *(Ujj+1) -= *Lj * (*Uj1);
                if (pivot && j+2<N) *(Uj1+1) -= *Lj * (*Uj2); // (if !pivot, *Uj2 == 0)
            }
            // j == N-1
            P[j] = N-1;
        }
    };

    // algo 71: Check whether m is tri-diagonal or not
    template <ptrdiff_t cs, ptrdiff_t rs, class M>
    struct BandLUDecompose_Helper<81,cs,rs,M>
    {
        static void call(M& m, ptrdiff_t* P)
        {
#ifdef PRINTALGO_BandLU
            std::cout<<"BandLUDecompose algo 71: cs,rs = "<<cs<<','<<rs<<std::endl;
#endif
            const ptrdiff_t lo = m._nlo==Unknown ? m.nlo() : m._nlo;
            const ptrdiff_t hi = m._nhi==Unknown ? m.nhi() : m._nhi;
            // Note: hi is the original lo+hi, so if the original was tridiagonal, then hi==2.
            if (lo == 1 && hi == 2) 
                BandLUDecompose_Helper<11,cs,rs,Mcm>::call(m,P);
            else
                BandLUDecompose_Helper<21,cs,rs,Mcm>::call(m,P);
        }
    };

    // algo 81: Copy to colmajor
    template <ptrdiff_t cs, ptrdiff_t rs, class M>
    struct BandLUDecompose_Helper<81,cs,rs,M>
    {
        static void call(M& m, ptrdiff_t* P)
        {
#ifdef PRINTALGO_BandLU
            std::cout<<"BandLUDecompose algo 81: cs,rs = "<<cs<<','<<rs<<std::endl;
#endif
            typedef typename M::value_type T;
            typedef typename MCopyHelper<T,Band,cs,rs,ColMajor>::type Mcm;
            Mcm mcm = m;
            BandLUDecompose_Helper<-2,cs,rs,Mcm>::call(mcm,P);
            m.noAlias() = mcm;
        }
    };

    // algo 82: Copy to diagmajor
    template <ptrdiff_t cs, ptrdiff_t rs, class M>
    struct BandLUDecompose_Helper<82,cs,rs,M>
    {
        static void call(M& m, ptrdiff_t* P)
        {
#ifdef PRINTALGO_BandLU
            std::cout<<"BandLUDecompose algo 82: cs,rs = "<<cs<<','<<rs<<std::endl;
#endif
            typedef typename M::value_type T;
            typedef typename MCopyHelper<T,Band,cs,rs,DiagMajor>::type Mdm;
            Mdm mdm = m;
            BandLUDecompose_Helper<-2,cs,rs,Mdm>::call(mdm,P);
            m.noAlias() = mdm;
        }
    };

    // algo 90: call InstBandLU_Decompose
    template <ptrdiff_t cs, ptrdiff_t rs, class M>
    struct BandLUDecompose_Helper<90,cs,rs,M>
    {
        static TMV_INLINE void call(M& m, ptrdiff_t* P)
        { InstBandLU_Decompose(m.xView(),P); }
    };

    // algo 97: Conjugate
    template <ptrdiff_t cs, ptrdiff_t rs, class M>
    struct BandLUDecompose_Helper<97,cs,rs,M>
    {
        static TMV_INLINE void call(M& m, ptrdiff_t* P)
        {
            typedef typename M::conjugate_type Mc;
            Mc mc = m.conjugate();
            BandLUDecompose_Helper<-2,cs,rs,Mc>::call(mc,P);
        }
    };

    // algo -4: No copies or branches -- CM algorithm
    template <ptrdiff_t cs, ptrdiff_t rs, class M1>
    struct BandLUDecompose_Helper<-4,cs,rs,M1>
    {
        typedef typename M1::value_type T;
        static TMV_INLINE void call(M1& m, ptrdiff_t* P)
        {
            const int algo = 
                cs == 0 || rs == 0 || cs == 1 ? 0 :
                TMV_OPT == 0 ? 11 :
                (M::_nlo == 1 && M::_nhi == 2) ? 21 :
                11;
#ifdef PRINTALGO_BandLU
            std::cout<<"Inline BandLUDecompose: \n";
            std::cout<<"m = "<<TMV_Text(m)<<std::endl;
            std::cout<<"cs = "<<cs<<"  rs = "<<rs<<std::endl;
            std::cout<<"sizes = "<<m.colsize()<<"  "<<m.rowsize()<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
#endif
#ifdef XDEBUG_BandLU
            typedef typename M3::real_type RT;
            Matrix<T> mi = m;
#endif
            BandLUDecompose_Helper<algo,cs,rs,M1>::call(m,P);
#ifdef XDEBUG_BandLU
            const ptrdiff_t N = m.colsize();
            Matrix<T> L(N,N,T(0));
            L.diag().setAllTo(T(1));
            for(ptrdiff_t i=0;i<N;++i) {
                Swap(L.row(i,0,i),L.row(P[i],0,i));
                ptrdiff_t end = TMV_MIN(i+nlo+1,N);
                L.col(i,i+1,end) = m.col(i,i+1,end);
            }
            Matrix<T> U = BandMatrixViewOf(m,0,m.nhi());
            Matrix<T> lu = L*U;
            Matrix<T> plu = lu;
            plu.reversePermuteRows(P);
            if (Norm(plu-mi) > 1.e-3*Norm(mi)) {
                std::cout<<"m => "<<m<<std::endl;
                std::cout<<"L = "<<L<<std::endl;
                std::cout<<"U = "<<U<<std::endl;
                std::cout<<"L*U = "<<lu<<std::endl;
                std::cout<<"P*L*U = "<<plu<<std::endl;
                abort();
            }
#endif
        }
    };

    // algo -3: Determine which algorithm to use
    template <ptrdiff_t cs, ptrdiff_t rs, class M1>
    struct BandLUDecompose_Helper<-3,cs,rs,M1>
    {
        static TMV_INLINE void call(M1& m, ptrdiff_t* P)
        {
            const int cmalgo = (
                ( cs != Unknown && rs != Unknown && cs <= 16 && rs <= 16 ) ? -4 :
                !M1::_colmajor ? 81 :
                -4 );
            const int dmalgo = (
                ( cs != Unknown && rs != Unknown && cs <= 16 && rs <= 16 ) ? -4 :
                !M1::_diagmajor ? 82 :
                -4 );

            const int algo = 
                cs == 0 || rs == 0 || cs == 1 ? 0 :
                TMV_OPT == 0 ? 11 :
                (M::_nlo == 1 && M::_nhi == 2) ? dmalgo :
                (N::_nlo != 1 && M::_nlo != Unknown) ? cmalgo :
                (N::_nhi != 2 && M::_nhi != Unknown) ? cmalgo :
                (M::_nlo == Unknown || M::_nhi == Unknown) ? 71 :
                cmalgo;
#ifdef PRINTALGO_BandLU
            const ptrdiff_t M = cs==Unknown ? m.colsize() : cs;
            const ptrdiff_t N = rs==Unknown ? m.rowsize() : rs;
            std::cout<<"BandLUDecompose algo -3: M,N,cs,rs = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<std::endl;
#endif
            BandLUDecompose_Helper<algo,cs,rs,M1>::call(m,P);
        }
    };

    // algo -2: Check for inst
    template <ptrdiff_t cs, ptrdiff_t rs, class M>
    struct BandLUDecompose_Helper<-2,cs,rs,M>
    {
        static TMV_INLINE void call(M& m, ptrdiff_t* P)
        {
            typedef typename M::value_type T;
            const bool inst = 
                (cs == Unknown || cs > 16) &&
                (rs == Unknown || rs > 16) &&
                Traits<T>::isinst;
            const int algo = 
                cs == 0 || rs == 0 || cs == 1 ? 0 :
                M::_conj ? 97 :
                inst ? 90 :
                -3;
            BandLUDecompose_Helper<algo,cs,rs,M>::call(m,P);
        }
    };

    template <ptrdiff_t cs, ptrdiff_t rs, class M>
    struct BandLUDecompose_Helper<-1,cs,rs,M>
    {
        static TMV_INLINE void call(M& m, ptrdiff_t* P)
        { BandLUDecompose_Helper<-2,cs,rs,M>::call(m,P); }
    };

    template <class M>
    inline void InlineBandLU_Decompose(
        BaseMatrix_Band_Mutable<M>& m, ptrdiff_t* P)
    {
        const ptrdiff_t cs = M::_colsize;
        const ptrdiff_t rs = M::_rowsize;
        typedef typename M::cview_type Mv;
        TMV_MAYBE_REF(M,Mv) mv = m.cView();
        BandLUDecompose_Helper<-3,cs,rs,Mv>::call(mv,P);
    }

    template <class M>
    inline void BandLU_Decompose(
        BaseMatrix_Band_Mutable<M>& m, ptrdiff_t* P)
    {
        const ptrdiff_t cs = M::_colsize;
        const ptrdiff_t rs = M::_rowsize;
        typedef typename M::cview_type Mv;
        TMV_MAYBE_REF(M,Mv) mv = m.cView();
        BandLUDecompose_Helper<-2,cs,rs,Mv>::call(mv,P);
    }

    // This function is a friend of Permutation class.
    template <class M>
    inline void BandLU_Decompose(BaseMatrix_Band_Mutable<M>& m, Permutation& P)
    {
        TMVAssert(P.size() == m.colsize());
        P.allocateMem();
        BandLU_Decompose(m,P.getMem());
        P.isinv = true;
    }

    // Allow views as an argument by value (for convenience)
    template <class T, int A>
    TMV_INLINE void BandLU_Decompose(BandMatrixView<T,A> m, Permutation& P)
    {
        typedef MatrixView<T,A> M;
        BandLU_Decompose(static_cast<BaseMatrix_Band_Mutable<M>&>(m),P); 
    }
    template <class T, ptrdiff_t M, ptrdiff_t N, ptrdiff_t lo, ptrdiff_t hi, ptrdiff_t Si, ptrdiff_t Sj, int A>
    TMV_INLINE void BandLU_Decompose(SmallBandMatrixView<T,M,N,lo,hi,Si,Sj,A> m, Permutation& P)
    {
        typedef SmallMatrixView<T,M,N,lo,hi,Si,Sj,A> MM;
        BandLU_Decompose(static_cast<BaseMatrix_Band_Mutable<MM>&>(m),P); 
    }

} // namespace tmv

#undef TMV_BandLU_INLINE_VVM

#endif

