

#ifndef TMV_InvertU_H
#define TMV_InvertU_H

#include "TMV_BaseMatrix_Tri.h"
#include "TMV_Scaling.h"
#include "TMV_MultUV.h"
#include "TMV_MultUM.h"
#include "TMV_InvertD.h"

#ifdef PRINTALGO_InvU
#include <iostream>
#include "tmv/TMV_TriMatrixIO.h"
#include "tmv/TMV_MatrixIO.h"
#endif

// The maximum nops to unroll
#if TMV_OPT >= 3
#define TMV_INVU_UNROLL 200 
#elif TMV_OPT >= 2
#define TMV_INVU_UNROLL 25
#elif TMV_OPT >= 1
#define TMV_INVU_UNROLL 9
#else
#define TMV_INVU_UNROLL 0
#endif

// Inline the MV (MultUV) calls.
#if TMV_OPT >= 1
#define TMV_INVU_INLINE_MV
#endif

// The size to stop recursing.
#if TMV_OPT >= 3
#define TMV_INVU_RECURSE 8
#else
#define TMV_INVU_RECURSE 1
#endif

namespace tmv {

    // Defined in TMV_InvertU.cpp
    template <class T>
    void InstInvertSelf(UpperTriMatrixView<T> m);
    template <class T>
    void InstInvertSelf(LowerTriMatrixView<T> m);

    //
    // U = U^-1
    //

    template <int algo, ptrdiff_t s, class M>
    struct InvertU_Helper;

    // algo 0: nothing to do (either s==0, or invert diag of unitdiag m.)
    template <ptrdiff_t s, class M>
    struct InvertU_Helper<0,s,M>
    { static TMV_INLINE void call(M& ) {} };

    // algo 1: s == 1, so simplifies to a scalar quotient
    template <ptrdiff_t s, class M>
    struct InvertU_Helper<1,s,M>
    {
        // m.diag is already inverted, so nothing to do.
        static TMV_INLINE void call(M& ) {}
    };

    // algo 2: transpose
    template <ptrdiff_t s, class M>
    struct InvertU_Helper<2,s,M>
    {
        static TMV_INLINE void call(M& m)
        {
#ifdef PRINTALGO_InvU
            ptrdiff_t N = (s == Unknown ? m.size() : s);
            std::cout<<"InvU algo 2: N,s,x = "<<N<<','<<s<<std::endl;
#endif
            typedef typename M::transpose_type Mt;
            Mt mt = m.transpose();
            InvertU_Helper<-4,s,Mt>::call(mt);
        }
    };

    // algo 11: column major loop 
    template <ptrdiff_t s, class M>
    struct InvertU_Helper<11,s,M>
    {
        static void call(M& m)
        {
            ptrdiff_t N = (s == Unknown ? m.size() : s);
#ifdef PRINTALGO_InvU
            std::cout<<"InvU algo 11: N,s,x = "<<N<<','<<s<<std::endl;
#endif
            const bool u = M::_unit;
            typedef typename M::col_sub_type Mc;
            typedef typename M::const_col_sub_type Mcc;
            typedef typename M::const_subtrimatrix_type Mst;
            const int ix = u?-1:0;
            typedef typename M::value_type T;
            typedef typename M::real_type RT;
            typedef typename TypeSelect<u,RT,T>::type XT;
            const ptrdiff_t xx = Unknown;
#ifdef TMV_INVU_INLINE_MV
            const int algo2 = -4;
#else
            const int algo2 = -2;
#endif

            for(ptrdiff_t j=0;j<N;++j) {
                // m.col(j,0,j) = -m.subTraiMatrix(0,j)^-1 *
                //                   m.col(j,0,j) * m(j,j)^-1
                Mc mj = m.get_col(j,0,j);
                Mst mst = m.cSubTriMatrix(0,j); 
                Scaling<ix,XT> invAjj(-Maybe<u>::real(m.cref(j,j)));
                MultUV_Helper<algo2,xx,false,ix,XT,Mst,Mcc,Mc>::call(
                    invAjj,mst,mj,mj);
            }
        }
    };

    // algo 12: row major loop
    template <ptrdiff_t s, class M>
    struct InvertU_Helper<12,s,M>
    {
        static void call(M& m)
        {
            ptrdiff_t N = (s == Unknown ? m.size() : s);
#ifdef PRINTALGO_InvU
            std::cout<<"InvU algo 11: N,s,x = "<<N<<','<<s<<std::endl;
#endif
            const bool u = M::_unit;
            typedef typename M::row_sub_type Mr;
            typedef typename M::const_row_sub_type Mcr;
            typedef typename M::const_subtrimatrix_type Mst;
            typedef typename Mst::const_transpose_type Mstct;
            const int ix = u?-1:0;
            typedef typename M::value_type T;
            typedef typename M::real_type RT;
            typedef typename TypeSelect<u,RT,T>::type XT;
            const ptrdiff_t xx = Unknown;
#ifdef TMV_INVU_INLINE_MV
            const int algo2 = -4;
#else
            const int algo2 = -2;
#endif

            for(ptrdiff_t i=N-1;i>=0;--i) {
                // m.row(i,i+1,N) = -(m(i,i)^-1) * m.row(i,i+1,N) * 
                //                      m.subTraiMatrix(i+1,N)^-1
                Mr mi = m.get_row(i,i+1,N);
                Mstct mstt = m.cSubTriMatrix(i+1,N).transpose();
                Scaling<ix,XT> invAii(-Maybe<u>::real(m.cref(i,i)));
                MultUV_Helper<algo2,xx,false,ix,XT,Mstct,Mcr,Mr>::call(
                    invAii,mstt,mi,mi);
            }
        }
    };

    // algo 16: UpperTri: Unroll small case
    template <ptrdiff_t s, class M>
    struct InvertU_Helper<16,s,M>
    {
        template <ptrdiff_t I, ptrdiff_t N>
        struct Unroller
        {
            static inline void unroll(M& m)
            {
                // [ B00 B01 ] * [ A00 A01 ] = [ 1 0 ]
                // [  0  B11 ]   [  0  A11 ]   [ 0 1 ]
                // B00 A00 = 1
                // B00 A01 + B01 A11 = 0
                // B11 A11 = 1
                //
                // B11 = A11^-1
                // B00 = A00^-1
                // B01 = -B00 A01 B11

                const ptrdiff_t Nx = N/2;
                const ptrdiff_t Ny = N-Nx;
                const ptrdiff_t I1 = I+Nx;
                const ptrdiff_t I2 = I+N;
                typedef typename M::submatrix_type Msm;
                typedef typename M::subtrimatrix_type Mst;
                typedef typename M::const_submatrix_type Msmc;
                typedef typename M::const_subtrimatrix_type Mstc;
                typedef typename Msm::transpose_type Msmt;
                typedef typename Mst::const_transpose_type Mstct;
                typedef typename Msm::const_transpose_type Msmct;
                typedef typename M::real_type RT;

                Mst A00 = m.cSubTriMatrix(I,I1);
                Msm A01 = m.cSubMatrix(I,I1,I1,I2);
                Msmt A01t = A01.transpose();
                Msmct A01ct = A01.transpose();
                Mstct A11t = m.cSubTriMatrix(I1,I2).transpose();

                // B00 = A00^-1
                Unroller<I,Nx>::unroll(m);

                // B11 = A11^-1
                Unroller<I1,Ny>::unroll(m);

                // B01 = -B00 A01
                Scaling<-1,RT> mone;
                MultUM_Helper<-4,Nx,Ny,false,-1,RT,Mstc,Msmc,Msm>::call(
                    mone,A00,A01,A01);

                // B01 = B01 B11
                Scaling<1,RT> one;
                MultUM_Helper<-4,Ny,Nx,false,1,RT,Mstct,Msmct,Msmt>::call(
                    one,A11t,A01ct,A01t);
            }
        };
        template <ptrdiff_t I>
        struct Unroller<I,1> // diagonal is already done, so nothing to do.
        { static TMV_INLINE void unroll(M&) {} };
        template <ptrdiff_t I>
        struct Unroller<I,0>
        { static TMV_INLINE void unroll(M&) {} };
        static inline void call(M& m)
        {
#ifdef PRINTALGO_InvU
            const ptrdiff_t N = m.size();
            std::cout<<"InvU algo 16: N,s,x = "<<N<<','<<s<<std::endl;
#endif
            Unroller<0,s>::unroll(m); 
        }
    };

    // algo 17: Split U into 3 sections and recurse 
    // TODO: Combine these the way I do for MultUU, etc.
    template <ptrdiff_t s, class M>
    struct InvertU_Helper<17,s,M>
    {
        static void call(M& m)
        {
            const ptrdiff_t N = m.size();
#ifdef PRINTALGO_InvU
            std::cout<<"InvU algo 17: N,s,x = "<<N<<','<<s<<std::endl;
#endif
            const ptrdiff_t sp1 = IntTraits<s>::Sp1;
            const ptrdiff_t sp2 = IntTraits<sp1>::Sp1;
            // nops = 1/6 n(n+1)(n+2)
            const ptrdiff_t nops =
                IntTraits2<IntTraits2<s,sp1>::safeprod,sp2>::safeprod/6;
            const bool unroll =
                s > 10 ? false :
                s == Unknown ? false :
                nops <= TMV_INVU_UNROLL;
            const int algo2 =
                s == 0 ? 0 :
                s == 1 ? ( M::_unit ? 0 : 1 ) :
                unroll ? 16 :
                // For known s, always recurse down to unroll size
                s != Unknown ? 0 :
                (TMV_INVU_RECURSE == 1) ? ( M::_unit ? 0 : 1 ) :
                M::_rowmajor ? 12 : 11;
            const int algo3 =  // The algorithm for N > INVU_RECURSE
                unroll || s == 1 ? 0 : 17;
            const int algo4 =  // The algorithm for MultUM
                unroll || s == 1 ? 0 : -2;
#ifdef PRINTALGO_InvU
            std::cout<<"algo2,3,4 = "<<algo2<<"  "<<algo3<<"  "<<algo4<<std::endl;
#endif

            if (s==Unknown ? N > TMV_INVU_RECURSE : (s > 1 && !unroll)) {
                // [ B00 B01 ] * [ A00 A01 ] = [ 1 0 ]
                // [  0  B11 ]   [  0  A11 ]   [ 0 1 ]
                // B00 A00 = 1
                // B00 A01 + B01 A11 = 0
                // B11 A11 = 1
                //
                // B11 = A11^-1
                // B00 = A00^-1
                // B01 = -B00 A01 B11

                const ptrdiff_t Nx = N > 16 ? ((((N-1)>>5)+1)<<4) : (N>>1);
                // (If N > 16, round N/2 up to a multiple of 16.)
                const ptrdiff_t sx = IntTraits<s>::half_roundup;
                const ptrdiff_t sy = IntTraits2<s,sx>::diff;

                typedef typename M::submatrix_type Msm;
                typedef typename M::subtrimatrix_type Mst;
                typedef typename Msm::transpose_type Msmt;
                typedef typename M::const_submatrix_type Msmc;
                typedef typename M::const_subtrimatrix_type Mstc;
                typedef typename Mst::const_transpose_type Mstct;
                typedef typename Msm::const_transpose_type Msmct;
                typedef typename M::real_type RT;

                Mst A00 = m.cSubTriMatrix(0,Nx);
                Msm A01 = m.cSubMatrix(0,Nx,Nx,N);
                Msmt A01t = A01.transpose();
                Mst A11 = m.cSubTriMatrix(Nx,N);
                Mstc A00c = A00;
                Msmct A01ct = A01.transpose();
                Mstct A11ct = A11.transpose();

                // B00 = A00^-1
                // F(n)
                InvertU_Helper<algo3,sx,Mst>::call(A00);

                // B11 = A11^-1
                // F(n)
                InvertU_Helper<algo3,sy,Mst>::call(A11);

                // B01 = -B00 A01
                // n*n*(n+1)/2
                Scaling<-1,RT> mone;
                MultUM_Helper<algo4,sx,sy,false,-1,RT,Mstc,Msmc,Msm>::call(
                    mone,A00c,A01,A01);

                // B01 = B01 B11
                // n*n*(n+1)/2
                Scaling<1,RT> one;
                MultUM_Helper<algo4,sy,sx,false,1,RT,Mstct,Msmct,Msmt>::call(
                    one,A11ct,A01ct,A01t);
            } else {
                InvertU_Helper<algo2,s,M>::call(m);
            }
        }
    };

    // algo 50: Invert the diagonal.
    // I invert the diagonal first, since it is more efficient to 
    // use the sse commands if possible.
    template <ptrdiff_t s, class M>
    struct InvertU_Helper<50,s,M>
    {
        static inline void call(M& m)
        {
            TMVStaticAssert(M::_upper);
            TMVStaticAssert(!M::_unit);
            typedef typename M::diag_type Md;
            Md md = m.diag();
            ElemInvert_Helper<-4,s,Md>::call(md);
        }
    };

    // algo 51: Maybe invert the diagonal.  (m is UnknownDiag)
    template <ptrdiff_t s, class M>
    struct InvertU_Helper<51,s,M>
    {
        static inline void call(M& m)
        {
            if (!m.isunit()) 
                InvertU_Helper<50,s,M>::call(m); 
        }
    };

    // algo 90: call InstInvertSelf
    template <ptrdiff_t s, class M>
    struct InvertU_Helper<90,s,M>
    {
        static TMV_INLINE void call(M& m)
        { InstInvertSelf(m.xView()); }
    };

    // algo 97: Conjugate
    template <ptrdiff_t s, class M>
    struct InvertU_Helper<97,s,M>
    {
        static TMV_INLINE void call(M& m)
        {
            typedef typename M::conjugate_type Mc;
            Mc mc = m.conjugate();
            InvertU_Helper<-1,s,Mc>::call(mc);
        }
    };

    // algo -4: No branches or copies
    template <ptrdiff_t s, class M>
    struct InvertU_Helper<-4,s,M>
    {
        static TMV_INLINE void call(M& m)
        {
            TMVStaticAssert(M::_upper);
            const ptrdiff_t s2 = s > 20 ? Unknown : s;
            const ptrdiff_t s2p1 = IntTraits<s2>::Sp1;
            const ptrdiff_t s2p2 = IntTraits<s2p1>::Sp1;
            // nops = 1/6 n(n+1)(n+2)
            const ptrdiff_t nops =
                IntTraits2<IntTraits2<s2,s2p1>::safeprod,s2p2>::safeprod / 6;
            const bool unroll = 
                s > 10 ? false :
                s == Unknown ? false :
                nops <= TMV_INVU_UNROLL;
            const int algo1 = 
                M::_unit ? 0 : M::_unknowndiag ? 51 : 50;
            const int algo2 = 
                s == 0 ? 0 :
                s == 1 ? ( M::_unit ? 0 : 1 ) :
                unroll ? 16 :
                TMV_OPT == 0 ? ( M::_rowmajor ? 12 : 11 ) :
                17;
#ifdef PRINTALGO_InvU
            std::cout<<"InlineInvertU (algo -4): \n";
            std::cout<<"m "<<TMV_Text(m)<<std::endl;
            std::cout<<"s,algo = "<<s<<"  "<<algo1<<"  "<<algo2<<std::endl;
#endif
            InvertU_Helper<algo1,s,M>::call(m);
            InvertU_Helper<algo2,s,M>::call(m);
        }
    };

    // algo -3: Determine which algorithm to use
    template <ptrdiff_t s, class M>
    struct InvertU_Helper<-3,s,M>
    {
        static TMV_INLINE void call(M& m)
        {
            // Possible algorithms to choose from:
            //
            // Trivial:
            //  0 = s == 0, so nothing to do
            //  1 = s == 1: reduces to trivial scalar quotient
            //  2 = transpose
            //
            // UpperTri:
            // 11 = column major
            // 12 = row major
            // 16 = unroll
            // 17 = recurse

            const int algo = 
                ( s == 0 ) ? 0 :
                M::_lower ? 2 :
                -4;
#ifdef PRINTALGO_InvU
            std::cout<<"InlineInvertU: \n";
            std::cout<<"m "<<TMV_Text(m)<<std::endl;
            std::cout<<"s,algo = "<<s<<"  "<<algo<<std::endl;
#endif
            if (m.isSingular()) ThrowSingular("TriMatrix");
            InvertU_Helper<algo,s,M>::call(m);
        }
    };

    // algo -2: Check for inst
    template <ptrdiff_t s, class M>
    struct InvertU_Helper<-2,s,M>
    {
        static TMV_INLINE void call(M& m)
        {
            typedef typename M::value_type T;
            const bool inst = 
                (s == Unknown || s > 16) &&
                Traits<T>::isinst;
            const int algo = 
                ( s == 0 ) ? 0 :
                M::_conj ? 97 :
                inst ? 90 : 
                -3;
            InvertU_Helper<algo,s,M>::call(m);
        }
    };

    template <ptrdiff_t s, class M>
    struct InvertU_Helper<-1,s,M>
    {
        static TMV_INLINE void call(M& m)
        { InvertU_Helper<-2,s,M>::call(m); }
    };

    template <class M>
    inline void InvertSelf(BaseMatrix_Tri_Mutable<M>& m)
    {
        typedef typename M::cview_type Mv;
        TMV_MAYBE_REF(M,Mv) mv = m.cView();
        InvertU_Helper<-2,M::_size,Mv>::call(mv);
    }

    template <class M>
    inline void InlineInvertSelf(BaseMatrix_Tri_Mutable<M>& m)
    {
        typedef typename M::cview_type Mv;
        TMV_MAYBE_REF(M,Mv) mv = m.cView();
        InvertU_Helper<-3,M::_size,Mv>::call(mv);
    }

#undef TMV_INVU_UNROLL
#undef TMV_INVU_RECURSE

} // namespace tmv

#endif 
