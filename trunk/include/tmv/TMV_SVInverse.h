

#ifndef TMV_SVInverse_H
#define TMV_SVInverse_H

#include "TMV_BaseMatrix_Rec.h"
#include "TMV_BaseMatrix_Diag.h"

#ifdef PRINTALGO_SVD
#include <iostream>
#include "TMV_MatrixIO.h"
#include "TMV_TriMatrixIO.h"
#endif

namespace tmv {

    // Defined in TMV_SVInverse.cpp
    template <class T1, int C1, class RT1, class T2>
    void InstSV_Inverse(
        const ConstMatrixView<T1,C1>& U, const ConstDiagMatrixView<RT1>& S,
        const ConstMatrixView<T1,C1>& V, MatrixView<T2> minv);
    template <class T1, int C1, class RT1, class T2>
    void InstSV_InverseATA(
        const ConstMatrixView<T1,C1>& U, const ConstDiagMatrixView<RT1>& S,
        const ConstMatrixView<T1,C1>& V, MatrixView<T2> ata);


    // Note: cs,rs refer to M1, not M2 (which is reverse of M1)
    template <int algo, int cs, int rs, class M1u, class M1s, class M1v, class M2>
    struct SV_Inverse_Helper;

    // algo 0: Trivial, nothing to do (M == 0 or N == 0)
    template <int cs, int rs, class M1u, class M1s, class M1v, class M2>
    struct SV_Inverse_Helper<0,cs,rs,M1u,M1s,M1v,M2>
    {
        static TMV_INLINE void call(
            const M1u& , const M1s& , const M1v& , int, M2& ) {} 
    };

    // algo 11: Normal case
    template <int cs, int rs, class M1u, class M1s, class M1v, class M2>
    struct SV_Inverse_Helper<11,cs,rs,M1u,M1s,M1v,M2>
    {
        static void call(
            const M1u& U, const M1s& S, const M1v& V, M2& minv)
        {
#ifdef PRINTALGO_SVD
            std::cout<<"SVInverse algo 11: cs,rs = "<<cs<<','<<rs<<std::endl;
#endif
            // m2 = (USV)^-1
            //    = Vt S^-1 Ut
            typedef typename M1u::value_type T;
            const int xx = TMV_UNKNOWN;
            typename MCopyHelper<T,Rec,xx,cs>::type SinvUt = U.adjoint() / S;
            minv = V.adjoint() * SinvUt;
        }
    };

    // algo 90: call InstSV_Inverse
    template <int cs, int rs, class M1u, class M1s, class M1v, class M2>
    struct SV_Inverse_Helper<90,cs,rs,M1u,M1s,M1v,M2>
    {
        static TMV_INLINE void call(
            const M1u& U, const M1s& S, const M1v& V, M2& minv)
        { InstSV_Inverse(U.xView(),S.xView(),V.xView(),minv.xView()); }
    };

    // algo 97: Conjugate
    template <int cs, int rs, class M1u, class M1s, class M1v, class M2>
    struct SV_Inverse_Helper<97,cs,rs,M1u,M1s,M1v,M2>
    {
        static TMV_INLINE void call(
            const M1u& U, const M1s& S, const M1v& V, M2& minv)
        {
            typedef typename M1u::const_conjugate_type M1uc;
            typedef typename M1v::const_conjugate_type M1vc;
            typedef typename M2::conjugate_type M2c;
            M1uc Uc = U.conjugate();
            M1vc Vc = V.conjugate();
            M2c minvc = minv.conjugate();
            SV_Inverse_Helper<-2,cs,rs,M1uc,M1s,M1vc,M2c>::call(
                Uc,S,Vc,minvc);
        }
    };

    // algo -3: Determine which algorithm to use
    template <int cs, int rs, class M1u, class M1s, class M1v, class M2>
    struct SV_Inverse_Helper<-3,cs,rs,M1u,M1s,M1v,M2>
    {
        static TMV_INLINE void call(
            const M1u& U, const M1s& S, const M1v& V, M2& minv)
        {
            const int algo = 
                cs == 0 || rs == 0 ? 0 : 
                11;
#ifdef PRINTALGO_SVD
            std::cout<<"Inline SVInverse\n";
            std::cout<<"cs,rs = "<<cs<<','<<rs<<std::endl;
            std::cout<<"U = "<<TMV_Text(U)<<std::endl;
            std::cout<<"S = "<<TMV_Text(S)<<std::endl;
            std::cout<<"V = "<<TMV_Text(V)<<std::endl;
            std::cout<<"minv = "<<TMV_Text(minv)<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
#endif
            SV_Inverse_Helper<algo,cs,rs,M1u,M1s,M1v,M2>::call(U,S,V,minv);
        }
    };

    // algo -2: Check for inst
    template <int cs, int rs, class M1u, class M1s, class M1v, class M2>
    struct SV_Inverse_Helper<-2,cs,rs,M1u,M1s,M1v,M2>
    {
        static TMV_INLINE void call(
            const M1u& U, const M1s& S, const M1v& V, M2& minv)
        {
            typedef typename M1u::value_type T1;
            typedef typename M2::value_type T2;
            const bool inst = 
                (cs == TMV_UNKNOWN || cs > 16) &&
                (rs == TMV_UNKNOWN || rs > 16) &&
#ifdef TMV_INST_MIX
                Traits2<T1,T2>::samebase &&
#else
                Traits2<T1,T2>::sametype &&
#endif
                Traits<T2>::isinst;
            const int algo = 
                cs == 0 || rs == 0 ? 0 : 
                M2::_conj ? 97 :
                inst ? 90 :
                -3;
            SV_Inverse_Helper<algo,cs,rs,M1u,M1s,M1v,M2>::call(U,S,V,minv);
        }
    };

    template <int cs, int rs, class M1u, class M1s, class M1v, class M2>
    struct SV_Inverse_Helper<-1,cs,rs,M1u,M1s,M1v,M2>
    {
        static TMV_INLINE void call(
            const M1u& U, const M1s& S, const M1v& V, M2& minv)
        { SV_Inverse_Helper<-2,cs,rs,M1u,M1s,M1v,M2>::call(U,S,V,minv); }
    };

    template <class M1u, class M1s, class M1v, class M2>
    inline void InlineSV_Inverse(
        const BaseMatrix_Rec<M1u>& U, const BaseMatrix_Diag<M1s>& S,
        const BaseMatrix_Rec<M1v>& V, BaseMatrix_Rec_Mutable<M2>& minv)
    {
        TMVStaticAssert((Sizes<M1u::_rowsize,M1s::_size>::same));
        TMVStaticAssert((Sizes<M1u::_rowsize,M1v::_colsize>::same));
        TMVStaticAssert((Sizes<M1u::_rowsize,M1v::_rowsize>::same));
        TMVStaticAssert((Sizes<M1u::_colsize,M2::_rowsize>::same));
        TMVStaticAssert((Sizes<M1u::_rowsize,M2::_colsize>::same));
        TMVAssert(U.rowsize() == S.size());
        TMVAssert(U.rowsize() == V.colsize());
        TMVAssert(U.colsize() == minv.rowsize());
        TMVAssert(V.rowsize() == minv.colsize());

        const int cs = Sizes<M2::_rowsize,M1u::_colsize>::size;
        const int rs = Sizes<M2::_colsize,M1v::_rowsize>::size;
        typedef typename M1u::const_cview_type M1uv;
        typedef typename M1s::const_cview_type M1sv;
        typedef typename M1v::const_cview_type M1vv;
        typedef typename M2::cview_type M2v;
        TMV_MAYBE_CREF(M1u,M1uv) Uv = U.cView();
        TMV_MAYBE_CREF(M1s,M1sv) Sv = S.cView();
        TMV_MAYBE_CREF(M1v,M1vv) Vv = V.cView();
        TMV_MAYBE_REF(M2,M2v) minvv = minv.cView();
        SV_Inverse_Helper<-3,cs,rs,M1uv,M1sv,M1vv,M2v>::call(Uv,Sv,Vv,minvv);
    }

    template <class M1u, class M1s, class M1v, class M2>
    inline void SV_Inverse(
        const BaseMatrix_Rec<M1u>& U, const BaseMatrix_Diag<M1s>& S,
        const BaseMatrix_Rec<M1v>& V, BaseMatrix_Rec_Mutable<M2>& minv)
    {
        TMVStaticAssert((Sizes<M1u::_rowsize,M1s::_size>::same));
        TMVStaticAssert((Sizes<M1u::_rowsize,M1v::_colsize>::same));
        TMVStaticAssert((Sizes<M1u::_rowsize,M1v::_rowsize>::same));
        TMVStaticAssert((Sizes<M1u::_colsize,M2::_rowsize>::same));
        TMVStaticAssert((Sizes<M1u::_rowsize,M2::_colsize>::same));
        TMVAssert(U.rowsize() == S.size());
        TMVAssert(U.rowsize() == V.colsize());
        TMVAssert(U.colsize() == minv.rowsize());
        TMVAssert(V.rowsize() == minv.colsize());

        const int cs = Sizes<M2::_rowsize,M1u::_colsize>::size;
        const int rs = Sizes<M2::_colsize,M1u::_rowsize>::size;
        typedef typename M1u::const_cview_type M1uv;
        typedef typename M1s::const_cview_type M1sv;
        typedef typename M1v::const_cview_type M1vv;
        typedef typename M2::cview_type M2v;
        TMV_MAYBE_CREF(M1u,M1uv) Uv = U.cView();
        TMV_MAYBE_CREF(M1s,M1sv) Sv = S.cView();
        TMV_MAYBE_CREF(M1v,M1vv) Vv = V.cView();
        TMV_MAYBE_REF(M2,M2v) minvv = minv.cView();
        SV_Inverse_Helper<-2,cs,rs,M1uv,M1sv,M1vv,M2v>::call(Uv,Sv,Vv,minvv);
    }

    template <int algo, int cs, int rs, class M1u, class M1s, class M1v, class M2>
    struct SV_InverseATA_Helper;

    // algo 0: Trivial, nothing to do (M == 0 or N == 0)
    template <int cs, int rs, class M1u, class M1s, class M1v, class M2>
    struct SV_InverseATA_Helper<0,cs,rs,M1u,M1s,M1v,M2>
    {
        static TMV_INLINE void call(
            const M1u& , const M1s& , const M1v& , int, M2& ) {} 
    };

    // algo 11: Normal case
    template <int cs, int rs, class M1u, class M1s, class M1v, class M2>
    struct SV_InverseATA_Helper<11,cs,rs,M1u,M1s,M1v,M2>
    {
        static void call(
            const M1u& , const M1s& S, const M1v& V, M2& ata)
        {
#ifdef PRINTALGO_SVD
            std::cout<<"SVInverseATA algo 11: cs,rs = "<<cs<<','<<rs<<std::endl;
#endif
            // A = USV
            // AtA = Vt S^2 V
            //     = Vt S^-2 V
            typedef typename M1u::value_type T;
            const int xx = TMV_UNKNOWN;
            typename MCopyHelper<T,Rec,xx,rs>::type SinvV = V / S;
            ata = SinvV.adjoint() * SinvV;
        }
    };

    // algo 90: call InstSV_InverseATA
    template <int cs, int rs, class M1u, class M1s, class M1v, class M2>
    struct SV_InverseATA_Helper<90,cs,rs,M1u,M1s,M1v,M2>
    {
        static TMV_INLINE void call(
            const M1u& U, const M1s& S, const M1v& V, M2& ata)
        { InstSV_InverseATA(U.xView(),S.xView(),V.xView(),ata.xView()); }
    };

    // algo 97: Conjugate
    template <int cs, int rs, class M1u, class M1s, class M1v, class M2>
    struct SV_InverseATA_Helper<97,cs,rs,M1u,M1s,M1v,M2>
    {
        static TMV_INLINE void call(
            const M1u& U, const M1s& S, const M1v& V, M2& ata)
        {
            typedef typename M1u::const_conjugate_type M1uc;
            typedef typename M1v::const_conjugate_type M1vc;
            typedef typename M2::conjugate_type M2c;
            M1uc Uc = U.conjugate();
            M1vc Vc = V.conjugate();
            M2c atac = ata.conjugate();
            SV_InverseATA_Helper<-2,cs,rs,M1uc,M1s,M1vc,M2c>::call(
                Uc,S,Vc,atac);
        }
    };

    // algo -3: Determine which algorithm to use
    template <int cs, int rs, class M1u, class M1s, class M1v, class M2>
    struct SV_InverseATA_Helper<-3,cs,rs,M1u,M1s,M1v,M2>
    {
        static TMV_INLINE void call(
            const M1u& U, const M1s& S, const M1v& V, M2& ata)
        {
            const int algo = 
                cs == 0 || rs == 0 ? 0 : 
                11;
#ifdef PRINTALGO_SVD
            std::cout<<"Inline SVInverseATA\n";
            std::cout<<"cs,rs = "<<cs<<','<<rs<<std::endl;
            std::cout<<"U = "<<TMV_Text(U)<<std::endl;
            std::cout<<"S = "<<TMV_Text(S)<<std::endl;
            std::cout<<"V = "<<TMV_Text(V)<<std::endl;
            std::cout<<"ata = "<<TMV_Text(ata)<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
#endif
            SV_InverseATA_Helper<algo,cs,rs,M1u,M1s,M1v,M2>::call(
                U,S,V,ata);
        }
    };

    // algo -2: Check for inst
    template <int cs, int rs, class M1u, class M1s, class M1v, class M2>
    struct SV_InverseATA_Helper<-2,cs,rs,M1u,M1s,M1v,M2>
    {
        static TMV_INLINE void call(
            const M1u& U, const M1s& S, const M1v& V, M2& ata)
        {
            typedef typename M1u::value_type T1;
            typedef typename M2::value_type T2;
            const bool inst = 
                (cs == TMV_UNKNOWN || cs > 16) &&
                (rs == TMV_UNKNOWN || rs > 16) &&
#ifdef TMV_INST_MIX
                Traits2<T1,T2>::samebase &&
#else
                Traits2<T1,T2>::sametype &&
#endif
                Traits<T2>::isinst;
            const int algo = 
                cs == 0 || rs == 0 ? 0 : 
                M2::_conj ? 97 :
                inst ? 90 :
                -3;
            SV_InverseATA_Helper<algo,cs,rs,M1u,M1s,M1v,M2>::call(
                U,S,V,ata);
        }
    };

    template <int cs, int rs, class M1u, class M1s, class M1v, class M2>
    struct SV_InverseATA_Helper<-1,cs,rs,M1u,M1s,M1v,M2>
    {
        static TMV_INLINE void call(
            const M1u& U, const M1s& S, const M1v& V, M2& ata)
        { 
            SV_InverseATA_Helper<-2,cs,rs,M1u,M1s,M1v,M2>::call(
                U,S,V,ata); 
        }
    }; 

    template <class M1u, class M1s, class M1v, class M2>
    inline void InlineSV_InverseATA(
        const BaseMatrix_Rec<M1u>& U, const BaseMatrix_Diag<M1s>& S,
        const BaseMatrix_Rec<M1v>& V, BaseMatrix_Rec_Mutable<M2>& ata)
    {
        TMVStaticAssert((Sizes<M1u::_rowsize,M1s::_size>::same));
        TMVStaticAssert((Sizes<M1u::_rowsize,M1v::_colsize>::same));
        TMVStaticAssert((Sizes<M1u::_rowsize,M1v::_rowsize>::same));
        TMVStaticAssert((Sizes<M1u::_rowsize,M2::_rowsize>::same));
        TMVStaticAssert((Sizes<M1u::_rowsize,M2::_colsize>::same));
        TMVAssert(U.rowsize() == S.size());
        TMVAssert(U.rowsize() == V.colsize());
        TMVAssert(V.rowsize() == ata.rowsize());
        TMVAssert(V.rowsize() == ata.colsize());

        const int cs = M1u::_colsize;
        const int rs = Sizes<M2::_colsize,
              Sizes<M2::_rowsize,M1u::_rowsize>::size>::size;
        typedef typename M1u::const_cview_type M1uv;
        typedef typename M1s::const_cview_type M1sv;
        typedef typename M1v::const_cview_type M1vv;
        typedef typename M2::cview_type M2v;
        TMV_MAYBE_CREF(M1u,M1uv) Uv = U.cView();
        TMV_MAYBE_CREF(M1s,M1sv) Sv = S.cView();
        TMV_MAYBE_CREF(M1v,M1vv) Vv = V.cView();
        TMV_MAYBE_REF(M2,M2v) atav = ata.cView();
        SV_InverseATA_Helper<-3,cs,rs,M1uv,M1sv,M1vv,M2v>::call(
            Uv,Sv,Vv,atav);
    }

    template <class M1u, class M1s, class M1v, class M2>
    inline void SV_InverseATA(
        const BaseMatrix_Rec<M1u>& U, const BaseMatrix_Diag<M1s>& S,
        const BaseMatrix_Rec<M1v>& V, BaseMatrix_Rec_Mutable<M2>& ata)
    {
        TMVStaticAssert((Sizes<M1u::_rowsize,M1s::_size>::same));
        TMVStaticAssert((Sizes<M1u::_rowsize,M1v::_colsize>::same));
        TMVStaticAssert((Sizes<M1v::_rowsize,M2::_rowsize>::same));
        TMVStaticAssert((Sizes<M1v::_rowsize,M2::_colsize>::same));
        TMVAssert(U.rowsize() == S.size());
        TMVAssert(U.rowsize() == V.colsize());
        TMVAssert(V.rowsize() == ata.rowsize());
        TMVAssert(V.rowsize() == ata.colsize());

        const int cs = M1u::_colsize;
        const int rs = Sizes<M2::_colsize,
              Sizes<M2::_rowsize,M1u::_rowsize>::size>::size;
        typedef typename M1u::const_cview_type M1uv;
        typedef typename M1s::const_cview_type M1sv;
        typedef typename M1v::const_cview_type M1vv;
        typedef typename M2::cview_type M2v;
        TMV_MAYBE_CREF(M1u,M1uv) Uv = U.cView();
        TMV_MAYBE_CREF(M1s,M1sv) Sv = S.cView();
        TMV_MAYBE_CREF(M1v,M1vv) Vv = V.cView();
        TMV_MAYBE_REF(M2,M2v) atav = ata.cView();
        SV_InverseATA_Helper<-2,cs,rs,M1uv,M1sv,M1vv,M2v>::call(
            Uv,Sv,Vv,atav);
    }

} // namespace tmv

#endif

