

#ifndef TMV_SVDiv_H
#define TMV_SVDiv_H

#include "TMV_BaseMatrix_Rec.h"
#include "TMV_BaseMatrix_Tri.h"
#include "TMV_DivVM_Funcs.h"
#include "TMV_DivMM_Funcs.h"
#include "TMV_MultMM_Funcs.h"
#include "TMV_Permutation.h"
#include "TMV_PackedQ.h"

#ifdef PRINTALGO_SVD
#include <iostream>
#include "TMV_ProdMV.h"
#include "TMV_ProdMM.h"
#include "TMV_MatrixIO.h"
#include "TMV_VectorIO.h"
#include "TMV_AddMM.h"
#include "TMV_AddVV.h"
#include "TMV_MultMM.h"
#include "TMV_MultMV.h"
#include "TMV_MultUV.h"
#endif

namespace tmv {

    // Defined in TMV_SVDiv.cpp
    template <class T1, int C1, class RT1, class T2, int C2, class T3>
    void InstSV_Solve(
        const ConstMatrixView<T1,C1>& U, const ConstDiagMatrixView<RT1>& S,
        const ConstMatrixView<T1,C1>& V, 
        const ConstMatrixView<T2,C2>& m2, MatrixView<T3> m3);
    template <class T1, int C1, class RT1, class T2, int C2, class T3>
    void InstSV_Solve(
        const ConstMatrixView<T1,C1>& U, const ConstDiagMatrixView<RT1>& S,
        const ConstMatrixView<T1,C1>& V, 
        const ConstVectorView<T2,C2>& v2, VectorView<T3> v3);

    template <int algo, int cs, int rs, class M1u, class M1s, class M1v, class M2, class M3>
    struct SV_Solve_Helper;

    // algo 0: Trivial, nothing to do (M == 0 or N == 0)
    // Also used for invalid real/complex combination from the virtual calls.
    template <int cs, int rs, class M1u, class M1s, class M1v, class M2, class M3>
    struct SV_Solve_Helper<0,cs,rs,M1u,M1s,M1v,M2,M3>
    {
        static TMV_INLINE void call(
            const M1u& , const M1s& , const M1v& , int , const M2& , M3& ) {}
    };

    // algo 11: Normal case
    template <int cs, int rs, class M1u, class M1s, class M1v, class M2, class M3>
    struct SV_Solve_Helper<11,cs,rs,M1u,M1s,M1v,M2,M3>
    {
        static void call(
            const M1u& U, const M1s& S, const M1v& V, const M2& m2, M3& m3)
        {
            typedef typename M1u::value_type T1;
            typedef typename M2::value_type T2;
            typedef typename Traits2<T1,T2>::type T12;
            const int xx = TMV_UNKNOWN;
            const int xs = Sizes<M2::_rowsize,M3::_rowsize>::size;
            typedef typename MCopyHelper<T12,Rec,xx,xs>::type M2c;

#ifdef PRINTALGO_SVD
            std::cout<<"SVSolve algo 11: cs,rs = "<<cs<<','<<rs<<std::endl;
            std::cout<<"U = "<<U<<std::endl;
            std::cout<<"S = "<<S<<std::endl;
            std::cout<<"V = "<<V<<std::endl;
            std::cout<<"m2 = "<<m2<<std::endl;
            std::cout<<"m3 = "<<m3<<std::endl;
#endif

            // m3 = (USV)^-1 m2
            //    = Vt S^-1 Ut m2
            M2c temp = U.adjoint() * m2;
            temp /= S;
            m3 = V.adjoint()*temp;
        }
    };

    // algo 12: Normal case - M2, M3 are vectors
    template <int cs, int rs, class M1u, class M1s, class M1v, class M2, class M3>
    struct SV_Solve_Helper<12,cs,rs,M1u,M1s,M1v,M2,M3>
    {
        static void call(
            const M1u& U, const M1s& S, const M1v& V, const M2& v2, M3& v3)
        {
            typedef typename M1u::value_type T1;
            typedef typename M2::value_type T2;
            typedef typename Traits2<T1,T2>::type T12;
            const int xx = TMV_UNKNOWN;
            typedef typename VCopyHelper<T12,xx>::type V2c;

#ifdef PRINTALGO_SVD
            std::cout<<"SVSolve algo 12: cs,rs = "<<cs<<','<<rs<<std::endl;
            std::cout<<"U = "<<U<<std::endl;
            std::cout<<"S = "<<S<<std::endl;
            std::cout<<"V = "<<V<<std::endl;
            std::cout<<"v2 = "<<v2<<std::endl;
            std::cout<<"v3 = "<<v3<<std::endl;
#endif
            // v3 = (USV)^-1 v2
            //    = Vt S^-1 Ut v2
            V2c temp = U.adjoint() * v2;
            temp /= S;
            v3 = V.adjoint()*temp;
        }
    };

    // algo 90: call InstSV_Solve
    template <int cs, int rs, class M1u, class M1s, class M1v, class M2, class M3>
    struct SV_Solve_Helper<90,cs,rs,M1u,M1s,M1v,M2,M3>
    {
        static TMV_INLINE void call(
            const M1u& U, const M1s& S, const M1v& V, const M2& m2, M3& m3)
        {
            InstSV_Solve(
                U.xView(),S.xView(),V.xView(),m2.xView(),m3.xView()); 
        }
    };

    // algo 97: Conjugate
    template <int cs, int rs, class M1u, class M1s, class M1v, class M2, class M3>
    struct SV_Solve_Helper<97,cs,rs,M1u,M1s,M1v,M2,M3>
    {
        static TMV_INLINE void call(
            const M1u& U, const M1s& S, const M1v& V, const M2& m2, M3& m3)
        {
            typedef typename M1u::const_conjugate_type M1uc;
            typedef typename M1v::const_conjugate_type M1vc;
            typedef typename M2::const_conjugate_type M2c;
            typedef typename M3::conjugate_type M3c;
            M1uc Uc = U.conjugate();
            M1vc Vc = V.conjugate();
            M2c m2c = m2.conjugate();
            M3c m3c = m3.conjugate();
            SV_Solve_Helper<-2,cs,rs,M1uc,M1s,M1vc,M2c,M3c>::call(
                Uc,S,Vc,m2c,m3c);
        }
    };

    // algo -3: Determine which algorithm to use
    template <int cs, int rs, class M1u, class M1s, class M1v, class M2, class M3>
    struct SV_Solve_Helper<-3,cs,rs,M1u,M1s,M1v,M2,M3>
    {
        static TMV_INLINE void call(
            const M1u& U, const M1s& S, const M1v& V, const M2& m2, M3& m3)
        {
            TMVStaticAssert((
                    ShapeTraits<M2::_shape>::vector == 
                    ShapeTraits<M3::_shape>::vector));

            const bool invalid =
                (M1u::iscomplex || M2::iscomplex) && M3::isreal;
            const int algo = 
                cs == 0 || rs == 0 || invalid ? 0 : 
                ShapeTraits<M3::_shape>::vector ? 12 : 11;
#ifdef PRINTALGO_SVD
            std::cout<<"Inline SVSolve\n";
            std::cout<<"cs,rs = "<<cs<<','<<rs<<std::endl;
            std::cout<<"U = "<<TMV_Text(U)<<std::endl;
            std::cout<<"S = "<<TMV_Text(S)<<std::endl;
            std::cout<<"V = "<<TMV_Text(V)<<std::endl;
            std::cout<<"m2 = "<<TMV_Text(m2)<<std::endl;
            std::cout<<"m3 = "<<TMV_Text(m3)<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
#endif
            SV_Solve_Helper<algo,cs,rs,M1u,M1s,M1v,M2,M3>::call(
                U,S,V,m2,m3);
        }
    };

    // algo -2: Check for inst
    template <int cs, int rs, class M1u, class M1s, class M1v, class M2, class M3>
    struct SV_Solve_Helper<-2,cs,rs,M1u,M1s,M1v,M2,M3>
    {
        static TMV_INLINE void call(
            const M1u& U, const M1s& S, const M1v& V, const M2& m2, M3& m3)
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
                Traits<T1>::isinst;
            const bool invalid =
                (M1u::iscomplex || M2::iscomplex) && M3::isreal;
            const int algo = 
                cs == 0 || rs == 0 || invalid ? 0 : 
                M3::_conj ? 97 :
                inst ? 90 :
                -3;
            SV_Solve_Helper<algo,cs,rs,M1u,M1s,M1v,M2,M3>::call(
                U,S,V,m2,m3);
        }
    };

    template <int cs, int rs, class M1u, class M1s, class M1v, class M2, class M3>
    struct SV_Solve_Helper<-1,cs,rs,M1u,M1s,M1v,M2,M3>
    {
        static TMV_INLINE void call(
            const M1u& U, const M1s& S, const M1v& V, const M2& m2, M3& m3)
        {
            SV_Solve_Helper<-2,cs,rs,M1u,M1s,M1v,M2,M3>::call(
                U,S,V,m2,m3); 
        }
    };

    template <class M1u, class M1s, class M1v, class M2, class M3>
    static inline void InlineSV_Solve(
        const BaseMatrix_Rec<M1u>& U, const BaseMatrix_Diag<M1s>& S,
        const BaseMatrix_Rec<M1v>& V, 
        const BaseMatrix_Rec<M2>& m2, BaseMatrix_Rec_Mutable<M3>& m3)
    {
        TMVStaticAssert((Sizes<M1u::_rowsize,M1s::_size>::same));
        TMVStaticAssert((Sizes<M1u::_rowsize,M1v::_colsize>::same));
        TMVStaticAssert((Sizes<M1u::_rowsize,M1v::_rowsize>::same));
        TMVStaticAssert((Sizes<M1u::_colsize,M2::_colsize>::same));
        TMVStaticAssert((Sizes<M1u::_rowsize,M3::_colsize>::same));
        TMVStaticAssert((Sizes<M2::_rowsize,M3::_rowsize>::same));
        TMVAssert(U.rowsize() == S.size());
        TMVAssert(U.rowsize() == V.colsize());
        TMVAssert(U.colsize() == m2.colsize());
        TMVAssert(V.rowsize() == m3.colsize());
        TMVAssert(m2.rowsize() == m3.rowsize());

        const int cs = Sizes<M3::_colsize,M1v::_rowsize>::size;
        const int rs = Sizes<M2::_colsize,M1u::_colsize>::size;
        typedef typename M1u::const_cview_type M1uv;
        typedef typename M1s::const_cview_type M1sv;
        typedef typename M1v::const_cview_type M1vv;
        typedef typename M2::const_cview_type M2v;
        typedef typename M3::cview_type M3v;
        TMV_MAYBE_CREF(M1u,M1uv) Uv = U.cView();
        TMV_MAYBE_CREF(M1s,M1sv) Sv = S.cView();
        TMV_MAYBE_CREF(M1v,M1vv) Vv = V.cView();
        TMV_MAYBE_CREF(M2,M2v) m2v = m2.cView();
        TMV_MAYBE_REF(M3,M3v) m3v = m3.cView();
        SV_Solve_Helper<-3,cs,rs,M1uv,M1sv,M1vv,M2v,M3v>::call(
            Uv,Sv,Vv,m2v,m3v);
    }

    template <class M1u, class M1s, class M1v, class M2, class M3>
    static inline void SV_Solve(
        const BaseMatrix_Rec<M1u>& U, const BaseMatrix_Diag<M1s>& S,
        const BaseMatrix_Rec<M1v>& V, 
        const BaseMatrix_Rec<M2>& m2, BaseMatrix_Rec_Mutable<M3>& m3)
    {
        TMVStaticAssert((Sizes<M1u::_rowsize,M1s::_size>::same));
        TMVStaticAssert((Sizes<M1u::_rowsize,M1v::_colsize>::same));
        TMVStaticAssert((Sizes<M1u::_rowsize,M1v::_rowsize>::same));
        TMVStaticAssert((Sizes<M1u::_colsize,M2::_colsize>::same));
        TMVStaticAssert((Sizes<M1u::_rowsize,M3::_colsize>::same));
        TMVStaticAssert((Sizes<M2::_rowsize,M3::_rowsize>::same));
        TMVAssert(U.rowsize() == S.size());
        TMVAssert(U.rowsize() == V.colsize());
        TMVAssert(U.colsize() == m2.colsize());
        TMVAssert(V.rowsize() == m3.colsize());
        TMVAssert(m2.rowsize() == m3.rowsize());

        const int cs = Sizes<M3::_colsize,M1u::_rowsize>::size;
        const int rs = Sizes<M2::_colsize,M1u::_colsize>::size;
        typedef typename M1u::const_cview_type M1uv;
        typedef typename M1s::const_cview_type M1sv;
        typedef typename M1v::const_cview_type M1vv;
        typedef typename M2::const_cview_type M2v;
        typedef typename M3::cview_type M3v;
        TMV_MAYBE_CREF(M1u,M1uv) Uv = U.cView();
        TMV_MAYBE_CREF(M1s,M1sv) Sv = S.cView();
        TMV_MAYBE_CREF(M1v,M1vv) Vv = V.cView();
        TMV_MAYBE_CREF(M2,M2v) m2v = m2.cView();
        TMV_MAYBE_REF(M3,M3v) m3v = m3.cView();
        SV_Solve_Helper<-2,cs,rs,M1uv,M1sv,M1vv,M2v,M3v>::call(
            Uv,Sv,Vv,m2v,m3v);
    }

    template <class M1u, class M1s, class M1v, class V2, class V3>
    static inline void InlineSV_Solve(
        const BaseMatrix_Rec<M1u>& U, const BaseMatrix_Diag<M1s>& S,
        const BaseMatrix_Rec<M1v>& V, 
        const BaseVector_Calc<V2>& v2, BaseVector_Mutable<V3>& v3)
    {
        TMVStaticAssert((Sizes<M1u::_rowsize,M1s::_size>::same));
        TMVStaticAssert((Sizes<M1u::_rowsize,M1v::_colsize>::same));
        TMVStaticAssert((Sizes<M1u::_rowsize,M1v::_rowsize>::same));
        TMVStaticAssert((Sizes<M1u::_colsize,V2::_size>::same));
        TMVStaticAssert((Sizes<M1u::_rowsize,V3::_size>::same));
        TMVAssert(U.rowsize() == S.size());
        TMVAssert(U.rowsize() == V.colsize());
        TMVAssert(U.colsize() == v2.size());
        TMVAssert(V.rowsize() == v3.size());

        const int cs = Sizes<V3::_size,M1u::_rowsize>::size;
        const int rs = Sizes<V2::_size,M1u::_colsize>::size;
        typedef typename M1u::const_cview_type M1uv;
        typedef typename M1s::const_cview_type M1sv;
        typedef typename M1v::const_cview_type M1vv;
        typedef typename V2::const_cview_type V2v;
        typedef typename V3::cview_type V3v;
        TMV_MAYBE_CREF(M1u,M1uv) Uv = U.cView();
        TMV_MAYBE_CREF(M1s,M1sv) Sv = S.cView();
        TMV_MAYBE_CREF(M1v,M1vv) Vv = V.cView();
        TMV_MAYBE_CREF(V2,V2v) v2v = v2.cView();
        TMV_MAYBE_REF(V3,V3v) v3v = v3.cView();
        SV_Solve_Helper<-3,cs,rs,M1uv,M1sv,M1vv,V2v,V3v>::call(
            Uv,Sv,Vv,v2v,v3v);
    }

    template <class M1u, class M1s, class M1v, class V2, class V3>
    static inline void SV_Solve(
        const BaseMatrix_Rec<M1u>& U, const BaseMatrix_Diag<M1s>& S,
        const BaseMatrix_Rec<M1v>& V,
        const BaseVector_Calc<V2>& v2, BaseVector_Mutable<V3>& v3)
    {
        TMVStaticAssert((Sizes<M1u::_rowsize,M1s::_size>::same));
        TMVStaticAssert((Sizes<M1u::_rowsize,M1v::_colsize>::same));
        TMVStaticAssert((Sizes<M1u::_rowsize,M1v::_rowsize>::same));
        TMVStaticAssert((Sizes<M1u::_colsize,V2::_size>::same));
        TMVStaticAssert((Sizes<M1u::_rowsize,V3::_size>::same));
        TMVAssert(U.rowsize() == S.size());
        TMVAssert(U.rowsize() == V.colsize());
        TMVAssert(U.colsize() == v2.size());
        TMVAssert(V.rowsize() == v3.size());

        const int cs = Sizes<V3::_size,M1u::_rowsize>::size;
        const int rs = Sizes<V2::_size,M1u::_colsize>::size;
        typedef typename M1u::const_cview_type M1uv;
        typedef typename M1s::const_cview_type M1sv;
        typedef typename M1v::const_cview_type M1vv;
        typedef typename V2::const_cview_type V2v;
        typedef typename V3::cview_type V3v;
        TMV_MAYBE_CREF(M1u,M1uv) Uv = U.cView();
        TMV_MAYBE_CREF(M1s,M1sv) Sv = S.cView();
        TMV_MAYBE_CREF(M1v,M1vv) Vv = V.cView();
        TMV_MAYBE_CREF(V2,V2v) v2v = v2.cView();
        TMV_MAYBE_REF(V3,V3v) v3v = v3.cView();
        SV_Solve_Helper<-2,cs,rs,M1uv,M1sv,M1vv,V2v,V3v>::call(
            Uv,Sv,Vv,v2v,v3v);
    }

} // namespace tmv

#endif

