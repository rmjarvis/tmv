

#ifndef TMV_LUDiv_H
#define TMV_LUDiv_H

#include "TMV_BaseMatrix_Rec.h"
#include "TMV_BaseMatrix_Tri.h"
#include "TMV_Permutation.h"
#include "TMV_DivVM_Funcs.h"
#include "TMV_DivMM_Funcs.h"
#include "TMV_MultMM_Funcs.h"

#ifdef PRINTALGO_LU
#include <iostream>
#include "TMV_ProdMV.h"
#include "TMV_ProdMM.h"
#endif

namespace tmv {

    // Defined in TMV_LUDiv.cpp
    template <class T1, int C1, class T2>
    void InstLU_SolveInPlace(
        const ConstMatrixView<T1,C1>& m1, const Permutation& P,
        MatrixView<T2> m2);
    template <class T1, int C1, class T2>
    void InstLU_SolveTransposeInPlace(
        const ConstMatrixView<T1,C1>& m1, const Permutation& P,
        MatrixView<T2> m2);
    template <class T1, int C1, class T2>
    void InstLU_SolveInPlace(
        const ConstMatrixView<T1,C1>& m1, const Permutation& P,
        VectorView<T2> v2);
    template <class T1, int C1, class T2>
    void InstLU_SolveTransposeInPlace(
        const ConstMatrixView<T1,C1>& m1, const Permutation& P,
        VectorView<T2> v2);

    template <int algo, bool trans, int cs, int rs, class M1, class M2>
    struct LU_Solve_Helper;

    // algo 0: Trivial, nothing to do (M == 0 or N == 0)
    // Also used for invalid real/complex combination from the virtual calls.
    template <bool trans, int cs, int rs, class M1, class M2>
    struct LU_Solve_Helper<0,trans,cs,rs,M1,M2>
    { static TMV_INLINE void call(const M1& , const Permutation& , M2& ) {} };

    // algo 11: Normal case
    template <int cs, int rs, class M1, class M2>
    struct LU_Solve_Helper<11,false,cs,rs,M1,M2>
    {
        static void call(const M1& m1, const Permutation& P, M2& m2)
        {
#ifdef PRINTALGO_LU
            std::cout<<"LUSolve algo 11: trans,cs,rs = "<<
                false<<','<<cs<<','<<rs<<std::endl;
            //std::cout<<"m1 = "<<m1<<std::endl;
            //std::cout<<"m2 = "<<m2<<std::endl;
#endif
            // m2 = (PLU)^-1 m2
            //    = U^-1 L^-1 Pt m2
            P.inverse().applyOnLeft(m2);
#ifdef PRINTALGO_LU
            //std::cout<<"after permute m2 => "<<m2<<std::endl;
#endif
            NoAliasTriLDivEq(m2,m1.unitLowerTri());
#ifdef PRINTALGO_LU
            //std::cout<<"after L^-1 m2 => "<<m2<<std::endl;
            //std::cout<<"L*m2 = "<<m1.unitLowerTri()*m2<<std::endl;
#endif
            NoAliasTriLDivEq(m2,m1.upperTri());
#ifdef PRINTALGO_LU
            //std::cout<<"after U^-1 m2 => "<<m2<<std::endl;
            //std::cout<<"U*m2 = "<<m1.upperTri()*m2<<std::endl;
#endif
        }
    };
    template <int cs, int rs, class M1, class M2>
    struct LU_Solve_Helper<11,true,cs,rs,M1,M2>
    {
        static void call(const M1& m1, const Permutation& P, M2& m2)
        {
#ifdef PRINTALGO_LU
            std::cout<<"LUSolve algo 11: trans,cs,rs = "<<
                true<<','<<cs<<','<<rs<<std::endl;
            //std::cout<<"m1 = "<<m1<<std::endl;
            //std::cout<<"m2 = "<<m2<<std::endl;
#endif
            // m2 = (PLU)^-1t m2
            //    = (U^-1 L^-1 Pt)t m2
            //    = P L^-1t U^-1t m2
            NoAliasTriLDivEq(m2,m1.upperTri().transpose());
#ifdef PRINTALGO_LU
            //std::cout<<"after UT^-1 m2 => "<<m2<<std::endl;
            //std::cout<<"Ut*m2 = "<<m1.upperTri().transpose()*m2<<std::endl;
#endif
            NoAliasTriLDivEq(m2,m1.unitLowerTri().transpose());
#ifdef PRINTALGO_LU
            //std::cout<<"after LT^-1 m2 => "<<m2<<std::endl;
            //std::cout<<"LT*m2 = "<<m1.unitLowerTri().transpose()*m2<<std::endl;
#endif
            P.applyOnLeft(m2);
#ifdef PRINTALGO_LU
            //std::cout<<"after reverse permute m2 => "<<m2<<std::endl;
#endif
        }
    };

    // algo 90: call InstLU_Solve
    template <int cs, int rs, class M1, class M2>
    struct LU_Solve_Helper<90,false,cs,rs,M1,M2>
    {
        static TMV_INLINE void call(const M1& m1, const Permutation& P, M2& m2)
        { InstLU_SolveInPlace(m1.xView(),P,m2.xView()); }
    };
    template <int cs, int rs, class M1, class M2>
    struct LU_Solve_Helper<90,true,cs,rs,M1,M2>
    {
        static TMV_INLINE void call(const M1& m1, const Permutation& P, M2& m2)
        { InstLU_SolveTransposeInPlace(m1.xView(),P,m2.xView()); }
    };

    // algo 95: Turn m1,m3 into vector
    template <bool trans, int cs, int rs, class M1, class M2>
    struct LU_Solve_Helper<95,trans,cs,rs,M1,M2>
    {
        static TMV_INLINE void call(const M1& m1, const Permutation& P, M2& m2)
        {
            TMVStaticAssert(!ShapeTraits<M2::_shape>::vector);
            typedef typename M2::col_type M2c;
            M2c m2c = m2.col(0);
            LU_Solve_Helper<-2,trans,cs,rs,M1,M2c>::call(m1,P,m2c);
        }
    };

    // algo 97: Conjugate
    template <bool trans, int cs, int rs, class M1, class M2>
    struct LU_Solve_Helper<97,trans,cs,rs,M1,M2>
    {
        static TMV_INLINE void call(const M1& m1, const Permutation& P, M2& m2)
        {
            typedef typename M1::const_conjugate_type M1c;
            typedef typename M2::conjugate_type M2c;
            M1c m1c = m1.conjugate();
            M2c m2c = m2.conjugate();
            LU_Solve_Helper<-2,trans,cs,rs,M1c,M2c>::call(m1c,P,m2c);
        }
    };

    // algo -3: Determine which algorithm to use
    template <bool trans, int cs, int rs, class M1, class M2>
    struct LU_Solve_Helper<-3,trans,cs,rs,M1,M2>
    {
        static TMV_INLINE void call(const M1& m1, const Permutation& P, M2& m2)
        {
            const bool invalid =
                M1::iscomplex && M2::isreal;
            const int algo = 
                cs == 0 || rs == 0 || invalid ? 0 : 
                11;
#ifdef PRINTALGO_LU
            std::cout<<"Inline LUSolve\n";
            std::cout<<"trans,cs,rs = "<<trans<<','<<cs<<','<<rs<<std::endl;
            std::cout<<"m1 = "<<TMV_Text(m1)<<std::endl;
            std::cout<<"m2 = "<<TMV_Text(m2)<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
            //std::cout<<"m1 = "<<m1<<std::endl;
            //std::cout<<"m2 = "<<m2<<std::endl;
#endif
            LU_Solve_Helper<algo,trans,cs,rs,M1,M2>::call(m1,P,m2);
#ifdef PRINTALGO_LU
            //std::cout<<"m2 => "<<m2<<std::endl;
#endif
        }
    };

    // algo -2: Check for inst
    template <bool trans, int cs, int rs, class M1, class M2>
    struct LU_Solve_Helper<-2,trans,cs,rs,M1,M2>
    {
        static TMV_INLINE void call(const M1& m1, const Permutation& P, M2& m2)
        {
            typedef typename M1::value_type T1;
            typedef typename M2::value_type T2;
            const bool inst = 
                (cs == TMV_UNKNOWN || cs > 16) &&
                (rs == TMV_UNKNOWN || rs > 16 || rs == 1) &&
#ifdef TMV_INST_MIX
                Traits2<T1,T2>::samebase &&
#else
                Traits2<T1,T2>::sametype &&
#endif
                M1::_colmajor && 
                Traits<T1>::isinst;
            const bool makevector =
                rs == 1 && !ShapeTraits<M2::_shape>::vector;
            const bool invalid =
                M1::iscomplex && M2::isreal;
            const int algo = 
                cs == 0 || rs == 0 || invalid ? 0 : 
                makevector ? 95 :
                M2::_conj ? 97 :
                inst ? 90 :
                -3;
            LU_Solve_Helper<algo,trans,cs,rs,M1,M2>::call(m1,P,m2);
        }
    };

    template <bool trans, int cs, int rs, class M1, class M2>
    struct LU_Solve_Helper<-1,trans,cs,rs,M1,M2>
    {
        static TMV_INLINE void call(const M1& m1, const Permutation& P, M2& m2)
        { LU_Solve_Helper<-2,trans,cs,rs,M1,M2>::call(m1,P,m2); }
    };

    template <class M1, class M2>
    static inline void InlineLU_SolveInPlace(
        const BaseMatrix_Rec<M1>& m1, const Permutation& P,
        BaseMatrix_Rec_Mutable<M2>& m2)
    {
        const int cs = M2::_colsize;
        const int rs = M2::_rowsize;
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::cview_type M2v;
        TMV_MAYBE_CREF(M1,M1v) m1v = m1.cView();
        TMV_MAYBE_REF(M2,M2v) m2v = m2.cView();
        LU_Solve_Helper<-3,false,cs,rs,M1v,M2v>::call(m1v,P,m2v);
    }

    template <class M1, class M2>
    static inline void LU_SolveInPlace(
        const BaseMatrix_Rec<M1>& m1, const Permutation& P,
        BaseMatrix_Rec_Mutable<M2>& m2)
    {
        const int cs = M2::_colsize;
        const int rs = M2::_rowsize;
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::cview_type M2v;
        TMV_MAYBE_CREF(M1,M1v) m1v = m1.cView();
        TMV_MAYBE_REF(M2,M2v) m2v = m2.cView();
        LU_Solve_Helper<-2,false,cs,rs,M1v,M2v>::call(m1v,P,m2v);
    }

    template <class M1, class V2>
    static inline void InlineLU_SolveInPlace(
        const BaseMatrix_Rec<M1>& m1, const Permutation& P,
        BaseVector_Mutable<V2>& v2)
    {
        const int cs = V2::_size;
        typedef typename M1::const_cview_type M1v;
        typedef typename V2::cview_type V2v;
        TMV_MAYBE_CREF(M1,M1v) m1v = m1.cView();
        TMV_MAYBE_REF(V2,V2v) v2v = v2.cView();
        LU_Solve_Helper<-3,false,cs,1,M1v,V2v>::call(m1v,P,v2v);
    }

    template <class M1, class V2>
    static inline void LU_SolveInPlace(
        const BaseMatrix_Rec<M1>& m1, const Permutation& P,
        BaseVector_Mutable<V2>& v2)
    {
        const int cs = V2::_size;
        typedef typename M1::const_cview_type M1v;
        typedef typename V2::cview_type V2v;
        TMV_MAYBE_CREF(M1,M1v) m1v = m1.cView();
        TMV_MAYBE_REF(V2,V2v) v2v = v2.cView();
        LU_Solve_Helper<-1,false,cs,1,M1v,V2v>::call(m1v,P,v2v);
    }

    template <class M1, class M2>
    static inline void InlineLU_SolveTransposeInPlace(
        const BaseMatrix_Rec<M1>& m1, const Permutation& P,
        BaseMatrix_Rec_Mutable<M2>& m2)
    {
        const int cs = M2::_colsize;
        const int rs = M2::_rowsize;
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::cview_type M2v;
        TMV_MAYBE_CREF(M1,M1v) m1v = m1.cView();
        TMV_MAYBE_REF(M2,M2v) m2v = m2.cView();
        LU_Solve_Helper<-3,true,cs,rs,M1v,M2v>::call(m1v,P,m2v);
    }

    template <class M1, class M2>
    static inline void LU_SolveTransposeInPlace(
        const BaseMatrix_Rec<M1>& m1, const Permutation& P,
        BaseMatrix_Rec_Mutable<M2>& m2)
    {
        const int cs = M2::_colsize;
        const int rs = M2::_rowsize;
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::cview_type M2v;
        TMV_MAYBE_CREF(M1,M1v) m1v = m1.cView();
        TMV_MAYBE_REF(M2,M2v) m2v = m2.cView();
        LU_Solve_Helper<-1,true,cs,rs,M1v,M2v>::call(m1v,P,m2v);
    }

    template <class M1, class V2>
    static inline void InlineLU_SolveTransposeInPlace(
        const BaseMatrix_Rec<M1>& m1, const Permutation& P,
        BaseVector_Mutable<V2>& v2)
    {
        const int cs = V2::_size;
        typedef typename M1::const_cview_type M1v;
        typedef typename V2::cview_type V2v;
        TMV_MAYBE_CREF(M1,M1v) m1v = m1.cView();
        TMV_MAYBE_REF(V2,V2v) v2v = v2.cView();
        LU_Solve_Helper<-3,true,cs,1,M1v,V2v>::call(m1v,P,v2v);
    }

    template <class M1, class V2>
    static inline void LU_SolveTransposeInPlace(
        const BaseMatrix_Rec<M1>& m1, const Permutation& P,
        BaseVector_Mutable<V2>& v2)
    {
        const int cs = V2::_size;
        typedef typename M1::const_cview_type M1v;
        typedef typename V2::cview_type V2v;
        TMV_MAYBE_CREF(M1,M1v) m1v = m1.cView();
        TMV_MAYBE_REF(V2,V2v) v2v = v2.cView();
        LU_Solve_Helper<-1,true,cs,1,M1v,V2v>::call(m1v,P,v2v);
    }

} // namespace tmv

#endif

