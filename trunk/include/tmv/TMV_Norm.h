
#ifndef TMV_Norm_H
#define TMV_Norm_H

#include "TMV_BaseVector.h"
#include "TMV_BaseMatrix_Rec.h"
#include "TMV_BaseMatrix_Tri.h"
#include "TMV_SVD.h"

namespace tmv {

    //
    // The algorithm selection for Norm2(v) and NormF(m) are 
    // basically the same, so we have them all here in one file
    // to use the same Helper structure.
    //

    // Defined in TMV_Vector.cpp
    template <class T>
    typename ConstVectorView<T>::float_type InstNorm2(
        const ConstVectorView<T>& v);

    // Defined in TMV_Matrix.cpp
    template <class T>
    typename ConstMatrixView<T>::float_type InstNormF(
        const ConstMatrixView<T>& m);

    // Defined in TMV_TriMatrix.cpp
    template <class T>
    typename ConstUpperTriMatrixView<T>::float_type InstNormF(
        const ConstUpperTriMatrixView<T>& m);
    

    // This helper struct works for either Vector or Matrix "V"
    template <int algo, class V>
    struct Norm_Helper;

    // algo 11: simple: sqrt(NormSq(v))
    template <class V>
    struct Norm_Helper<11,V>
    {
        typedef typename V::float_type RT;
        static inline RT call(const V& v)
        { return TMV_SQRT(v.normSq()); }
    };

    // algo 12: Robust algorithm with checks for overflow and underflow.
    // This one always calls MaxAbsElement and then NormSq.
    // This is inefficient if there are no problems.
    // Since no problems is the usual case, I switched to the below
    // version (algo 13) that calls NormSq first and then redoes it if there
    // are problems.
    template <class V>
    struct Norm_Helper<12,V>
    {
        typedef typename V::float_type RT;
        static RT call(const V& v)
        {
            const RT eps = TMV_Epsilon<RT>();

            // Start with the maximum |v(i)|.  It will tell us how (and if)
            // we need to use a scaling for NormSq().
            RT vmax = v.maxAbs2Element();

            if (vmax == RT(0)) {
                // If vmax = 0, then norm2 = 0:
                return RT(0);
            } else if (TMV_Underfloat(vmax * vmax)) {
                // If vmax^2 underflows but vmax != 0, then a naive NormSq()
                // will produce underflow rounding errors.  Find a better 
                // scaling.  eps is a pure power of 2, so no rounding errors 
                // from rescaling by a power of eps.
                const RT inveps = RT(1)/eps;
                RT scale = inveps;
                vmax *= scale;
                const RT eps2 = eps*eps;
                while (vmax < eps2) { scale *= inveps; vmax *= inveps; }
                return TMV_SQRT(v.normSq(scale))/scale;
            } else if (RT(1)/vmax == RT(0)) {
                // If 1/vmax == 0, then vmax is already inf, so no hope of
                // making it more accurate.  (And need to check, since otherwise
                // the next section would cause an infinite loop.)
                return vmax;
            } else if (TMV_Underflow(RT(1)/(v.nElements()*vmax*vmax))) {
                // If 1/(n*vmax^2) underflows, then a naive NormSq() will 
                // produce overflow./ Find a better scaling.
                const RT inveps = RT(1)/eps;
                RT scale = eps;
                vmax *= scale;
                while (vmax > inveps) { scale *= eps; vmax *= eps; }
                return TMV_SQRT(v.normSq(scale))/scale;
            } else {
                // No problems with overflow or underflow.
                return TMV_SQRT(v.normSq());
            }
        }
    };

    // algo 13: Robust algorithm with checks for overflow and underflow.
    // This version is slower if there is a problem, but since
    // there usually isn't a problem, it is generally faster.
    template <class V>
    struct Norm_Helper<13,V>
    {
        typedef typename V::float_type RT;
        static RT call(const V& v)
        {
            const RT eps = TMV_Epsilon<RT>();
            const RT vnormsq = v.normSq();

            if (TMV_Underflow(vnormsq)) {
                // Possible underflow errors:

                // If vmax = 0, then norm2 = 0:
                RT vmax = v.maxAbs2Element();
                if (vmax == RT(0)) {
                    return RT(0);
                } else if (TMV_Underflow(vmax * vmax)) {
                    // If vmax^2 underflows, but vmax != 0, then vnormsq has
                    // underflow rounding errors.  Find a better scaling.
                    // eps is a pure power of 2, so no rounding errors from
                    // rescaling by a power of eps.
                    const RT inveps = RT(1)/eps;
                    RT scale = inveps;
                    vmax *= scale;
                    RT eps2 = eps*eps;
                    while (vmax < eps2) { scale *= inveps; vmax *= inveps; }
                    return TMV_SQRT(v.normSq(scale))/scale;
                } else {
                    return TMV_SQRT(vnormsq);
                }
            } else if (TMV_Underflow(RT(1)/vnormsq)) {
                // Possible overflow errors:

                // If 1/vmax == 0, then vmax is already inf, so no hope of
                // making it more accurate.  (And need to check, since 
                // otherwise the next section would cause an infinite loop.)
                RT vmax = v.maxAbs2Element();
                if (RT(1)/vmax == RT(0)) {
                    return vmax;
                } else if (TMV_Underflow(RT(1)/(v.nElements()*vmax*vmax))) {
                    // If 1/(vmax^2) underflows, then vnormsq has overflow 
                    // errors.  Find a better scaling.
                    RT scale = eps;
                    vmax *= scale;
                    while (vmax > RT(1)) { scale *= eps; vmax *= eps; }
                    return TMV_SQRT(v.normSq(scale))/scale;
                } else {
                    return TMV_SQRT(vnormsq);
                }
            } else {
                // No problems with overflow or underflow.
                return TMV_SQRT(vnormsq);
            }
        }
    };

    // algo 90: Call inst
    template <class V>
    struct Norm_Helper<90,V>
    {
        typedef typename V::float_type RT;
        template <class V2>
        static TMV_INLINE RT call(const BaseVector<V2>& v)
        { return InstNorm2(v.vec().xView()); }
        template <class M2>
        static TMV_INLINE RT call(const BaseMatrix_Rec<M2>& m)
        { return InstNormF(m.mat().xView()); }
        template <class M2>
        static TMV_INLINE RT call(const BaseMatrix_Tri<M2>& m)
        { return InstNormF(m.mat().xView()); }
    };

    // algo 95: Conjugate Matrix
    template <class M>
    struct Norm_Helper<95,M>
    {
        typedef typename M::float_type RT;
        static TMV_INLINE RT call(const M& m)
        {
            typedef typename M::const_nonconj_type Mnc;
            Mnc mnc = m.nonConj();
            return Norm_Helper<-2,Mnc>::call(mnc);
        }
    };

    // algo 96: Transpose
    template <class V>
    struct Norm_Helper<96,V>
    {
        typedef typename V::float_type RT;
        static TMV_INLINE RT call(const V& v)
        {
            typedef typename V::const_transpose_type Vt;
            Vt vt = v.transpose();
            return Norm_Helper<-2,Vt>::call(vt);
        }
    };

    // algo 97: Conjugate Vector
    template <class V>
    struct Norm_Helper<97,V>
    {
        typedef typename V::float_type RT;
        static TMV_INLINE RT call(const V& v)
        {
            typedef typename V::const_nonconj_type Vnc;
            Vnc vnc = v.nonConj();
            return Norm_Helper<-1,Vnc>::call(vnc);
        }
    };

    // algo -3: Select algorithm
    template <class V>
    struct Norm_Helper<-3,V>
    {
        typedef typename V::float_type RT;
        static TMV_INLINE RT call(const V& v)
        {
            typedef typename V::value_type T;
            const int algo = 
                TMV_OPT == 0 ? 11 :
                Traits<T>::isinteger ? 11 :
                13;
            return Norm_Helper<algo,V>::call(v);
        }
    };

    // algo -2: Check for inst - Matrix
    template <class M>
    struct Norm_Helper<-2,M>
    {
        typedef typename M::float_type RT;
        static TMV_INLINE RT call(const M& m)
        {
            typedef typename M::value_type VT;
            const bool inst = 
                (M::_colsize == TMV_UNKNOWN || M::_colsize > 16) &&
                (M::_rowsize == TMV_UNKNOWN || M::_rowsize > 16) &&
                Traits<VT>::isinst;
            const bool up = ShapeTraits<M::_shape>::upper;
            const bool lo = ShapeTraits<M::_shape>::lower;
            const int algo = 
                (lo && !up) ? 96 :
                M::_conj ? 95 :
                inst ? 90 : 
                -3;
            return Norm_Helper<algo,M>::call(m);
        }
    };

    // algo -1: Check for inst - Vector
    template <class V>
    struct Norm_Helper<-1,V>
    {
        typedef typename V::float_type RT;
        static TMV_INLINE RT call(const V& v)
        {
            typedef typename V::value_type VT;
            const bool inst = 
                (V::_size == TMV_UNKNOWN || V::_size > 16) &&
                Traits<VT>::isinst;
            const int algo = 
                V::_conj ? 97 :
                inst ? 90 : 
                -3;
            return Norm_Helper<algo,V>::call(v);
        }
    };



    template <class V>
    static inline typename V::float_type InlineNorm2(
        const BaseVector_Calc<V>& v)
    {
        typedef typename V::const_cview_type Vv;
        TMV_MAYBE_CREF(V,Vv) vv = v.cView();
        return Norm_Helper<-3,Vv>::call(vv);
    }

    template <class V>
    static inline typename V::float_type DoNorm2(const BaseVector_Calc<V>& v)
    {
        typedef typename V::const_cview_type Vv;
        TMV_MAYBE_CREF(V,Vv) vv = v.cView();
        return Norm_Helper<-1,Vv>::call(vv);
    }

    template <class M>
    static inline typename M::float_type InlineNormF(const BaseMatrix_Rec<M>& m)
    {   
        typedef typename M::const_cview_type Mv;
        TMV_MAYBE_CREF(M,Mv) mv = m.cView();
        return Norm_Helper<-3,Mv>::call(mv);
    }   

    template <class M>
    static inline typename M::float_type DoNormF(const BaseMatrix_Rec<M>& m)
    {
        typedef typename M::const_cview_type Mv;
        TMV_MAYBE_CREF(M,Mv) mv = m.cView();
        return Norm_Helper<-2,Mv>::call(mv);
    }


    template <class M>
    static inline typename M::float_type InlineNormF(const BaseMatrix_Tri<M>& m)
    {   
        typedef typename M::const_cview_type Mv;
        TMV_MAYBE_CREF(M,Mv) mv = m.cView();
        return Norm_Helper<-3,Mv>::call(mv);
    }   

    template <class M>
    static inline typename M::float_type DoNormF(const BaseMatrix_Tri<M>& m)
    {
        typedef typename M::const_cview_type Mv;
        TMV_MAYBE_CREF(M,Mv) mv = m.cView();
        return Norm_Helper<-2,Mv>::call(mv);
    }


    //
    // Norm2
    //

    template <int algo, class M>
    struct Norm2M_Helper;

    // algo 1: Transpose m
    template <class M>
    struct Norm2M_Helper<1,M>
    {
        typedef typename M::float_type FT;
        static inline FT call(const M& m)
        {
#ifdef PRINTALGO_NormM
            const int N = m.rowsize();
            std::cout<<"Norm2 algo 1: N = "<<N<<std::endl;
#endif
            typedef typename M::const_transpose_type Mt;
            Mt mt = m.transpose();
            return Norm2M_Helper<-4,Mt>::call(mt);
        }
    };

    // algo 11: Calculate SV decomposition on the spot.
    template <class M>
    struct Norm2M_Helper<11,M>
    {
        typedef typename M::float_type FT;
        static inline FT call(const M& m)
        {
#ifdef PRINTALGO_NormM
            const int N = m.rowsize();
            std::cout<<"Norm2 algo 11: N = "<<N<<std::endl;
#endif
            typedef typename M::zfloat_type ZFT;
            const int cs = M::_colsize;
            const int rs = M::_rowsize;
            const bool rm = M::_rowmajor;
            typedef typename MCopyHelper<ZFT,Rec,cs,rs,rm>::type Mc;
            Mc mc = m;
            DiagMatrix<FT> S(m.rowsize());
            SV_Decompose(mc,S,false);
            return S.cref(0);
        }
    };

    // algo 12: Use Divider
    template <class M>
    struct Norm2M_Helper<12,M>
    {
        typedef typename M::float_type FT;
        static FT call(const M& m)
        {
#ifdef PRINTALGO_NormM
            const int N = m.rowsize();
            std::cout<<"Norm2 algo 12: N = "<<N<<std::endl;
#endif
            m.setDiv();
            // The Traits<> thing is to get rid of any reference.
            typedef typename Traits<typename M::svd_type>::type svd_type;
            TMVAssert(dynamic_cast<const svd_type*>(m.getDiv()));
            const svd_type* div = static_cast<const svd_type*>(m.getDiv());
            FT norm2 = div->norm2();
            m.doneDiv();
            return norm2;
        }
    };

    // algo 31: Use Divider if it is SV
    template <class M>
    struct Norm2M_Helper<31,M>
    {
        typedef typename M::float_type FT;
        static FT call(const M& m)
        {
#ifdef PRINTALGO_NormM
            const int N = m.rowsize();
            std::cout<<"Norm2 algo 31: N = "<<N<<std::endl;
#endif
            if (m.getDivType() == tmv::SV) 
                return Norm2M_Helper<12,M>::call(m);
            else 
                return Norm2M_Helper<-4,M>::call(m);
        }
    };

    // algo 32: Check if we need to transpose
    template <class M>
    struct Norm2M_Helper<32,M>
    {
        typedef typename M::float_type FT;
        static FT call(const M& m)
        {
#ifdef PRINTALGO_NormM
            const int N = m.rowsize();
            std::cout<<"Norm2 algo 32: N = "<<N<<std::endl;
#endif
            if (m.colsize() < m.rowsize()) 
                return Norm2M_Helper<1,M>::call(m);
            else
                return Norm2M_Helper<11,M>::call(m);
        }
    };

    // algo -4: Don't use a divider object
    template <class M>
    struct Norm2M_Helper<-4,M>
    {
        typedef typename M::float_type FT;
        static TMV_INLINE FT call(const M& m)
        {
            typedef typename M::real_type RT;
            const bool up = ShapeTraits<M::_shape>::upper;
            const bool lo = ShapeTraits<M::_shape>::lower;
            const int cs = M::_colsize;
            const int rs = M::_rowsize;
            const int algo = 
                cs != TMV_UNKNOWN && rs != TMV_UNKNOWN ? (
                    cs < rs ? 1 : 11 ) :
                ShapeTraits<M::_shape>::square ? 11 :
                32;
#ifdef PRINTALGO_NormM
            std::cout<<"Inline Norm2 (Non-divider) \n";
            std::cout<<"m = "<<TMV_Text(m)<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
#endif
            return Norm2M_Helper<algo,M>::call(m);
        }
    };

    // algo -3: Determine which algorithm to use
    template <class M>
    struct Norm2M_Helper<-3,M>
    {
        typedef typename M::float_type FT;
        static TMV_INLINE FT call(const M& m)
        {
            typedef typename M::real_type RT;
            const bool up = ShapeTraits<M::_shape>::upper;
            const bool lo = ShapeTraits<M::_shape>::lower;
            const int algo = 
                M::_hasdivider ? 31 :
                -4;
#ifdef PRINTALGO_NormM
            const int N = m.rowsize();
            std::cout<<"Inline Norm2 \n";
            std::cout<<"m = "<<TMV_Text(m)<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
#endif
            return Norm2M_Helper<algo,M>::call(m);
        }
    };

    template <class M>
    static inline typename M::float_type DoNorm2(const BaseMatrix<M>& m)
    {
        // Don't make a view, since we want to make sure we keep 
        // a divider object if one is present.
        return Norm2M_Helper<-3,M>::call(m.mat());
    }



    //
    // Condition
    //

    template <int algo, class M>
    struct ConditionM_Helper;

    // algo 1: Transpose m
    template <class M>
    struct ConditionM_Helper<1,M>
    {
        typedef typename M::float_type FT;
        static inline FT call(const M& m)
        {
#ifdef PRINTALGO_NormM
            const int N = m.rowsize();
            std::cout<<"Condition algo 1: N = "<<N<<std::endl;
#endif
            typedef typename M::const_transpose_type Mt;
            Mt mt = m.transpose();
            return ConditionM_Helper<-4,Mt>::call(mt);
        }
    };

    // TODO: Is the 2-condition really the one we want as the default?
    // algo 11: Calculate SV decomposition on the spot.
    template <class M>
    struct ConditionM_Helper<11,M>
    {
        typedef typename M::float_type FT;
        static inline FT call(const M& m)
        {
#ifdef PRINTALGO_NormM
            const int N = m.rowsize();
            std::cout<<"Condition algo 11: N = "<<N<<std::endl;
#endif
            typedef typename M::zfloat_type ZFT;
            const int cs = M::_colsize;
            const int rs = M::_rowsize;
            const bool rm = M::_rowmajor;
            typedef typename MCopyHelper<ZFT,Rec,cs,rs,rm>::type Mc;
            Mc mc = m;
            DiagMatrix<FT> S(m.rowsize());
            SV_Decompose(mc,S,false);
            return S.cref(0) / S.cref(S.size()-1);
        }
    };

    // algo 12: Use Divider
    template <class M>
    struct ConditionM_Helper<12,M>
    {
        typedef typename M::float_type FT;
        static FT call(const M& m)
        {
#ifdef PRINTALGO_NormM
            const int N = m.rowsize();
            std::cout<<"Condition algo 12: N = "<<N<<std::endl;
#endif
            m.setDiv();
            // SV doesn't use the norminf argument, so don't calculate it.
            FT norminf = m.getDivType() == tmv::SV ? FT(0) : m.normInf();
            FT cond = m.getDiv()->condition(norminf);
            m.doneDiv();
            return cond;
        }
    };

    // algo 32: Check if we need to transpose
    template <class M>
    struct ConditionM_Helper<32,M>
    {
        typedef typename M::float_type FT;
        static FT call(const M& m)
        {
#ifdef PRINTALGO_NormM
            const int N = m.rowsize();
            std::cout<<"Condition algo 32: N = "<<N<<std::endl;
#endif
            if (m.colsize() < m.rowsize()) 
                return ConditionM_Helper<1,M>::call(m);
            else
                return ConditionM_Helper<11,M>::call(m);
        }
    };

    // algo -4: Don't use a divider object
    template <class M>
    struct ConditionM_Helper<-4,M>
    {
        typedef typename M::float_type FT;
        static TMV_INLINE FT call(const M& m)
        {
            typedef typename M::real_type RT;
            const bool up = ShapeTraits<M::_shape>::upper;
            const bool lo = ShapeTraits<M::_shape>::lower;
            const int cs = M::_colsize;
            const int rs = M::_rowsize;
            const int algo = 
                cs != TMV_UNKNOWN && rs != TMV_UNKNOWN ? (
                    cs < rs ? 1 : 11 ) :
                ShapeTraits<M::_shape>::square ? 11 :
                32;
#ifdef PRINTALGO_NormM
            std::cout<<"Inline Condition (Non-divider) \n";
            std::cout<<"m = "<<TMV_Text(m)<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
#endif
            return ConditionM_Helper<algo,M>::call(m);
        }
    };

    // algo -3: Determine which algorithm to use
    template <class M>
    struct ConditionM_Helper<-3,M>
    {
        typedef typename M::float_type FT;
        static TMV_INLINE FT call(const M& m)
        {
            typedef typename M::real_type RT;
            const bool up = ShapeTraits<M::_shape>::upper;
            const bool lo = ShapeTraits<M::_shape>::lower;
            const int algo = 
                M::_hasdivider ? 12 :
                -4;
#ifdef PRINTALGO_NormM
            const int N = m.rowsize();
            std::cout<<"Inline Condition \n";
            std::cout<<"m = "<<TMV_Text(m)<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
#endif
            return ConditionM_Helper<algo,M>::call(m);
        }
    };

    template <class M>
    static inline typename M::float_type DoCondition(const BaseMatrix<M>& m)
    {
        // Don't make a view, since we want to make sure we keep 
        // a divider object if one is present.
        return ConditionM_Helper<-3,M>::call(m.mat());
    }



} // namespace tmv

#endif
