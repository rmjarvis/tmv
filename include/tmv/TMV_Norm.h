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

#ifndef TMV_Norm_H
#define TMV_Norm_H

#include "TMV_BaseVector.h"
#include "TMV_BaseMatrix_Rec.h"
#include "TMV_BaseMatrix_Tri.h"

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
        static RT call(const V& v)
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

            // If vmax = 0, then norm2 = 0:
            if (vmax == RT(0)) return RT(0);

            // If vmax^2 underflows but vmax != 0, then a naive NormSq()
            // will produce underflow rounding errors.  Find a better scaling.
            // eps is a pure power of 2, so no rounding errors from
            // rescaling by a power of eps.
            else if (TMV_Underfloat(vmax * vmax)) {
                const RT inveps = RT(1)/eps;
                RT scale = inveps;
                vmax *= scale;
                const RT eps2 = eps*eps;
                while (vmax < eps2) { scale *= inveps; vmax *= inveps; }
                return TMV_SQRT(v.normSq(scale))/scale;
            }

            // If 1/vmax == 0, then vmax is already inf, so no hope of
            // making it more accurate.  (And need to check, since otherwise
            // the next section would cause an infinite loop.)
            else if (RT(1)/vmax == RT(0)) {
                return vmax;
            }

            // If 1/(n*vmax^2) underflows, then a naive NormSq() will produce 
            // overflow.  Find a better scaling.
            else if (TMV_Underflow(RT(1)/(v.nElements()*vmax*vmax))) {
                const RT inveps = RT(1)/eps;
                RT scale = eps;
                vmax *= scale;
                while (vmax > inveps) { scale *= eps; vmax *= eps; }
                return TMV_SQRT(v.normSq(scale))/scale;
            }

            // No problems with overflow or underflow.
            else return TMV_SQRT(v.normSq());
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
                if (vmax == RT(0)) return RT(0);

                // If vmax^2 underflows, but vmax != 0, then vnormsq has
                // underflow rounding errors.  Find a better scaling.
                // eps is a pure power of 2, so no rounding errors from
                // rescaling by a power of eps.
                else if (TMV_Underflow(vmax * vmax)) {
                    const RT inveps = RT(1)/eps;
                    RT scale = inveps;
                    vmax *= scale;
                    RT eps2 = eps*eps;
                    while (vmax < eps2) { scale *= inveps; vmax *= inveps; }
                    return TMV_SQRT(v.normSq(scale))/scale;
                }

                else return TMV_SQRT(vnormsq);
            }

            else if (TMV_Underflow(RT(1)/vnormsq)) {
                // Possible overflow errors:

                // If 1/vmax == 0, then vmax is already inf, so no hope of
                // making it more accurate.  (And need to check, since 
                // otherwise the next section would cause an infinite loop.)
                RT vmax = v.maxAbs2Element();
                if (RT(1)/vmax == RT(0)) {
                    return vmax;
                }

                // If 1/(vmax^2) underflows, then vnormsq has overflow errors. 
                // Find a better scaling.
                else if (TMV_Underflow(RT(1)/(v.nElements()*vmax*vmax))) {
                    RT scale = eps;
                    vmax *= scale;
                    while (vmax > RT(1)) { scale *= eps; vmax *= eps; }
                    return TMV_SQRT(v.normSq(scale))/scale;
                }

                else return TMV_SQRT(vnormsq);
            }

            // No problems with overflow or underflow.
            else return TMV_SQRT(vnormsq);
        }
    };

    // algo 90: Call inst
    template <class V>
    struct Norm_Helper<90,V>
    {
        typedef typename V::float_type RT;
        template <class V2>
        static RT call(const BaseVector<V2>& v)
        { return InstNorm2(v.vec().xView()); }
        template <class M2>
        static RT call(const BaseMatrix_Rec<M2>& m)
        { return InstNormF(m.mat().xView()); }
        template <class M2>
        static RT call(const BaseMatrix_Tri<M2>& m)
        { return InstNormF(m.mat().xView()); }
    };

    // algo 95: Conjugate Matrix
    template <class M>
    struct Norm_Helper<95,M>
    {
        typedef typename M::float_type RT;
        static RT call(const M& m)
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
        static RT call(const V& v)
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
        static RT call(const V& v)
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
        static RT call(const V& v)
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
        static RT call(const M& m)
        {
            typedef typename M::value_type VT;
            const bool inst = 
                (M::_colsize == UNKNOWN || M::_colsize > 16) &&
                (M::_rowsize == UNKNOWN || M::_rowsize > 16) &&
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
        static RT call(const V& v)
        {
            typedef typename V::value_type VT;
            const bool inst = 
                (V::_size == UNKNOWN || V::_size > 16) &&
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

} // namespace tmv

#endif
