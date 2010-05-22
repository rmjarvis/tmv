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

#ifndef TMV_Det_H
#define TMV_Det_H

#include "TMV_BaseVector.h"
#include "TMV_BaseMatrix.h"
#include "TMV_Scaling.h"
#include "TMV_MinMax.h"

#ifdef PRINTALGO_Det
#include <iostream>
#endif

namespace tmv {

    // Defined below:
    template <class M>
    inline typename M::value_type Det(const BaseMatrix_Calc<M>& m);
    template <class M>
    inline typename M::value_type InlineDet(const BaseVector_Calc<M>& m);
    template <class M>
    inline typename M::real_type LogDet(
        const BaseMatrix_Calc<M>& m, typename M::value_type* sign);
    template <class M>
    inline typename M::real_type InlineLogDet(
        const BaseVector_Calc<M>& m, typename M::value_type* sign);
    template <class M>
    inline bool IsSingular(const BaseMatrix_Calc<M>& m);
    template <class M>
    inline bool InlineIsSingular(const BaseVector_Calc<M>& m);

    template <class V>
    inline typename V::value_type ProdElements(const BaseVector_Calc<V>& v);
    template <class V>
    inline typename V::value_type InlineProdElements(
        const BaseVector_Calc<V>& v);

    template <class V>
    inline typename V::real_type LogProdElements(
        const BaseVector_Calc<V>& v, typename V::value_type* sign);
    template <class V>
    inline typename V::real_type InlineLogProdElements(
        const BaseVector_Calc<V>& v, typename V::value_type* sign);

    template <class V>
    inline bool HasZeroElement(const BaseVector_Calc<V>& v);
    template <class V>
    inline bool InlineHasZeroElement(const BaseVector_Calc<V>& v);


    // Defined in TMV_Det.cpp:
    template <class T>
    T InstProdElements(const ConstVectorView<T>& v);
    template <class T>
    typename Traits<T>::real_type InstLogProdElements(
        const ConstVectorView<T>& v, T* sign);
    template <class T>
    bool InstHasZeroElement(const ConstVectorView<T>& v);


    //
    // Det
    //

    template <int algo, int size, class M> 
    struct DetM_Helper;

    // algo 0: size == 0, det = 1 (by definition)
    // Also unit-diag triangle matrices.
    template <int size, class M>
    struct DetM_Helper<0,size,M> 
    {
        typedef typename M::value_type T;
        static inline T call(const M& m)
        {
#ifdef PRINTALGO_Det
            const int N = m.rowsize();
            std::cout<<"Det algo 0: N,s = "<<N<<','<<size<<std::endl;
#endif
            return T(1); 
        }
    };

    // algo 1: size == 1
    template <class M>
    struct DetM_Helper<1,1,M> 
    {
        typedef typename M::value_type T;
        static inline T call(const M& m)
        { 
#ifdef PRINTALGO_Det
            const int N = m.rowsize();
            std::cout<<"Det algo 1: N,s = "<<N<<','<<1<<std::endl;
#endif
            return m.cref(0,0); 
        }
    };

    // algo 2: size == 2
    template <class M>
    struct DetM_Helper<2,2,M> 
    {
        typedef typename M::value_type T;
        static inline T call(const M& m)
        {
#ifdef PRINTALGO_Det
            const int N = m.rowsize();
            std::cout<<"Det algo 2: N,s = "<<N<<','<<2<<std::endl;
#endif
            return m.cref(0,0)*m.cref(1,1) - m.cref(0,1)*m.cref(1,0); 
        }
    };

    // algo 3: size == 3
    template <class M>
    struct DetM_Helper<3,3,M> 
    {
        typedef typename M::value_type T;
        static inline T call(const M& m)
        { 
#ifdef PRINTALGO_Det
            const int N = m.rowsize();
            std::cout<<"Det algo 3: N,s = "<<N<<','<<3<<std::endl;
#endif
            return 
                m.cref(0,0)*(m.cref(1,1)*m.cref(2,2)-m.cref(1,2)*m.cref(2,1)) -
                m.cref(0,1)*(m.cref(1,0)*m.cref(2,2)-m.cref(1,2)*m.cref(2,0)) +
                m.cref(0,2)*(m.cref(1,0)*m.cref(2,1)-m.cref(1,1)*m.cref(2,0));
        }
    };

    // TODO: There is a fast 4x4 calculation.  Probably worth putting in.
    // Also would benefit from SSE commands.  Especially for float.

    // algo 11: Direct product of diagonal:
    template <int size, class M>
    struct DetM_Helper<11,size,M> 
    {
        typedef typename M::value_type T;
        static inline T call(const M& m)
        {
#ifdef PRINTALGO_Det
            const int N = m.rowsize();
            std::cout<<"Det algo 11: N,s = "<<N<<','<<size<<std::endl;
            std::cout<<"m = "<<m<<std::endl;
            T det = m.diag().prodElements(); 
            std::cout<<"det = "<<det<<std::endl;
            return det;
#else
            return m.diag().prodElements(); 
#endif
        }
    };

    // algo 12: Use Divider
    template <int size, class M>
    struct DetM_Helper<12,size,M> 
    {
        typedef typename M::value_type T;
        static inline T call(const M& m)
        {
#ifdef PRINTALGO_Det
            const int N = m.rowsize();
            std::cout<<"Det algo 12: N,s = "<<N<<','<<size<<std::endl;
            std::cout<<"m = "<<m<<std::endl;
#endif
            m.setDiv();
            T det = m.getDiv()->det();
#ifdef PRINTALGO_Det
            std::cout<<"det = "<<det<<std::endl;
#endif
            m.doneDiv();
            return det;
        }
    };

    // algo 13: Calculate LU decomposition on the spot.
    template <int size, class M>
    struct DetM_Helper<13,size,M> 
    {
        typedef typename M::value_type T;
        static inline T call(const M& m)
        {
#ifdef PRINTALGO_Det
            const int N = m.rowsize();
            std::cout<<"Det algo 13: N,s = "<<N<<','<<size<<std::endl;
#endif
            return m.lud().det(); 
        }
    };

    // algo -3: Determine which algorithm to use
    template <int size, class M>
    struct DetM_Helper<-3,size,M> 
    {
        typedef typename M::value_type T;
        static inline T call(const M& m)
        {
            const int algo = 
                size == 0 ? 0 :
                ShapeTraits<M::_shape>::unit ? 0 :
                size == 1 ? 1 :
                !ShapeTraits<M::_shape>::upper ? 11 :
                !ShapeTraits<M::_shape>::lower ? 11 :
                size == 2 ? 2 : 
                size == 3 ? 3 :
                M::_hasdivider ? 12 :
                13;
#ifdef PRINTALGO_Det
            const int N = m.rowsize();
            std::cout<<"Inline Det N,s = "<<N<<','<<size<<std::endl;
            std::cout<<"m = "<<TMV_Text(m)<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
#endif
            return DetM_Helper<algo,size,M>::call(m);
        }
    };

    template <class M>
    inline typename M::value_type Det(const BaseMatrix_Calc<M>& m)
    {
        TMVStaticAssert((Sizes<M::_colsize,M::_rowsize>::same));
        TMVAssert(m.colsize() == m.rowsize());
        const int size = Sizes<M::_colsize,M::_rowsize>::size;
        // Don't make a view, since we want to make sure we keep 
        // a divider object if one is present.
        return DetM_Helper<-3,size,M>::call(m.mat());
    }


    //
    // LogDet
    //

    template <int algo, int size, class M>
    struct LogDetM_Helper;

    // algo 0: Det = 1.
    template <int size, class M>
    struct LogDetM_Helper<0,size,M> 
    {
        typedef typename M::real_type RT;
        typedef typename M::value_type T;
        static inline RT call(const M& m, T* sign)
        { 
            if (sign) *sign = T(1);
            return 0;
        }
    };

    // algo 1: Log of direct det calculation.
    template <int size, class M>
    struct LogDetM_Helper<1,size,M> 
    {
        typedef typename M::real_type RT;
        typedef typename M::value_type T;
        static inline RT call(const M& m, T* sign)
        { 
            T det = Det(m);
            RT absdet = TMV_ABS(det);
            RT logdet = TMV_LOG(absdet);
            if (sign) *sign = TMV_SIGN(det,absdet);
            return logdet;
        }
    };

    // algo 11: Direct log product of diagonal:
    template <int size, class M>
    struct LogDetM_Helper<11,size,M> 
    {
        typedef typename M::real_type RT;
        typedef typename M::value_type T;
        static inline RT call(const M& m, T* sign)
        { return m.diag().logProdElements(sign); }
    };

    // algo 12: Use Divider 
    template <int size, class M>
    struct LogDetM_Helper<12,size,M> 
    {
        typedef typename M::real_type RT;
        typedef typename M::value_type T;
        static inline RT call(const M& m, T* sign)
        { 
            m.setDiv();
            RT logdet = m.getDiv()->logDet(sign);
            m.doneDiv();
            return logdet;
        }
    };

    // algo 13: Calculate LU decomposition on the spot.
    template <int size, class M>
    struct LogDetM_Helper<13,size,M> 
    {
        typedef typename M::real_type RT;
        typedef typename M::value_type T;
        static inline RT call(const M& m, T* sign)
        { return m.lud().logDet(sign); }
    };

    // algo -3: LogDetermine which algorithm to use
    template <int size, class M>
    struct LogDetM_Helper<-3,size,M> 
    {
        typedef typename M::real_type RT;
        typedef typename M::value_type T;
        static inline RT call(const M& m, T* sign)
        {
            const int algo = 
                size == 0 ? 0 :
                ShapeTraits<M::_shape>::unit ? 0 :
                size == 1 ? 1 :
                !ShapeTraits<M::_shape>::upper ? 11 :
                !ShapeTraits<M::_shape>::lower ? 11 :
                size != UNKNOWN && size <= 3 ? 1 :
                M::_hasdivider ? 12 :
                13;
#ifdef PRINTALGO_Det
            const int N = m.rowsize();
            std::cout<<"Inline LogDet N,s = "<<N<<','<<size<<std::endl;
            std::cout<<"m = "<<TMV_Text(m)<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
#endif
            return LogDetM_Helper<algo,size,M>::call(m,sign);
        }
    };

    template <class M>
    inline typename M::real_type LogDet(
        const BaseMatrix_Calc<M>& m, typename M::value_type* sign)
    {
        TMVStaticAssert((Sizes<M::_colsize,M::_rowsize>::same));
        TMVAssert(m.colsize() == m.rowsize());
        const int size = Sizes<M::_colsize,M::_rowsize>::size;
        // Don't make a view, since we want to make sure we keep 
        // a divider object if one is present.
        return LogDetM_Helper<-3,size,M>::call(m.mat(),sign);
    }


    //
    // IsSingular
    //

    template <int algo, int size, class M>
    struct IsSingularM_Helper;

    // algo 0: Trivially non-singular
    template <int size, class M>
    struct IsSingularM_Helper<0,size,M> 
    {
        static inline bool call(const M& m)
        { return false; }
    };

    // algo 1: Check if Det == 0
    template <int size, class M>
    struct IsSingularM_Helper<1,size,M> 
    {
        static inline bool call(const M& m)
        { return Det(m) == typename M::real_type(0); }
    };

    // algo 11: Look for a zero on the diagonal
    template <int size, class M>
    struct IsSingularM_Helper<11,size,M> 
    {
        static inline bool call(const M& m)
        { return m.diag().hasZeroElement(); }
    };

    // algo 12: Use Divider
    template <int size, class M>
    struct IsSingularM_Helper<12,size,M> 
    {
        static inline bool call(const M& m)
        { 
            m.setDiv();
            bool sing = m.getDiv()->isSingular();
            m.doneDiv();
            return sing;
        }
    };

    // algo 13: Calculate LU decomposition on the spot.
    template <int size, class M>
    struct IsSingularM_Helper<13,size,M> 
    {
        static inline bool call(const M& m)
        { return m.lud().isSingular(); }
    };

    // algo 21: TriMatrix -- need to check for UnknownDiag
    template <int size, class M>
    struct IsSingularM_Helper<21,size,M> 
    {
        static inline bool call(const M& m)
        { 
            const int algo2 = 
                M::_unit ? 0 :
                M::_unknowndiag ? 22 :
                11;
            return IsSingularM_Helper<algo2,size,M>::call(m); 
        }
    };

    // algo 22: UnknownDiag TriMatrix, so might be trivial.
    template <int size, class M>
    struct IsSingularM_Helper<22,size,M> 
    {
        static inline bool call(const M& m)
        { 
            if (m.isunit()) return false;
            else return IsSingularM_Helper<11,size,M>::call(m); 
        }
    };

    // algo -3: Determine which algorithm to use
    template <int size, class M>
    struct IsSingularM_Helper<-3,size,M> 
    {
        static inline bool call(const M& m)
        {
            const bool up = ShapeTraits<M::_shape>::upper;
            const bool lo = ShapeTraits<M::_shape>::lower;
            const int algo = 
                size == 0 ? 0 :
                ShapeTraits<M::_shape>::unit ? 0 :
                size == 1 ? 1 :
                !up && !lo ? 11 :
                !up || !lo ? 21 :
                size != UNKNOWN && size <= 3 ? 1 :
                M::_hasdivider ? 12 :
                13;
#ifdef PRINTALGO_Det
            const int N = m.rowsize();
            std::cout<<"Inline IsSingular N,s = "<<N<<','<<size<<std::endl;
            std::cout<<"m = "<<TMV_Text(m)<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
#endif
            return IsSingularM_Helper<algo,size,M>::call(m);
        }
    };

    template <class M>
    inline bool IsSingular(const BaseMatrix_Calc<M>& m)
    {
        TMVStaticAssert((Sizes<M::_colsize,M::_rowsize>::same));
        TMVAssert(m.colsize() == m.rowsize());
        const int size = Sizes<M::_colsize,M::_rowsize>::size;
        // Don't make a view, since we want to make sure we keep 
        // a divider object if one is present.
        return IsSingularM_Helper<-3,size,M>::call(m.mat());
    }



    // 
    // ProdElements
    //

    // TODO: Write SSE algo's.
    template <int algo, int size, class V>
    struct ProdElementsV_Helper;

    // algo 11: simple for loop
    template <int size, class V>
    struct ProdElementsV_Helper<11,size,V> 
    {
        typedef typename V::value_type T;
        static inline T call(const V& v)
        {
            const int n = size == UNKNOWN ? int(v.size()) : size;
            T prod = T(1);
            for(int i=0;i<n;++i) prod *= v.cref(i);
            return prod;
        }
    };

    // algo 12: 2 at a time
    template <int size, class V>
    struct ProdElementsV_Helper<12,size,V> 
    {
        typedef typename V::value_type T;
        static inline T call(const V& v)
        {
            const int n = size == UNKNOWN ? int(v.size()) : size;
            T prod0(1), prod1(1);
            typename V::const_iterator it = v.begin();
            int n_2 = (n>>1);
            const int nb = n-(n_2<<1);

            if (n_2) {
                do {
                    prod0 *= it[0];
                    prod1 *= it[1]; it += 2;
                } while (--n_2);
                prod0 *= prod1;
            }
            if (nb) {
                prod0 *= *it;
            }
            return prod0;
        }
    };

    // algo 13: 4 at a time
    template <int size, class V>
    struct ProdElementsV_Helper<13,size,V> 
    {
        typedef typename V::value_type T;
        static inline T call(const V& v)
        {
            const int n = size == UNKNOWN ? int(v.size()) : size;
            T prod0(1), prod1(1);
            typename V::const_iterator it = v.begin();
            int n_4 = (n>>2);
            int nb = n-(n_4<<2);

            if (n_4) {
                do {
                    prod0 *= it[0];
                    prod1 *= it[1];
                    prod0 *= it[2];
                    prod1 *= it[3]; it += 4;
                } while (--n_4);
                prod0 *= prod1;
            }
            if (nb) do {
                prod0 *= *it++;
            } while (--nb);
            return prod0;
        }
    };

    // algo 14: 8 at a time
    template <int size, class V>
    struct ProdElementsV_Helper<14,size,V> 
    {
        typedef typename V::value_type T;
        static inline T call(const V& v)
        {
            const int n = size == UNKNOWN ? int(v.size()) : size;
            T prod0(1), prod1(1), prod2(1), prod3(1);
            typename V::const_iterator it = v.begin();
            int n_8 = (n>>3);
            int nb = n-(n_8<<3);

            if (n_8) {
                do {
                    prod0 *= it[0];
                    prod1 *= it[1];
                    prod2 *= it[2];
                    prod3 *= it[3];
                    prod0 *= it[4];
                    prod1 *= it[5];
                    prod2 *= it[6];
                    prod3 *= it[7]; it += 8;
                } while (--n_8);
                prod0 *= prod2;
                prod1 *= prod3;
                prod0 *= prod1;
            }
            if (nb) do {
                prod0 *= *it++;
            } while (--nb);
            return prod0;
        }
    };

    // algo 15: fully unroll
    template <int size, class V>
    struct ProdElementsV_Helper<15,size,V>
    {
        typedef typename V::value_type T;

        template <int I, int N>
        struct Unroller
        {
            static inline T unroll(const V& v)
            {
                return (
                    Unroller<I,N/2>::unroll(v) *
                    Unroller<I+N/2,N-N/2>::unroll(v));
            }
        };
        template <int I>
        struct Unroller<I,1>
        { static inline T unroll(const V& v) { return v.cref(I); } };
        template <int I>
        struct Unroller<I,0>
        { static inline T unroll(const V& v) { return T(1); } };
        static inline T call(const V& v)
        { return Unroller<0,size>::unroll(v); }
    };

    // algo 31: If the direct product might cause an overflow, use
    // LogDet instead.
    template <int size, class V>
    struct ProdElementsV_Helper<31,size,V> 
    {
        typedef typename V::value_type T;
        typedef typename V::real_type RT;
        static inline T call(const V& v)
        {
            RT max = v.maxAbs2Element();
            RT min = v.minAbs2Element();
            if (max > RT(1) && min < RT(1)) {
                // Then it's possible for a direct product to overflow,
                // but the actual product to be calculable.
                const int n = size == UNKNOWN ? int(v.size()) : size;
                // TODO: there is probably a more efficient way to 
                // do this.  This requires 2*log(n) multiplies, which isn't
                // large compared to what will actually be done in the 
                // ProdElements or logProdElements loop, but still I 
                // suspect there is something faster than this.
                for(int k=1;k<n;k*=2) {
                    max *= max;
                    min *= min;
                }
                const RT eps = TMV_Epsilon<RT>();
                if (eps*min == RT(0) || eps/max == RT(0)) {
                    // Need to do LogDet version
                    T sign;
                    RT logdet = LogProdElements(v,&sign);
                    if (sign == T(0)) return T(0);
                    else return sign * TMV_EXP(logdet);
                } // else direct calc should be ok.  Drop through.
            }
            return ProdElementsV_Helper<-4,size,V>::call(v);
        }
    };

    // algo -4: No branches or copies
    template <int size, class V>
    struct ProdElementsV_Helper<-4,size,V> 
    {
        typedef typename V::value_type T;
        static inline T call(const V& v)
        {
#if TMV_OPT == 0
            const int algo = 11;
#else
            typedef typename V::real_type RT;
            const bool unit = V::_step == 1;
            const int algo = 
                size == 0 ? 0 :
                size == 1 ? 1 :
                (size != UNKNOWN && size <= int(128/sizeof(T))) ? 15 :
                (unit && sizeof(RT) == 8) ? (V::iscomplex ? 11 : 13) :
                (unit && sizeof(RT) == 4) ? (V::iscomplex ? 12 : 14) :
                11;
#endif
            return ProdElementsV_Helper<algo,size,V>::call(v);
        }
    };

    // algo -3: Determine which algorithm to use
    template <int size, class V>
    struct ProdElementsV_Helper<-3,size,V> 
    {
        typedef typename V::value_type T;
        static inline T call(const V& v)
        {
            const int algo = 
#if TMV_OPT >= 2
                // 20 here is a bit arbitrary.  
                // For N <= 20, the product is not very likely to overflow.
                (size == UNKNOWN || size > 20) ? 31 :
#endif
                -4;
            return ProdElementsV_Helper<algo,size,V>::call(v);
        }
    };

    // algo 97: Conjugate
    template <int size, class V>
    struct ProdElementsV_Helper<97,size,V> 
    {
        typedef typename V::value_type T;
        static inline T call(const V& v)
        { 
            typedef typename V::const_conjugate_type Vc;
            Vc vc = v.conjugate();
            return ProdElementsV_Helper<-2,size,Vc>::call(vc);
        }
    };

    // algo 98: Call inst
    template <int size, class V>
    struct ProdElementsV_Helper<98,size,V> 
    {
        typedef typename V::value_type T;
        static inline T call(const V& v)
        { return InstProdElements(v.xView()); }
    };

    // algo -2: Check for inst
    template <int size, class V>
    struct ProdElementsV_Helper<-2,size,V> 
    {
        typedef typename V::value_type T;
        static inline T call(const V& v)
        {
            const bool inst =
                V::unknownsizes &&
                Traits<T>::isinst;
            const int algo =
                V::_conj ? 97 : 
                inst ? 98 : 
                -3;
            return ProdElementsV_Helper<algo,size,V>::call(v); 
        }
    };

    // algo -1: Check for aliases? No.
    template <int size, class V>
    struct ProdElementsV_Helper<-1,size,V> 
    {
        typedef typename V::value_type T;
        static inline T call(const V& v)
        { return ProdElementsV_Helper<-1,size,V>::call(v); }
    };

    template <class V>
    inline typename V::value_type InlineProdElements(
        const BaseVector_Calc<V>& v)
    {
        typedef typename V::const_cview_type Vv;
        Vv vv = v.cView();
        return ProdElementsV_Helper<-3,V::_size,Vv>::call(vv);
    }

    template <class V>
    inline typename V::value_type ProdElements(const BaseVector_Calc<V>& v)
    {
        typedef typename V::const_cview_type Vv;
        Vv vv = v.cView();
        return ProdElementsV_Helper<-2,V::_size,Vv>::call(vv);
    }


    // 
    // LogProdElements
    //

    template <int algo, int size, class V>
    struct LogProdElementsV_Helper;

    // algo 1: If no sign given, call version without sign.
    template <int size, class V>
    struct LogProdElementsV_Helper<1,size,V> 
    {
        typedef typename V::real_type RT;
        typedef typename V::value_type T;
        static inline RT call(const V& v, T* sign)
        { 
            if (sign) {
                return LogProdElementsV_Helper<-4,size,V>::call(v,sign);
            } else {
                return LogProdElementsV_Helper<-4,size,V>::call(v);
            }
        }
    };

    // algo 11: simple for loop, real
    template <int size, class V>
    struct LogProdElementsV_Helper<11,size,V> 
    {
        typedef typename V::real_type RT;
        typedef typename V::value_type T;
        static inline RT call(const V& v, T* sign)
        {
            TMVStaticAssert(Traits<T>::isreal);
            const int n = size == UNKNOWN ? int(v.size()) : size;
            RT sum(0);
            T v0;
            typename V::const_iterator it = v.begin();
            *sign = T(1);

            for(int i=0;i<n;++i) {
                v0 = *it++;
                sum += TMV_LOG(TMV_ABS(v0));
                if (v0 > T(0)) continue;
                else if (v0 < T(0)) *sign *= T(-1);
                else { *sign = T(0); return sum; }
            }
            return sum;
        }
        static inline RT call(const V& v)
        {
            TMVStaticAssert(Traits<T>::isreal);
            const int n = size == UNKNOWN ? int(v.size()) : size;
            RT sum(0);
            T v0;
            typename V::const_iterator it = v.begin();
            for(int i=0;i<n;++i) {
                v0 = *it++;
                sum += TMV_LOG(TMV_ABS(v0));
            }
            return sum;
        }
    };

    // algo 12: 2 at a time, real
    template <int size, class V>
    struct LogProdElementsV_Helper<12,size,V> 
    {
        typedef typename V::real_type RT;
        typedef typename V::value_type T;
        static inline RT call(const V& v, T* sign)
        {
            const int n = size == UNKNOWN ? int(v.size()) : size;
            RT sum0(0), sum1(0);
            T v0, v1;
            int n_2 = (n>>1);
            const int nb = n-(n_2<<1);
            typename V::const_iterator it = v.begin();
            *sign = T(1);

            if (n_2) {
                do {
                    v0 = it[0]; v1 = it[1]; it += 2;
                    sum0 += TMV_LOG(TMV_ABS(v0));
                    sum1 += TMV_LOG(TMV_ABS(v1));
                    T temp = v0*v1;
                    if (temp > T(0)) continue;
                    else if (temp < T(0)) *sign *= T(-1);
                    else { *sign = T(0); return sum0+sum1; }
                } while (--n_2);
                sum0 += sum1;
            }
            if (nb) {
                v0 = *it;
                sum0 += TMV_LOG(TMV_ABS(v0));
                if (v0 > T(0)) {}
                else if (v0 < T(0)) *sign *= T(-1);
                else { *sign = T(0); return sum0+sum1; }
            }
            return sum0;
        }
        static inline RT call(const V& v)
        {
            const int n = size == UNKNOWN ? int(v.size()) : size;
            RT sum0(1), sum1(1);
            T v0, v1;
            typename V::const_iterator it = v.begin();
            int n_2 = (n>>1);
            const int nb = n-(n_2<<1);

            if (n_2) {
                do {
                    v0 = it[0]; v1 = it[1]; it += 2;
                    sum0 += TMV_LOG(TMV_ABS(v0));
                    sum1 += TMV_LOG(TMV_ABS(v1));
                } while (--n_2);
                sum0 += sum1;
            }
            if (nb) {
                v0 = *it;
                sum0 += TMV_LOG(TMV_ABS(v0));
            }
            return sum0;
        }
    };

    // algo 16: simple for loop, complex
    template <int size, class V>
    struct LogProdElementsV_Helper<16,size,V> 
    {
        typedef typename V::real_type RT;
        typedef typename V::value_type T;
        static inline RT call(const V& v, T* sign)
        {
            TMVStaticAssert(Traits<T>::iscomplex);
            const int n = size == UNKNOWN ? int(v.size()) : size;
            RT sum(0);
            RT absv0;
            T v0, v1;
            typename V::const_iterator it = v.begin();
            *sign = T(1);

            for(int i=0;i<n;++i) {
                v0 = *it++;
                absv0 = std::abs(v0);
                sum += std::log(absv0);
                if (v0 != T(0)) *sign *= v0/absv0; 
                else { *sign = T(0); return sum; }
            }
            return sum;
        }
        static inline RT call(const V& v)
        {
            TMVStaticAssert(Traits<T>::iscomplex);
            const int n = size == UNKNOWN ? int(v.size()) : size;
            RT sum(0);
            T v0;
            typename V::const_iterator it = v.begin();
            for(int i=0;i<n;++i) {
                v0 = *it++;
                sum += TMV_LOG(TMV_ABS(v0));
            }
            return sum;
        }
    };

    // algo -4: No branches or copies
    template <int size, class V>
    struct LogProdElementsV_Helper<-4,size,V> 
    {
        typedef typename V::real_type RT;
        typedef typename V::value_type T;
#if TMV_OPT == 0
        enum { algo = V::isreal ? 11 : 16 };
#else
        enum { algo = (
                V::iscomplex ? 16 :
                (V::_step == 1 && sizeof(RT) == 4) ? 12 : 11 ) };
#endif
        static inline RT call(const V& v, T* sign)
        { return LogProdElementsV_Helper<algo,size,V>::call(v,sign); }
        static inline RT call(const V& v)
        { return LogProdElementsV_Helper<algo,size,V>::call(v); }
    };

    // algo -3: Determine which algorithm to use
    template <int size, class V>
    struct LogProdElementsV_Helper<-3,size,V> 
    {
        typedef typename V::real_type RT;
        typedef typename V::value_type T;
        static inline RT call(const V& v, T* sign)
        { return LogProdElementsV_Helper<1,size,V>::call(v,sign); }
    };

    // algo 97: Conjugate
    template <int size, class V>
    struct LogProdElementsV_Helper<97,size,V> 
    {
        typedef typename V::real_type RT;
        typedef typename V::value_type T;
        static inline RT call(const V& v, T* sign)
        { 
            typedef typename V::const_conjugate_type Vc;
            Vc vc = v.conjugate();
            RT ret = LogProdElementsV_Helper<-2,size,Vc>::call(vc,sign);
            if (V::iscomplex && sign) *sign = TMV_CONJ(*sign);
        }
    };

    // algo 98: Call inst
    template <int size, class V>
    struct LogProdElementsV_Helper<98,size,V> 
    {
        typedef typename V::real_type RT;
        typedef typename V::value_type T;
        static inline RT call(const V& v, T* sign)
        { return InstLogProdElements(v.xView(),sign); }
    };

    // algo -2: Check for inst
    template <int size, class V>
    struct LogProdElementsV_Helper<-2,size,V> 
    {
        typedef typename V::real_type RT;
        typedef typename V::value_type T;
        static inline RT call(const V& v, T* sign)
        {
            const bool inst =
                V::unknownsizes &&
                Traits<T>::isinst;
            const int algo =
                V::_conj ? 97 : 
                inst ? 98 : 
                -3;
            return LogProdElementsV_Helper<algo,size,V>::call(v,sign); 
        }
    };

    // algo -1: Check for aliases? No.
    template <int size, class V>
    struct LogProdElementsV_Helper<-1,size,V> 
    {
        typedef typename V::real_type RT;
        typedef typename V::value_type T;
        static inline RT call(const V& v, T* sign)
        { return LogProdElementsV_Helper<-1,size,V>::call(v,sign); }
    };

    template <class V>
    inline typename V::real_type InlineLogProdElements(
        const BaseVector_Calc<V>& v, typename V::value_type* sign)
    {
        typedef typename V::const_cview_type Vv;
        Vv vv = v.cView();
        return LogProdElementsV_Helper<-3,V::_size,Vv>::call(vv,sign);
    }

    template <class V>
    inline typename V::real_type LogProdElements(
        const BaseVector_Calc<V>& v, typename V::value_type* sign)
    {
        typedef typename V::const_cview_type Vv;
        Vv vv = v.cView();
        return LogProdElementsV_Helper<-2,V::_size,Vv>::call(vv,sign);
    }


    // 
    // HasZeroElement
    //

    template <int algo, int size, class V>
    struct HasZeroElementV_Helper;

    // algo 11: simple for loop
    template <int size, class V>
    struct HasZeroElementV_Helper<11,size,V> 
    {
        typedef typename V::value_type T;
        static inline bool call(const V& v)
        {
            const int n = size == UNKNOWN ? int(v.size()) : size;
            for(int i=0;i<n;++i) {
                if (v.cref(i) != T(0)) continue;
                else return true;
            }
            return false;
        }
    };

    // algo 12: 2 at a time
    template <int size, class V>
    struct HasZeroElementV_Helper<12,size,V> 
    {
        typedef typename V::value_type T;
        static inline bool call(const V& v)
        {
            const int n = size == UNKNOWN ? int(v.size()) : size;
            typename V::const_iterator it = v.begin();
            int n_2 = (n>>1);
            const int nb = n-(n_2<<1);

            if (n_2) do {
                if (it[0] != T(0) && it[1] != T(0)) it+=2;
                else return true;
            } while (--n_2);
            if (nb) return *it == T(0);
            else return false;
        }
    };

    // algo 13: 4 at a time
    template <int size, class V>
    struct HasZeroElementV_Helper<13,size,V> 
    {
        typedef typename V::value_type T;
        static inline bool call(const V& v)
        {
            const int n = size == UNKNOWN ? int(v.size()) : size;
            typename V::const_iterator it = v.begin();
            int n_4 = (n>>2);
            int nb = n-(n_4<<2);

            if (n_4) do {
                if (it[0] != T(0) && it[1] != T(0) &&
                    it[2] != T(0) && it[3] != T(0)) it+=4;
                else return true;
            } while (--n_4);
            if (nb) do {
                if (*it++ != T(0)) continue;
                else return true;
            } while (--nb);
            return false;
        }
    };

    // algo 14: 8 at a time
    template <int size, class V>
    struct HasZeroElementV_Helper<14,size,V> 
    {
        typedef typename V::value_type T;
        static inline bool call(const V& v)
        {
            const int n = size == UNKNOWN ? int(v.size()) : size;
            typename V::const_iterator it = v.begin();
            int n_8 = (n>>3);
            int nb = n-(n_8<<3);

            if (n_8) do {
                if (it[0] != T(0) && it[1] != T(0) &&
                    it[2] != T(0) && it[3] != T(0) &&
                    it[4] != T(0) && it[5] != T(0) &&
                    it[6] != T(0) && it[7] != T(0)) it+=8;
                else return true;
            } while (--n_8);
            if (nb) do {
                if (*it++ != T(0)) continue;
                else return true;
            } while (--nb);
            return false;
        }
    };

    // algo 15: fully unroll
    template <int size, class V>
    struct HasZeroElementV_Helper<15,size,V>
    {
        typedef typename V::value_type T;

        template <int I, int N>
        struct Unroller
        {
            static inline bool unroll(const V& v)
            {
                return (
                    Unroller<I,N/2>::unroll(v) ||
                    Unroller<I+N/2,N-N/2>::unroll(v));
            }
        };
        template <int I>
        struct Unroller<I,1>
        { static inline bool unroll(const V& v) { return v.cref(I) == T(0); } };
        template <int I>
        struct Unroller<I,0>
        { static inline bool unroll(const V& v) { return false; } };
        static inline bool call(const V& v)
        { return Unroller<0,size>::unroll(v); }
    };

    // algo -4: No branches or copies
    template <int size, class V>
    struct HasZeroElementV_Helper<-4,size,V> 
    {
        typedef typename V::value_type T;
        static inline bool call(const V& v)
        {
#if TMV_OPT == 0
            const int algo = 11;
#else
            const bool unit = V::_step == 1;
            typedef typename V::real_type RT;
            const int algo = 
                (size != UNKNOWN && size <= int(128/sizeof(T))) ? 15 :
                (unit && sizeof(RT) == 8) ? (V::iscomplex ? 12 : 13) :
                (unit && sizeof(RT) == 4) ? (V::iscomplex ? 13 : 14) :
                11;
#endif
            return HasZeroElementV_Helper<algo,size,V>::call(v);
        }
    };

    // algo -3: Determine which algorithm to use
    template <int size, class V>
    struct HasZeroElementV_Helper<-3,size,V> 
    {
        typedef typename V::value_type T;
        static inline bool call(const V& v)
        { return HasZeroElementV_Helper<-4,size,V>::call(v); }
    };

    // algo 97: Conjugate
    template <int size, class V>
    struct HasZeroElementV_Helper<97,size,V> 
    {
        typedef typename V::value_type T;
        static inline bool call(const V& v)
        { 
            typedef typename V::const_conjugate_type Vc;
            Vc vc = v.conjugate();
            return HasZeroElementV_Helper<-2,size,Vc>::call(vc);
        }
    };

    // algo 98: Call inst
    template <int size, class V>
    struct HasZeroElementV_Helper<98,size,V> 
    {
        typedef typename V::value_type T;
        static inline bool call(const V& v)
        { return InstHasZeroElement(v.xView()); }
    };

    // algo -2: Check for inst
    template <int size, class V>
    struct HasZeroElementV_Helper<-2,size,V> 
    {
        typedef typename V::value_type T;
        static inline bool call(const V& v)
        {
            const bool inst =
                V::unknownsizes &&
                Traits<T>::isinst;
            const int algo =
                V::_conj ? 97 : 
                inst ? 98 : 
                -3;
            return HasZeroElementV_Helper<algo,size,V>::call(v); 
        }
    };

    // algo -1: Check for aliases? No.
    template <int size, class V>
    struct HasZeroElementV_Helper<-1,size,V> 
    {
        typedef typename V::value_type T;
        static inline bool call(const V& v)
        { return HasZeroElementV_Helper<-1,size,V>::call(v); }
    };

    template <class V>
    inline bool InlineHasZeroElement(const BaseVector_Calc<V>& v)
    {
        typedef typename V::const_cview_type Vv;
        Vv vv = v.cView();
        return HasZeroElementV_Helper<-4,V::_size,Vv>::call(vv);
    }

    template <class V>
    inline bool HasZeroElement(const BaseVector_Calc<V>& v)
    {
        typedef typename V::const_cview_type Vv;
        Vv vv = v.cView();
        return HasZeroElementV_Helper<-2,V::_size,Vv>::call(vv);
    }


} // namespace tmv

#endif
