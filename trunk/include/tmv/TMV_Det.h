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
#include "TMV_BaseMatrix_Rec.h"
#include "TMV_Scaling.h"
#include "TMV_MinMax.h"
#include "TMV_IntegerDet.h"

#ifdef PRINTALGO_Det
#include <iostream>
#include "TMV_VectorIO.h"
#endif

namespace tmv {

#define DETM_MAX_DIRECT_SIZE 3

    // Defined below:
    template <class M>
    static typename M::value_type Det(const BaseMatrix_Calc<M>& m);
    template <class M>
    static typename M::value_type InlineDet(const BaseVector_Calc<M>& m);
    template <class M>
    static typename M::float_type LogDet(
        const BaseMatrix_Calc<M>& m, typename M::zfloat_type* sign);
    template <class M>
    static typename M::float_type InlineLogDet(
        const BaseVector_Calc<M>& m, typename M::zfloat_type* sign);
    template <class M>
    static bool IsSingular(const BaseMatrix_Calc<M>& m);
    template <class M>
    static bool InlineIsSingular(const BaseVector_Calc<M>& m);

    template <class V>
    static typename V::value_type ProdElements(const BaseVector_Calc<V>& v);
    template <class V>
    static typename V::value_type InlineProdElements(
        const BaseVector_Calc<V>& v);

    template <class V>
    static typename V::float_type LogProdElements(
        const BaseVector_Calc<V>& v, typename V::zfloat_type* sign);
    template <class V>
    static typename V::float_type InlineLogProdElements(
        const BaseVector_Calc<V>& v, typename V::zfloat_type* sign);

    template <class V>
    static bool HasZeroElement(const BaseVector_Calc<V>& v);
    template <class V>
    static bool InlineHasZeroElement(const BaseVector_Calc<V>& v);


    // Defined in TMV_Det.cpp:
    template <class T>
    T InstProdElements(const ConstVectorView<T>& v);
    template <class T>
    typename ConstVectorView<T>::float_type InstLogProdElements(
        const ConstVectorView<T>& v, 
        typename ConstVectorView<T>::zfloat_type* sign);
    template <class T>
    bool InstHasZeroElement(const ConstVectorView<T>& v);


    //
    // Det
    //

    template <int algo, int s, class M> 
    struct DetM_Helper;

    // algo 0: s == 0, det = 1 (by definition)
    // Also unit-diag triangle matrices.
    template <int s, class M>
    struct DetM_Helper<0,s,M> 
    {
        typedef typename M::value_type T;
        static T call(const M& m)
        {
#ifdef PRINTALGO_Det
            const int N = m.rowsize();
            std::cout<<"Det algo 0: N,s = "<<N<<','<<s<<std::endl;
#endif
            return T(1); 
        }
    };

    // algo 1: s == 1
    template <class M>
    struct DetM_Helper<1,1,M> 
    {
        typedef typename M::value_type T;
        static T call(const M& m)
        { 
#ifdef PRINTALGO_Det
            const int N = m.rowsize();
            std::cout<<"Det algo 1: N,s = "<<N<<','<<1<<std::endl;
#endif
            return m.cref(0,0); 
        }
    };

    // algo 2: s == 2
    template <class M>
    struct DetM_Helper<2,2,M> 
    {
        typedef typename M::value_type T;
        static T call(const M& m)
        {
#ifdef PRINTALGO_Det
            const int N = m.rowsize();
            std::cout<<"Det algo 2: N,s = "<<N<<','<<2<<std::endl;
#endif
            return m.cref(0,0)*m.cref(1,1) - m.cref(0,1)*m.cref(1,0); 
        }
    };

    // algo 3: s == 3
    template <class M>
    struct DetM_Helper<3,3,M> 
    {
        typedef typename M::value_type T;
        static T call(const M& m)
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
    template <int s, class M>
    struct DetM_Helper<11,s,M> 
    {
        typedef typename M::value_type T;
        static T call(const M& m)
        {
#ifdef PRINTALGO_Det
            const int N = m.rowsize();
            std::cout<<"Det algo 11: N,s = "<<N<<','<<s<<std::endl;
#endif
            return m.diag().prodElements(); 
        }
    };

    // algo 12: Use Divider
    template <int s, class M>
    struct DetM_Helper<12,s,M> 
    {
        typedef typename M::value_type T;
        static T call(const M& m)
        {
#ifdef PRINTALGO_Det
            const int N = m.rowsize();
            std::cout<<"Det algo 12: N,s = "<<N<<','<<s<<std::endl;
#endif
            m.setDiv();
            T det = m.getDiv()->det();
            m.doneDiv();
            return det;
        }
    };

    // algo 13: Calculate LU decomposition on the spot.
    template <int s, class M>
    struct DetM_Helper<13,s,M> 
    {
        typedef typename M::value_type T;
        static T call(const M& m)
        {
#ifdef PRINTALGO_Det
            const int N = m.rowsize();
            std::cout<<"Det algo 13: N,s = "<<N<<','<<s<<std::endl;
#endif
            TMVStaticAssert(!Traits<typename M::real_type>::isinteger);
            LUD<M> lud(m,false);
            return lud.det();
        }
    };

    // algo 21: TriMatrix -- need to check for UnknownDiag
    template <int s, class M>
    struct DetM_Helper<21,s,M> 
    {
        typedef typename M::value_type T;
        static T call(const M& m)
        { 
            const int algo2 = 
                M::_unit ? 0 :
                M::_unknowndiag ? 22 :
                11;
#ifdef PRINTALGO_Det
            const int N = m.rowsize();
            std::cout<<"Det algo 21: N,s = "<<N<<','<<s<<std::endl;
            std::cout<<" -> algo "<<algo2<<std::endl;
#endif
            return DetM_Helper<algo2,s,M>::call(m); 
        }
    };

    // algo 22: UnknownDiag TriMatrix, so might be trivial.
    template <int s, class M>
    struct DetM_Helper<22,s,M> 
    {
        typedef typename M::value_type T;
        static T call(const M& m)
        { 
#ifdef PRINTALGO_Det
            const int N = m.rowsize();
            std::cout<<"Det algo 22: N,s = "<<N<<','<<s<<std::endl;
            std::cout<<"m.isunit() = "<<m.isunit()<<std::endl;
#endif
            if (m.isunit()) return T(1);
            else return DetM_Helper<11,s,M>::call(m); 
        }
    };
 
    // algo 31: For integer, unknown size
    template <int s, class M>
    struct DetM_Helper<31,s,M> 
    {
        typedef typename M::value_type T;
        static T call(const M& m)
        {
            const int N = m.rowsize();
#ifdef PRINTALGO_Det
            std::cout<<"Det algo 31: N,s = "<<N<<','<<s<<std::endl;
#endif
            if (N == 0)
                return DetM_Helper<0,0,M>::call(m);
            else if (N == 1)
                return DetM_Helper<1,1,M>::call(m);
            else if (N == 2)
                return DetM_Helper<2,2,M>::call(m);
            else if (N == 3)
                return DetM_Helper<3,3,M>::call(m);
            else
                return DetM_Helper<32,s,M>::call(m);
        }
    };

    // algo 32: Use 1x1 Bareiss algorithm.
    template <int s, class M>
    struct DetM_Helper<32,s,M> 
    {
        typedef typename M::value_type T;
        static T call(const M& m)
        {
#ifdef PRINTALGO_Det
            const int N = m.rowsize();
            std::cout<<"Det algo 32: N,s = "<<N<<','<<s<<std::endl;
#endif
            return IntegerDet(m);
        }
    };

    // algo -3: Determine which algorithm to use
    template <int s, class M>
    struct DetM_Helper<-3,s,M> 
    {
        typedef typename M::value_type T;
        static T call(const M& m)
        {
            const bool up = ShapeTraits<M::_shape>::upper;
            const bool lo = ShapeTraits<M::_shape>::lower;
            const int algo = 
                s == 0 ? 0 :
                ShapeTraits<M::_shape>::unit ? 0 :
                s == 1 ? 1 :
                !up && !lo ? 11 :
                !up || !lo ? 21 :
                s == 2 ? 2 : 
                s == 3 ? 3 :
                Traits<T>::isinteger ? 31 :
                M::_hasdivider ? 12 :
                13;

#ifdef PRINTALGO_Det
            const int N = m.rowsize();
            std::cout<<"Inline Det N,s = "<<N<<','<<s<<std::endl;
            std::cout<<"m = "<<TMV_Text(m)<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
#endif
            return DetM_Helper<algo,s,M>::call(m);
        }
    };

    template <class M>
    static typename M::value_type Det(const BaseMatrix_Calc<M>& m)
    {
        TMVStaticAssert((Sizes<M::_colsize,M::_rowsize>::same));
        TMVAssert(m.colsize() == m.rowsize());
        const int s = Sizes<M::_colsize,M::_rowsize>::size;
        // Don't make a view, since we want to make sure we keep 
        // a divider object if one is present.
        return DetM_Helper<-3,s,M>::call(m.mat());
    }


    //
    // LogDet
    //

    template <int algo, int s, class M>
    struct LogDetM_Helper;

    // algo 0: Det = 1.
    template <int s, class M>
    struct LogDetM_Helper<0,s,M> 
    {
        typedef typename M::float_type RT;
        typedef typename M::zfloat_type T;
        static RT call(const M& m, T* sign)
        { 
#ifdef PRINTALGO_Det
            const int N = m.rowsize();
            std::cout<<"LogDet algo 0: N,s = "<<N<<','<<s<<std::endl;
#endif
            if (sign) *sign = T(1);
            return RT(0);
        }
    };

    // algo 1: Log of direct det calculation.
    template <int s, class M>
    struct LogDetM_Helper<1,s,M> 
    {
        typedef typename M::float_type RT;
        typedef typename M::zfloat_type T;
        static RT call(const M& m, T* sign)
        { 
            T det = Traits<T>::convert(Det(m));
            RT absdet = TMV_ABS(det);
            RT logdet = TMV_LOG(absdet);
            if (sign) *sign = TMV_SIGN(det,absdet);
#ifdef PRINTALGO_Det
            const int N = m.rowsize();
            std::cout<<"LogDet algo 1: N,s = "<<N<<','<<s<<std::endl;
            std::cout<<"det = "<<det<<std::endl;
            std::cout<<"logdet = "<<logdet<<std::endl;
            if (sign) std::cout<<"sign = "<<*sign<<std::endl;
#endif
            return logdet;
        }
    };

    // algo 11: Direct log product of diagonal:
    template <int s, class M>
    struct LogDetM_Helper<11,s,M> 
    {
        typedef typename M::float_type RT;
        typedef typename M::zfloat_type T;
        static RT call(const M& m, T* sign)
        {
#ifdef PRINTALGO_Det
            const int N = m.rowsize();
            std::cout<<"LogDet algo 11: N,s = "<<N<<','<<s<<std::endl;
#endif
            return m.diag().logProdElements(sign); 
        }
    };

    // algo 12: Use Divider 
    template <int s, class M>
    struct LogDetM_Helper<12,s,M> 
    {
        typedef typename M::float_type RT;
        typedef typename M::zfloat_type T;
        static RT call(const M& m, T* sign)
        { 
#ifdef PRINTALGO_Det
            const int N = m.rowsize();
            std::cout<<"LogDet algo 12: N,s = "<<N<<','<<s<<std::endl;
            std::cout<<"divisset? "<<m.divIsSet()<<std::endl;
#endif
            m.setDiv();
            RT logdet = m.getDiv()->logDet(sign);
            m.doneDiv();
            return logdet;
        }
    };

    // algo 13: Calculate LU decomposition on the spot.
    template <int s, class M>
    struct LogDetM_Helper<13,s,M> 
    {
        typedef typename M::float_type RT;
        typedef typename M::zfloat_type T;
        static RT call(const M& m, T* sign)
        {
#ifdef PRINTALGO_Det
            const int N = m.rowsize();
            std::cout<<"LogDet algo 13: N,s = "<<N<<','<<s<<std::endl;
#endif
            TMVStaticAssert(!Traits<typename M::real_type>::isinteger);
            LUD<M> lud(m,false);
            return lud.logDet(sign);
        }
    };

    // algo 21: TriMatrix -- need to check for UnknownDiag
    template <int s, class M>
    struct LogDetM_Helper<21,s,M> 
    {
        typedef typename M::float_type RT;
        typedef typename M::zfloat_type T;
        static RT call(const M& m, T* sign)
        { 
            const int algo2 = 
                M::_unit ? 0 :
                M::_unknowndiag ? 22 :
                11;
#ifdef PRINTALGO_Det
            const int N = m.rowsize();
            std::cout<<"LogDet algo 21: N,s = "<<N<<','<<s<<std::endl;
            std::cout<<" -> algo "<<algo2<<std::endl;
#endif
            return LogDetM_Helper<algo2,s,M>::call(m,sign); 
        }
    };

    // algo 22: UnknownDiag TriMatrix, so might be trivial.
    template <int s, class M>
    struct LogDetM_Helper<22,s,M> 
    {
        typedef typename M::float_type RT;
        typedef typename M::zfloat_type T;
        static RT call(const M& m, T* sign)
        { 
#ifdef PRINTALGO_Det
            const int N = m.rowsize();
            std::cout<<"LogDet algo 22: N,s = "<<N<<','<<s<<std::endl;
            std::cout<<"m.isunit() = "<<m.isunit()<<std::endl;
#endif
            if (m.isunit()) return LogDetM_Helper<0,s,M>::call(m,sign); 
            else return LogDetM_Helper<11,s,M>::call(m,sign); 
        }
    };
 
    // algo -3: Determine which algorithm to use
    template <int s, class M>
    struct LogDetM_Helper<-3,s,M> 
    {
        typedef typename M::float_type RT;
        typedef typename M::zfloat_type T;
        static RT call(const M& m, T* sign)
        {
            typedef typename M::real_type MRT;
            const bool up = ShapeTraits<M::_shape>::upper;
            const bool lo = ShapeTraits<M::_shape>::lower;
            const int algo = 
                s == 0 ? 0 :
                ShapeTraits<M::_shape>::unit ? 0 :
                s == 1 ? 1 :
                Traits<MRT>::isinteger ? 1 :
                !up && !lo ? 11 :
                !up || !lo ? 21 :
                s != UNKNOWN && s <= 3 ? 1 :
                M::_hasdivider ? 12 :
                13;
#ifdef PRINTALGO_Det
            const int N = m.rowsize();
            std::cout<<"Inline LogDet N,s = "<<N<<','<<s<<std::endl;
            std::cout<<"m = "<<TMV_Text(m)<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
#endif
            return LogDetM_Helper<algo,s,M>::call(m,sign);
        }
    };

    template <class M>
    static typename M::float_type LogDet(
        const BaseMatrix_Calc<M>& m, typename M::zfloat_type* sign)
    {
        TMVStaticAssert((Sizes<M::_colsize,M::_rowsize>::same));
        TMVAssert(m.colsize() == m.rowsize());
        const int s = Sizes<M::_colsize,M::_rowsize>::size;
        // Don't make a view, since we want to make sure we keep 
        // a divider object if one is present.
        return LogDetM_Helper<-3,s,M>::call(m.mat(),sign);
    }


    //
    // IsSingular
    //

    template <int algo, int s, class M>
    struct IsSingularM_Helper;

    // algo 0: Trivially non-singular
    template <int s, class M>
    struct IsSingularM_Helper<0,s,M> 
    {
        static bool call(const M& m)
        { 
#ifdef PRINTALGO_Det
            const int N = m.rowsize();
            std::cout<<"IsSingular algo 0: N,s = "<<N<<','<<s<<std::endl;
#endif
            return false; 
        }
    };

    // algo 1: Check if Det == 0
    template <int s, class M>
    struct IsSingularM_Helper<1,s,M> 
    {
        static bool call(const M& m)
        {
            typedef typename M::real_type RT;
            RT det = Det(m);
#ifdef PRINTALGO_Det
            const int N = m.rowsize();
            std::cout<<"IsSingular algo 1: N,s = "<<N<<','<<s<<std::endl;
            std::cout<<"det = "<<det<<std::endl;
#endif
            return det == typename M::real_type(0); 
        }
    };

    // algo 11: Look for a zero on the diagonal
    template <int s, class M>
    struct IsSingularM_Helper<11,s,M> 
    {
        static bool call(const M& m)
        {
#ifdef PRINTALGO_Det
            const int N = m.rowsize();
            std::cout<<"IsSingular algo 11: N,s = "<<N<<','<<s<<std::endl;
#endif
            return m.diag().hasZeroElement(); 
        }
    };

    // algo 12: Use Divider
    template <int s, class M>
    struct IsSingularM_Helper<12,s,M> 
    {
        static bool call(const M& m)
        { 
#ifdef PRINTALGO_Det
            const int N = m.rowsize();
            std::cout<<"IsSingular algo 12: N,s = "<<N<<','<<s<<std::endl;
            std::cout<<"divisset? "<<m.divIsSet()<<std::endl;
#endif
            m.setDiv();
            bool sing = m.getDiv()->isSingular();
            m.doneDiv();
            return sing;
        }
    };

    // algo 13: Calculate LU decomposition on the spot.
    template <int s, class M>
    struct IsSingularM_Helper<13,s,M> 
    {
        static bool call(const M& m)
        {
#ifdef PRINTALGO_Det
            const int N = m.rowsize();
            std::cout<<"IsSingular algo 12: N,s = "<<N<<','<<s<<std::endl;
#endif
            TMVStaticAssert(!Traits<typename M::real_type>::isinteger);
            LUD<M> lud(m,false);
            return lud.isSingular(); 
        }
    };

    // algo 21: TriMatrix -- need to check for UnknownDiag
    template <int s, class M>
    struct IsSingularM_Helper<21,s,M> 
    {
        static bool call(const M& m)
        { 
            const int algo2 = 
                M::_unit ? 0 :
                M::_unknowndiag ? 22 :
                11;
#ifdef PRINTALGO_Det
            const int N = m.rowsize();
            std::cout<<"IsSingular algo 21: N,s = "<<N<<','<<s<<std::endl;
            std::cout<<" -> algo ="<<algo2<<std::endl;
#endif
            return IsSingularM_Helper<algo2,s,M>::call(m); 
        }
    };

    // algo 22: UnknownDiag TriMatrix, so might be trivial.
    template <int s, class M>
    struct IsSingularM_Helper<22,s,M> 
    {
        static bool call(const M& m)
        { 
#ifdef PRINTALGO_Det
            const int N = m.rowsize();
            std::cout<<"IsSingular algo 21: N,s = "<<N<<','<<s<<std::endl;
            std::cout<<"m.isunit() = "<<m.isunit()<<std::endl;
#endif
            if (m.isunit()) return false;
            else return IsSingularM_Helper<11,s,M>::call(m); 
        }
    };

    // algo -3: Determine which algorithm to use
    template <int s, class M>
    struct IsSingularM_Helper<-3,s,M> 
    {
        static bool call(const M& m)
        {
            typedef typename M::value_type T;
            const bool up = ShapeTraits<M::_shape>::upper;
            const bool lo = ShapeTraits<M::_shape>::lower;
            const int algo = 
                s == 0 ? 0 :
                ShapeTraits<M::_shape>::unit ? 0 :
                s == 1 ? 1 :
                Traits<T>::isinteger ? 1 :
                !up && !lo ? 11 :
                !up || !lo ? 21 :
                s != UNKNOWN && s <= 3 ? 1 :
                M::_hasdivider ? 12 :
                13;
#ifdef PRINTALGO_Det
            const int N = m.rowsize();
            std::cout<<"Inline IsSingular N,s = "<<N<<','<<s<<std::endl;
            std::cout<<"m = "<<TMV_Text(m)<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
#endif
            return IsSingularM_Helper<algo,s,M>::call(m);
        }
    };

    template <class M>
    static bool IsSingular(const BaseMatrix_Calc<M>& m)
    {
        TMVStaticAssert((Sizes<M::_colsize,M::_rowsize>::same));
        TMVAssert(m.colsize() == m.rowsize());
        const int s = Sizes<M::_colsize,M::_rowsize>::size;
        // Don't make a view, since we want to make sure we keep 
        // a divider object if one is present.
        return IsSingularM_Helper<-3,s,M>::call(m.mat());
    }



    // 
    // ProdElements
    //

    // Defined below, but used here.
    template <int algo, int s, class V>
    struct HasZeroElementV_Helper;
    template <int algo, int s, class V>
    struct LogProdElementsV_Helper;

    // TODO: Write SSE algo's.
    template <int algo, int s, class V>
    struct ProdElementsV_Helper;

    // algo 0: s == 0, define prod = 1
    template <int s, class V>
    struct ProdElementsV_Helper<0,s,V> 
    {
        typedef typename V::value_type T;
        static T call(const V& v)
        {
#ifdef PRINTALGO_Det
            const int n = s == UNKNOWN ? int(v.size()) : s;
            std::cout<<"Prod Elements algo 0: n,s = "<<n<<','<<s<<std::endl;
#endif
            return T(1);
        }
    };
    // algo 1: s == 1
    template <int s, class V>
    struct ProdElementsV_Helper<1,s,V> 
    {
        typedef typename V::value_type T;
        static T call(const V& v)
        {
#ifdef PRINTALGO_Det
            const int n = s == UNKNOWN ? int(v.size()) : s;
            std::cout<<"Prod Elements algo 1: n,s = "<<n<<','<<s<<std::endl;
#endif
            return v.cref(0);
        }
    };

    // algo 11: simple for loop
    template <int s, class V>
    struct ProdElementsV_Helper<11,s,V> 
    {
        typedef typename V::value_type T;
        static T call(const V& v)
        {
            const int n = s == UNKNOWN ? int(v.size()) : s;
#ifdef PRINTALGO_Det
            std::cout<<"Prod Elements algo 11: n,s = "<<n<<','<<s<<std::endl;
#endif
            if (n > 0) {
                T prod = v.cref(0);
                for(int i=1;i<n;++i) 
                prod = ZProd<false,false>::prod(prod,v.cref(i));
                return prod;
            } else return T(1);
        }
    };

    // algo 12: 2 at a time
    template <int s, class V>
    struct ProdElementsV_Helper<12,s,V> 
    {
        typedef typename V::value_type T;
        static T call(const V& v)
        {
            const int n = s == UNKNOWN ? int(v.size()) : s;
#ifdef PRINTALGO_Det
            std::cout<<"Prod Elements algo 12: n,s = "<<n<<','<<s<<std::endl;
#endif
            T prod0(1), prod1(1);
            typename V::const_iterator it = v.begin();
            int n_2 = (n>>1);
            const int nb = n-(n_2<<1);

            if (n_2) {
                do {
                    prod0 = ZProd<false,false>::prod(prod0,it[0]);
                    prod1 = ZProd<false,false>::prod(prod1,it[1]); it += 2;
                } while (--n_2);
            }
            if (nb) prod0 = ZProd<false,false>::prod(prod0,*it);
#ifdef PRINTALGO_Det
            std::cout<<"prod0, prod1 = "<<prod0<<"  "<<prod1<<std::endl;
            std::cout<<"product = "<<ZProd<false,false>::prod(prod0,prod1)<<std::endl;
#endif
            return ZProd<false,false>::prod(prod0,prod1);
        }
    };

    // algo 13: 4 at a time
    template <int s, class V>
    struct ProdElementsV_Helper<13,s,V> 
    {
        typedef typename V::value_type T;
        static T call(const V& v)
        {
            const int n = s == UNKNOWN ? int(v.size()) : s;
#ifdef PRINTALGO_Det
            std::cout<<"Prod Elements algo 13: n,s = "<<n<<','<<s<<std::endl;
#endif
            T prod0(1), prod1(1);
            typename V::const_iterator it = v.begin();
            int n_4 = (n>>2);
            int nb = n-(n_4<<2);

            if (n_4) {
                do {
                    prod0 = ZProd<false,false>::prod(prod0,it[0]);
                    prod1 = ZProd<false,false>::prod(prod1,it[1]); 
                    prod0 = ZProd<false,false>::prod(prod0,it[2]);
                    prod1 = ZProd<false,false>::prod(prod1,it[3]); it += 4;
                } while (--n_4);
            }
            if (nb) do {
                prod0 = ZProd<false,false>::prod(prod0,*it++);
            } while (--nb);
            return ZProd<false,false>::prod(prod0,prod1);
        }
    };

    // algo 14: 8 at a time
    template <int s, class V>
    struct ProdElementsV_Helper<14,s,V> 
    {
        typedef typename V::value_type T;
        static T call(const V& v)
        {
            const int n = s == UNKNOWN ? int(v.size()) : s;
#ifdef PRINTALGO_Det
            std::cout<<"Prod Elements algo 14: n,s = "<<n<<','<<s<<std::endl;
#endif
            T prod0(1), prod1(1), prod2(1), prod3(1);
            typename V::const_iterator it = v.begin();
            int n_8 = (n>>3);
            int nb = n-(n_8<<3);

            if (n_8) {
                do {
                    prod0 = ZProd<false,false>::prod(prod0,it[0]);
                    prod1 = ZProd<false,false>::prod(prod1,it[1]); 
                    prod2 = ZProd<false,false>::prod(prod2,it[2]);
                    prod3 = ZProd<false,false>::prod(prod3,it[3]);
                    prod0 = ZProd<false,false>::prod(prod0,it[4]);
                    prod1 = ZProd<false,false>::prod(prod1,it[5]); 
                    prod2 = ZProd<false,false>::prod(prod2,it[6]);
                    prod3 = ZProd<false,false>::prod(prod3,it[7]); it += 8;
                } while (--n_8);
                prod0 = ZProd<false,false>::prod(prod0,prod2);
                prod1 = ZProd<false,false>::prod(prod1,prod3);
            }
            if (nb) do {
                prod0 = ZProd<false,false>::prod(prod0,*it++);
            } while (--nb);
            return ZProd<false,false>::prod(prod0,prod1);
        }
    };

    // algo 15: fully unroll
    template <int s, class V>
    struct ProdElementsV_Helper<15,s,V>
    {
        typedef typename V::value_type T;

        template <int I, int N>
        struct Unroller
        {
            static T unroll(const V& v)
            {
                return ZProd<false,false>::prod(
                    Unroller<I,N/2>::unroll(v),
                    Unroller<I+N/2,N-N/2>::unroll(v));
            }
        };
        template <int I>
        struct Unroller<I,1>
        { static T unroll(const V& v) { return v.cref(I); } };
        template <int I>
        struct Unroller<I,0>
        { static T unroll(const V& v) { return T(1); } };
        static T call(const V& v)
        { 
#ifdef PRINTALGO_Det
            const int n = v.size();
            std::cout<<"Prod Elements algo 15: n,s = "<<n<<','<<s<<std::endl;
#endif
            return Unroller<0,s>::unroll(v); 
        }
    };

#ifdef __SSE2__
    // algo 31: double precision SSE2: real
    template <int s, class V>
    struct ProdElementsV_Helper<31,s,V> 
    {
        typedef typename V::const_iterator IT;
        static double call(const V& v)
        {
            int n = s == UNKNOWN ? int(v.size()) : s;
#ifdef PRINTALGO_Det
            std::cout<<"Prod Elements algo 31: n,s = "<<n<<','<<s<<std::endl;
#endif
            if (n) {
                IT it = v.begin();
                const bool unit = V::_step == 1;

                double prod0(1), prod1;

                if (unit) {
                    while (n && !TMV_Aligned(it.get()) ) {
                        prod0 *= *it++;
                        --n;
                    }
                }

                int n_2 = (n>>1);
                int nb = n-(n_2<<1);

                if (n_2) {
                    IT it1 = it+1;

                    union { __m128d xm; double xd[2]; } xprod;
                    xprod.xm = _mm_set1_pd(1.);
                    __m128d x1;
                    do {
                        Maybe<unit>::sse_load(x1,it.get(),it1.get());
                        it+=2; it1+=2;
                        xprod.xm = _mm_mul_pd(xprod.xm,x1);
                    } while (--n_2);
                    prod1 = xprod.xd[0] * xprod.xd[1];
                } else { prod1 = 1.; }

                if (nb) do {
                    prod0 *= *it++;
                } while (--nb);
                return prod0 * prod1;
            } else return double(1);
        }
    };

    // algo 32: double precision SSE2: complex
    template <int s, class V>
    struct ProdElementsV_Helper<32,s,V> 
    {
        typedef typename V::const_iterator IT;
        static std::complex<double> call(const V& v)
        {
            int n = s == UNKNOWN ? int(v.size()) : s;
#ifdef PRINTALGO_Det
            std::cout<<"Prod Elements algo 32: n,s = "<<n<<','<<s<<std::endl;
#endif
            if (n) {
                IT it = v.begin();

                union { __m128d xm; double xd[2]; } xprod;
                xprod.xm = _mm_set_pd(0.,1.); // == 1, since args are backwards!
                __m128d x0,x1,x2,x3,xr,xi;
                if (TMV_Aligned(it.get())) {
                    do {
                        // r = v1r * v2r - v1i * v2i
                        // i = v1r * v2i + v1i * v2r
                        xr = _mm_set_pd(xprod.xd[0],xprod.xd[0]);
                        xi = _mm_set_pd(xprod.xd[1],-xprod.xd[1]);
                        Maybe<true>::sse_load(x0,it.get()); ++it;
                        x1 = _mm_shuffle_pd(x0,x0,_MM_SHUFFLE2(0,1));
                        x2 = _mm_mul_pd(xr,x0);
                        x3 = _mm_mul_pd(xi,x1);
                        xprod.xm = _mm_add_pd(x2,x3);
                    } while (--n);
                } else {
                    do {
                        // r = v1r * v2r - v1i * v2i
                        // i = v1r * v2i + v1i * v2r
                        xr = _mm_set_pd(xprod.xd[0],xprod.xd[0]);
                        xi = _mm_set_pd(-xprod.xd[1],xprod.xd[1]);
                        Maybe<true>::sse_loadu(x0,it.get()); ++it;
                        x1 = _mm_shuffle_pd(x0,x0,_MM_SHUFFLE2(0,1));
                        x2 = _mm_mul_pd(xr,x0);
                        x3 = _mm_mul_pd(xi,x1);
                        xprod.xm = _mm_add_pd(x2,x3);
                    } while (--n);
                }
                return std::complex<double>(xprod.xd[0],xprod.xd[1]);
            } else return std::complex<double>(1);
        }
    };
#endif

    // algo 41: If the direct product might cause an overflow, use
    // LogDet instead.
    template <int s, class V>
    struct ProdElementsV_Helper<41,s,V> 
    {
        typedef typename V::value_type T;
        typedef typename V::real_type RT;
        static T call(const V& v)
        {
#ifdef PRINTALGO_Det
            const int n = v.size();
            std::cout<<"Prod Elements algo 41: n,s = "<<n<<','<<s<<std::endl;
#endif
            RT max = v.maxAbs2Element();
            RT min = v.minAbs2Element();
            if (max > RT(1) && min < RT(1)) {
                // Then it's possible for a direct product to overflow,
                // but the actual product to be calculable.
                const int n = s == UNKNOWN ? int(v.size()) : s;
                // There is probably a more efficient way to do this.
                // This requires 2*log(n) multiplies, which isn't
                // large compared to what will actually be done in the 
                // ProdElements or logProdElements loop, but still I 
                // suspect there is something faster than this.
                for(int k=1;k<n;k*=2) {
                    max *= max;
                    min *= min;
                }
                const RT eps = TMV_Epsilon<RT>();
                if (TMV_Underfloat(eps*min) || TMV_Underfloat(eps/max)) {
                    // Need to do LogDet version
                    T sign;
                    RT logdet = LogProdElements(v,&sign);
#ifdef PRINTALGO_Det
                    std::cout<<"Underflow or overflow found:\n";
                    std::cout<<"logdet, sign = "<<logdet<<" , "<<sign<<std::endl;
                    std::cout<<"exp(logdet) = "<<TMV_EXP(logdet)<<std::endl;
                    std::cout<<"sign*exp(logdet) = "<<sign*TMV_EXP(logdet)<<std::endl;
#endif
                    if (sign == T(0)) return T(0);
                    else return sign * TMV_EXP(logdet);
                } // else direct calc should be ok.  Drop through.
            }
            return ProdElementsV_Helper<-4,s,V>::call(v);
        }
    };

    // algo 42: Similar to 41, but try direct product first, and only 
    // check min,max if overflow or underflow is found.
    template <int s, class V>
    struct ProdElementsV_Helper<42,s,V> 
    {
        typedef typename V::value_type T;
        typedef typename V::real_type RT;
        static T call(const V& v)
        {
#ifdef PRINTALGO_Det
            const int n = v.size();
            std::cout<<"Prod Elements algo 42: n,s = "<<n<<','<<s<<std::endl;
#endif
            T det1 = ProdElementsV_Helper<-4,s,V>::call(v);

            // If det1 isn't 0 or inf or nan, then it's ok to return it.
            if (TMV_ABS(det1) > 0 && RT(1)/TMV_ABS(det1) > 0) return det1;

            // Check if exactly zero and not from underflow.
            if (det1 == T(0) &&
                HasZeroElementV_Helper<-3,s,V>::call(v)) return T(0);

            // See if a better value might be possible with logDet:
            RT max = v.maxAbs2Element();
            RT min = v.minAbs2Element();
            if (max > RT(1) && min < RT(1)) {
                // Then it's possible for a direct product to overflow,
                // but the actual product to be calculable using the 
                // LogDet version.
                T sign;
                RT logdet = LogProdElementsV_Helper<-4,s,V>::call(v,&sign);
#ifdef PRINTALGO_Det
                std::cout<<"Underflow or overflow found:\n";
                std::cout<<"logdet, sign = "<<logdet<<" , "<<sign<<std::endl;
                std::cout<<"exp(logdet) = "<<TMV_EXP(logdet)<<std::endl;
                std::cout<<"sign*exp(logdet) = "<<sign*TMV_EXP(logdet)<<std::endl;
#endif
                if (sign == T(0)) return T(0);
                else return sign * TMV_EXP(logdet);
            } else {
#ifdef PRINTALGO_Det
                std::cout<<"Underflow or overflow found:\n";
                std::cout<<"But can't do better.  det1 = "<<det1<<std::endl;
                std::cout<<"max,min = "<<max<<" , "<<min<<std::endl;
#endif
                // direct calc was the best we can do, so just return det1
                return det1;
            }
        }
    };

    // algo -4: No branches or copies
    template <int s, class V>
    struct ProdElementsV_Helper<-4,s,V> 
    {
        typedef typename V::value_type T;
        static T call(const V& v)
        {
#if TMV_OPT == 0
            const int algo = 11;
#else
            typedef typename V::real_type RT;
            const bool vdouble = Traits2<RT,double>::sametype;
            const bool vreal = V::isreal;
            const bool vcomplex = V::iscomplex;
            const int algo = 
                s == 0 ? 0 :
                s == 1 ? 1 :
                (s != UNKNOWN && s <= int(128/sizeof(T))) ? 15 :
#ifdef __SSE2__
                (vdouble && vreal) ? 31 :
                (vdouble && vcomplex) ? 32 :
#endif
                (sizeof(RT) == 8) ? (V::iscomplex ? 11 : 13) :
                (sizeof(RT) == 4) ? (V::iscomplex ? 12 : 14) :
                11;
#endif
#ifdef PRINTALGO_Det
            std::cout<<"No branch Prod Elements\n";
            std::cout<<"v = "<<TMV_Text(v)<<std::endl;
            std::cout<<"s = "<<s<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
#endif
            return ProdElementsV_Helper<algo,s,V>::call(v);
        }
    };

    // algo -3: Determine which algorithm to use
    template <int s, class V>
    struct ProdElementsV_Helper<-3,s,V> 
    {
        typedef typename V::value_type T;
        static T call(const V& v)
        {
            const int algo = 
#if TMV_OPT >= 2
                // 20 here is a bit arbitrary.  
                // For N <= 20, the product is not very likely to overflow.
                (!Traits<T>::isinteger && (s == UNKNOWN || s > 20)) ? 42 :
#endif
                -4;
#ifdef PRINTALGO_Det
            std::cout<<"Inline ProdElements\n";
            std::cout<<"v = "<<TMV_Text(v)<<std::endl;
            std::cout<<"v = "<<v<<std::endl;
            std::cout<<"s = "<<s<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
#endif
            return ProdElementsV_Helper<algo,s,V>::call(v);
        }
    };

    // algo 97: Conjugate
    template <int s, class V>
    struct ProdElementsV_Helper<97,s,V> 
    {
        typedef typename V::value_type T;
        static T call(const V& v)
        { 
            typedef typename V::const_conjugate_type Vc;
            Vc vc = v.conjugate();
            T ret = ProdElementsV_Helper<-1,s,Vc>::call(vc);
            return TMV_CONJ(ret);
        }
    };

    // algo 98: Call inst
    template <int s, class V>
    struct ProdElementsV_Helper<98,s,V> 
    {
        typedef typename V::value_type T;
        static T call(const V& v)
        { return InstProdElements(v.xView()); }
    };

    // algo -1: Check for inst
    template <int s, class V>
    struct ProdElementsV_Helper<-1,s,V> 
    {
        typedef typename V::value_type T;
        static T call(const V& v)
        {
            const bool inst =
                (s == UNKNOWN || s > 16) &&
                Traits<T>::isinst;
            const int algo =
                V::_conj ? 97 : 
                inst ? 98 : 
                -3;
            return ProdElementsV_Helper<algo,s,V>::call(v); 
        }
    };

    template <class V>
    static typename V::value_type InlineProdElements(
        const BaseVector_Calc<V>& v)
    {
        typedef typename V::const_cview_type Vv;
        Vv vv = v.cView();
        return ProdElementsV_Helper<-3,V::_size,Vv>::call(vv);
    }

    template <class V>
    static typename V::value_type ProdElements(const BaseVector_Calc<V>& v)
    {
        typedef typename V::const_cview_type Vv;
        Vv vv = v.cView();
        return ProdElementsV_Helper<-1,V::_size,Vv>::call(vv);
    }


    // 
    // LogProdElements
    //

    // algo 1: If no sign given, call version without sign.
    template <int s, class V>
    struct LogProdElementsV_Helper<1,s,V> 
    {
        typedef typename V::float_type RT;
        typedef typename V::zfloat_type T;
        static RT call(const V& v, T* sign)
        { 
#ifdef PRINTALGO_Det
            const int n = v.size();
            std::cout<<"LogProd Elements algo 1: n,s = "<<n<<','<<s<<std::endl;
#endif
            if (sign) {
                return LogProdElementsV_Helper<-4,s,V>::call(v,sign);
            } else {
                return LogProdElementsV_Helper<-4,s,V>::call(v);
            }
        }
    };

    // algo 11: simple for loop, real
    template <int s, class V>
    struct LogProdElementsV_Helper<11,s,V> 
    {
        typedef typename V::float_type RT;
        typedef typename V::zfloat_type T;
        static RT call(const V& v, T* sign)
        {
            TMVStaticAssert(Traits<T>::isreal);
            const int n = s == UNKNOWN ? int(v.size()) : s;
#ifdef PRINTALGO_Det
            std::cout<<"LogProd Elements algo 11: n,s = "<<n<<','<<s<<std::endl;
#endif
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
        static RT call(const V& v)
        {
            TMVStaticAssert(Traits<T>::isreal);
            const int n = s == UNKNOWN ? int(v.size()) : s;
#ifdef PRINTALGO_Det
            std::cout<<"LogProd Elements algo 11: n,s = "<<n<<','<<s<<std::endl;
#endif
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
    template <int s, class V>
    struct LogProdElementsV_Helper<12,s,V> 
    {
        typedef typename V::float_type RT;
        typedef typename V::zfloat_type T;
        static RT call(const V& v, T* sign)
        {
            const int n = s == UNKNOWN ? int(v.size()) : s;
#ifdef PRINTALGO_Det
            std::cout<<"LogProd Elements algo 12: n,s = "<<n<<','<<s<<std::endl;
#endif
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
        static RT call(const V& v)
        {
            const int n = s == UNKNOWN ? int(v.size()) : s;
#ifdef PRINTALGO_Det
            std::cout<<"LogProd Elements algo 12: n,s = "<<n<<','<<s<<std::endl;
#endif
            RT sum0(0), sum1(0);
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
    template <int s, class V>
    struct LogProdElementsV_Helper<16,s,V> 
    {
        typedef typename V::float_type RT;
        typedef typename V::zfloat_type T;
        static RT call(const V& v, T* sign)
        {
            TMVStaticAssert(Traits<T>::iscomplex);
            const int n = s == UNKNOWN ? int(v.size()) : s;
#ifdef PRINTALGO_Det
            std::cout<<"LogProd Elements algo 16: n,s = "<<n<<','<<s<<std::endl;
#endif
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
        static RT call(const V& v)
        {
            TMVStaticAssert(Traits<T>::iscomplex);
            const int n = s == UNKNOWN ? int(v.size()) : s;
#ifdef PRINTALGO_Det
            std::cout<<"LogProd Elements algo 16: n,s = "<<n<<','<<s<<std::endl;
#endif
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
    template <int s, class V>
    struct LogProdElementsV_Helper<-4,s,V> 
    {
        typedef typename V::float_type RT;
        typedef typename V::zfloat_type T;
#if TMV_OPT == 0
        enum { algo = V::isreal ? 11 : 16 };
#else
        enum { algo = (
                V::iscomplex ? 16 :
                (V::_step == 1 && sizeof(RT) == 4) ? 12 : 11 ) };
#endif
        static RT call(const V& v, T* sign)
        { 
#ifdef PRINTALGO_Det
            std::cout<<"Inline LogProdElements with sign\n";
            std::cout<<"v = "<<TMV_Text(v)<<std::endl;
            std::cout<<"s = "<<s<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
#endif
            return LogProdElementsV_Helper<algo,s,V>::call(v,sign); 
        }
        static RT call(const V& v)
        { 
#ifdef PRINTALGO_Det
            std::cout<<"Inline LogProdElements no sign\n";
            std::cout<<"v = "<<TMV_Text(v)<<std::endl;
            std::cout<<"s = "<<s<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
#endif
            return LogProdElementsV_Helper<algo,s,V>::call(v); 
        }
    };

    // algo -3: Determine which algorithm to use
    // Just check if sign is 0 and call correct version of algo -4.
    template <int s, class V>
    struct LogProdElementsV_Helper<-3,s,V> 
    {
        typedef typename V::float_type RT;
        typedef typename V::zfloat_type T;
        static RT call(const V& v, T* sign)
        { return LogProdElementsV_Helper<1,s,V>::call(v,sign); }
    };

    // algo 97: Conjugate
    template <int s, class V>
    struct LogProdElementsV_Helper<97,s,V> 
    {
        typedef typename V::float_type RT;
        typedef typename V::zfloat_type T;
        static RT call(const V& v, T* sign)
        { 
            typedef typename V::const_conjugate_type Vc;
            Vc vc = v.conjugate();
            RT ret = LogProdElementsV_Helper<-1,s,Vc>::call(vc,sign);
            if (sign) *sign = TMV_CONJ(*sign);
            return ret;
        }
    };

    // algo 98: Call inst
    template <int s, class V>
    struct LogProdElementsV_Helper<98,s,V> 
    {
        typedef typename V::float_type RT;
        typedef typename V::zfloat_type T;
        static RT call(const V& v, T* sign)
        { return InstLogProdElements(v.xView(),sign); }
    };

    // algo -1: Check for inst
    template <int s, class V>
    struct LogProdElementsV_Helper<-1,s,V> 
    {
        typedef typename V::float_type RT;
        typedef typename V::zfloat_type T;
        static RT call(const V& v, T* sign)
        {
            const bool inst =
                (s == UNKNOWN || s > 16) &&
                Traits<T>::isinst;
            const int algo =
                V::_conj ? 97 : 
                inst ? 98 : 
                -3;
            return LogProdElementsV_Helper<algo,s,V>::call(v,sign); 
        }
    };

    template <class V>
    static typename V::float_type InlineLogProdElements(
        const BaseVector_Calc<V>& v, typename V::zfloat_type* sign)
    {
        typedef typename V::const_cview_type Vv;
        Vv vv = v.cView();
        return LogProdElementsV_Helper<-3,V::_size,Vv>::call(vv,sign);
    }

    template <class V>
    static typename V::float_type LogProdElements(
        const BaseVector_Calc<V>& v, typename V::zfloat_type* sign)
    {
        typedef typename V::const_cview_type Vv;
        Vv vv = v.cView();
        return LogProdElementsV_Helper<-1,V::_size,Vv>::call(vv,sign);
    }


    // 
    // HasZeroElement
    //

    // algo 11: simple for loop
    template <int s, class V>
    struct HasZeroElementV_Helper<11,s,V> 
    {
        typedef typename V::value_type T;
        static bool call(const V& v)
        {
            const int n = s == UNKNOWN ? int(v.size()) : s;
#ifdef PRINTALGO_Det
            std::cout<<"HaxZeroElement algo 11: n,s = "<<n<<','<<s<<std::endl;
#endif
            for(int i=0;i<n;++i) {
                if (v.cref(i) != T(0)) continue;
                else return true;
            }
            return false;
        }
    };

    // algo 12: 2 at a time
    template <int s, class V>
    struct HasZeroElementV_Helper<12,s,V> 
    {
        typedef typename V::value_type T;
        static bool call(const V& v)
        {
            const int n = s == UNKNOWN ? int(v.size()) : s;
#ifdef PRINTALGO_Det
            std::cout<<"HaxZeroElement algo 12: n,s = "<<n<<','<<s<<std::endl;
#endif
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
    template <int s, class V>
    struct HasZeroElementV_Helper<13,s,V> 
    {
        typedef typename V::value_type T;
        static bool call(const V& v)
        {
            const int n = s == UNKNOWN ? int(v.size()) : s;
#ifdef PRINTALGO_Det
            std::cout<<"HaxZeroElement algo 13: n,s = "<<n<<','<<s<<std::endl;
#endif
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
    template <int s, class V>
    struct HasZeroElementV_Helper<14,s,V> 
    {
        typedef typename V::value_type T;
        static bool call(const V& v)
        {
            const int n = s == UNKNOWN ? int(v.size()) : s;
#ifdef PRINTALGO_Det
            std::cout<<"HaxZeroElement algo 14: n,s = "<<n<<','<<s<<std::endl;
#endif
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
    template <int s, class V>
    struct HasZeroElementV_Helper<15,s,V>
    {
        typedef typename V::value_type T;

        template <int I, int N>
        struct Unroller
        {
            static bool unroll(const V& v)
            {
                return (
                    Unroller<I,N/2>::unroll(v) ||
                    Unroller<I+N/2,N-N/2>::unroll(v));
            }
        };
        template <int I>
        struct Unroller<I,1>
        { static bool unroll(const V& v) { return v.cref(I) == T(0); } };
        template <int I>
        struct Unroller<I,0>
        { static bool unroll(const V& v) { return false; } };
        static bool call(const V& v)
        {
#ifdef PRINTALGO_Det
            std::cout<<"HaxZeroElement algo 15: n,s = "<<v.size()<<','<<s<<std::endl;
#endif
            return Unroller<0,s>::unroll(v); 
        }
    };

    // algo -3: Determine which algorithm to use
    template <int s, class V>
    struct HasZeroElementV_Helper<-3,s,V> 
    {
        typedef typename V::value_type T;
        static bool call(const V& v)
        {
#if TMV_OPT == 0
            const int algo = 11;
#else
            typedef typename V::real_type RT;
            const int algo = 
                (s != UNKNOWN && s <= int(128/sizeof(T))) ? 15 :
                (sizeof(RT) == 8) ? (V::iscomplex ? 12 : 13) :
                (sizeof(RT) == 4) ? (V::iscomplex ? 13 : 14) :
                11;
#endif
#ifdef PRINTALGO_Det
            std::cout<<"Inline HaxZeroElement\n";
            std::cout<<"v = "<<TMV_Text(v)<<std::endl;
            std::cout<<"s = "<<s<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
#endif
            return HasZeroElementV_Helper<algo,s,V>::call(v);
        }
    };

    // algo 97: Conjugate
    template <int s, class V>
    struct HasZeroElementV_Helper<97,s,V> 
    {
        typedef typename V::value_type T;
        static bool call(const V& v)
        { 
            typedef typename V::const_conjugate_type Vc;
            Vc vc = v.conjugate();
            return HasZeroElementV_Helper<-1,s,Vc>::call(vc);
        }
    };

    // algo 98: Call inst
    template <int s, class V>
    struct HasZeroElementV_Helper<98,s,V> 
    {
        typedef typename V::value_type T;
        static bool call(const V& v)
        { return InstHasZeroElement(v.xView()); }
    };

    // algo -1: Check for inst
    template <int s, class V>
    struct HasZeroElementV_Helper<-1,s,V> 
    {
        typedef typename V::value_type T;
        static bool call(const V& v)
        {
            const bool inst =
                (s == UNKNOWN || s > 16) &&
                Traits<T>::isinst;
            const int algo =
                V::_conj ? 97 : 
                inst ? 98 : 
                -3;
            return HasZeroElementV_Helper<algo,s,V>::call(v); 
        }
    };

    template <class V>
    static bool InlineHasZeroElement(const BaseVector_Calc<V>& v)
    {
        typedef typename V::const_cview_type Vv;
        Vv vv = v.cView();
        return HasZeroElementV_Helper<-3,V::_size,Vv>::call(vv);
    }

    template <class V>
    static bool HasZeroElement(const BaseVector_Calc<V>& v)
    {
        typedef typename V::const_cview_type Vv;
        Vv vv = v.cView();
        return HasZeroElementV_Helper<-1,V::_size,Vv>::call(vv);
    }


} // namespace tmv

#endif
