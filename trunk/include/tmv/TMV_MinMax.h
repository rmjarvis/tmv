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

#ifndef TMV_MinMax_H
#define TMV_MinMax_H

namespace tmv {

    // Defined in TMV_Vector.cpp
    template <class T>
    T InstMaxElement(const ConstVectorView<T>& v, int* imax);

    template <class T>
    typename ConstVectorView<T>::float_type InstMaxAbsElement(
        const ConstVectorView<T>& v, int* imax); 

    template <class T>
    typename Traits<T>::real_type InstMaxAbs2Element(
        const ConstVectorView<T>& v, int* imax);

    template <class T>
    T InstMinElement(const ConstVectorView<T>& v, int* imin);

    template <class T>
    typename ConstVectorView<T>::float_type InstMinAbsElement(
        const ConstVectorView<T>& v, int* imin);

    template <class T>
    typename Traits<T>::real_type InstMinAbs2Element(
        const ConstVectorView<T>& v, int* imin);

 

    // 
    // I combine all of the above into a single Helper class
    // that has the component and min/max as template arguments,
    // since the underlying algorithms are basically the same.
    //

    template <int algo, CompType comp, bool max, class V>
    struct MinMaxElement_Helper;

    // algo 0: size == 0
    template <CompType comp, bool max, class V>
    struct MinMaxElement_Helper<0,comp,max,V>
    {
        typedef typename Component<comp,typename V::value_type>::ret_type ret;
        static TMV_INLINE ret call(const V& , int* ibest)
        {
            if (ibest) *ibest = -1;
            return ret(0);
        }
    };

    // algo 11: simple for loop
    template <CompType comp, bool max, class V>
    struct MinMaxElement_Helper<11,comp,max,V>
    {
        typedef typename Component<comp,typename V::value_type>::ret_type ret;
        static ret call(const V& v, int* ibest)
        {
            typedef typename V::value_type VT;
            typedef typename V::const_iterator IT;
            typedef typename TypeSelect<comp==AbsComp ,
                    typename V::zfloat_type , typename V::value_type >::type ZT;
            ZT value;
            ret best;
            int n = V::_size == UNKNOWN ? v.size() : V::_size;
            if (n == 0) 
                return MinMaxElement_Helper<0,comp,max,V>::call(v,ibest);
            IT it = v.begin();
            best = Component<comp,VT>::f(*it);
            IT bestit = it;

            if (n) do {
                value = Traits<ZT>::convert(*it++);
                Component<comp,ZT>::applyf(value);

                // It turns out that doing the comparison this way makes 
                // a huge difference in speed.
                // The reason is basically that the processor anticipates
                // if statements to be true and then backtracks if it turned 
                // out to be false.
                // So you always want to test for the most likely choice
                // first.  In this case, it is more likely that value
                // is not the new best.
                if (Maybe<max>::less(TMV_REAL(value),TMV_REAL(best))) continue;
                else { best = Component<comp,ZT>::get(value); bestit = it-1; }
#if 0
                // This is the slower but more intuitive way:
                if (Maybe<max>::less(TMV_REAL(best),TMV_REAL(value)))
                { best = Component<comp,ZT>::get(value); bestit = it-1; }
#endif

            } while (--n);
            if (ibest) *ibest = bestit - v.begin();
            return best;
        }
    };

    // algo 12: 2 at a time
    template <CompType comp, bool max, class V>
    struct MinMaxElement_Helper<12,comp,max,V>
    {
        typedef typename Component<comp,typename V::value_type>::ret_type ret;
        static ret call(const V& v, int* ibest)
        {
            typedef typename V::value_type VT;
            typedef typename V::const_iterator IT;
            typedef typename TypeSelect<comp==AbsComp ,
                    typename V::zfloat_type , typename V::value_type >::type ZT;
            ret best;
            ZT value, value1;
            IT it = v.begin();
            int n = V::_size == UNKNOWN ? v.size() : V::_size;
            if (n == 0) 
                return MinMaxElement_Helper<0,comp,max,V>::call(v,ibest);
            best = Component<comp,VT>::f(*it);
            IT bestit = it;

            int n_2 = ((n-1)>>1);
            const int nb = n-1-(n_2<<1);

            if (nb) {
                value = Traits<ZT>::convert(it[1]);
                it += 2;
                Component<comp,ZT>::applyf(value);
                if (Maybe<max>::less(TMV_REAL(best),TMV_REAL(value))) 
                { best = Component<comp,ZT>::get(value); ++bestit; }
            }
            else ++it;

            if (n_2) do {
                value = Traits<ZT>::convert(it[0]);
                value1 = Traits<ZT>::convert(it[1]);
                it += 2;
                Component<comp,ZT>::applyf(value);
                Component<comp,ZT>::applyf(value1);

                if (Maybe<max>::less(TMV_REAL(value),TMV_REAL(best)) &&
                    Maybe<max>::less(TMV_REAL(value1),TMV_REAL(best))) continue;
                else if (Maybe<max>::less(TMV_REAL(value1),TMV_REAL(value)))
                { best = Component<comp,ZT>::get(value); bestit = it-2; }
                else { best = Component<comp,ZT>::get(value1); bestit = it-1; }
            } while (--n_2);
            if (ibest) *ibest = bestit - v.begin();
            return best;
        }
    };

    // algo 13: 3 at a time
    template <CompType comp, bool max, class V>
    struct MinMaxElement_Helper<13,comp,max,V>
    {
        typedef typename Component<comp,typename V::value_type>::ret_type ret;
        static ret call(const V& v, int* ibest)
        {
            typedef typename V::value_type VT;
            typedef typename V::const_iterator IT;
            typedef typename TypeSelect<comp==AbsComp ,
                    typename V::zfloat_type , typename V::value_type >::type ZT;
            ret best;
            ZT value, value1, value2;
            IT it = v.begin();
            int n = V::_size == UNKNOWN ? v.size() : V::_size;
            if (n == 0) 
                return MinMaxElement_Helper<0,comp,max,V>::call(v,ibest);
            const IT end = it + n;
            best = Component<comp,VT>::f(*it);
            IT bestit = it;

            int n_3 = n/3;

            if (n % 3 != 1) {
                value = Traits<ZT>::convert(it[1]);
                Component<comp,ZT>::applyf(value);
                if (Maybe<max>::less(TMV_REAL(best),TMV_REAL(value))) 
                { best = Component<comp,ZT>::get(value); ++bestit; }
                if (n % 3 == 0) {
                    value = Traits<ZT>::convert(it[2]);
                    Component<comp,ZT>::applyf(value);
                    if (Maybe<max>::less(TMV_REAL(best),TMV_REAL(value))) 
                    { best = Component<comp,ZT>::get(value); bestit = it + 2; }
                    it += 3;
                }
                else it += 2;
            }
            else ++it;

            if (n_3) do {
                value = Traits<ZT>::convert(it[0]);
                value1 = Traits<ZT>::convert(it[1]);
                value2 = Traits<ZT>::convert(it[2]);
                it += 3;
                Component<comp,ZT>::applyf(value);
                Component<comp,ZT>::applyf(value1);
                Component<comp,ZT>::applyf(value2);

                if (Maybe<max>::less(TMV_REAL(value),TMV_REAL(best)) &&
                    Maybe<max>::less(TMV_REAL(value1),TMV_REAL(best)) &&
                    Maybe<max>::less(TMV_REAL(value2),TMV_REAL(best))) continue;
                else if (
                    Maybe<max>::less(TMV_REAL(value1),TMV_REAL(value)) &&
                    Maybe<max>::less(TMV_REAL(value2),TMV_REAL(value)))
                { best = Component<comp,ZT>::get(value); bestit = it-3; }
                else if (Maybe<max>::less(TMV_REAL(value2),TMV_REAL(value1)))
                { best = Component<comp,ZT>::get(value1); bestit = it-2; }
                else { best = Component<comp,ZT>::get(value2); bestit = it-1; }
            } while (--n_3);
            if (ibest) *ibest = bestit - v.begin();
            return best;
        }
    };

    // algo 14: 4 at a time
    template <CompType comp, bool max, class V>
    struct MinMaxElement_Helper<14,comp,max,V>
    {
        typedef typename Component<comp,typename V::value_type>::ret_type ret;
        static ret call(const V& v, int* ibest)
        {
            typedef typename V::value_type VT;
            typedef typename V::const_iterator IT;
            typedef typename TypeSelect<comp==AbsComp ,
                    typename V::zfloat_type , typename V::value_type >::type ZT;
            ret best;
            ZT value, value1, value2, value3;
            IT it = v.begin();
            int n = V::_size == UNKNOWN ? v.size() : V::_size;
            if (n == 0) 
                return MinMaxElement_Helper<0,comp,max,V>::call(v,ibest);
            int n_4 = ((n-1)>>2);
            best = Component<comp,VT>::f(*it);
            IT bestit = it;

            if (n % 4 != 1) {
                value = Traits<ZT>::convert(it[1]);
                Component<comp,ZT>::applyf(value);
                if (Maybe<max>::less(TMV_REAL(best),TMV_REAL(value))) 
                { best = Component<comp,ZT>::get(value); ++bestit; }
                if (n % 4 != 2) {
                    value = Traits<ZT>::convert(it[2]);
                    Component<comp,ZT>::applyf(value);
                    if (Maybe<max>::less(TMV_REAL(best),TMV_REAL(value))) 
                    { best = Component<comp,ZT>::get(value); bestit = it + 2; }
                    if (n % 4 == 0) {
                        value = Traits<ZT>::convert(it[3]);
                        Component<comp,ZT>::applyf(value);
                        if (Maybe<max>::less(TMV_REAL(best),TMV_REAL(value))) 
                        { best = Component<comp,ZT>::get(value); bestit = it + 3; }
                        it += 4;
                    }
                    else it += 3;
                }
                else it += 2;
            }
            else ++it;

            if (n_4) do {
                value = Traits<ZT>::convert(it[0]);
                value1 = Traits<ZT>::convert(it[1]);
                value2 = Traits<ZT>::convert(it[2]);
                value3 = Traits<ZT>::convert(it[3]);
                Component<comp,ZT>::applyf(value);
                Component<comp,ZT>::applyf(value1);
                Component<comp,ZT>::applyf(value2);
                Component<comp,ZT>::applyf(value3);
                it += 4;

                if (Maybe<max>::less(TMV_REAL(value),TMV_REAL(best)) &&
                    Maybe<max>::less(TMV_REAL(value1),TMV_REAL(best)) &&
                    Maybe<max>::less(TMV_REAL(value2),TMV_REAL(best)) &&
                    Maybe<max>::less(TMV_REAL(value3),TMV_REAL(best))) continue;
                else if (
                    Maybe<max>::less(TMV_REAL(value1),TMV_REAL(value)) &&
                    Maybe<max>::less(TMV_REAL(value2),TMV_REAL(value)) &&
                    Maybe<max>::less(TMV_REAL(value3),TMV_REAL(value)))
                { best = Component<comp,ZT>::get(value); bestit = it-4; }
                else if (
                    Maybe<max>::less(TMV_REAL(value2),TMV_REAL(value1)) &&
                    Maybe<max>::less(TMV_REAL(value3),TMV_REAL(value1)))
                { best = Component<comp,ZT>::get(value1); bestit = it-3; }
                else if (Maybe<max>::less(TMV_REAL(value3),TMV_REAL(value2)))
                { best = Component<comp,ZT>::get(value2); bestit = it-2; }
                else { best = Component<comp,ZT>::get(value3); bestit = it-1; }
            } while (--n_4);
            if (ibest) *ibest = bestit - v.begin();
            return best;
        }
    };

    // algo 15: fully unroll
    template <CompType comp, bool max, class V>
    struct MinMaxElement_Helper<15,comp,max,V>
    {
        typedef typename Component<comp,typename V::value_type>::ret_type ret;
        typedef typename V::value_type VT;
        typedef typename V::const_iterator IT;

        template <int n, int x>
        struct Helper1
        {
            static inline IT unroll1(const IT& v0)
            {
                const IT temp1 = Helper1<n-1,x>::unroll1(v0);
                const IT temp2 = v0 + n-1;
                return Maybe<max>::less(
                    Component<comp,VT>::f(*temp1),
                    Component<comp,VT>::f(*temp2)) ? temp2 : temp1;
            }
            static inline ret unroll2(const V& v)
            {
                const ret temp1 = Helper1<n-1,x>::unroll2(v);
                const ret temp2 = Component<comp,VT>::f(v.cref(n-1));
                return Maybe<max>::less(TMV_REAL(temp1),TMV_REAL(temp2)) ?
                    temp2 : temp1;
            }
            static inline IT unroll3(const IT& v0, ret& best)
            {
                const IT temp1 = Helper1<n-1,x>::unroll3(v0,best);
                const IT temp2 = v0 + n-1;
                const ret f2 = Component<comp,VT>::f(*temp2);
                if (Maybe<max>::less(TMV_REAL(f2),TMV_REAL(best))) 
                    return temp1;
                else { best = f2; return temp2; }
            }
        };
        template <int x>
        struct Helper1<1,x>
        {
            static inline IT unroll1(const IT& v0)
            { return v0; }
            static inline ret unroll2(const V& v)
            { return Component<comp,VT>::f(v.cref(0)); }
            static inline IT unroll3(const IT& v0, ret& best)
            { best = Component<comp,VT>::f(*v0); return v0; }
        };
        template <int x>
        struct Helper1<0,x>; // not defined so give a compile error

        static inline ret call(const V& v, int* ibest)
        {
            if (ibest) {
                IT bestit = Helper1<V::_size,1>::unroll1(v.begin());
                *ibest = bestit - v.begin();
                return Component<comp,VT>::f(*bestit);
            } else {
                return Helper1<V::_size,1>::unroll2(v);
            }
        }
    };

    // algo 16: fully unroll with storage of best along the way
    template <CompType comp, bool max, class V>
    struct MinMaxElement_Helper<16,comp,max,V>
    {
        typedef typename Component<comp,typename V::value_type>::ret_type ret;
        typedef typename V::value_type VT;
        typedef typename V::const_iterator IT;

        template <int n, int x>
        struct Helper1
        {
            static inline void unroll(ret& best, IT& bestit, const IT& v0)
            {
                Helper1<n-1,x>::unroll(best,bestit,v0);
                const ret temp = Component<comp,VT>::f(v0[n-1]);
                if (Maybe<max>::less(TMV_REAL(temp),TMV_REAL(best))) return;
                else { best = temp; bestit = v0+n-1; }
            }
        };

        template <int x>
        struct Helper1<1,x>
        { static inline void unroll(ret& best, IT& bestit, const IT& v0) {} };

        template <int x>
        struct Helper1<0,x>; 

        static inline ret call(const V& v, int* ibest)
        {
            IT it = v.begin();
            ret best = Component<comp,VT>::f(*it);
            IT bestit = it;
            Helper1<V::_size,1>::unroll(best,bestit,it);
            if (ibest) *ibest = bestit - v.begin();
            return best;
        }
    };

    // algo 21: special algorithm for complex abs(), max
    // Try to avoid as many abs calls as possible by checking the components
    // to see if it's even possible for the abs value to be greater.
    template <class V>
    struct MinMaxElement_Helper<21,AbsComp,true,V>
    {
        typedef typename V::float_type ret;
        static ret call(const V& v, int* ibest)
        {
            typedef typename V::value_type VT;
            typedef typename V::const_iterator IT;
            typedef typename V::zfloat_type ZT;
            ZT value;
            ret best;
            int n = V::_size == UNKNOWN ? v.size() : V::_size;
            if (n == 0) 
                return MinMaxElement_Helper<0,AbsComp,true,V>::call(v,ibest);
            IT it = v.begin();
            best = TMV_ABS(*it);
            IT bestit = it;

            if (n) do {
                value = Traits<ZT>::convert(*it++);
                // Only if sum of abs(real) + abs(imag) > best
                // should we even both doing the complex abs()
                if (TMV_ABS(real(value)) + TMV_ABS(imag(value)) < best) 
                    continue;
                else {
                    Component<AbsComp,ZT>::applyf(value);
                    if (TMV_REAL(value) < best) continue;
                    else { best = TMV_REAL(value); bestit = it-1; }
                }
            } while (--n);
            if (ibest) *ibest = bestit - v.begin();
            return best;
        }
    };
    // algo 21: min version uses a similar principle
    template <class V>
    struct MinMaxElement_Helper<21,AbsComp,false,V>
    {
        typedef typename V::float_type ret;
        static ret call(const V& v, int* ibest)
        {
            typedef typename V::value_type VT;
            typedef typename V::const_iterator IT;
            typedef typename V::zfloat_type ZT;
            ZT value;
            ret best;
            int n = V::_size == UNKNOWN ? v.size() : V::_size;
            if (n == 0) 
                return MinMaxElement_Helper<0,AbsComp,false,V>::call(v,ibest);
            IT it = v.begin();
            best = TMV_ABS(*it);
            IT bestit = it;

            if (n) do {
                value = Traits<ZT>::convert(*it++);
                // Only if both abs(real) and abs(imag) < best
                // should we even both doing the complex abs()
                if (TMV_ABS(real(value)) >= best && TMV_ABS(imag(value)) >= best) 
                    continue;
                else {
                    Component<AbsComp,ZT>::applyf(value);
                    if (TMV_REAL(value) > best) continue;
                    else { best = TMV_REAL(value); bestit = it-1; }
                }
            } while (--n);
            if (ibest) *ibest = bestit - v.begin();
            return best;
        }
    };

    // algo 24: 4 at a time, but compared in two batches
    // Based on ATLAS algorithm iamax_abs2p24_x1
    template <CompType comp, bool max, class V>
    struct MinMaxElement_Helper<24,comp,max,V>
    {
        typedef typename Component<comp,typename V::value_type>::ret_type ret;
        static ret call(const V& v, int* ibest)
        {
            typedef typename V::value_type VT;
            typedef typename V::const_iterator IT;
            typedef typename TypeSelect<comp==AbsComp ,
                    typename V::zfloat_type , typename V::value_type >::type ZT;
            ret best;
            ZT value, value1;
            IT it = v.begin();
            int n = V::_size == UNKNOWN ? v.size() : V::_size;
            if (n == 0) 
                return MinMaxElement_Helper<0,comp,max,V>::call(v,ibest);
            int n_4 = ((n-1)>>2);
            best = Component<comp,VT>::f(*it);
            IT bestit = it;

            if (n % 4 != 1) {
                value = Traits<ZT>::convert(it[1]);
                Component<comp,ZT>::applyf(value);
                if (Maybe<max>::less(TMV_REAL(best),TMV_REAL(value))) 
                { best = Component<comp,ZT>::get(value); ++bestit; }
                if (n % 4 != 2) {
                    value = Traits<ZT>::convert(it[2]);
                    Component<comp,ZT>::applyf(value);
                    if (Maybe<max>::less(TMV_REAL(best),TMV_REAL(value))) 
                    { best = Component<comp,ZT>::get(value); bestit = it + 2; }
                    if (n % 4 == 0) {
                        value = Traits<ZT>::convert(it[3]);
                        Component<comp,ZT>::applyf(value);
                        if (Maybe<max>::less(TMV_REAL(best),TMV_REAL(value))) 
                        { best = Component<comp,ZT>::get(value); bestit = it + 3; }
                        it += 4;
                    }
                    else it += 3;
                }
                else it += 2;
            }
            else ++it;

            if (n_4) do {
                value = Traits<ZT>::convert(it[0]);
                value1 = Traits<ZT>::convert(it[1]);
                Component<comp,ZT>::applyf(value);
                Component<comp,ZT>::applyf(value1);

                if (Maybe<max>::less(TMV_REAL(value),TMV_REAL(best)) &&
                    Maybe<max>::less(TMV_REAL(value1),TMV_REAL(best))) goto L1;
                else if (Maybe<max>::less(TMV_REAL(value1),TMV_REAL(value)))
                { best = Component<comp,ZT>::get(value); bestit = it; }
                else { best = Component<comp,ZT>::get(value1); bestit = it+1; }

L1:
                value = Traits<ZT>::convert(it[2]);
                value1 = Traits<ZT>::convert(it[3]);
                Component<comp,ZT>::applyf(value);
                Component<comp,ZT>::applyf(value1);

                if (Maybe<max>::less(TMV_REAL(value),TMV_REAL(best)) &&
                    Maybe<max>::less(TMV_REAL(value1),TMV_REAL(best))) goto L2;
                else if (Maybe<max>::less(TMV_REAL(value1),TMV_REAL(value)))
                { best = Component<comp,ZT>::get(value); bestit = it+2; }
                else { best = Component<comp,ZT>::get(value1); bestit = it+3; }
L2:
                it += 4;
            } while (--n_4);
            if (ibest) *ibest = bestit - v.begin();
            return best;
        }
    };

    // algo 32: 11 or 12 depending on size
    template <CompType comp, bool max, class V>
    struct MinMaxElement_Helper<32,comp,max,V>
    {
        typedef typename Component<comp,typename V::value_type>::ret_type ret;
        static ret call(const V& v, int* ibest)
        {
#if TMV_OPT >= 2
            if (!ibest && v.size() < 250) 
                return MinMaxElement_Helper<11,comp,max,V>::call(v,ibest);
            else 
#endif
                return MinMaxElement_Helper<12,comp,max,V>::call(v,ibest);
        }
    };

    // algo 33: 11 or 13 depending on size
    template <CompType comp, bool max, class V>
    struct MinMaxElement_Helper<33,comp,max,V>
    {
        typedef typename Component<comp,typename V::value_type>::ret_type ret;
        static ret call(const V& v, int* ibest)
        {
#if TMV_OPT >= 2
            if (!ibest && v.size() < 250) 
                return MinMaxElement_Helper<11,comp,max,V>::call(v,ibest);
            else 
#endif
                return MinMaxElement_Helper<13,comp,max,V>::call(v,ibest);
        }
    };

    // algo 34: 11 or 14 depending on size
    template <CompType comp, bool max, class V>
    struct MinMaxElement_Helper<34,comp,max,V>
    {
        typedef typename Component<comp,typename V::value_type>::ret_type ret;
        static ret call(const V& v, int* ibest)
        {
#if TMV_OPT >= 2
            if (!ibest && v.size() < 250) 
                return MinMaxElement_Helper<11,comp,max,V>::call(v,ibest);
            else 
#endif
                return MinMaxElement_Helper<14,comp,max,V>::call(v,ibest);
        }
    };

    // algo 35: 11 or 24 depending on size
    template <CompType comp, bool max, class V>
    struct MinMaxElement_Helper<35,comp,max,V>
    {
        typedef typename Component<comp,typename V::value_type>::ret_type ret;
        static ret call(const V& v, int* ibest)
        {
#if TMV_OPT >= 2
            if (!ibest && v.size() < 250) 
                return MinMaxElement_Helper<11,comp,max,V>::call(v,ibest);
            else 
#endif
                return MinMaxElement_Helper<24,comp,max,V>::call(v,ibest);
        }
    };

    // algo 90: Call inst
    template <class V>
    struct MinMaxElement_Helper<90,ValueComp,true,V>
    {
        typedef typename V::value_type ret;
        static TMV_INLINE ret call(const V& v, int* ibest)
        { return InstMaxElement(v.xView(),ibest); }
    };
    template <class V>
    struct MinMaxElement_Helper<90,AbsComp,true,V>
    {
        typedef typename V::float_type ret;
        static TMV_INLINE ret call(const V& v, int* ibest)
        { return InstMaxAbsElement(v.xView(),ibest); }
    };
    template <class V>
    struct MinMaxElement_Helper<90,Abs2Comp,true,V>
    {
        typedef typename V::real_type ret;
        static TMV_INLINE ret call(const V& v, int* ibest)
        { return InstMaxAbs2Element(v.xView(),ibest); }
    };
    template <class V>
    struct MinMaxElement_Helper<90,ValueComp,false,V>
    {
        typedef typename V::value_type ret;
        static TMV_INLINE ret call(const V& v, int* ibest)
        { return InstMinElement(v.xView(),ibest); }
    };
    template <class V>
    struct MinMaxElement_Helper<90,AbsComp,false,V>
    {
        typedef typename V::float_type ret;
        static TMV_INLINE ret call(const V& v, int* ibest)
        { return InstMinAbsElement(v.xView(),ibest); }
    };
    template <class V>
    struct MinMaxElement_Helper<90,Abs2Comp,false,V>
    {
        typedef typename V::real_type ret;
        static TMV_INLINE ret call(const V& v, int* ibest)
        { return InstMinAbs2Element(v.xView(),ibest); }
    };
 
    // algo 97: Conjugate
    template <CompType comp, bool max, class V>
    struct MinMaxElement_Helper<97,comp,max,V>
    {
        typedef typename Component<comp,typename V::value_type>::ret_type ret;
        static TMV_INLINE ret call(const V& v, int* ibest)
        {
            typedef typename V::const_nonconj_type Vnc;
            Vnc vnc = v.nonConj();
            return MinMaxElement_Helper<-2,comp,max,Vnc>::call(vnc,ibest);
        }
    };
    template <bool max, class V>
    struct MinMaxElement_Helper<97,ValueComp,max,V>
    {
        typedef typename V::value_type ret;
        static TMV_INLINE ret call(const V& v, int* ibest)
        {
            typedef typename V::const_conjugate_type Vc;
            Vc vc = v.conjugate();
            return TMV_CONJ(
                MinMaxElement_Helper<-2,ValueComp,max,Vc>::call(vc));
        }
    };
 
    // algo -3: Determine which algorithm to use
    template <class V>
    struct MinMaxElement_Helper<-3,ValueComp,true,V>
    {
        typedef typename V::value_type ret;
        static TMV_INLINE ret call(const V& v, int* ibest)
        {
#if TMV_OPT >= 2
            const int algo1 = 
                V::_size == 0 ? 0 :
                V::_size != UNKNOWN ? (
                    V::_size <= 160 ? 15 :
                    (V::_step == 1 && V::_size > 250) ? (
                        ( V::iscomplex ? 12 : 24 ) ) :
                    11 ) :
                V::_step == 1 ? ( V::iscomplex ? 32 : 34 ) :
                11;
#endif
            const int algo2 = 
                V::_size == 0 ? 0 :
                TMV_OPT == 0 ? 11 :
                V::_size != UNKNOWN ? (
                    V::_size <= 12 ? 15 :
                    V::_size <= 40 ? 16 :
                    (V::_step == 1 && V::_size > 250) ? (
                        ( V::iscomplex ? 11 : 24 ) ) :
                    11 ) :
                V::_step == 1 ? ( V::iscomplex ? 11 : 34 ) :
                11;
#ifdef PRINTALGO
            std::cout<<"InlineMaxElement: algo = "<<algo1<<" "<<algo2<<std::endl;
#endif

#if TMV_OPT >= 2
            if (!ibest) 
                return MinMaxElement_Helper<algo1,ValueComp,true,V>::call(
                    v.vec(),ibest);
            else 
#endif
                return MinMaxElement_Helper<algo2,ValueComp,true,V>::call(
                    v.vec(),ibest);
        }
    };
    template <class V>
    struct MinMaxElement_Helper<-3,AbsComp,true,V>
    {
        typedef typename V::float_type ret;
        static TMV_INLINE ret call(const V& v, int* ibest)
        {
#if TMV_OPT >= 2
            const int algo1 = 
                V::_size == 0 ? 0 :
                V::_size != UNKNOWN ? (
                    V::iscomplex ? (
                        ( V::_size <= 16 ? 15 : V::_step == 1 ? 21 : 11 ) ) :
                    V::_size <= 160 ? 15 :
                    (V::_step == 1 && V::_size > 250) ?  24 : 11 ) :
                V::_step == 1 ? ( V::iscomplex ? 21 : 35 ) :
                11;
#endif
            const int algo2 = 
                V::_size == 0 ? 0 :
                TMV_OPT == 0 ? 11 :
                V::_size != UNKNOWN ? (
                    V::iscomplex ? (
                        ( V::_step == 1 ? 21 : 11 ) ) :
                    V::_size <= 12 ? 15 :
                    V::_size <= 40 ? 16 :
                    (V::_step == 1 && V::_size > 250) ? 24 : 11 ) :
                V::_step == 1 ? ( V::iscomplex ? 21 : 35 ) :
                11;
#ifdef PRINTALGO
            std::cout<<"InlineMaxAbsElement: algo = "<<algo1<<
                " "<<algo2<<std::endl;
#endif

#if TMV_OPT >= 2
            if (!ibest) 
                return MinMaxElement_Helper<algo1,AbsComp,true,V>::call(
                    v.vec(),ibest);
            else 
#endif
                return MinMaxElement_Helper<algo2,AbsComp,true,V>::call(
                    v.vec(),ibest);
        }
    };
    template <class V>
    struct MinMaxElement_Helper<-3,Abs2Comp,true,V>
    {
        typedef typename V::real_type ret;
        static TMV_INLINE ret call(const V& v, int* ibest)
        {
#if TMV_OPT >= 2
            const int algo1 = 
                V::_size == 0 ? 0 :
                V::_size != UNKNOWN ? (
                    V::_size <= 160 ? 15 :
                    (V::_step == 1 && V::_size > 250) ? (
                        ( V::iscomplex ? 12 : 24 ) ) :
                    11 ) :
                V::_step == 1 ? ( V::iscomplex ? 32 : 35 ) :
                11;
#endif
            const int algo2 = 
                V::_size == 0 ? 0 :
                TMV_OPT == 0 ? 11 :
                V::_size != UNKNOWN ? (
                    V::_size <= 12 ? 15 :
                    V::_size <= 40 ? 16 :
                    (V::_step == 1 && V::_size > 250) ? (
                        ( V::iscomplex ? 11 : 24 ) ) :
                    11 ) :
                V::_step == 1 ? ( V::iscomplex ? 11 : 35 ) :
                11;
#ifdef PRINTALGO
            std::cout<<"InlineMaxAbs2Element: algo = "<<algo1<<
                " "<<algo2<<std::endl;
#endif

#if TMV_OPT >= 2
            if (!ibest) 
                return MinMaxElement_Helper<algo1,Abs2Comp,true,V>::call(
                    v.vec(),ibest);
            else 
#endif
                return MinMaxElement_Helper<algo2,Abs2Comp,true,V>::call(
                    v.vec(),ibest);
        }
    };
    template <class V>
    struct MinMaxElement_Helper<-3,ValueComp,false,V>
    {
        typedef typename V::value_type ret;
        static TMV_INLINE ret call(const V& v, int* ibest)
        {
#if TMV_OPT >= 2
            const int algo1 = 
                V::_size == 0 ? 0 :
                V::_size != UNKNOWN ? (
                    V::_size <= 160 ? 15 :
                    (V::_step == 1 && V::_size > 250) ? (
                        ( V::iscomplex ? 12 : 24 ) ) :
                    11 ) :
                V::_step == 1 ? ( V::iscomplex ? 32 : 34 ) :
                11;
#endif
            const int algo2 = 
                V::_size == 0 ? 0 :
                TMV_OPT == 0 ? 11 :
                V::_size != UNKNOWN ? (
                    V::_size <= 12 ? 15 :
                    V::_size <= 40 ? 16 :
                    (V::_step == 1 && V::_size > 250) ? (
                        ( V::iscomplex ? 11 : 24 ) ) :
                    11 ) :
                V::_step == 1 ? ( V::iscomplex ? 11 : 34 ) :
                11;
#ifdef PRINTALGO
            std::cout<<"InlineMinElement: algo = "<<algo1<<" "<<algo2<<std::endl;
#endif

#if TMV_OPT >= 2
            if (!ibest) 
                return MinMaxElement_Helper<algo1,ValueComp,false,V>::call(
                    v.vec(),ibest);
            else 
#endif
                return MinMaxElement_Helper<algo2,ValueComp,false,V>::call(
                    v.vec(),ibest);
        }
    };
    template <class V>
    struct MinMaxElement_Helper<-3,AbsComp,false,V>
    {
        typedef typename V::float_type ret;
        static TMV_INLINE ret call(const V& v, int* ibest)
        
        {
#if TMV_OPT >= 2
            const int algo1 = 
                V::_size == 0 ? 0 :
                V::_size != UNKNOWN ? (
                    V::iscomplex ? (
                        ( V::_size <= 40 ? 15 : V::_step == 1 ? 21 : 11 ) ) :
                    V::_size <= 160 ? 15 :
                    (V::_step == 1 && V::_size > 250) ?  24 : 11 ) :
                V::_step == 1 ? ( V::iscomplex ? 21 : 35 ) :
                11;
#endif
            const int algo2 = 
                V::_size == 0 ? 0 :
                TMV_OPT == 0 ? 11 :
                V::_size != UNKNOWN ? (
                    V::iscomplex ? (
                        ( V::_step == 1 ? 
                          ( V::_size <= 20 ? 12 : 21 ) : 
                          11 ) ) :
                    V::_size <= 12 ? 15 :
                    V::_size <= 40 ? 16 :
                    (V::_step == 1 && V::_size > 250) ? 24 : 11 ) :
                V::_step == 1 ? ( V::iscomplex ? 21 : 35 ) :
                11;
#ifdef PRINTALGO
            std::cout<<"InlineMinAbsElement: algo = "<<algo1<<
                " "<<algo2<<std::endl;
#endif

#if TMV_OPT >= 2
            if (!ibest) 
                return MinMaxElement_Helper<algo1,AbsComp,false,V>::call(
                    v.vec(),ibest);
            else 
#endif
                return MinMaxElement_Helper<algo2,AbsComp,false,V>::call(
                    v.vec(),ibest);
        }
    };
    template <class V>
    struct MinMaxElement_Helper<-3,Abs2Comp,false,V>
    {
        typedef typename V::real_type ret;
        static TMV_INLINE ret call(const V& v, int* ibest)
        {
#if TMV_OPT >= 2
            const int algo1 = 
                V::_size == 0 ? 0 :
                V::_size != UNKNOWN ? (
                    V::_size <= 160 ? 15 :
                    (V::_step == 1 && V::_size > 250) ? (
                        ( V::iscomplex ? 12 : 24 ) ) :
                    11 ) :
                V::_step == 1 ? ( V::iscomplex ? 32 : 35 ) :
                11;
#endif
            const int algo2 = 
                V::_size == 0 ? 0 :
                TMV_OPT == 0 ? 11 :
                V::_size != UNKNOWN ? (
                    V::_size <= 12 ? 15 :
                    V::_size <= 40 ? 16 :
                    (V::_step == 1 && V::_size > 250) ? (
                        ( V::iscomplex ? 11 : 24 ) ) :
                    11 ) :
                V::_step == 1 ? ( V::iscomplex ? 11 : 35 ) :
                11;
#ifdef PRINTALGO
            std::cout<<"InlineMinAbs2Element: algo = "<<algo1<<
                " "<<algo2<<std::endl;
#endif

#if TMV_OPT >= 2
            if (!ibest) 
                return MinMaxElement_Helper<algo1,Abs2Comp,false,V>::call(
                    v.vec(),ibest);
            else 
#endif
                return MinMaxElement_Helper<algo2,Abs2Comp,false,V>::call(
                    v.vec(),ibest);
        }
    };

    // algo -2: Check for inst
    template <CompType comp, bool max, class V>
    struct MinMaxElement_Helper<-2,comp,max,V>
    {
        typedef typename Component<comp,typename V::value_type>::ret_type ret;
        static TMV_INLINE ret call(const V& v, int* ibest)
        {
            typedef typename V::value_type VT;
            const bool inst = 
                (V::_size == UNKNOWN || V::_size > 16) &&
                Traits<VT>::isinst;
            const int algo =
                V::_conj ? 97 :
                inst ? 90 :
                -3;
            return MinMaxElement_Helper<algo,comp,max,V>::call(v,ibest);
        }
    };

    template <CompType comp, bool max, class V>
    struct MinMaxElement_Helper<-1,comp,max,V>
    {
        typedef typename Component<comp,typename V::value_type>::ret_type ret;
        static TMV_INLINE ret call(const V& v, int* ibest)
        { return MinMaxElement_Helper<-2,comp,max,V>::call(v,ibest); }
    };

    template <class V>
    static inline typename V::value_type DoMaxElement(
        const BaseVector_Calc<V>& v, int* imax)
    {
        typedef typename V::const_cview_type Vv;
        TMV_MAYBE_CREF(V,Vv) vv = v.cView();
        return MinMaxElement_Helper<-2,ValueComp,true,Vv>::call(vv,imax);
    }

    template <class V>
    static inline typename V::float_type DoMaxAbsElement(
        const BaseVector_Calc<V>& v, int* imax)
    {
        typedef typename V::const_cview_type Vv;
        TMV_MAYBE_CREF(V,Vv) vv = v.cView();
        return MinMaxElement_Helper<-2,AbsComp,true,Vv>::call(vv,imax);
    }

    template <class V>
    static inline typename V::real_type DoMaxAbs2Element(
        const BaseVector_Calc<V>& v, int* imax)
    {
        typedef typename V::const_cview_type Vv;
        TMV_MAYBE_CREF(V,Vv) vv = v.cView();
        return MinMaxElement_Helper<-2,Abs2Comp,true,Vv>::call(vv,imax);
    }

    template <class V>
    static inline typename V::value_type DoMinElement(
        const BaseVector_Calc<V>& v, int* imin)
    {
        typedef typename V::const_cview_type Vv;
        TMV_MAYBE_CREF(V,Vv) vv = v.cView();
        return MinMaxElement_Helper<-2,ValueComp,false,Vv>::call(vv,imin);
    }

    template <class V>
    static inline typename V::float_type DoMinAbsElement(
        const BaseVector_Calc<V>& v, int* imin)
    {
        typedef typename V::const_cview_type Vv;
        TMV_MAYBE_CREF(V,Vv) vv = v.cView();
        return MinMaxElement_Helper<-2,AbsComp,false,Vv>::call(vv,imin);
    }

    template <class V>
    static inline typename V::real_type DoMinAbs2Element(
        const BaseVector_Calc<V>& v, int* imin)
    {
        typedef typename V::const_cview_type Vv;
        TMV_MAYBE_CREF(V,Vv) vv = v.cView();
        return MinMaxElement_Helper<-2,Abs2Comp,false,Vv>::call(vv,imin);
    }

    template <class V>
    static inline typename V::value_type InlineMaxElement(
        const BaseVector_Calc<V>& v, int* imax)
    {
        typedef typename V::const_cview_type Vv;
        TMV_MAYBE_CREF(V,Vv) vv = v.cView();
        return MinMaxElement_Helper<-3,ValueComp,true,Vv>::call(vv,imax);
    }

    template <class V>
    static inline typename V::float_type InlineMaxAbsElement(
        const BaseVector_Calc<V>& v, int* imax)
    {
        typedef typename V::const_cview_type Vv;
        TMV_MAYBE_CREF(V,Vv) vv = v.cView();
        return MinMaxElement_Helper<-3,AbsComp,true,Vv>::call(vv,imax);
    }

    template <class V>
    static inline typename V::real_type InlineMaxAbs2Element(
        const BaseVector_Calc<V>& v, int* imax)
    {
        typedef typename V::const_cview_type Vv;
        TMV_MAYBE_CREF(V,Vv) vv = v.cView();
        return MinMaxElement_Helper<-3,Abs2Comp,true,Vv>::call(vv,imax);
    }

    template <class V>
    static inline typename V::value_type InlineMinElement(
        const BaseVector_Calc<V>& v, int* imin)
    {
        typedef typename V::const_cview_type Vv;
        TMV_MAYBE_CREF(V,Vv) vv = v.cView();
        return MinMaxElement_Helper<-3,ValueComp,false,Vv>::call(vv,imin);
    }

    template <class V>
    static inline typename V::float_type InlineMinAbsElement(
        const BaseVector_Calc<V>& v, int* imin)
    {
        typedef typename V::const_cview_type Vv;
        TMV_MAYBE_CREF(V,Vv) vv = v.cView();
        return MinMaxElement_Helper<-3,AbsComp,false,Vv>::call(vv,imin);
    }

    template <class V>
    static inline typename V::real_type InlineMinAbs2Element(
        const BaseVector_Calc<V>& v, int* imin)
    {
        typedef typename V::const_cview_type Vv;
        TMV_MAYBE_CREF(V,Vv) vv = v.cView();
        return MinMaxElement_Helper<-3,Abs2Comp,false,Vv>::call(vv,imin);
    }


} // namespace tmv

#endif
