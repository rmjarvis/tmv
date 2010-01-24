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

    // 
    // MaxElement
    //

    template <int algo, CompType comp, bool max, class V>
    struct MinMaxElement_Helper;

    // algo 1: simple for loop
    template <CompType comp, bool max, class V>
    struct MinMaxElement_Helper<1,comp,max,V>
    {
        typedef typename V::value_type VT;
        typedef typename V::real_type RT;
        typedef typename V::const_iterator IT;

        static inline RT call(const V& v, int*const ibest)
        {
            VT value;
            RT best;
            int n = v.size();
            IT it = v.begin();
            best = Component<comp,VT>::f(*it);
            IT bestit = it;

            if (n) do {
                value = *it++;
                Component<comp,VT>::applyf(value);

                // It turns out that doing the comparison this way makes 
                // a huge difference in speed.
                // The reason is basically that the processor anticipates
                // if statements to be true and then backtracks if it turned 
                // out to be false.
                // So you always want to test for the most likely choice
                // first.  In this case, it is more likely that value
                // is not the new best.
#if 1
                if (Maybe<max>::less(TMV_REAL(value),best)) continue;
                else { best = TMV_REAL(value); bestit = it-1; }
#else
                // This is the slower but more intuitive way:
                if (Maybe<max>::less(best,TMV_REAL(value)))
                { best = TMV_REAL(value); bestit = it-1; }
#endif

            } while (--n);
            if (ibest) *ibest = bestit - v.begin();
            return best;
        }
    };

    // algo 2: 2 at a time
    template <CompType comp, bool max, class V>
    struct MinMaxElement_Helper<2,comp,max,V>
    {
        typedef typename V::value_type VT;
        typedef typename V::real_type RT;
        typedef typename V::const_iterator IT;

        static inline RT call(const V& v, int*const ibest)
        {
            RT best;
            VT value, value1;
            IT it = v.begin();
            const int n = v.size();
            best = Component<comp,VT>::f(*it);
            IT bestit = it;

            int n_2 = ((n-1)>>1);
            const int nb = n-1-(n_2<<1);

            if (nb) {
                value = it[1];
                it += 2;
                Component<comp,VT>::applyf(value);
                if (Maybe<max>::less(best,TMV_REAL(value))) 
                { best = TMV_REAL(value); ++bestit; }
            }
            else ++it;

            if (n_2) do {
                value = it[0];
                value1 = it[1];
                it += 2;
                Component<comp,VT>::applyf(value);
                Component<comp,VT>::applyf(value1);

                if (Maybe<max>::less(TMV_REAL(value),best) &&
                    Maybe<max>::less(TMV_REAL(value1),best)) continue;
                else if (Maybe<max>::less(TMV_REAL(value1),TMV_REAL(value)))
                { best = TMV_REAL(value); bestit = it-2; }
                else { best = TMV_REAL(value1); bestit = it-1; }
            } while (--n_2);
            if (ibest) *ibest = bestit - v.begin();
            return best;
        }
    };

    // algo 3: 3 at a time
    template <CompType comp, bool max, class V>
    struct MinMaxElement_Helper<3,comp,max,V>
    {
        typedef typename V::value_type VT;
        typedef typename V::real_type RT;
        typedef typename V::const_iterator IT;

        static inline RT call(const V& v, int*const ibest)
        {
            RT best;
            VT value, value1, value2;
            IT it = v.begin();
            const int n = v.size();
            const IT end = it + n;
            best = Component<comp,VT>::f(*it);
            IT bestit = it;

            int n_3 = n/3;

            if (n % 3 != 1) {
                value = it[1];
                Component<comp,VT>::applyf(value);
                if (Maybe<max>::less(best,TMV_REAL(value))) 
                { best = TMV_REAL(value); ++bestit; }
                if (n % 3 == 0) {
                    value = it[2];
                    Component<comp,VT>::applyf(value);
                    if (Maybe<max>::less(best,TMV_REAL(value))) 
                    { best = TMV_REAL(value); bestit = it + 2; }
                    it += 3;
                }
                else it += 2;
            }
            else ++it;

            if (n_3) do {
                value = it[0];
                value1 = it[1];
                value2 = it[2];
                it += 3;
                Component<comp,VT>::applyf(value);
                Component<comp,VT>::applyf(value1);
                Component<comp,VT>::applyf(value2);

                if (Maybe<max>::less(TMV_REAL(value),best) &&
                    Maybe<max>::less(TMV_REAL(value1),best) &&
                    Maybe<max>::less(TMV_REAL(value2),best)) continue;
                else if (
                    Maybe<max>::less(TMV_REAL(value1),TMV_REAL(value)) &&
                    Maybe<max>::less(TMV_REAL(value2),TMV_REAL(value)))
                { best = TMV_REAL(value); bestit = it-3; }
                else if (Maybe<max>::less(TMV_REAL(value2),TMV_REAL(value1)))
                { best = TMV_REAL(value1); bestit = it-2; }
                else { best = TMV_REAL(value2); bestit = it-1; }
            } while (--n_3);
            if (ibest) *ibest = bestit - v.begin();
            return best;
        }
    };

    // algo 4: 4 at a time
    template <CompType comp, bool max, class V>
    struct MinMaxElement_Helper<4,comp,max,V>
    {
        typedef typename V::value_type VT;
        typedef typename V::real_type RT;
        typedef typename V::const_iterator IT;

        static inline RT call(const V& v, int*const ibest)
        {
            RT best;
            VT value, value1, value2, value3;
            IT it = v.begin();
            const int n = v.size();
            int n_4 = ((n-1)>>2);
            best = Component<comp,VT>::f(*it);
            IT bestit = it;

            if (n % 4 != 1) {
                value = it[1];
                Component<comp,VT>::applyf(value);
                if (Maybe<max>::less(best,TMV_REAL(value))) 
                { best = TMV_REAL(value); ++bestit; }
                if (n % 4 != 2) {
                    value = it[2];
                    Component<comp,VT>::applyf(value);
                    if (Maybe<max>::less(best,TMV_REAL(value))) 
                    { best = TMV_REAL(value); bestit = it + 2; }
                    if (n % 4 == 0) {
                        value = it[3];
                        Component<comp,VT>::applyf(value);
                        if (Maybe<max>::less(best,TMV_REAL(value))) 
                        { best = TMV_REAL(value); bestit = it + 3; }
                        it += 4;
                    }
                    else it += 3;
                }
                else it += 2;
            }
            else ++it;

            if (n_4) do {
                value = it[0];
                value1 = it[1];
                value2 = it[2];
                value3 = it[3];
                Component<comp,VT>::applyf(value);
                Component<comp,VT>::applyf(value1);
                Component<comp,VT>::applyf(value2);
                Component<comp,VT>::applyf(value3);
                it += 4;

                if (Maybe<max>::less(TMV_REAL(value),best) &&
                    Maybe<max>::less(TMV_REAL(value1),best) &&
                    Maybe<max>::less(TMV_REAL(value2),best) &&
                    Maybe<max>::less(TMV_REAL(value3),best)) continue;
                else if (
                    Maybe<max>::less(TMV_REAL(value1),TMV_REAL(value)) &&
                    Maybe<max>::less(TMV_REAL(value2),TMV_REAL(value)) &&
                    Maybe<max>::less(TMV_REAL(value3),TMV_REAL(value)))
                { best = TMV_REAL(value); bestit = it-4; }
                else if (
                    Maybe<max>::less(TMV_REAL(value2),TMV_REAL(value1)) &&
                    Maybe<max>::less(TMV_REAL(value3),TMV_REAL(value1)))
                { best = TMV_REAL(value1); bestit = it-3; }
                else if (Maybe<max>::less(TMV_REAL(value3),TMV_REAL(value2)))
                { best = TMV_REAL(value2); bestit = it-2; }
                else { best = TMV_REAL(value3); bestit = it-1; }
            } while (--n_4);
            if (ibest) *ibest = bestit - v.begin();
            return best;
        }
    };

    // algo 5: 4 at a time, but compared in two batches
    // Based on ATLAS algorithm iamax_abs2p24_x1
    template <CompType comp, bool max, class V>
    struct MinMaxElement_Helper<5,comp,max,V>
    {
        typedef typename V::value_type VT;
        typedef typename V::real_type RT;
        typedef typename V::const_iterator IT;

        static inline RT call(const V& v, int*const ibest)
        {
            RT best;
            VT value, value1;
            IT it = v.begin();
            const int n = v.size();
            int n_4 = ((n-1)>>2);
            best = Component<comp,VT>::f(*it);
            IT bestit = it;

            if (n % 4 != 1) {
                value = it[1];
                Component<comp,VT>::applyf(value);
                if (Maybe<max>::less(best,TMV_REAL(value))) 
                { best = TMV_REAL(value); ++bestit; }
                if (n % 4 != 2) {
                    value = it[2];
                    Component<comp,VT>::applyf(value);
                    if (Maybe<max>::less(best,TMV_REAL(value))) 
                    { best = TMV_REAL(value); bestit = it + 2; }
                    if (n % 4 == 0) {
                        value = it[3];
                        Component<comp,VT>::applyf(value);
                        if (Maybe<max>::less(best,TMV_REAL(value))) 
                        { best = TMV_REAL(value); bestit = it + 3; }
                        it += 4;
                    }
                    else it += 3;
                }
                else it += 2;
            }
            else ++it;

            if (n_4) do {
                value = it[0];
                value1 = it[1];
                Component<comp,VT>::applyf(value);
                Component<comp,VT>::applyf(value1);

                if (Maybe<max>::less(TMV_REAL(value),best) &&
                    Maybe<max>::less(TMV_REAL(value1),best)) goto L1;
                else if (Maybe<max>::less(TMV_REAL(value1),TMV_REAL(value)))
                { best = TMV_REAL(value); bestit = it; }
                else { best = TMV_REAL(value1); bestit = it+1; }

L1:
                value = it[2];
                value1 = it[3];
                Component<comp,VT>::applyf(value);
                Component<comp,VT>::applyf(value1);

                if (Maybe<max>::less(TMV_REAL(value),best) &&
                    Maybe<max>::less(TMV_REAL(value1),best)) goto L2;
                else if (Maybe<max>::less(TMV_REAL(value1),TMV_REAL(value)))
                { best = TMV_REAL(value); bestit = it+2; }
                else { best = TMV_REAL(value1); bestit = it+3; }
L2:
                it += 4;
            } while (--n_4);
            if (ibest) *ibest = bestit - v.begin();
            return best;
        }
    };

    // algo 6: fully unroll
    template <CompType comp, bool max, class V>
    struct MinMaxElement_Helper<6,comp,max,V>
    {
        typedef typename V::value_type VT;
        typedef typename V::real_type RT;
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
            static inline RT unroll2(const V& v)
            {
                const RT temp1 = Helper1<n-1,x>::unroll2(v);
                const RT temp2 = Component<comp,VT>::f(v.cref(n-1));
                return Maybe<max>::less(temp1,temp2) ? temp2 : temp1;
            }
            static inline IT unroll3(const IT& v0, RT& best)
            {
                const IT temp1 = Helper1<n-1,x>::unroll3(v0,best);
                const IT temp2 = v0 + n-1;
                const RT f2 = Component<comp,VT>::f(*temp2);
                if (Maybe<max>::less(f2,best)) return temp1;
                else { best = f2; return temp2; }
            }
        };
        template <int x>
        struct Helper1<1,x>
        {
            static inline IT unroll1(const IT& v0)
            { return v0; }
            static inline RT unroll2(const V& v)
            { return Component<comp,VT>::f(v.cref(0)); }
            static inline IT unroll3(const IT& v0, RT& best)
            { best = Component<comp,VT>::f(*v0); return v0; }
        };
        template <int x>
        struct Helper1<0,x>; // not defined so give a compile error

        static inline RT call(const V& v, int*const ibest)
        {
            if (ibest) {
                IT bestit = Helper1<V::vsize,1>::unroll1(v.begin());
                *ibest = bestit - v.begin();
                return Component<comp,VT>::f(*bestit);
            } else {
                return Helper1<V::vsize,1>::unroll2(v);
            }
        }
    };

    // algo 7: fully unroll with storage of best along the way
    template <CompType comp, bool max, class V>
    struct MinMaxElement_Helper<7,comp,max,V>
#if 1
    {
        typedef typename V::value_type VT;
        typedef typename V::real_type RT;
        typedef typename V::const_iterator IT;

        template <int n, int x>
        struct Helper1
        {
            static inline void unroll(RT& best, IT& bestit, const IT& v0)
            {
                Helper1<n-1,x>::unroll(best,bestit,v0);
                const RT temp = Component<comp,VT>::f(v0[n-1]);
                if (Maybe<max>::less(temp,best)) return;
                else { best = temp; bestit = v0+n-1; }
            }
        };

        template <int x>
        struct Helper1<1,x>
        { static inline void unroll(RT& best, IT& bestit, const IT& v0) {} };

        template <int x>
        struct Helper1<0,x>; 

        static inline RT call(const V& v, int*const ibest)
        {
            IT it = v.begin();
            RT best = Component<comp,VT>::f(*it);
            IT bestit = it;
            Helper1<V::vsize,1>::unroll(best,bestit,it);
            if (ibest) *ibest = bestit - v.begin();
            return best;
        }
    };
#else
    {
        typedef typename V::value_type VT;
        typedef typename V::real_type RT;

        template <int n, int x>
        struct Helper1
        {
            static inline RT unroll(const V& v, int& ibest)
            {
                const RT temp1 = Helper1<n-1,x>::unroll(v,ibest);
                const RT temp2 = Component<comp,VT>::f(v.cref(n-1));
                if (Maybe<max>::less(temp2,temp1)) return temp1;
                else { ibest = n-1; return temp2; }
            }
        };

        template <int x>
        struct Helper1<1,x>
        {
            static inline RT unroll(const V& v, int& ibest)
            { ibest = 0; return Component<comp,VT>::f(v.cref(0)); }
        };

        template <int x>
        struct Helper1<0,x>; 

        static inline RT call(const V& v, int*const ibest)
        { return Helper1<V::vsize,1>::unroll(v,*ibest); }
    };
#endif

    // algo 8: special algorithm for complex abs(), max
    // Try to avoid as many abs calls as possible by checking the components
    // to see if it's even possible for the abs value to be greater.
    template <class V>
    struct MinMaxElement_Helper<8,AbsComp,true,V>
    {
        typedef typename V::value_type VT;
        typedef typename V::real_type RT;
        typedef typename V::const_iterator IT;

        static inline RT call(const V& v, int*const ibest)
        {
            VT value;
            RT best;
            int n = v.size();
            IT it = v.begin();
            best = TMV_ABS(*it);
            IT bestit = it;

            if (n) do {
                value = *it++;
                // Only if sum of abs(real) + abs(imag) > best
                // should we even both doing the complex abs()
                if (TMV_ABS(real(value)) + TMV_ABS(imag(value)) < best) 
                    continue;
                else { 
                    Component<AbsComp,VT>::applyf(value);
                    if (TMV_REAL(value) < best) continue;
                    else { best = TMV_REAL(value); bestit = it-1; }
                }
            } while (--n);
            if (ibest) *ibest = bestit - v.begin();
            return best;
        }
    };
    // algo 8: min version uses a similar principle
    template <class V>
    struct MinMaxElement_Helper<8,AbsComp,false,V>
    {
        typedef typename V::value_type VT;
        typedef typename V::real_type RT;
        typedef typename V::const_iterator IT;

        static inline RT call(const V& v, int*const ibest)
        {
            VT value;
            RT best;
            int n = v.size();
            IT it = v.begin();
            best = TMV_ABS(*it);
            IT bestit = it;

            if (n) do {
                value = *it++;
                // Only if both abs(real) and abs(imag) < best
                // should we even both doing the complex abs()
                if (TMV_ABS(real(value)) > best && TMV_ABS(imag(value)) > best) 
                    continue;
                else { 
                    Component<AbsComp,VT>::applyf(value);
                    if (TMV_REAL(value) > best) continue;
                    else { best = TMV_REAL(value); bestit = it-1; }
                }
            } while (--n);
            if (ibest) *ibest = bestit - v.begin();
            return best;
        }
    };

    // algo 12: 1 or 2 depending on size
    template <CompType comp, bool max, class V>
    struct MinMaxElement_Helper<12,comp,max,V>
    {
        typedef typename V::value_type VT;
        typedef typename V::real_type RT;
        typedef typename V::const_iterator IT;

        static inline RT call(const V& v, int*const ibest)
        {
#if TMV_OPT >= 2
            if (!ibest && v.size() < 250) 
                return MinMaxElement_Helper<1,comp,max,V>::call(v,ibest);
            else 
#endif
                return MinMaxElement_Helper<2,comp,max,V>::call(v,ibest);
        }
    };

    // algo 13: 1 or 3 depending on size
    template <CompType comp, bool max, class V>
    struct MinMaxElement_Helper<13,comp,max,V>
    {
        typedef typename V::real_type RT;
        static inline RT call(const V& v, int*const ibest)
        {
#if TMV_OPT >= 2
            if (!ibest && v.size() < 250) 
                return MinMaxElement_Helper<1,comp,max,V>::call(v,ibest);
            else 
#endif
                return MinMaxElement_Helper<3,comp,max,V>::call(v,ibest);
        }
    };

    // algo 14: 1 or 4 depending on size
    template <CompType comp, bool max, class V>
    struct MinMaxElement_Helper<14,comp,max,V>
    {
        typedef typename V::real_type RT;
        static inline RT call(const V& v, int*const ibest)
        {
#if TMV_OPT >= 2
            if (!ibest && v.size() < 250) 
                return MinMaxElement_Helper<1,comp,max,V>::call(v,ibest);
            else 
#endif
                return MinMaxElement_Helper<4,comp,max,V>::call(v,ibest);
        }
    };

    // algo 15: 1 or 5 depending on size
    template <CompType comp, bool max, class V>
    struct MinMaxElement_Helper<15,comp,max,V>
    {
        typedef typename V::real_type RT;
        static inline RT call(const V& v, int*const ibest)
        {
#if TMV_OPT >= 2
            if (!ibest && v.size() < 250) 
                return MinMaxElement_Helper<1,comp,max,V>::call(v,ibest);
            else 
#endif
                return MinMaxElement_Helper<5,comp,max,V>::call(v,ibest);
        }
    };

    template <class V>
    inline typename V::value_type InlineMaxElement(
        const BaseVector_Calc<V>& v, int*const imax)
    {
#if TMV_OPT >= 2
        const int algo1 = 
            V::vsize != UNKNOWN ? (
                V::vsize <= 160 ? 6 :
                (V::vstep == 1 && V::vsize > 250) ? (
                    ( V::viscomplex ? 2 : 5 ) ) :
                1 ) :
            V::vstep == 1 ? ( V::viscomplex ? 2 : 14 ) :
            1;
#endif
        const int algo2 = 
#if TMV_OPT >= 1
            V::vsize != UNKNOWN ? (
                V::vsize <= 12 ? 6 :
                V::vsize <= 40 ? 7 :
                (V::vstep == 1 && V::vsize > 250) ? (
                    ( V::viscomplex ? 1 : 5 ) ) :
                1 ) :
            V::vstep == 1 ? ( V::viscomplex ? 1 : 14 ) :
#endif
            1;
#ifdef PRINTALGO
        std::cout<<"InlineMaxElement: algo = "<<algo1<<" "<<algo2<<std::endl;
#endif

#if TMV_OPT >= 2
        if (!imax) 
            return MinMaxElement_Helper<algo1,RealComp,true,V>::call(
                v.vec(),imax);
        else 
#endif
            return MinMaxElement_Helper<algo2,RealComp,true,V>::call(
                v.vec(),imax);
    }

    // Defined in TMV_Vector.cpp
    template <class T>
    T InstMaxElement(const ConstVectorView<T>& v, int*const imax);

    template <bool conj, bool inst, class V>
    struct CallMaxElementv // conj = true
    {
        static inline typename V::value_type call(const V& v, int*const imax)
        { 
            typedef typename V::const_conjugate_type Vc;
            return TMV_CONJ(CallMaxElementv<false,inst,Vc>::call(
                    v.conjugate(),imax));
        }
    };
    template <class V>
    struct CallMaxElementv<false,false,V> // inst = false
    {
        static inline typename V::value_type call(const V& v, int*const imax)
        { return InlineMaxElement(v,imax); }
    };
    template <class V>
    struct CallMaxElementv<false,true,V> // inst = true
    {
        static inline typename V::value_type call(const V& v, int*const imax)
        { return InstMaxElement(v.xView(),imax); }
    };

    template <class V>
    inline typename V::value_type MaxElement(
        const BaseVector_Calc<V>& v, int*const imax)
    {
        TMVAssert(v.size() > 0);
        typedef typename V::value_type T;
        const bool inst = 
            Traits<T>::isinst &&
            V::vsize == UNKNOWN;
        return CallMaxElementv<V::vconj,inst,V>::call(v.vec(),imax);
    }


    // 
    // MaxAbsElement
    //

    template <class V>
    inline typename V::real_type InlineMaxAbsElement(
        const BaseVector_Calc<V>& v, int*const imax)
    {
#if TMV_OPT >= 2
        const int algo1 = 
            V::vsize != UNKNOWN ? (
                V::viscomplex ? (
                    ( V::vsize <= 16 ? 6 : V::vstep == 1 ? 8 : 1 ) ) :
                V::vsize <= 160 ? 6 :
                (V::vstep == 1 && V::vsize > 250) ?  5 : 1 ) :
            V::vstep == 1 ? ( V::viscomplex ? 8 : 15 ) :
            1;
#endif
        const int algo2 = 
#if TMV_OPT >= 1
            V::vsize != UNKNOWN ? (
                V::viscomplex ? (
                    ( V::vstep == 1 ? 8 : 1 ) ) :
                V::vsize <= 12 ? 6 :
                V::vsize <= 40 ? 7 :
                (V::vstep == 1 && V::vsize > 250) ? 5 : 1 ) :
            V::vstep == 1 ? ( V::viscomplex ? 8 : 15 ) :
#endif
            1;
#ifdef PRINTALGO
        std::cout<<"InlineMaxAbsElement: algo = "<<algo1<<
            " "<<algo2<<std::endl;
#endif

#if TMV_OPT >= 2
        if (!imax) 
            return MinMaxElement_Helper<algo1,AbsComp,true,V>::call(
                v.vec(),imax);
        else 
#endif
            return MinMaxElement_Helper<algo2,AbsComp,true,V>::call(
                v.vec(),imax);
    }

    // Defined in TMV_Vector.cpp
    template <class T>
    typename Traits<T>::real_type InstMaxAbsElement(
        const ConstVectorView<T>& v, int*const imax); 

    template <bool inst, class V>
    struct CallMaxAbsElementv // inst = false
    {
        static inline typename V::real_type call(const V& v, int*const imax)
        { return InlineMaxAbsElement(v,imax); }
    };
    template <class V>
    struct CallMaxAbsElementv<true,V>
    {
        static inline typename V::real_type call(const V& v, int*const imax)
        { return InstMaxAbsElement(v.xView(),imax); }
    };

    template <class V>
    inline typename V::real_type MaxAbsElement(
        const BaseVector_Calc<V>& v, int*const imax)
    {
        TMVAssert(v.size() > 0);
        typedef typename V::value_type T;
        typedef typename V::const_nonconj_type Vn;
        const bool inst = 
            Traits<T>::isinst &&
            V::vsize == UNKNOWN;
        return CallMaxAbsElementv<inst,Vn>::call(v.nonConj(),imax);
    }


    //
    // MaxAbs2Element
    //

    template <class V>
    inline typename V::real_type InlineMaxAbs2Element(
        const BaseVector_Calc<V>& v, int*const imax)
    {
#if TMV_OPT >= 2
        const int algo1 = 
            V::vsize != UNKNOWN ? (
                V::vsize <= 160 ? 6 :
                (V::vstep == 1 && V::vsize > 250) ? (
                    ( V::viscomplex ? 2 : 5 ) ) :
                1 ) :
            V::vstep == 1 ? ( V::viscomplex ? 2 : 15 ) :
            1;
#endif
        const int algo2 = 
#if TMV_OPT >= 1
            V::vsize != UNKNOWN ? (
                V::vsize <= 12 ? 6 :
                V::vsize <= 40 ? 7 :
                (V::vstep == 1 && V::vsize > 250) ? (
                    ( V::viscomplex ? 1 : 5 ) ) :
                1 ) :
            V::vstep == 1 ? ( V::viscomplex ? 1 : 15 ) :
#endif
            1;
#ifdef PRINTALGO
        std::cout<<"InlineMaxAbs2Element: algo = "<<algo1<<
            " "<<algo2<<std::endl;
#endif

#if TMV_OPT >= 2
        if (!imax) 
            return MinMaxElement_Helper<algo1,Abs2Comp,true,V>::call(
                v.vec(),imax);
        else 
#endif
            return MinMaxElement_Helper<algo2,Abs2Comp,true,V>::call(
                v.vec(),imax);
    }

    // Defined in TMV_Vector.cpp
    template <class T>
    typename Traits<T>::real_type InstMaxAbs2Element(
        const ConstVectorView<T>& v, int*const imax);

    template <class T>
    inline typename Traits<T>::real_type InstMaxAbs2Element(
        const ConstVectorView<T,UNKNOWN,true>& v, int*const imax)
    { return InstMaxAbs2Element(v.conjugate(),imax); }

    template <bool inst, class V>
    struct CallMaxAbs2Elementv // inst = false
    {
        static inline typename V::real_type call(const V& v, int*const imax)
        { return InlineMaxAbs2Element(v,imax); }
    };
    template <class V>
    struct CallMaxAbs2Elementv<true,V>
    {
        static inline typename V::real_type call(const V& v, int*const imax)
        { return InstMaxAbs2Element(v.xView(),imax); }
    };

    template <class V>
    inline typename V::real_type MaxAbs2Element(
        const BaseVector_Calc<V>& v, int*const imax)
    {
        TMVAssert(v.size() > 0);
        typedef typename V::value_type T;
        typedef typename V::const_nonconj_type Vn;
        const bool inst = 
            Traits<T>::isinst &&
            V::vsize == UNKNOWN;
        return CallMaxAbs2Elementv<inst,Vn>::call(v.nonConj(),imax);
    }


    // 
    // MinElement
    //

    template <class V>
    inline typename V::value_type InlineMinElement(
        const BaseVector_Calc<V>& v, int*const imin)
    {
#if TMV_OPT >= 2
        const int algo1 = 
            V::vsize != UNKNOWN ? (
                V::vsize <= 160 ? 6 :
                (V::vstep == 1 && V::vsize > 250) ? (
                    ( V::viscomplex ? 2 : 5 ) ) :
                1 ) :
            V::vstep == 1 ? ( V::viscomplex ? 2 : 14 ) :
            1;
#endif
        const int algo2 = 
#if TMV_OPT >= 1
            V::vsize != UNKNOWN ? (
                V::vsize <= 12 ? 6 :
                V::vsize <= 40 ? 7 :
                (V::vstep == 1 && V::vsize > 250) ? (
                    ( V::viscomplex ? 1 : 5 ) ) :
                1 ) :
            V::vstep == 1 ? ( V::viscomplex ? 1 : 14 ) :
#endif
            1;
#ifdef PRINTALGO
        std::cout<<"InlineMinElement: algo = "<<algo1<<" "<<algo2<<std::endl;
#endif

#if TMV_OPT >= 2
        if (!imin) 
            return MinMaxElement_Helper<algo1,RealComp,false,V>::call(
                v.vec(),imin);
        else 
#endif
            return MinMaxElement_Helper<algo2,RealComp,false,V>::call(
                v.vec(),imin);
    }

    // Defined in TMV_Vector.cpp
    template <class T>
    T InstMinElement(const ConstVectorView<T>& v, int*const imin);

    template <bool conj, bool inst, class V>
    struct CallMinElementv // conj = true
    {
        static inline typename V::value_type call(const V& v, int*const imin)
        { 
            typedef typename V::const_conjugate_type Vc;
            return TMV_CONJ(CallMinElementv<false,inst,Vc>::call(
                    v.conjugate(),imin));
        }
    };
    template <class V>
    struct CallMinElementv<false,false,V> // inst = false
    {
        static inline typename V::value_type call(const V& v, int*const imin)
        { return InlineMinElement(v,imin); }
    };
    template <class V>
    struct CallMinElementv<false,true,V> // inst = true
    {
        static inline typename V::value_type call(const V& v, int*const imin)
        { return InstMinElement(v.xView(),imin); }
    };

    template <class V>
    inline typename V::value_type MinElement(
        const BaseVector_Calc<V>& v, int*const imin)
    {
        TMVAssert(v.size() > 0);
        typedef typename V::value_type T;
        const bool inst = 
            Traits<T>::isinst &&
            V::vsize == UNKNOWN;
        return CallMinElementv<V::vconj,inst,V>::call(v.vec(),imin);
    }


    // 
    // MinAbsElement
    //

    template <class V>
    inline typename V::real_type InlineMinAbsElement(
        const BaseVector_Calc<V>& v, int*const imin)
    {
#if TMV_OPT >= 2
        const int algo1 = 
            V::vsize != UNKNOWN ? (
                V::viscomplex ? (
                    ( V::vsize <= 40 ? 6 : V::vstep == 1 ? 8 : 1 ) ) :
                V::vsize <= 160 ? 6 :
                (V::vstep == 1 && V::vsize > 250) ?  5 : 1 ) :
            V::vstep == 1 ? ( V::viscomplex ? 8 : 15 ) :
            1;
#endif
        const int algo2 = 
#if TMV_OPT >= 1
            V::vsize != UNKNOWN ? (
                V::viscomplex ? (
                    ( V::vstep == 1 ? ( V::vsize <= 20 ? 2 : 8 ) : 1 ) ) :
                V::vsize <= 12 ? 6 :
                V::vsize <= 40 ? 7 :
                (V::vstep == 1 && V::vsize > 250) ? 5 : 1 ) :
            V::vstep == 1 ? ( V::viscomplex ? 8 : 15 ) :
#endif
            1;
#ifdef PRINTALGO
        std::cout<<"InlineMinAbsElement: algo = "<<algo1<<
            " "<<algo2<<std::endl;
#endif

#if TMV_OPT >= 2
        if (!imin) 
            return MinMaxElement_Helper<algo1,AbsComp,false,V>::call(
                v.vec(),imin);
        else 
#endif
            return MinMaxElement_Helper<algo2,AbsComp,false,V>::call(
                v.vec(),imin);
    }

    // Defined in TMV_Vector.cpp
    template <class T>
    typename Traits<T>::real_type InstMinAbsElement(
        const ConstVectorView<T>& v, int*const imin);

    template <class T>
    inline typename Traits<T>::real_type InstMinAbsElement(
        const ConstVectorView<T,UNKNOWN,true>& v, int*const imin)
    { return InstMinAbsElement(v.conjugate(),imin); }

    template <bool inst, class V>
    struct CallMinAbsElementv // inst = false
    {
        static inline typename V::real_type call(const V& v, int*const imin)
        { return InlineMinAbsElement(v,imin); }
    };
    template <class V>
    struct CallMinAbsElementv<true,V>
    {
        static inline typename V::real_type call(const V& v, int*const imin)
        { return InstMinAbsElement(v.xView(),imin); }
    };

    template <class V>
    inline typename V::real_type MinAbsElement(
        const BaseVector_Calc<V>& v, int*const imin)
    {
        TMVAssert(v.size() > 0);
        typedef typename V::value_type T;
        typedef typename V::const_nonconj_type Vn;
        const bool inst = 
            Traits<T>::isinst &&
            V::vsize == UNKNOWN;
        return CallMinAbsElementv<inst,Vn>::call(v.nonConj(),imin);
    }

    // 
    // MinAbs2Element
    //

    template <class V>
    inline typename V::real_type InlineMinAbs2Element(
        const BaseVector_Calc<V>& v, int*const imin)
    {
#if TMV_OPT >= 2
        const int algo1 = 
            V::vsize != UNKNOWN ? (
                V::vsize <= 160 ? 6 :
                (V::vstep == 1 && V::vsize > 250) ? (
                    ( V::viscomplex ? 2 : 5 ) ) :
                1 ) :
            V::vstep == 1 ? ( V::viscomplex ? 2 : 15 ) :
            1;
#endif
        const int algo2 = 
#if TMV_OPT >= 1
            V::vsize != UNKNOWN ? (
                V::vsize <= 12 ? 6 :
                V::vsize <= 40 ? 7 :
                (V::vstep == 1 && V::vsize > 250) ? (
                    ( V::viscomplex ? 1 : 5 ) ) :
                1 ) :
            V::vstep == 1 ? ( V::viscomplex ? 1 : 15 ) :
#endif
            1;
#ifdef PRINTALGO
        std::cout<<"InlineMinAbs2Element: algo = "<<algo1<<
            " "<<algo2<<std::endl;
#endif

#if TMV_OPT >= 2
        if (!imin) 
            return MinMaxElement_Helper<algo1,Abs2Comp,false,V>::call(
                v.vec(),imin);
        else 
#endif
            return MinMaxElement_Helper<algo2,Abs2Comp,false,V>::call(
                v.vec(),imin);
    }

    // Defined in TMV_Vector.cpp
    template <class T>
    typename Traits<T>::real_type InstMinAbs2Element(
        const ConstVectorView<T>& v, int*const imin);

    template <class T>
    inline typename Traits<T>::real_type InstMinAbs2Element(
        const ConstVectorView<T,UNKNOWN,true>& v, int*const imin)
    { return InstMinAbs2Element(v.conjugate(),imin); }

    template <bool inst, class V>
    struct CallMinAbs2Elementv // inst = false
    {
        static inline typename V::real_type call(const V& v, int*const imin)
        { return InlineMinAbs2Element(v,imin); }
    };
    template <class V>
    struct CallMinAbs2Elementv<true,V>
    {
        static inline typename V::real_type call(const V& v, int*const imin)
        { return InstMinAbs2Element(v.xView(),imin); }
    };

    template <class V>
    inline typename V::real_type MinAbs2Element(
        const BaseVector_Calc<V>& v, int*const imin)
    {
        TMVAssert(v.size() > 0);
        typedef typename V::value_type T;
        typedef typename V::const_nonconj_type Vn;
        const bool inst = 
            Traits<T>::isinst &&
            V::vsize == UNKNOWN;
        return CallMinAbs2Elementv<inst,Vn>::call(v.nonConj(),imin);
    }

} // namespace tmv

#endif
