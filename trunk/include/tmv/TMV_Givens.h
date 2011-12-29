
//
// This file contains the code to implement Givens rotations
//
// A Givens rotation is a 2x2 matrix of the form:
//
// G = [  c  s  ]
//     [ -s* c* ]
//
// where |c|^2 + |s|^2 = 1.
// Normally, people consider real Givens matrices where the *'s
// (conjugations) are not necessary.  But they generalize easily
// to the complex case as above.
//
// Givens matrices are used to zero out a single element in a matrix by
// rotating a 2 element vector so that its second element is 0:
//
// [  c  s  ] [ x ] = [ r ]
// [ -s* c* ] [ y ]   [ 0 ]
// 
// [ Note: often (eg. in Golub and Van Loan) Givens matrices are defined
//   as the transpose (or adjoint) of this.  I prefer this definition. ]
//
// (1)  c x + s y = r
// (2)  c* y - s* x = 0
//
// From (2), we get  s* = yc*/x, so 
// |x|^2 |s|^2 = |y|^2 |c|^2 = |y|^2 (1 - |s|^2)
//
// |s| = |y|/sqrt(|x|^2+|y|^2)  and |c| = |x|/sqrt(|x|^2+|y|^2)
//
// From (1), r = cx+sy = cx+yy*c/x* = c/x* (|x|^2+|y|^2)
// |r| = sqrt(|x|^2+|y|^2)
//
// All that remain is to determine the phases of c,s
//
// Let x = |x| e^ip
//     y = |y| e^iq
//     r = |r| e^it
//     c = |x|/|r| e^ia
//     s = |y|/|r| e^ib
//
// Plugging these into (1), we find that  t = a+p = b+q
//
// Since a priori, we should consider all three of (a,b,t)
// to be unknown, we are left with one free parameter.
// There are (at least) 3 reasonable options:
//
// Choice 1: a = -p, b = -q, t = 0
//
// This has the advantage that r is real and positive.
// Also, the equations for c,s look much like the usual
// forms for real Givens matrices:
//
// c = x*/|r|  s = y*/|r|   r = |r|
//
// Choice 2: a = 0, b = p-q, t = p
//
// This has the advantage that c is real, so the Givens
// matrix only needs to store 3 real numbers, not 4.
// Multiplication later is also faster by a factor
// of 4/3.
//
// c = |x|/|r|  s = (x/|x|)y*/|r|  r = (x/|x|)|r|
//
// Choice 3: a = q-p, b = 0, t = q
//
// Similar to choice 2, but s is real:
// 
// c = (y/|y|)x*/|r|  s = |y|/|r|  r = (y/|y|)|r|
//
// 
// The speed of calculation is usually more important than the 
// convenience of a rotating to a real value, so choice 2 or 3
// is preferred.  
// Choice 2 is a bit better when |y|<|x|, since we prefer to 
// calculate sqrt(1+|y/x|^2), and the calculations for 2 work 
// better in this case.
// Since this is more common than |x|<|y|, we always use choice 2.
//

#ifndef TMV_Givens_H
#define TMV_Givens_H

#include "tmv/TMV_BaseVector.h"
#include "tmv/TMV_BaseMatrix_Rec.h"

namespace tmv {

    
    template <class T, class RT>
    void InstGivensRotate(T& x, T& y, RT& c, T& s);
    template <class T, class RT, class T2>
    void InstGivensMultV(RT c, T s, VectorView<T2> v1, VectorView<T2> v2);

    // Use a Givens matrix G to rotate the vector so y = 0:
    // G [ x ] = [ r ]
    //   [ y ]   [ 0 ]
    // Also, return the Givens matrix used to do so.
    template <int algo, class T>
    struct GivensRotate_Helper;

    template <class T>
    struct GivensRotate_Helper<11,T>
    {
        typedef typename Traits<T>::real_type RT;
        static void call(T& x, T& y, RT& c, T& s)
        {
            // Use a Givens matrix G to rotate the vector so y = 0:
            // G [ x ] = [ r ]
            //   [ y ]   [ 0 ]
            // Also, return the Givens matrix used to do so.

            const RT sqrteps = TMV_SQRT(TMV_Epsilon<T>());

            RT absx = TMV_ABS2(x);
            RT absy = TMV_ABS2(y);
            if (absy == RT(0)) {
                c = RT(1);
                s = T(0);
                y = RT(0);
#ifdef PRINTALGO_GIVENS
                std::cout<<"Simple1: return 1,0\n";
#endif
            } else if (absx == RT(0)) {
                x = RT(0);
                RT absy = TMV_ABS(y);
                c = RT(0);
                s = TMV_SIGN(TMV_CONJ(y),absy);
                x = absy; y = RT(0);
#ifdef PRINTALGO_GIVENS
                std::cout<<"Simple2: return 0,"<<s<<std::endl;
#endif
            } else if (absx > absy) {
                if (absy <= sqrteps*absx) {
                    // Then everything simplifies:
                    // c = 1
                    // s = (y/x)* 
                    // r = f
                    c = RT(1);
                    s = TMV_CONJ(y)/TMV_CONJ(x);
                    y = RT(0);
#ifdef PRINTALGO_GIVENS
                    std::cout<<"Simple3: return 1,"<<s<<std::endl;
#endif
                } else {
                    // c = 1/sqrt(1+|y/x|^2)
                    // s = (y/x)*/sqrt(1+|y/x|^2)
                    // r = x sqrt(1+|y/x|^2)
                    // We get a slightly more accurate calculation of r
                    // if we calculate r-x and add this to x:
                    // r-x = x (sqrt(1+|y/x|^2)-1) = x |y/x|^2 / (1 + sqrt(1+|y/x|^2))
                    T yoverx = y/x;
                    RT n = TMV_NORM(yoverx);
                    RT sqrtfactor = TMV_SQRT(RT(1)+n);
                    c = RT(1)/sqrtfactor;
                    s = TMV_CONJ(yoverx)*c;
                    x += x*(n/(RT(1)+sqrtfactor));
                    y = RT(0);
#ifdef PRINTALGO_GIVENS
                    std::cout<<"x>y: return "<<c<<','<<s<<std::endl;
#endif
                }
            } else {
                T xovery = x/y;
                RT n = TMV_NORM(xovery);
                RT absxovery = TMV_SQRT(n);
                if (n <= TMV_Epsilon<T>()) {
                    // c = |x/y|
                    // s = (x/y)/|x/y|
                    // r = x/|x/y|
                    c = absxovery;
                    s = TMV_SIGN(xovery,absxovery);
                    x = s*y;
                    y = RT(0);
#ifdef PRINTALGO_GIVENS
                    std::cout<<"Simple4: return "<<absxovery<<','<<s<<std::endl;
#endif
                } else {
                    // c = |x/y|/sqrt(1+|x/y|^2)
                    // s = (x/y)/|x/y|/sqrt(1+|x/y|^2)
                    // r = x/|x/y| sqrt(|x|^2+|y|^2)
                    RT sqrtfactor = TMV_SQRT(RT(1)+n);
                    RT invsqrtfactor = RT(1)/sqrtfactor;
                    T signxovery = TMV_SIGN(xovery,absxovery); // (x/y)/|x/y|
                    c = absxovery*invsqrtfactor;
                    s = signxovery*invsqrtfactor;
                    x = y * signxovery * sqrtfactor;
                    y = T(0);
#ifdef PRINTALGO_GIVENS
                    std::cout<<"x<y: return "<<c<<','<<s<<std::endl;
#endif
                }
            }
        }
    };

    // algo 90: Call inst
    template <class T>
    struct GivensRotate_Helper<90,T>
    {
        typedef typename Traits<T>::real_type RT;
        static void call(T& x, T& y, RT& c, T& s)
        { InstGivensRotate(x,y,c,s); }
    };

    // algo -3: Select algorithm
    template <class T>
    struct GivensRotate_Helper<-3,T>
    {
        typedef typename Traits<T>::real_type RT;
        static void call(T& x, T& y, RT& c, T& s)
        {
            const int algo = 11;
#ifdef PRINTALGO_GIVENS
            std::cout<<"Inline GivensRotate\n";
            std::cout<<"x,y = "<<x<<','<<y<<std::endl;
#endif
            GivensRotate_Helper<algo,T>::call(x,y,c,s);
#ifdef PRINTALGO_GIVENS
            std::cout<<"x,y => "<<x<<','<<y<<std::endl;
            std::cout<<"c,s = "<<s<<','<<s<<std::endl;
#endif
        }
    };

    // algo -2: Check for inst
    template <class T>
    struct GivensRotate_Helper<-2,T>
    {
        typedef typename Traits<T>::real_type RT;
        static void call(T& x, T& y, RT& c, T& s)
        {
            const bool inst = Traits<T>::isinst;
            const int algo = 
                inst ? 90 :
                -3;
            GivensRotate_Helper<algo,T>::call(x,y,c,s);
        }
    };

    template <class T>
    struct GivensRotate_Helper<-1,T>
    {
        typedef typename Traits<T>::real_type RT;
        static void call(T& x, T& y, RT& c, T& s)
        { GivensRotate_Helper<-2,T>::call(x,y,c,s); }
    };

    template <class T, class RT>
    inline void InlineGivensRotate(T& x, T& y, RT& c, T& s)
    { GivensRotate_Helper<-3,T>::call(x,y,c,s); }

    // Note: GivensRotate(x,y) without the "Do" only takes the x,y argument 
    // and returns a Givens<T> object.  This function is found below.
    template <class T, class RT>
    inline void DoGivensRotate(T& x, T& y, RT& c, T& s)
    { GivensRotate_Helper<-2,T>::call(x,y,c,s); }


    //
    // GivensMult of a single element:
    // [ x' ] = [  c  s ] [ x ] = [ cx+sy  ]
    // [ y' ]   [ -s* c ] [ y ]   [ cy-s*x ]
    //
    
    // Better to leave this one inline.  So no Helper struct here.
    template <class RTg, class Tg, class T> 
    inline void GivensMult(RTg c, Tg s, T& x, T& y)
    {
        T xx = c*x+s*y;
        y = c*y-TMV_CONJ(s)*x;
        x = xx;
    }
         
    template <class RTg, class Tg, class T> 
    TMV_INLINE void GivensMult(RTg c, Tg s, ConjRef<T> x, ConjRef<T> y)
    { GivensMult(c,TMV_CONJ(s),x.getRef(),y.getRef()); }

    //
    // GivensMult of a pair of vectors
    //
    
    template <int also, int size, class T, class V1, class V2>
    struct GivensMultV_Helper;

    // algo 0: size == 0, nothing to do
    template <int size, class T, class V1, class V2>
    struct GivensMultV_Helper<0,size,T,V1,V2>
    {
        typedef typename Traits<T>::real_type RT;
        static TMV_INLINE void call(RT , T , V1& , V2& ) {}
    };

    // algo 11: simple for loop
    template <int size, class T, class V1, class V2>
    struct GivensMultV_Helper<11,size,T,V1,V2>
    {
        typedef typename Traits<T>::real_type RT;
        static void call(RT c, T s, V1& v1, V2& v2)
        {
            TMVAssert(v1.size() == v2.size());
            TMVAssert(!SameStorage(v1,v2));
            TMVStaticAssert(!V1::_conj);
            TMVStaticAssert(!V2::_conj);
            int N = size == TMV_UNKNOWN ? v2.size() : size;

            typedef typename V1::iterator IT1;
            typedef typename V2::iterator IT2;
            typedef typename V1::value_type T2;

            IT1 it1 = v1.begin();
            IT2 it2 = v2.begin();

            T2 x,y;

            if (N) do {
                x = *it1;
                y = *it2;
                *it1++ = c*x+s*y;
                *it2++ = c*y-TMV_CONJ(s)*x;
            } while (--N);
        }
    };

    // algo 12: 2 at a time
    template <int size, class T, class V1, class V2>
    struct GivensMultV_Helper<12,size,T,V1,V2>
    {
        typedef typename Traits<T>::real_type RT;
        static void call(RT c, T s, V1& v1, V2& v2)
        {
            TMVAssert(v1.size() == v2.size());
            TMVAssert(!SameStorage(v1,v2));
            TMVStaticAssert(!V1::_conj);
            TMVStaticAssert(!V2::_conj);
            const int N = size == TMV_UNKNOWN ? v2.size() : size;
            int N_2 = (N >> 1);
            const int N_x = N - (N_2 << 1);

            typedef typename V1::iterator IT1;
            typedef typename V2::iterator IT2;
            typedef typename V1::value_type T2;

            IT1 it1 = v1.begin();
            IT2 it2 = v2.begin();

            T2 x0,y0,x1,y1;

            if (N_2) do {
                x0 = it1[0];  x1 = it1[1];
                y0 = it2[0];  y1 = it2[1];
                it1[0] = c*x0+s*y0;
                it1[1] = c*x1+s*y1;
                it2[0] = c*y0-TMV_CONJ(s)*x0;
                it2[1] = c*y1-TMV_CONJ(s)*x1;
                it1 += 2; it2 += 2;
            } while (--N_2);
            if (N_x) {
                x0 = *it1;
                y0 = *it2;
                *it1 = c*x0+s*y0;
                *it2 = c*y0-TMV_CONJ(s)*x0;
            }
        }
    };

    // algo 13: 4 at a time
    template <int size, class T, class V1, class V2>
    struct GivensMultV_Helper<13,size,T,V1,V2>
    {
        typedef typename Traits<T>::real_type RT;
        static void call(RT c, T s, V1& v1, V2& v2)
        {
            TMVAssert(v1.size() == v2.size());
            TMVAssert(!SameStorage(v1,v2));
            TMVStaticAssert(!V1::_conj);
            TMVStaticAssert(!V2::_conj);
            const int N = size == TMV_UNKNOWN ? v2.size() : size;
            int N_4 = (N >> 2);
            int N_x = N - (N_4 << 2);

            typedef typename V1::iterator IT1;
            typedef typename V2::iterator IT2;
            typedef typename V1::value_type T2;

            IT1 it1 = v1.begin();
            IT2 it2 = v2.begin();

            T2 x0,y0,x1,y1,x2,y2,x3,y3;

            if (N_4) do {
                x0 = it1[0];  x1 = it1[1]; x2 = it1[2];  x3 = it1[3];
                y0 = it2[0];  y1 = it2[1]; y2 = it2[2];  y3 = it2[3];
                it1[0] = c*x0+s*y0;
                it1[1] = c*x1+s*y1;
                it1[2] = c*x2+s*y2;
                it1[3] = c*x3+s*y3;
                it2[0] = c*y0-TMV_CONJ(s)*x0;
                it2[1] = c*y1-TMV_CONJ(s)*x1;
                it2[2] = c*y2-TMV_CONJ(s)*x2;
                it2[3] = c*y3-TMV_CONJ(s)*x3;
                it1 += 4; it2 += 4;
            } while (--N_4);
            if (N_x) do {
                x0 = *it1;
                y0 = *it2;
                *it1++ = c*x0+s*y0;
                *it2++ = c*y0-TMV_CONJ(s)*x0;
            } while (--N_x);
        }
    };

    // TODO: SSE versions

    // algo 90: Call inst
    template <int size, class T, class V1, class V2>
    struct GivensMultV_Helper<90,size,T,V1,V2>
    {
        typedef typename Traits<T>::real_type RT;
        static void call(RT c, T s, V1& v1, V2& v2)
        { InstGivensMultV(c,s,v1.xView(),v2.xView()); }
    };

    // algo 97: Conjugate
    template <int size, class T, class V1, class V2>
    struct GivensMultV_Helper<97,size,T,V1,V2>
    {
        typedef typename Traits<T>::real_type RT;
        static void call(RT c, T s, V1& v1, V2& v2)
        {
            typedef typename V1::conjugate_type V1c;
            typedef typename V2::conjugate_type V2c;
            V1c v1c = v1.conjugate();
            V2c v2c = v2.conjugate();
            GivensMultV_Helper<-2,size,T,V1c,V2c>::call(c,TMV_CONJ(s),v1c,v2c);
        }
    };

    // algo -3: Determine which algorithm to use
    template <int size, class T, class V1, class V2>
    struct GivensMultV_Helper<-3,size,T,V1,V2>
    {
        typedef typename Traits<T>::real_type RT;
        static void call(RT c, T s, V1& v1, V2& v2)
        {
            TMVAssert(v1.size() == v2.size());
            TMVAssert(!SameStorage(v1,v2));
            TMVStaticAssert(!V1::_conj);
            TMVStaticAssert(!V2::_conj);
            
            typedef typename V1::value_type T1;
            typedef typename V2::value_type T2;
            TMVStaticAssert((Traits2<T1,T2>::sametype));

            const int allunit = V1::_step == 1 && V2::_step == 1;
            const bool sreal = Traits<T>::isreal;
            const bool vreal = Traits<T1>::isreal;
            const int algo = 
                size == 0 ? 0 :
                TMV_OPT == 0 ? 11 :
                allunit && sreal && vreal && sizeof(T1) == 4 ? 13 :
                allunit && sreal && vreal && sizeof(T1) == 8 ? 12 :
                11;
#ifdef PRINTALGO_GIVENS
            std::cout<<"Inline GivensMultMV: c,s = "<<c<<','<<s<<std::endl;
            std::cout<<"v1 = "<<TMV_Text(v1)<<std::endl;
            std::cout<<"v2 = "<<TMV_Text(v2)<<std::endl;
            std::cout<<"size = "<<size<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
#endif
            GivensMultV_Helper<algo,size,T,V1,V2>::call(c,s,v1,v2);
        }
    };

    // algo -2: Check for inst
    template <int size, class T, class V1, class V2>
    struct GivensMultV_Helper<-2,size,T,V1,V2>
    {
        typedef typename Traits<T>::real_type RT;
        static void call(RT c, T s, V1& v1, V2& v2)
        {
            TMVAssert(v1.size() == v2.size());
            TMVAssert(!SameStorage(v1,v2));
            TMVStaticAssert(V1::_conj == V2::_conj);
            
            typedef typename V1::value_type T1;
            typedef typename V2::value_type T2;
            TMVStaticAssert((Traits2<T1,T2>::sametype));

            const bool inst = 
                (size == TMV_UNKNOWN || size > 16) &&
                Traits<T1>::isinst;
            const int algo = 
                size == 0 ? 0 :
                V2::_conj ? 97 :
                inst ? 90 :
                -3;
            GivensMultV_Helper<algo,size,T,V1,V2>::call(c,s,v1,v2);
        }
    };

    template <int size, class T, class V1, class V2>
    struct GivensMultV_Helper<-1,size,T,V1,V2>
    {
        typedef typename Traits<T>::real_type RT;
        static void call(RT c, T s, V1& v1, V2& v2)
        { GivensMultV_Helper<-2,size,T,V1,V2>::call(c,s,v1,v2); }
    };

    template <class T, class V1, class V2>
    inline void InlineGivensMultV(
        typename Traits<T>::real_type c, T s, 
        BaseVector_Mutable<V1>& v1, BaseVector_Mutable<V2>& v2)
    {
        typedef typename V1::value_type T1;
        typedef typename V2::value_type T2;
        TMVStaticAssert((Traits2<T1,T2>::sametype));
        TMVStaticAssert(V1::_conj == V2::_conj);
        TMVAssert(!SameStorage(v1,v2));

        TMVStaticAssert((Sizes<V1::_size,V2::_size>::same));
        TMVAssert(v1.size() == v2.size());
        const int size = Sizes<V1::_size,V2::_size>::size;

        typedef typename V1::cview_type V1v;
        typedef typename V2::cview_type V2v;
        TMV_MAYBE_REF(V1,V1v) v1v = v1.cView();
        TMV_MAYBE_REF(V2,V2v) v2v = v2.cView();
        GivensMultV_Helper<-3,size,T,V1v,V2v>::call(c,s,v1v,v2v);
    }

    template <class T, class V1, class V2>
    inline void GivensMultV(
        typename Traits<T>::real_type c, T s, 
        BaseVector_Mutable<V1>& v1, BaseVector_Mutable<V2>& v2)
    {
        typedef typename V1::value_type T1;
        typedef typename V2::value_type T2;
        TMVStaticAssert((Traits2<T1,T2>::sametype));
        TMVStaticAssert(V1::_conj == V2::_conj);
        TMVAssert(!SameStorage(v1,v2));

        TMVStaticAssert((Sizes<V1::_size,V2::_size>::same));
        TMVAssert(v1.size() == v2.size());
        const int size = Sizes<V1::_size,V2::_size>::size;

        typedef typename V1::cview_type V1v;
        typedef typename V2::cview_type V2v;
        TMV_MAYBE_REF(V1,V1v) v1v = v1.cView();
        TMV_MAYBE_REF(V2,V2v) v2v = v2.cView();
        GivensMultV_Helper<-2,size,T,V1v,V2v>::call(c,s,v1v,v2v);
    }


    //
    // Symmetric Givens Mult
    //

    template <class RTg, class Tg, class T> 
    inline void GivensHermMult(RTg c, Tg s, T& d0, T& d1, T& e0)
    {
        // [ d0 e0* ] = [  c  s ] [ d0 e0* ] [ c  -s ]
        // [ e0 d1  ]   [ -s* c ] [ e0 d1  ] [ s*  c ]
        // =[c^2 d0 + 2c Re(s e0) + |s|^2 d1   cs(d1-d0) + c^2 e0* - s^2 e0   ]
        //  [cs*(d1-d0) + c^2 e0 - s*^2 e0*    c^2 d1 - 2c Re(s e0) + |s|^2 d0]
        // (using c^2 = 1-|s|^2:)
        // d0' = d0 + 2c Re(s e0) + |s|^2 (d1-d0)
        // e0' = e0 - 2s* Re(s e0) + c s* (d1-d0)
        // d1' = d1 - 2c Re(s e0) - |s|^2 (d1-d0)

        typedef typename Traits<T>::real_type RT;
        RTg s2 = TMV_NORM(s);
        RT Rese0 = TMV_REAL(s*e0);
        T d1md0 = d1-d0;
        T dd = RT(2)*c*Rese0 + s2*d1md0;
        d0 += dd;
        d1 -= dd;
        e0 += TMV_CONJ(s)*(c*d1md0 - RT(2)*Rese0);
    }
          
    template <class RTg, class Tg, class T> 
    inline void GivensSymMult(RTg c, Tg s, T& d0, T& d1, T& e0)
    {
        // [ d0 e0 ] = [  c  s ] [ d0 e0 ] [ c -s* ]
        // [ e0 d1 ]   [ -s* c ] [ e0 d1 ] [ s  c  ]
        // d0' = d0 + 2c (s e0) + s^2 (d1-d0)
        // e0' = e0 - 2|s|^2 e0 + c (s d1 - s* d0)
        // d1' = d1 - 2c (s e0) - s^2 (d1-d0)

        typedef typename Traits<T>::real_type RT;
        Tg s2 = s*s;
        T se0 = s*e0;
        T d1md0 = d1-d0;
        T dd = RT(2)*c*se0 + s2*d1md0;
        d0 += dd;
        d1 -= dd;
        if (Traits<T>::isreal)
            e0 += s*(c*d1md0 - RT(2)*se0);
        else
            e0 += c*(s*d1-TMV_CONJ(s)*d0) - RT(2)*TMV_CONJ(s)*se0;
    }


    //
    // Now a convenient helper class for calling the above functions.
    //
    
    template <class T> 
    class Givens 
    {
        typedef typename Traits<T>::real_type RT;

    public:

        // Constructors
        Givens(const RT& _c, const T& _s) : c(_c), s(_s) {}
        ~Givens() {}
        // Use default copy, op=

        Givens<T> transpose() const 
        { return Givens<T>(c,-TMV_CONJ(s)); }

        Givens<T> conjugate() const 
        { return Givens<T>(c,TMV_CONJ(s)); }

        Givens<T> adjoint() const 
        { return Givens<T>(c,-s); }

        template <class T2> 
        void mult(T2& x, T2& y) const 
        { GivensMult(c,s,x,y); }
        template <class T2> 
        void conjMult(T2& x, T2& y) const 
        { GivensMult(c,TMV_CONJ(s),x,y); }

        template <class T2> 
        void mult(ConjRef<T2> x, ConjRef<T2> y) const 
        { GivensMult(c,TMV_CONJ(s),x.getRef(),y.getRef()); }
        template <class T2> 
        void conjMult(ConjRef<T2> x, ConjRef<T2> y) const 
        { GivensMult(c,s,x.getRef(),y.getRef()); }

        template <class V2> 
        void mult(BaseVector_Mutable<V2>& v) const
        { TMVAssert(v.size()==2); GivensMult(c,s,v.ref(0),v.ref(1)); }
        template <class V2> 
        void conjMult(BaseVector_Mutable<V2>& v) const
        { TMVAssert(v.size()==2); GivensMult(c,TMV_CONJ(s),v.ref(0),v.ref(1)); }

        template <class M2> 
        void mult(BaseMatrix_Rec_Mutable<M2>& m)
        { 
            TMVAssert(m.colsize()==2);
            typedef typename M2::row_type M2r;
            M2r m0 = m.row(0);
            M2r m1 = m.row(1);
            GivensMultV(c,s,m0,m1);
        }
        template <class M2> 
        void conjMult(BaseMatrix_Rec_Mutable<M2>& m)
        { 
            TMVAssert(m.colsize()==2);
            typedef typename M2::row_type M2r;
            M2r m0 = m.row(0);
            M2r m1 = m.row(1);
            GivensMultV(c,TMV_CONJ(s),m0,m1);
        }

        template <class T2> 
        void hermMult(T2& d0, T2& d1, T2& e0) const 
        { GivensHermMult(c,s,d0,d1,e0); }
        template <class T2> 
        void symMult(T2& d0, T2& d1, T2& e0) const 
        { GivensSymMult(c,s,d0,d1,e0); }

    private:

        RT c;
        T s;
    };

    template <class T>
    inline Givens<T> GivensRotate(T& x, T& y)
    {
        typename Traits<T>::real_type c;
        T s;
        DoGivensRotate(x,y,c,s);
        return Givens<T>(c,s);
    }

} // namespace tmv

#endif
