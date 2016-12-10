

#ifndef TMV_Scaling_H
#define TMV_Scaling_H

#include "TMV_Base.h"

namespace tmv {

    // This structure implements a scaling factor with specialization
    // for x = 1 or -1
    // The ix is a flag integer that can be 1, -1 or 0.
    // When ix is 1 or -1 then x is 1 or -1 and isn't stored.
    // When ix is 0, then x is stored and acts normally.
    template <int ix, class T>
    struct Scaling;

    template <class T>
    struct Scaling<0,T>
    {
#ifdef TMV_EXTRA_DEBUG
        T x;
#else
        const T x;
#endif
        TMV_INLINE Scaling(const T _x) : x(_x) {}
        TMV_INLINE_ND ~Scaling() {
#ifdef TMV_EXTRA_DEBUG
            x = Traits<T>::destr_value();
#endif
        }
        TMV_INLINE operator T() const { return x; }
    };
    // We enforce that T is real for 1 and -1 to avoid accidentally 
    // having two different template paths produce functionally 
    // identical code.  
    // Always use real T for 1 and -1.
    template <class T>
    struct Scaling<1,T> // ix = 1
    {
        TMV_INLINE Scaling() { TMVStaticAssert(Traits<T>::isreal); }
        TMV_INLINE_ND Scaling(const T _x) 
        { TMVStaticAssert(Traits<T>::isreal); TMVAssert(_x == T(1));  }
        TMV_INLINE operator T() const { return T(1); }
    };
    template <class T>
    struct Scaling<-1,T> // ix = -1
    {
        TMV_INLINE Scaling() { TMVStaticAssert(Traits<T>::isreal); }
        TMV_INLINE_ND Scaling(const T _x) 
        { TMVStaticAssert(Traits<T>::isreal); TMVAssert(_x == T(-1));  }
        TMV_INLINE operator T() const { return T(-1); }
    };

    template <int ix, class T>
    TMV_INLINE Scaling<ix,T> TMV_CONJ(const Scaling<ix,T>& x)
    { return x; }
    template <class T>
    TMV_INLINE Scaling<0,T> TMV_CONJ(const Scaling<0,T>& x)
    { return Scaling<0,T>(TMV_CONJ(x.x)); }

    // Define how a Scaling object multiplies a regular number:
    template <class T1, class T>
    TMV_INLINE T1& operator*=(T1& y, const Scaling<0,T>& x)
    { y *= x.x; return y; }
    template <class T1, class T>
    TMV_INLINE T1& operator*=(T1& y, const Scaling<1,T>& x)
    { return y; }
    template <class T1, class T>
    TMV_INLINE T1& operator*=(T1& y, const Scaling<-1,T>& x)
    { y = -y; return y; }
    template <class T1, class T>
    TMV_INLINE T1& operator/=(T1& y, const Scaling<0,T>& x)
    { y /= x.x; return y; }
    template <class T1, class T>
    TMV_INLINE T1& operator/=(T1& y, const Scaling<1,T>& x)
    { return y; }
    template <class T1, class T>
    TMV_INLINE T1& operator/=(T1& y, const Scaling<-1,T>& x)
    { y = -y; return y; }

    template <class T1, class T>
    TMV_INLINE typename Traits2<T,T1>::type operator*(
        const Scaling<0,T>& x, const T1& y)
    { return ZProd<false,false>::prod(x.x,y); }
    template <class T1, class T>
    TMV_INLINE T1 operator*(const Scaling<1,T>& x, const T1& y)
    { return y; }
    template <class T1, class T>
    TMV_INLINE T1 operator*(const Scaling<-1,T>& x, const T1& y)
    { return -y; }

    template <class T1, class T>
    TMV_INLINE typename Traits2<T,T1>::type operator*(
        const T1& y, const Scaling<0,T>& x)
    { return ZProd<false,false>::prod(x.x,y); }
    template <class T1, class T>
    TMV_INLINE T1 operator*(const T1& y, const Scaling<1,T>& x)
    { return y; }
    template <class T1, class T>
    TMV_INLINE T1 operator*(const T1& y, const Scaling<-1,T>& x)
    { return -y; }

    template <class T1, class T>
    TMV_INLINE typename Traits2<T,T1>::type operator/(
        const Scaling<0,T>& x, const T1& y)
    { return x.x / y; }
    template <class T1, class T>
    TMV_INLINE T1 operator/(const Scaling<1,T>& x, const T1& y)
    { return typename Traits<T1>::real_type(1) / y; }
    template <class T1, class T>
    TMV_INLINE T1 operator/(const Scaling<-1,T>& x, const T1& y)
    { return typename Traits<T1>::real_type(-1) / y; }

    template <class T1, class T>
    TMV_INLINE typename Traits2<T,T1>::type operator/(
        const T1& y, const Scaling<0,T>& x)
    { return y / x.x; }
    template <class T1, class T>
    TMV_INLINE T1 operator/(const T1& y, const Scaling<1,T>& x)
    { return y; }
    template <class T1, class T>
    TMV_INLINE T1 operator/(const T1& y, const Scaling<-1,T>& x)
    { return -y; }

    template <int ix, class T>
    TMV_INLINE Scaling<-ix,T> operator-(const Scaling<ix,T>& x)
    { return Scaling<-ix,T>(-T(x)); }

    // Scaling * Scaling:
    template <class T1, class T2>
    TMV_INLINE typename Traits2<T1,T2>::type operator*(
        const Scaling<0,T1>& x, const Scaling<0,T2>& y)
    { return ZProd<false,false>::prod(x.x,y.x); }
    template <class T1, class T2>
    TMV_INLINE T2 operator*(const Scaling<1,T1>& x, const Scaling<0,T2>& y)
    { return y.x; }
    template <class T1, class T2>
    TMV_INLINE T2 operator*(const Scaling<-1,T1>& x, const Scaling<0,T2>& y)
    { return -y.x; }

    template <class T1, class T2>
    TMV_INLINE T1 operator*(const Scaling<0,T1>& x, const Scaling<1,T2>& y)
    { return x; }
    template <class T1, class T2>
    TMV_INLINE Scaling<1,T2> operator*(
        const Scaling<1,T1>& x, const Scaling<1,T2>& y)
    { return Scaling<1,T2>(); }
    template <class T1, class T2>
    TMV_INLINE Scaling<-1,T2> operator*(
        const Scaling<-1,T1>& x, const Scaling<1,T2>& y)
    { return Scaling<-1,T2>(); }

    template <class T1, class T2>
    TMV_INLINE T1 operator*(const Scaling<0,T1>& x, const Scaling<-1,T2>& y)
    { return -x.x; }
    template <class T1, class T2>
    TMV_INLINE Scaling<-1,T2> operator*(
        const Scaling<1,T1>& x, const Scaling<-1,T2>& y)
    { return Scaling<-1,T2>(); }
    template <class T1, class T2>
    TMV_INLINE Scaling<1,T2> operator*(
        const Scaling<-1,T1>& x, const Scaling<-1,T2>& y)
    { return Scaling<1,T2>(); }


    // Scaling / Scaling:
    template <class T1, class T2>
    TMV_INLINE typename Traits2<T1,T2>::type operator/(
        const Scaling<0,T1>& x, const Scaling<0,T2>& y)
    { return x.x/y.x; }
    template <class T1, class T2>
    TMV_INLINE T2 operator/(const Scaling<1,T1>& x, const Scaling<0,T2>& y)
    { return T1(1)/y.x; }
    template <class T1, class T2>
    TMV_INLINE T2 operator/(const Scaling<-1,T1>& x, const Scaling<0,T2>& y)
    { return T1(-1)/y.x; }

    template <class T1, class T2>
    TMV_INLINE T1 operator/(const Scaling<0,T1>& x, const Scaling<1,T2>& y)
    { return x; }
    template <class T1, class T2>
    TMV_INLINE Scaling<1,T2> operator/(
        const Scaling<1,T1>& x, const Scaling<1,T2>& y)
    { return Scaling<1,T2>(); }
    template <class T1, class T2>
    TMV_INLINE Scaling<-1,T2> operator/(
        const Scaling<-1,T1>& x, const Scaling<1,T2>& y)
    { return Scaling<-1,T2>(); }

    template <class T1, class T2>
    TMV_INLINE T1 operator/(const Scaling<0,T1>& x, const Scaling<-1,T2>& y)
    { return -x.x; }
    template <class T1, class T2>
    TMV_INLINE Scaling<-1,T2> operator/(
        const Scaling<1,T1>& x, const Scaling<-1,T2>& y)
    { return Scaling<-1,T2>(); }
    template <class T1, class T2>
    TMV_INLINE Scaling<1,T2> operator/(
        const Scaling<-1,T1>& x, const Scaling<-1,T2>& y)
    { return Scaling<1,T2>(); }


    // Define + and - too:
    template <class T1, int ix, class T2>
    TMV_INLINE T1& operator+=(T1& y, const Scaling<ix,T2>& x)
    { y += T2(x); return y; }
    template <class T1, int ix, class T2>
    TMV_INLINE T1& operator-=(T1& y, const Scaling<ix,T2>& x)
    { y -= T2(x); return y; }

    template <int ix, class T1, class T2>
    TMV_INLINE typename Traits2<T1,T2>::type operator+(
        const Scaling<ix,T1>& x, const T2& y)
    { return T1(x) + y; }
    template <int ix, class T1, class T2>
    TMV_INLINE typename Traits2<T1,T2>::type operator-(
        const Scaling<ix,T1>& x, const T2& y)
    { return T1(x) - y; }

    template <class T1, int ix, class T2>
    TMV_INLINE typename Traits2<T1,T2>::type operator+(
        const T1& x, const Scaling<ix,T2>& y)
    { return x + T2(y); }
    template <class T1, int ix, class T2>
    TMV_INLINE typename Traits2<T1,T2>::type operator-(
        const T1& x, const Scaling<ix,T2>& y)
    { return x - T2(y); }

    template <int ix1, class T1, int ix2, class T2>
    TMV_INLINE typename Traits2<T1,T2>::type operator+(
        const Scaling<ix1,T1>& x, const Scaling<ix2,T2>& y)
    { return T1(x) + T2(y); }
    template <int ix1, class T1, int ix2, class T2>
    TMV_INLINE typename Traits2<T1,T2>::type operator-(
        const Scaling<ix1,T1>& x, const Scaling<ix2,T2>& y)
    { return T1(x) - T2(y); }

}

#endif
