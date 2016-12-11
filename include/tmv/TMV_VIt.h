///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// The Template Matrix/Vector Library for C++ was created by Mike Jarvis     //
// Copyright (C) 1998 - 2016                                                 //
// All rights reserved                                                       //
//                                                                           //
// The project is hosted at https://code.google.com/p/tmv-cpp/               //
// where you can find the current version and current documention.           //
//                                                                           //
// For concerns or problems with the software, Mike may be contacted at      //
// mike_jarvis17 [at] gmail.                                                 //
//                                                                           //
// This software is licensed under a FreeBSD license.  The file              //
// TMV_LICENSE should have bee included with this distribution.              //
// It not, you can get a copy from https://code.google.com/p/tmv-cpp/.       //
//                                                                           //
// Essentially, you can use this software however you want provided that     //
// you include the TMV_LICENSE file in any distribution that uses it.        //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


//-----------------------------------------------------------------------------
//
// This file defines all the iterator classes, as well as the classes
// to deal with conjugate vectors/matrices.
//
// The vector iterators are similar to the standard library's
// iterator and const_iterator types.
//
// VIt and CVIt are iterators along a mutable vector (or view)
// and constant vector respectively.
//

#ifndef TMV_VIt_H
#define TMV_VIt_H

#include "tmv/TMV_Base.h"

namespace tmv {

    template <typename T>
    class ConjRef; // Undefined unless T is complex<T>

    template <typename T>
    class ConjRef<std::complex<T> >
    {
    public:

        typedef std::complex<T> CT;

        explicit ConjRef(CT& _val) : val(_val) {}
        ConjRef(const ConjRef<CT>& rhs) : val(rhs.getRef()) {}
        ~ConjRef() {}

        operator CT() const { return std::conj(val); }
        CT& getRef() { return val; }
        CT conj() const { return val; }
        T real() const { return val.real(); }
        T imag() const { return -val.imag(); }
        CT operator-() const { return -std::conj(val); }

        ConjRef<CT>& operator=(const ConjRef<CT>& rhs)
        { val = rhs.getRef(); return *this; }
        ConjRef<CT>& operator=(CT rhs)
        { val = std::conj(rhs); return *this; }
        ConjRef<CT>& operator=(T rhs)
        { val = rhs; return *this; }

        ConjRef<CT>& operator+=(const ConjRef<CT>& x2)
        { val += x2.val; return *this; }
        ConjRef<CT>& operator+=(CT x2)
        { val += std::conj(x2); return *this; }
        ConjRef<CT>& operator+=(T x2)
        { val += x2; return *this; }
        CT operator+(const ConjRef<CT>& x2)
        { return std::conj(val+x2.val); }
        friend CT operator+(const ConjRef<CT>& x1, CT x2)
        { return std::conj(x1.val)+x2; }
        friend CT operator+(const ConjRef<CT>& x1, T x2)
        { return std::conj(x1.val)+x2; }
        friend CT operator+(CT x1, const ConjRef<CT>& x2)
        { return x1+std::conj(x2.val); }
        friend CT operator+(T x1, const ConjRef<CT>& x2)
        { return x1+std::conj(x2.val); }
        //friend CT& operator+=(CT& x1, const ConjRef<CT>& x2)
        //{ return x1+=std::conj(x2.val); }

        ConjRef<CT>& operator-=(const ConjRef<CT>& x2)
        { val -= x2.val; return *this; }
        ConjRef<CT>& operator-=(CT x2)
        { val -= std::conj(x2); return *this; }
        ConjRef<CT>& operator-=(T x2)
        { val -= x2; return *this; }
        CT operator-(const ConjRef<CT>& x2)
        { return std::conj(val-x2.val); }
        friend CT operator-(const ConjRef<CT>& x1, CT x2)
        { return std::conj(x1.val)-x2; }
        friend CT operator-(const ConjRef<CT>& x1, T x2)
        { return std::conj(x1.val)-x2; }
        friend CT operator-(CT x1, const ConjRef<CT>& x2)
        { return x1-std::conj(x2.val); }
        friend CT operator-(T x1, const ConjRef<CT>& x2)
        { return x1-std::conj(x2.val); }
        //friend CT& operator-=(CT& x1, const ConjRef<CT>& x2)
        //{ return x1-=std::conj(x2.val); }

        ConjRef<CT>& operator*=(const ConjRef<CT>& x2)
        { val *= x2.val; return *this; }
        ConjRef<CT>& operator*=(CT x2)
        { val *= std::conj(x2); return *this; }
        ConjRef<CT>& operator*=(T x2)
        { val *= x2; return *this; }
        CT operator*(const ConjRef<CT> x2)
        { return std::conj(val*x2.val); }
        friend CT operator*(const ConjRef<CT>& x1, CT x2)
        { return std::conj(x1.val)*x2; }
        friend CT operator*(const ConjRef<CT>& x1, T x2)
        { return std::conj(x1.val)*x2; }
        friend CT operator*(CT x1, const ConjRef<CT>& x2)
        { return x1*std::conj(x2.val); }
        friend CT operator*(T x1, const ConjRef<CT>& x2)
        { return x1*std::conj(x2.val); }
        //friend CT& operator*=(CT& x1, const ConjRef<CT>& x2)
        //{ return x1*=std::conj(x2.val); }

        ConjRef<CT>& operator/=(const ConjRef<CT>& x2)
        { val /= x2.val; return *this; }
        ConjRef<CT>& operator/=(CT x2)
        { val /= std::conj(x2); return *this; }
        ConjRef<CT>& operator/=(T x2)
        { val /= x2; return *this; }
        CT operator/(const ConjRef<CT>& x2)
        { return std::conj(val/x2.val); }
        friend CT operator/(const ConjRef<CT>& x1, CT x2)
        { return std::conj(x1.val)/x2; }
        friend CT operator/(const ConjRef<CT>& x1, T x2)
        { return std::conj(x1.val)/x2; }
        friend CT operator/(CT x1, const ConjRef<CT>& x2)
        { return x1/std::conj(x2.val); }
        friend CT operator/(T x1, const ConjRef<CT>& x2)
        { return x1/std::conj(x2.val); }
        //friend CT& operator/=(CT& x1, const ConjRef<CT>& x2)
        //{ return x1/=std::conj(x2.val); }

        bool operator==(const ConjRef<CT>& x2) const
        { return val == x2.val; }
        bool operator==(CT x2) const
        { return std::conj(val) == x2; }
        bool operator==(T x2) const
        { return std::real(val) == x2 && std::imag(val) == T(0); }
        friend bool operator==(CT x1, const ConjRef<CT>& x2)
        { return x2==x1; }
        friend bool operator==(T x1, const ConjRef<CT>& x2)
        { return x2==x1; }
        bool operator!=(const ConjRef<CT>& x2) const
        { return !(operator==(x2)); }
        bool operator!=(CT x2) const
        { return !(operator==(x2)); }
        bool operator!=(T x2) const
        { return !(operator==(x2)); }
        friend bool operator!=(CT x1, const ConjRef<CT>& x2)
        { return !(x2==x1); }
        friend bool operator!=(T x1, const ConjRef<CT>& x2)
        { return !(x2==x1); }

        void swapWith(CT& x2)
        {
            TMVAssert(&val != &x2);
            CT temp = x2; x2 = std::conj(val); val = std::conj(temp);
        }
        void swapWith(ConjRef<CT> x2)
        {
            TMVAssert(&val != &x2);
            CT temp = x2.val; x2.val = val; val = temp;
        }

        friend std::ostream& operator<<(std::ostream& os, ConjRef<CT> x)
        { os << std::conj(x.val); return os; }
        friend std::istream& operator>>(std::istream& is, ConjRef<CT> x)
        { is >> x.val; x.val = std::conj(x.val); return is; }

    private:

        CT& getRef() const { return val; }
        CT& val;
    };

    // A helper structure to define the reference type correctly for VIt.
    // Also applies a conjugation if necessary.
    template <typename T, ConjType C>
    struct AuxRef // real T or C = NonConj
    {
        typedef T& reference;
        static T apply(const T& x) { return x; }
    };
    template <typename T>
    struct AuxRef<std::complex<T>,Conj>
    {
        typedef ConjRef<std::complex<T> > reference;
        static std::complex<T> apply(const std::complex<T>& x)
        { return std::conj(x); }
    };

    template <ConjType C, typename T>
    inline T DoConj(const T& x) { return AuxRef<T,C>::apply(x); }


    template <typename T, int S, ConjType C>
    class VIt
    {
    public :

        typedef VIt<T,S,C> type;
        typedef std::random_access_iterator_tag iterator_category;
        typedef T value_type;
        typedef ptrdiff_t difference_type;
        typedef T* pointer;
        typedef typename AuxRef<T,C>::reference reference;

        VIt(T* inp, ptrdiff_t step) : p(inp), s(step) {}
        explicit VIt(T* inp) : p(inp), s(S)
        { TMVAssert(S != Unknown); }
        VIt(const type& rhs) : p(rhs.get()), s(rhs.step()) {}

        template <int S2>
        VIt(const VIt<T,S2,C>& rhs) : p(rhs.get()), s(rhs.step()) {}

        type& operator=(const type& rhs)
        { TMVAssert(step()==rhs.step()); p=rhs.get(); return *this; }

        template <int S2>
        type& operator=(const VIt<T,S2,C>& rhs)
        { TMVAssert(step()==rhs.step()); p=rhs.get(); return *this; }

        ~VIt() {}

        T* get() const { return p; }
        ptrdiff_t step() const { return s; }

        bool operator==(const type& rhs) const
        { return p == rhs.get(); }
        bool operator!=(const type& rhs) const
        { return p != rhs.get(); }
        bool operator<(const type& rhs) const
        { return (step()>0 ? p < rhs.get() : p > rhs.get()); }

        reference operator*() const { return reference(*p); }

        type& operator++() { p+=step(); return *this; }
        type& operator--() { p-=step(); return *this; }
        type operator++(int)
        { type p2 = *this; p+=step(); return p2; }
        type operator--(int)
        { type p2 = *this; p-=step(); return p2; }

        type& operator+=(ptrdiff_t n) { p += n*step(); return *this; }
        type& operator-=(ptrdiff_t n) { p -= n*step(); return *this; }
        type operator+(ptrdiff_t n) const
        { return type(p+n*step(),step()); }
        type operator-(ptrdiff_t n) const
        { return type(p-n*step(),step()); }
        type& shiftP(ptrdiff_t n) { p += n; return *this; }

        ptrdiff_t operator-(const type& rhs) const
        { return (p-rhs.get())/step(); }

        reference operator[](ptrdiff_t n) const
        { return reference(p[n*step()]); }

        typedef VIt<T,S,NonConj> nonconj_type;
        nonconj_type nonConj() const
        { return nonconj_type(p,step()); }
        typedef typename Traits<T>::real_type real_type;
        typedef VIt<real_type,S,NonConj> flatten_type;
        flatten_type flatten() const
        { return flatten_type(reinterpret_cast<real_type*>(p),1); }

    private :

        T* p;
        const CheckedInt<S> s;
    };

    template <typename T, int S, ConjType C>
    class CVIt
    {
    public :

        typedef CVIt<T,S,C> type;
        typedef std::random_access_iterator_tag iterator_category;
        typedef T value_type;
        typedef ptrdiff_t difference_type;
        typedef const T* pointer;
        typedef const T& reference;

        CVIt(const T* inp, ptrdiff_t step) : p(inp), s(step) {}
        explicit CVIt(const T* inp) : p(inp), s(S)
        { TMVAssert(S != Unknown); }
        CVIt(const type& rhs) : p(rhs.get()), s(rhs.step()) {}

        template <int S2>
        CVIt(const CVIt<T,S2,C>& rhs) :
            p(rhs.get()), s(rhs.step()) {}

        template <int S2>
        CVIt(const VIt<T,S2,C>& rhs) :
            p(rhs.get()), s(rhs.step()) {}

        type& operator=(const type& rhs)
        { TMVAssert(step()==rhs.step()); p = rhs.get(); return *this; }

        template <int S2>
        type& operator=(const CVIt<T,S2,C>& rhs)
        { TMVAssert(step()==rhs.step()); p = rhs.get(); return *this; }

        template <int S2>
        type& operator=(const VIt<T,S2,C>& rhs)
        { TMVAssert(step()==rhs.step()); p = rhs.get(); return *this; }

        ~CVIt() {}

        const T* get() const { return p; }
        ptrdiff_t step() const { return s; }

        bool operator==(const type& rhs) const
        { return p == rhs.get(); }
        bool operator!=(const type& rhs) const
        { return p != rhs.get(); }
        bool operator<(const type& rhs) const
        { return (step()>0 ? p < rhs.get() : p > rhs.get()); }

        T operator*() const { return DoConj<C>(*p); }

        type& operator++() { p+=step(); return *this; }
        type& operator--() { p-=step(); return *this; }
        type operator++(int)
        { type p2 = *this; p+=step(); return p2; }
        type operator--(int)
        { type p2 = *this; p-=step(); return p2; }

        type& operator+=(ptrdiff_t n) { p += n*step(); return *this; }
        type& operator-=(ptrdiff_t n) { p -= n*step(); return *this; }
        type operator+(ptrdiff_t n) const
        { return type(p+n*step(),step()); }
        type operator-(ptrdiff_t n) const
        { return type(p-n*step(),step()); }
        type& shiftP(ptrdiff_t n) { p += n; return *this; }

        ptrdiff_t operator-(const type& rhs) const
        { return (p-rhs.get())/step(); }

        T operator[](ptrdiff_t n) const
        { return DoConj<C>(p[n*step()]); }

        typedef CVIt<T,S,NonConj> nonconj_type;
        nonconj_type nonConj() const
        { return nonconj_type(p,step()); }
        typedef typename Traits<T>::real_type real_type;
        typedef CVIt<real_type,S,NonConj> flatten_type;
        flatten_type flatten() const
        { return flatten_type(reinterpret_cast<const real_type*>(p),1); }

    private :

        const T* p;
        const CheckedInt<S> s;
    };

    template <typename T, int S, ConjType C>
    inline CVIt<T,S,C> operator+(ptrdiff_t i, const CVIt<T,S,C>& it)
    { return it + i; }
    template <typename T, int S, ConjType C>
    inline VIt<T,S,C> operator+(ptrdiff_t i, const VIt<T,S,C>& it)
    { return it + i; }


    // Overload some functions to work with ConjRef<T>
    template <typename T>
    inline T TMV_CONJ(const ConjRef<T>& x) { return x.conj(); }
    template <typename T>
    inline typename Traits<T>::real_type TMV_NORM(const ConjRef<T>& x)
    { return TMV_NORM(x.conj()); }
    template <typename T>
    inline typename Traits<T>::real_type TMV_ABS(const ConjRef<T>& x)
    { return TMV_ABS(x.conj()); }
    template <typename T>
    inline T TMV_SQR(const ConjRef<T>& x)
    { return TMV_SQR(x.conj()); }
    template <typename T>
    inline T TMV_SQRT(const ConjRef<T>& x)
    { return TMV_SQRT(x.conj()); }
    template <typename T>
    inline typename Traits<T>::real_type TMV_REAL(const ConjRef<T>& x)
    { return x.real(); }
    template <typename T>
    inline typename Traits<T>::real_type TMV_IMAG(const ConjRef<T>& x)
    { return x.imag(); }

    template <typename T>
    inline void TMV_SWAP(
        tmv::ConjRef<std::complex<T> > x1, tmv::ConjRef<std::complex<T> > x2)
    { return x1.swapWith(x2); }
    template <typename T>
    inline void TMV_SWAP(
        std::complex<T>& x1, tmv::ConjRef<std::complex<T> > x2)
    { return x2.swapWith(x1); }
    template <typename T>
    inline void TMV_SWAP(
        tmv::ConjRef<std::complex<T> > x1, std::complex<T>& x2)
    { return x1.swapWith(x2); }


    // Now some special classes to deal with complex<T> when there is
    // a possibility (unknown at compile time) of the reference
    // or iterator needing to be conjugated:

    template <typename T>
    class VarConjRef; // Undefined unless T is complex<T>
    template <typename T>
    class VarConjIter; // Undefined unless T is complex<T>
    template <typename T>
    class CVarConjIter; // Undefined unless T is complex<T>

    template <typename T>
    class VarConjRef<std::complex<T> >
    {
    public:

        typedef std::complex<T> CT;

        VarConjRef(CT& _val, ConjType _c) : val(_val), c(_c) {}
        VarConjRef(const VarConjRef<CT>& rhs) : val(rhs.val), c(rhs.c) {}
        ~VarConjRef() {}

        inline bool isconj() const { return c == Conj; }
        inline operator CT() const { return c==Conj?std::conj(val):val; }
        inline CT& getRef() { return val; }
        inline CT conj() const { return c==Conj?val:std::conj(val); }
        inline T real() const { return val.real(); }
        inline T imag() const
        { return c==Conj? -val.imag() : val.imag(); }
        inline CT operator-() const { return -CT(*this); }

        VarConjRef<CT>& operator=(const VarConjRef<CT>& rhs)
        { val = c==rhs.c ? rhs.val : std::conj(rhs.val); return *this; }
        inline VarConjRef<CT>& operator=(CT rhs)
        { val = c==Conj ? std::conj(rhs) : rhs; return *this; }
        inline VarConjRef<CT>& operator=(T rhs)
        { val = rhs; return *this; }

        inline VarConjRef<CT>& operator+=(const VarConjRef<CT>& x2)
        { val += c==x2.c ? x2.val : std::conj(x2.val); return *this; }
        inline VarConjRef<CT>& operator+=(CT x2)
        { val += c==Conj ? std::conj(x2) : x2; return *this; }
        inline VarConjRef<CT>& operator+=(T x2)
        { val += x2; return *this; }
        inline CT operator+(const VarConjRef<CT>& x2)
        { return CT(*this)+CT(x2); }
        inline friend CT operator+(const VarConjRef<CT>& x1, CT x2)
        { return CT(x1)+x2; }
        inline friend CT operator+(const VarConjRef<CT>& x1, T x2)
        { return CT(x1)+x2; }
        inline friend CT operator+(CT x1, const VarConjRef<CT>& x2)
        { return x1+CT(x2); }
        inline friend CT operator+(T x1, const VarConjRef<CT>& x2)
        { return x1+CT(x2); }
        //inline friend CT& operator+=(CT& x1, const VarConjRef<CT>& x2)
        //{ return x1 += (x2.c==Conj ? std::conj(x2.val) : x2.val); }

        inline VarConjRef<CT>& operator-=(const VarConjRef<CT>& x2)
        { val -= c==x2.c ? x2.val : std::conj(x2.val); return *this; }
        inline VarConjRef<CT>& operator-=(CT x2)
        { val -= c==Conj ? std::conj(x2) : x2; return *this; }
        inline VarConjRef<CT>& operator-=(T x2)
        { val -= x2; return *this; }
        inline CT operator-(const VarConjRef<CT>& x2)
        { return CT(*this)-CT(x2); }
        inline friend CT operator-(const VarConjRef<CT>& x1, CT x2)
        { return CT(x1)-x2; }
        inline friend CT operator-(const VarConjRef<CT>& x1, T x2)
        { return CT(x1)-x2; }
        inline friend CT operator-(CT x1, const VarConjRef<CT>& x2)
        { return x1-CT(x2); }
        inline friend CT operator-(T x1, const VarConjRef<CT>& x2)
        { return x1-CT(x2); }
        //inline friend CT& operator-=(CT& x1, const VarConjRef<CT>& x2)
        //{ return x1 -= (x2.c==Conj ? std::conj(x2.val) : x2.val); }

        inline VarConjRef<CT>& operator*=(const VarConjRef<CT>& x2)
        { val *= c==x2.c ? x2.val : std::conj(x2.val); return *this; }
        inline VarConjRef<CT>& operator*=(CT x2)
        { val *= c==Conj ? std::conj(x2) : x2; return *this; }
        inline VarConjRef<CT>& operator*=(T x2)
        { val *= x2; return *this; }
        inline CT operator*(const VarConjRef<CT> x2)
        { return CT(*this)*CT(x2); }
        inline friend CT operator*(const VarConjRef<T>& x1, CT x2)
        { return CT(x1)*x2; }
        inline friend CT operator*(const VarConjRef<CT>& x1, T x2)
        { return CT(x1)*x2; }
        inline friend CT operator*(CT x1, const VarConjRef<CT>& x2)
        { return x1*CT(x2); }
        inline friend CT operator*(T x1, const VarConjRef<CT>& x2)
        { return x1*CT(x2); }
        //inline friend CT& operator*=(CT& x1, const VarConjRef<CT>& x2)
        //{ return x1 *= (x2.c==Conj ? std::conj(x2.val) : x2.val); }

        inline VarConjRef<CT>& operator/=(const VarConjRef<CT>& x2)
        { val /= c==x2.c ? x2.val : std::conj(x2.val); return *this; }
        inline VarConjRef<CT>& operator/=(CT x2)
        { val /= c==Conj ? std::conj(x2) : x2; return *this; }
        inline VarConjRef<CT>& operator/=(T x2)
        { val /= x2; return *this; }
        inline CT operator/(const VarConjRef<CT>& x2)
        { return CT(*this)/CT(x2); }
        inline friend CT operator/(const VarConjRef<CT>& x1, CT x2)
        { return CT(x1)/x2; }
        inline friend CT operator/(const VarConjRef<CT>& x1, T x2)
        { return CT(x1)/x2; }
        inline friend CT operator/(CT x1, const VarConjRef<CT>& x2)
        { return x1/CT(x2); }
        inline friend CT operator/(T x1, const VarConjRef<CT>& x2)
        { return x1/CT(x2); }
        //inline friend CT& operator/=(CT& x1, const VarConjRef<CT>& x2)
        //{ return x1 /= (x2.c==Conj ? std::conj(x2.val) : x2.val); }

        inline bool operator==(const VarConjRef<CT>& x2) const
        { return val == (c==x2.c ? x2.val : std::conj(x2.val)); }
        inline bool operator==(CT x2) const
        { return val == (c==Conj ? std::conj(x2) : x2); }
        inline bool operator==(T x2) const
        { return std::real(val) == x2 && std::imag(val) == CT(0); }
        inline friend bool operator==(CT x1, const VarConjRef<CT>& x2)
        { return x2==x1; }
        inline friend bool operator==(T x1, const VarConjRef<CT>& x2)
        { return x2==x1; }
        inline bool operator!=(const VarConjRef<CT>& x2) const
        { return !(operator==(x2)); }
        inline bool operator!=(CT x2) const
        { return !(operator==(x2)); }
        inline bool operator!=(T x2) const
        { return !(operator==(x2)); }
        inline friend bool operator!=(CT x1, const VarConjRef<CT>& x2)
        { return !(x2==x1); }
        inline friend bool operator!=(T x1, const VarConjRef<CT>& x2)
        { return !(x2==x1); }

        inline void swapWith(CT& x2)
        {
            if (&val == &x2) {
                TMVAssert(c != Conj);
            } else {
                CT temp = x2; x2 = CT(*this); *this = temp;
            }
        }
        inline void swapWith(VarConjRef<CT> x2)
        {
            if (&val == &x2.val) {
                TMVAssert(c == x2.c);
            } else {
                if (c==x2.c) {
                    CT temp = x2.val; x2.val = val; val = temp;
                }
                else {
                    CT temp = x2; x2 = CT(*this); *this = temp;
                }
            }
        }

        inline friend std::ostream& operator<<(
            std::ostream& os, VarConjRef<CT> x)
        { os << (x.c==Conj ? std::conj(x.val) : x.val); return os; }
        inline friend std::istream& operator>>(
            std::istream& is, VarConjRef<CT> x)
        { is >> x.val; if(x.c==Conj) x.val = std::conj(x.val); return is; }

    private:

        CT& val;
        ConjType c;
    };

    template <typename T> inline T TMV_CONJ(const VarConjRef<T>& x)
    { return x.conj(); }
    template <typename T> inline TMV_RealType(T) TMV_NORM(const VarConjRef<T>& x)
    { return TMV_NORM(T(x)); }
    template <typename T> inline TMV_RealType(T) TMV_ABS(const VarConjRef<T>& x)
    { return TMV_ABS(T(x)); }
    template <typename T> inline T TMV_SQR(const VarConjRef<T>& x)
    { return TMV_SQR(T(x)); }
    template <typename T> inline T TMV_SQRT(const VarConjRef<T>& x)
    { return TMV_SQRT(T(x)); }
    template <typename T> inline TMV_RealType(T) TMV_REAL(const VarConjRef<T>& x)
    { return x.real(); }
    template <typename T> inline TMV_RealType(T) TMV_IMAG(const VarConjRef<T>& x)
    { return x.imag(); }
    template <typename T> inline void TMV_SWAP(
        tmv::VarConjRef<std::complex<T> > x1,
        tmv::VarConjRef<std::complex<T> > x2)
    { return x1.swapWith(x2); }
    template <typename T> inline void TMV_SWAP(
        std::complex<T>& x1, tmv::VarConjRef<std::complex<T> > x2)
    { return x2.swapWith(x1); }
    template <typename T> inline void TMV_SWAP(
        tmv::VarConjRef<std::complex<T> > x1, std::complex<T>& x2)
    { return x1.swapWith(x2); }

    template <typename T>
    class VarConjIter<std::complex<T> >
    {
        typedef std::complex<T> CT;
    public :

        VarConjIter() : p(0), s(0), c(NonConj) TMV_DEFFIRSTLAST(0,0) {}
        VarConjIter(CT* inp, ptrdiff_t instep, ConjType inc
              TMV_PARAMFIRSTLAST(CT) ) :
            p(inp), s(instep), c(inc) TMV_DEFFIRSTLAST(_first,_last) {}
        VarConjIter(const VarConjIter<CT >& rhs) :
            p(rhs.p), s(rhs.s), c(rhs.c) TMV_DEFFIRSTLAST(rhs._first,rhs._last)
        {}
        VarConjIter<CT >& operator=(const VarConjIter<CT >& rhs)
        {
            TMVAssert(s==rhs.s);
            TMVAssert(c==rhs.c);
            p=rhs.p;
            TMV_SETFIRSTLAST(rhs._first,rhs._last);
            return *this;
        }
        ~VarConjIter() {}

        CT* getP() const { return p; }
        ptrdiff_t step() const { return s; }
        ConjType getC() const { return c; }

        template <int S, ConjType C>
        inline operator VIt<CT,S,C>() const { return VIt<CT,S,C>(getP(), step()); }
        template <int S, ConjType C>
        inline operator CVIt<CT,S,C>() const { return CVIt<CT,S,C>(getP(), step()); }

        inline bool operator==(const VarConjIter<CT >& rhs) const
        { return p == rhs.p; }
        inline bool operator!=(const VarConjIter<CT >& rhs) const
        { return p != rhs.p; }
        inline bool operator<(const VarConjIter<CT >& rhs) const
        { return (s > 0 ? p < rhs.p : p > rhs.p); }

        inline VarConjRef<CT > operator*() const
        {
#ifdef TMVFLDEBUG
            if (!(p>=_first && p<_last)) {
                std::cerr<<"p = "<<p<<std::endl;
                std::cerr<<"first,last = "<<_first<<"  "<<_last<<std::endl;
            }
            TMVAssert(p>=_first);
            TMVAssert(p<_last);
#endif
            return VarConjRef<CT >(*p,c);
        }

        inline VarConjIter<CT >& operator++() { p += s; return *this; }
        inline VarConjIter<CT >& operator--() { p -= s; return *this; }
        inline VarConjIter<CT > operator++(int)
        { VarConjIter<CT > p2 = *this; p+=s; return p2; }
        inline VarConjIter<CT > operator--(int)
        { VarConjIter<CT > p2 = *this; p-=s; return p2; }

        inline VarConjIter<CT >& operator+=(ptrdiff_t n)
        { if(s==1) ++p; else p += n*s; return *this; }
        inline VarConjIter<CT >& operator-=(ptrdiff_t n)
        { if(s==1) --p; else p -= n*s; return *this; }
        inline VarConjIter<CT > operator+(ptrdiff_t n) const
        { return VarConjIter<CT >(s==1?p+n:p+n*s,s,c TMV_FIRSTLAST ); }
        inline VarConjIter<CT > operator-(ptrdiff_t n) const
        { return VarConjIter<CT >(s==1?p-n:p-n*s,s,c TMV_FIRSTLAST ); }

        inline ptrdiff_t operator-(const VarConjIter<CT >& rhs) const
        {
            TMVAssert(rhs.c==c);
            TMVAssert(rhs.s==s);
            return s==1 ? p-rhs.p : (p-rhs.p)/s;
        }

        inline VarConjRef<CT > operator[](ptrdiff_t n) const
        {
#ifdef TMVFLDEBUG
            CT* pn = s==1 ? p+n : p+n*s;
            if (!(pn>=_first && pn<_last)) {
                std::cerr<<"pn = "<<pn<<std::endl;
                std::cerr<<"first,last = "<<_first<<"  "<<_last<<std::endl;
            }
            TMVAssert(pn >= _first);
            TMVAssert(pn < _last);
            return VarConjRef<CT >(*pn,c);
#else
            return VarConjRef<CT >(*(s==1 ? p+n : p+n*s),c);
#endif
        }

        typedef std::random_access_iterator_tag   iterator_category;
        typedef CT                                value_type;
        typedef ptrdiff_t                         difference_type;
        typedef CT*                               pointer;
        typedef VarConjRef<CT >                   reference;

    private :

        CT* p;
        const ptrdiff_t s;
        ConjType c;

#ifdef TMVFLDEBUG
    public :
        const CT* _first;
        const CT* _last;
#endif
    };

    template <typename T> class CVarConjIter<std::complex<T> >
    {
        typedef std::complex<T> CT;
    public :

        CVarConjIter() : p(0), s(0), c(NonConj) {}
        CVarConjIter(const CT* inp, ptrdiff_t instep, ConjType inc) :
            p(inp), s(instep), c(inc) {}
        CVarConjIter(const CVarConjIter<CT >& rhs) :
            p(rhs.p), s(rhs.s), c(rhs.c) {}
        CVarConjIter(const VarConjIter<CT >& rhs) :
            p(rhs.getP()), s(rhs.step()), c(rhs.getC()) {}
        CVarConjIter<CT >& operator=(const CVarConjIter<CT >& rhs)
        {
            TMVAssert(s==rhs.s);
            TMVAssert(c==rhs.c);
            p=rhs.p;
            return *this;
        }
        CVarConjIter<CT >& operator=(const VarConjIter<CT >& rhs)
        {
            TMVAssert(s==rhs.step());
            TMVAssert(c==rhs.getC());
            p=rhs.getP();
            return *this;
        }
        ~CVarConjIter() {}

        const CT* getP() const { return p; }
        ptrdiff_t step() const { return s; }
        ConjType getC() const { return c; }

        template <int S, ConjType C>
        inline operator CVIt<CT,S,C>() const { return CVIt<CT,S,C>(getP(), step()); }

        inline bool operator==(const CVarConjIter<CT >& rhs) const
        { return p == rhs.p; }
        inline bool operator!=(const CVarConjIter<CT >& rhs) const
        { return p != rhs.p; }
        inline bool operator<(const CVarConjIter<CT >& rhs) const
        { return (s > 0 ? p < rhs.p : p > rhs.p); }

        inline CT operator*() const
        { return c==Conj ? std::conj(*p) : *p; }

        inline CVarConjIter<CT >& operator++() { p += s; return *this; }
        inline CVarConjIter<CT >& operator--() { p -= s; return *this; }
        inline CVarConjIter<CT > operator++(int)
        { CVarConjIter<CT > p2 = *this; p+=s; return p2; }
        inline CVarConjIter<CT > operator--(int)
        { CVarConjIter<CT > p2 = *this; p-=s; return p2; }

        inline CVarConjIter<CT >& operator+=(ptrdiff_t n)
        { if(s==1) ++p; else p += n*s; return *this; }
        inline CVarConjIter<CT >& operator-=(ptrdiff_t n)
        { if(s==1) --p; else p -= n*s; return *this; }
        inline CVarConjIter<CT > operator+(ptrdiff_t n) const
        { return CVarConjIter<CT >(s==1?p+n:p+n*s,s,c); }
        inline CVarConjIter<CT > operator-(ptrdiff_t n) const
        { return CVarConjIter<CT >(s==1?p-n:p-n*s,s,c); }

        inline ptrdiff_t operator-(const CVarConjIter<CT >& rhs) const
        {
            TMVAssert(rhs.c==c);
            TMVAssert(rhs.s==s);
            return s==1 ? p-rhs.p : (p-rhs.p)/s;
        }

        inline CT operator[](ptrdiff_t n) const
        { return c==Conj ? std::conj(*(p+(s==1?n:n*s))) : *(p+(s==1?n:n*s)); }

        typedef std::random_access_iterator_tag   iterator_category;
        typedef CT                                value_type;
        typedef ptrdiff_t                         difference_type;
        typedef const CT*                         pointer;
        typedef const VarConjRef<CT >             reference;

    private :

        const CT* p;
        const ptrdiff_t s;
        ConjType c;
    };

    template <typename T>
    struct RefHelper
    {
        typedef T& reference;
        typedef VIt<T,Unknown,NonConj> iterator;
        typedef CVIt<T,Unknown,NonConj> const_iterator;
        static reference makeRef(T* p, ConjType ) { return *p; }
        static iterator makeIter(T* p, ptrdiff_t s, ConjType )
        { return iterator(p,s); }
        static const_iterator makeIter(const T* p, ptrdiff_t s, ConjType )
        { return const_iterator(p,s); }
    };

    template <class RT>
    struct RefHelper<std::complex<RT> >
    {
        typedef std::complex<RT> T;
        typedef VarConjRef<T> reference;
        typedef VarConjIter<T>  iterator;
        typedef CVarConjIter<T> const_iterator;
        static reference makeRef(T* p, ConjType c) { return reference(*p,c); }
        static iterator makeIter(T* p, ptrdiff_t s, ConjType c)
        { return iterator(p,s,c); }
        static const_iterator makeIter(const T* p, ptrdiff_t s, ConjType c)
        { return const_iterator(p,s,c); }
    };

    template <typename T>
    inline std::string TMV_Text(ConjRef<T>)
    { return std::string("ConjRef<") + TMV_Text(T()) + ">"; }

    template <typename T, int S, ConjType C>
    inline std::string TMV_Text(VIt<T,S,C> it)
    {
        std::ostringstream s;
        s << "VIt<" << TMV_Text(T())<<",";
        if (S == Unknown) s << "Unknown ("<<it.step()<<")";
        else s << S;
        s << ","<< TMV_Text(C) << ">";
        return s.str();
    }

    template <typename T, int S, ConjType C>
    inline std::string TMV_Text(CVIt<T,S,C> it)
    {
        std::ostringstream s;
        s << "CVIt<" << TMV_Text(T())<<",";
        if (S == Unknown) s << "Unknown ("<<it.step()<<")";
        else s << S;
        s << ","<< TMV_Text(C) << ">";
        return s.str();
    }

    template <typename T>
    inline std::string TMV_Text(VarConjRef<T>)
    { return std::string("VarConjRef<") + TMV_Text(T()) + ">"; }

    template <typename T>
    inline std::string TMV_Text(VarConjIter<T> it)
    {
        return std::string("VarConjIter<") + TMV_Text(T()) + "," +
            it.step() + "," + TMV_Text(it.getC()) + ">";
    }

    template <typename T>
    inline std::string TMV_Text(CVarConjIter<T> it)
    {
        return std::string("CVarConjIter<") + TMV_Text(T()) + "," +
            it.step() + "," + TMV_Text(it.getC()) + ">";
    }

} // namespace tmv

#endif
