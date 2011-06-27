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

#include "TMV_Base.h"

namespace tmv {

    template <class T>
    class ConjRef; // Undefined unless T is complex<T>

    template <class T>
    class ConjRef<std::complex<T> >
    {
    public:

        typedef std::complex<T> CT;

        TMV_INLINE explicit ConjRef(CT& _val) : val(_val) {}
        TMV_INLINE ConjRef(const ConjRef<CT>& rhs) : val(rhs.getRef()) {}
        TMV_INLINE ~ConjRef() {}

        TMV_INLINE operator CT() const { return std::conj(val); }
        TMV_INLINE CT& getRef() { return val; }
        TMV_INLINE CT conj() const { return val; }
        TMV_INLINE T real() const { return val.real(); }
        TMV_INLINE T imag() const { return -val.imag(); }
        TMV_INLINE CT operator-() const { return -std::conj(val); }

        TMV_INLINE ConjRef<CT>& operator=(const ConjRef<CT>& rhs)
        { val = rhs.getRef(); return *this; }
        TMV_INLINE ConjRef<CT>& operator=(CT rhs)
        { val = std::conj(rhs); return *this; }
        TMV_INLINE ConjRef<CT>& operator=(T rhs)
        { val = rhs; return *this; }

        TMV_INLINE ConjRef<CT>& operator+=(const ConjRef<CT>& x2)
        { val += x2.val; return *this; }
        TMV_INLINE ConjRef<CT>& operator+=(CT x2)
        { val += std::conj(x2); return *this; }
        TMV_INLINE ConjRef<CT>& operator+=(T x2)
        { val += x2; return *this; }
        TMV_INLINE CT operator+(const ConjRef<CT>& x2)
        { return std::conj(val+x2.val); }
        TMV_INLINE friend CT operator+(const ConjRef<CT>& x1, CT x2)
        { return std::conj(x1.val)+x2; }
        TMV_INLINE friend CT operator+(const ConjRef<CT>& x1, T x2)
        { return std::conj(x1.val)+x2; }
        TMV_INLINE friend CT operator+(CT x1, const ConjRef<CT>& x2)
        { return x1+std::conj(x2.val); }
        TMV_INLINE friend CT operator+(T x1, const ConjRef<CT>& x2)
        { return x1+std::conj(x2.val); }
        //friend CT& operator+=(CT& x1, const ConjRef<CT>& x2)
        //{ return x1+=std::conj(x2.val); }

        TMV_INLINE ConjRef<CT>& operator-=(const ConjRef<CT>& x2) 
        { val -= x2.val; return *this; }
        TMV_INLINE ConjRef<CT>& operator-=(CT x2) 
        { val -= std::conj(x2); return *this; }
        TMV_INLINE ConjRef<CT>& operator-=(T x2) 
        { val -= x2; return *this; }
        TMV_INLINE CT operator-(const ConjRef<CT>& x2)
        { return std::conj(val-x2.val); }
        TMV_INLINE friend CT operator-(const ConjRef<CT>& x1, CT x2)
        { return std::conj(x1.val)-x2; }
        TMV_INLINE friend CT operator-(const ConjRef<CT>& x1, T x2)
        { return std::conj(x1.val)-x2; }
        TMV_INLINE friend CT operator-(CT x1, const ConjRef<CT>& x2)
        { return x1-std::conj(x2.val); }
        TMV_INLINE friend CT operator-(T x1, const ConjRef<CT>& x2)
        { return x1-std::conj(x2.val); }
        //friend CT& operator-=(CT& x1, const ConjRef<CT>& x2)
        //{ return x1-=std::conj(x2.val); }

        TMV_INLINE ConjRef<CT>& operator*=(const ConjRef<CT>& x2) 
        { val *= x2.val; return *this; }
        TMV_INLINE ConjRef<CT>& operator*=(CT x2) 
        { val *= std::conj(x2); return *this; }
        TMV_INLINE ConjRef<CT>& operator*=(T x2) 
        { val *= x2; return *this; }
        TMV_INLINE CT operator*(const ConjRef<CT> x2)
        { return std::conj(val*x2.val); }
        TMV_INLINE friend CT operator*(const ConjRef<CT>& x1, CT x2)
        { return std::conj(x1.val)*x2; }
        TMV_INLINE friend CT operator*(const ConjRef<CT>& x1, T x2)
        { return std::conj(x1.val)*x2; }
        TMV_INLINE friend CT operator*(CT x1, const ConjRef<CT>& x2)
        { return x1*std::conj(x2.val); }
        TMV_INLINE friend CT operator*(T x1, const ConjRef<CT>& x2)
        { return x1*std::conj(x2.val); }
        //friend CT& operator*=(CT& x1, const ConjRef<CT>& x2)
        //{ return x1*=std::conj(x2.val); }

        TMV_INLINE ConjRef<CT>& operator/=(const ConjRef<CT>& x2) 
        { val /= x2.val; return *this; }
        TMV_INLINE ConjRef<CT>& operator/=(CT x2) 
        { val /= std::conj(x2); return *this; }
        TMV_INLINE ConjRef<CT>& operator/=(T x2) 
        { val /= x2; return *this; }
        TMV_INLINE CT operator/(const ConjRef<CT>& x2)
        { return std::conj(val/x2.val); }
        TMV_INLINE friend CT operator/(const ConjRef<CT>& x1, CT x2)
        { return std::conj(x1.val)/x2; }
        TMV_INLINE friend CT operator/(const ConjRef<CT>& x1, T x2)
        { return std::conj(x1.val)/x2; }
        TMV_INLINE friend CT operator/(CT x1, const ConjRef<CT>& x2)
        { return x1/std::conj(x2.val); }
        TMV_INLINE friend CT operator/(T x1, const ConjRef<CT>& x2)
        { return x1/std::conj(x2.val); }
        //friend CT& operator/=(CT& x1, const ConjRef<CT>& x2)
        //{ return x1/=std::conj(x2.val); }

        TMV_INLINE bool operator==(const ConjRef<CT>& x2) const
        { return val == x2.val; }
        TMV_INLINE bool operator==(CT x2) const 
        { return std::conj(val) == x2; }
        TMV_INLINE bool operator==(T x2) const 
        { return std::real(val) == x2 && std::imag(val) == T(0); }
        TMV_INLINE friend bool operator==(CT x1, const ConjRef<CT>& x2)
        { return x2==x1; }
        TMV_INLINE friend bool operator==(T x1, const ConjRef<CT>& x2)
        { return x2==x1; }
        TMV_INLINE bool operator!=(const ConjRef<CT>& x2) const
        { return !(operator==(x2)); }
        TMV_INLINE bool operator!=(CT x2) const 
        { return !(operator==(x2)); }
        TMV_INLINE bool operator!=(T x2) const 
        { return !(operator==(x2)); }
        TMV_INLINE friend bool operator!=(CT x1, const ConjRef<CT>& x2)
        { return !(x2==x1); }
        TMV_INLINE friend bool operator!=(T x1, const ConjRef<CT>& x2)
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

        TMV_INLINE CT& getRef() const { return val; }
        CT& val;
    };

    // A helper struct to determine the reference type of VIt, VectorView
    template <class T, bool isconj>
    struct AuxRef
    { typedef T& reference; };
    template <class T>
    struct AuxRef<std::complex<T>,true>
    { typedef ConjRef<std::complex<T> > reference; };
    template <class T>
    struct AuxRef<std::complex<T>,false>
    { typedef std::complex<T>& reference; };

    template <class T, int S, bool C>
    class VIt 
    {
    public :

        typedef VIt<T,S,C> type;
        typedef std::random_access_iterator_tag iterator_category;
        typedef T value_type;
        typedef ptrdiff_t difference_type;
        typedef T* pointer;
        typedef typename AuxRef<T,C>::reference reference;

        TMV_INLINE VIt(T* inp, int step) : p(inp), s(step) {}
        TMV_INLINE explicit VIt(T* inp) : p(inp), s(S) 
        { TMVAssert(S != UNKNOWN); }
        TMV_INLINE VIt(const type& rhs) : p(rhs.get()), s(rhs.step()) {}

        template <int S2>
        TMV_INLINE VIt(const VIt<T,S2,C>& rhs) : p(rhs.get()), s(rhs.step()) {}

        TMV_INLINE type& operator=(const type& rhs) 
        { TMVAssert(step()==rhs.step()); p=rhs.get(); return *this; }

        template <int S2>
        TMV_INLINE type& operator=(const VIt<T,S2,C>& rhs) 
        { TMVAssert(step()==rhs.step()); p=rhs.get(); return *this; }

        TMV_INLINE ~VIt() {}

        TMV_INLINE T* get() const { return p; }
        TMV_INLINE int step() const { return s; }

        TMV_INLINE bool operator==(const type& rhs) const 
        { return p == rhs.get(); }
        TMV_INLINE bool operator!=(const type& rhs) const 
        { return p != rhs.get(); }
        TMV_INLINE bool operator<(const type& rhs) const 
        { return (step()>0 ? p < rhs.get() : p > rhs.get()); }

        TMV_INLINE reference operator*() const { return reference(*p); }

        TMV_INLINE type& operator++() { p+=step(); return *this; }
        TMV_INLINE type& operator--() { p-=step(); return *this; }
        TMV_INLINE type operator++(int) 
        { type p2 = *this; p+=step(); return p2; }
        TMV_INLINE type operator--(int) 
        { type p2 = *this; p-=step(); return p2; }

        TMV_INLINE type& operator+=(int n) { p += n*step(); return *this; }
        TMV_INLINE type& operator-=(int n) { p -= n*step(); return *this; }
        TMV_INLINE type operator+(int n) const 
        { return type(p+n*step(),step()); }
        TMV_INLINE type operator-(int n) const 
        { return type(p-n*step(),step()); }
        TMV_INLINE type& shiftP(int n) { p += n; return *this; }

        TMV_INLINE ptrdiff_t operator-(const type& rhs) const 
        { return (p-rhs.get())/step(); }

        TMV_INLINE reference operator[](int n) const 
        { return reference(p[n*step()]); }

        typedef VIt<T,S,false> nonconj_type;
        TMV_INLINE nonconj_type nonConj() const 
        { return nonconj_type(p,step()); }
        typedef typename Traits<T>::real_type real_type;
        typedef VIt<real_type,S,false> flatten_type;
        TMV_INLINE flatten_type flatten() const 
        { return flatten_type(reinterpret_cast<real_type*>(p),1); }

    private :

        T* p;
        const CheckedInt<S> s;
    };

    template <class T, int S, bool C>
    class CVIt
    {
    public :

        typedef CVIt<T,S,C> type;
        typedef std::random_access_iterator_tag iterator_category;
        typedef T value_type;
        typedef ptrdiff_t difference_type;
        typedef const T* pointer;
        typedef const T& reference;

        TMV_INLINE CVIt(const T* inp, int step) : p(inp), s(step) {}
        TMV_INLINE explicit CVIt(const T* inp) : p(inp), s(S) 
        { TMVAssert(S != UNKNOWN); }
        TMV_INLINE CVIt(const type& rhs) : p(rhs.get()), s(rhs.step()) {}

        template <int S2>
        TMV_INLINE CVIt(const CVIt<T,S2,C>& rhs) : 
            p(rhs.get()), s(rhs.step()) {}

        template <int S2>
        TMV_INLINE CVIt(const VIt<T,S2,C>& rhs) : 
            p(rhs.get()), s(rhs.step()) {}

        TMV_INLINE type& operator=(const type& rhs)
        { TMVAssert(step()==rhs.step()); p = rhs.get(); return *this; }

        template <int S2>
        TMV_INLINE type& operator=(const CVIt<T,S2,C>& rhs)
        { TMVAssert(step()==rhs.step()); p = rhs.get(); return *this; }

        template <int S2>
        TMV_INLINE type& operator=(const VIt<T,S2,C>& rhs)
        { TMVAssert(step()==rhs.step()); p = rhs.get(); return *this; }

        TMV_INLINE ~CVIt() {}

        TMV_INLINE const T* get() const { return p; }
        TMV_INLINE int step() const { return s; }

        TMV_INLINE bool operator==(const type& rhs) const 
        { return p == rhs.get(); }
        TMV_INLINE bool operator!=(const type& rhs) const 
        { return p != rhs.get(); }
        TMV_INLINE bool operator<(const type& rhs) const 
        { return (step()>0 ? p < rhs.get() : p > rhs.get()); }

        TMV_INLINE T operator*() const { return DoConj<C>(*p); }

        TMV_INLINE type& operator++() { p+=step(); return *this; }
        TMV_INLINE type& operator--() { p-=step(); return *this; }
        TMV_INLINE type operator++(int) 
        { type p2 = *this; p+=step(); return p2; }
        TMV_INLINE type operator--(int) 
        { type p2 = *this; p-=step(); return p2; }

        TMV_INLINE type& operator+=(int n) { p += n*step(); return *this; }
        TMV_INLINE type& operator-=(int n) { p -= n*step(); return *this; }
        TMV_INLINE type operator+(int n) const 
        { return type(p+n*step(),step()); }
        TMV_INLINE type operator-(int n) const 
        { return type(p-n*step(),step()); }
        TMV_INLINE type& shiftP(int n) { p += n; return *this; }

        TMV_INLINE ptrdiff_t operator-(const type& rhs) const 
        { return (p-rhs.get())/step(); }

        TMV_INLINE T operator[](int n) const 
        { return DoConj<C>(p[n*step()]); }

        typedef CVIt<T,S,false> nonconj_type;
        TMV_INLINE nonconj_type nonConj() const 
        { return nonconj_type(p,step()); }
        typedef typename Traits<T>::real_type real_type;
        typedef CVIt<real_type,S,false> flatten_type;
        TMV_INLINE flatten_type flatten() const 
        { return flatten_type(reinterpret_cast<const real_type*>(p),1); }

    private :

        const T* p;
        const CheckedInt<S> s;
    };

    template <class T, int S, bool C>
    static TMV_INLINE CVIt<T,S,C> operator+(int i, const CVIt<T,S,C>& it)
    { return it + i; }
    template <class T, int S, bool C>
    static TMV_INLINE VIt<T,S,C> operator+(int i, const VIt<T,S,C>& it)
    { return it + i; }


    // Overload some functions to work with ConjRef<T>
    template <class T>
    static TMV_INLINE T TMV_CONJ(const ConjRef<T>& x) { return x.conj(); }
    template <class T>
    static TMV_INLINE typename Traits<T>::real_type TMV_NORM(
        const ConjRef<T>& x) 
    { return TMV_NORM(x.conj()); }
    template <class T>
    static TMV_INLINE typename Traits<T>::real_type TMV_ABS(
        const ConjRef<T>& x) 
    { return TMV_ABS(x.conj()); }
    template <class T>
    static TMV_INLINE T TMV_SQR(const ConjRef<T>& x) 
    { return TMV_SQR(x.conj()); }
    template <class T>
    static TMV_INLINE T TMV_SQRT(const ConjRef<T>& x) 
    { return TMV_SQRT(x.conj()); }
    template <class T>
    static TMV_INLINE typename Traits<T>::real_type TMV_REAL(
        const ConjRef<T>& x) 
    { return x.real(); }
    template <class T>
    static TMV_INLINE typename Traits<T>::real_type TMV_IMAG(
        const ConjRef<T>& x) 
    { return x.imag(); }

    template <class T>
    static TMV_INLINE void TMV_SWAP(
        tmv::ConjRef<std::complex<T> > x1, tmv::ConjRef<std::complex<T> > x2)
    { return x1.swapWith(x2); }
    template <class T>
    static TMV_INLINE void TMV_SWAP(
        std::complex<T>& x1, tmv::ConjRef<std::complex<T> > x2)
    { return x2.swapWith(x1); }
    template <class T>
    static TMV_INLINE void TMV_SWAP(
        tmv::ConjRef<std::complex<T> > x1, std::complex<T>& x2)
    { return x1.swapWith(x2); }

    template <class T1, class T2>
    struct Traits2<T1,ConjRef<T2> >
    {
        enum { sametype = Traits2<T1,T2>::sametype };
        enum { samebase = Traits2<T1,T2>::samebase };
        typedef typename Traits2<T1,T2>::type type;
    };
    template <class T1, class T2>
    struct Traits2<ConjRef<T1>,T2>
    {
        enum { sametype = Traits2<T1,T2>::sametype };
        enum { samebase = Traits2<T1,T2>::samebase };
        typedef typename Traits2<T1,T2>::type type;
    };
    template <class T1, class T2>
    struct Traits2<ConjRef<T1>,ConjRef<T2> >
    {
        enum { sametype = Traits2<T1,T2>::sametype };
        enum { samebase = Traits2<T1,T2>::samebase };
        typedef typename Traits2<T1,T2>::type type;
    };

#ifdef TMV_TEXT
    template <class T>
    static inline std::string TMV_Text(ConjRef<T>)
    { return std::string("ConjRef<") + TMV_Text(T()) + ">"; }

    template <class T, int S, bool C>
    static inline std::string TMV_Text(VIt<T,S,C> it)
    {
        std::ostringstream s;
        s << "VIt<" << TMV_Text(T());
        s << ","<<IntTraits<S>::text();
        if (S == UNKNOWN) s << "("<<it.step()<<")";
        s << ","<< C << ">";
        return s.str();
    }

    template <class T, int S, bool C>
    static inline std::string TMV_Text(CVIt<T,S,C> it)
    {
        std::ostringstream s;
        s << "CVIt<" << TMV_Text(T());
        s << ","<<IntTraits<S>::text();
        if (S == UNKNOWN) s << "("<<it.step()<<")";
        s << ","<< C << ">";
        return s.str();
    }
#endif

} // namespace tmv

#endif
