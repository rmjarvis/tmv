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

#include "tmv/TMV_Base.h"

namespace tmv {

    // Start with real version so Conj isn't an issue.
    // Specialize for complex below.
    template <class T> class VIter
    { 
    public :

        VIter() : p(0), s(0) TMV_DEFFIRSTLAST(0,0) {}
        VIter(T* inp, int instep, ConjType 
              TMV_DEBUGPARAM(inc) TMV_PARAMFIRSTLAST(T) ) : 
            p(inp), s(instep) TMV_DEFFIRSTLAST(_first,_last) 
        { TMVAssert(inc==NonConj); }
        VIter(const VIter<T>& rhs) : 
            p(rhs.p), s(rhs.s) TMV_DEFFIRSTLAST(rhs._first,rhs._last) {}
        VIter<T>& operator=(const VIter<T>& rhs) 
        { 
            TMVAssert(s == rhs.s);
            p=rhs.p; 
            TMV_SETFIRSTLAST(rhs._first,rhs._last);
            return *this; 
        }
        ~VIter() {}

        T* getP() const { return p; }
        int step() const { return s; }
        ConjType getC() const { return NonConj; }

        inline bool operator==(const VIter<T>& rhs) const 
        { return p == rhs.p; }
        inline bool operator!=(const VIter<T>& rhs) const 
        { return p != rhs.p; }
        inline bool operator<(const VIter<T>& rhs) const 
        { return (s > 0 ? p < rhs.p : p > rhs.p); }

        inline T& operator*() const
        { 
#ifdef TMVFLDEBUG
            if (!(p>=_first && p<_last)) {
                std::cerr<<"p = "<<p<<std::endl;
                std::cerr<<"first,last = "<<_first<<"  "<<_last<<std::endl;
            }
            TMVAssert(p>=_first);
            TMVAssert(p<_last);
#endif
            return *p; 
        }

        inline VIter<T>& operator++() { p += s; return *this; }
        inline VIter<T>& operator--() { p -= s; return *this; }
        inline VIter<T> operator++(int) 
        { VIter<T> p2 = *this; p+=s; return p2; }
        inline VIter<T> operator--(int) 
        { VIter<T> p2 = *this; p-=s; return p2; }

        inline VIter<T>& operator+=(int n) { p += n*s; return *this; }
        inline VIter<T>& operator-=(int n) { p -= n*s; return *this; }
        inline VIter<T> operator+(int n) const 
        { return VIter<T>(p+n*s,s,NonConj TMV_FIRSTLAST ); }
        inline VIter<T> operator-(int n) const 
        { return VIter<T>(p-n*s,s,NonConj TMV_FIRSTLAST ); }

        inline ptrdiff_t operator-(const VIter<T>& rhs) const 
        { return (p-rhs.p)/s; }

        inline T& operator[](int n) const
        {
#ifdef TMVFLDEBUG
            T* pn = p+n*s;
            if (!(p>=_first && p<_last)) {
                std::cerr<<"p = "<<p<<std::endl;
                std::cerr<<"first,last = "<<_first<<"  "<<_last<<std::endl;
            }
            TMVAssert(pn >= _first);
            TMVAssert(pn < _last);
            return *pn; 
#else
            return *(p+n*s); 
#endif
        }

        typedef std::random_access_iterator_tag iterator_category;
        typedef T                               value_type;
        typedef ptrdiff_t                       difference_type;
        typedef T*                              pointer;
        typedef T&                              reference;

    private :

        T* p;
        const int s;

#ifdef TMVFLDEBUG
    public :
        const T* _first;
        const T* _last;
#endif

    };

    template <class T> class CVIter
    { 
    public :

        CVIter() : p(0), s(0) {}
        CVIter(const T* inp, int instep, ConjType TMV_DEBUGPARAM(inc)) : 
            p(inp), s(instep) 
        { TMVAssert(inc==NonConj); }
        CVIter(const CVIter<T>& rhs) : p(rhs.p), s(rhs.s) {}
        CVIter(const VIter<T>& rhs) : p(rhs.getP()), s(rhs.step()) {}
        CVIter<T>& operator=(const CVIter<T>& rhs) 
        { 
            TMVAssert(s == rhs.s);
            p=rhs.p; 
            return *this; 
        }
        ~CVIter() {}

        const T* getP() const { return p; }
        int step() const { return s; }
        ConjType getC() const { return NonConj; }

        inline bool operator==(const CVIter<T>& rhs) const 
        { return p == rhs.p; }
        inline bool operator!=(const CVIter<T>& rhs) const 
        { return p != rhs.p; }
        inline bool operator<(const CVIter<T>& rhs) const 
        { return (s > 0 ? p < rhs.p : p > rhs.p); }

        inline T operator*() const { return *p; }

        inline CVIter<T>& operator++() { p += s; return *this; }
        inline CVIter<T>& operator--() { p -= s; return *this; }
        inline CVIter<T> operator++(int) 
        { CVIter<T> p2 = *this; p+=s; return p2; }
        inline CVIter<T> operator--(int) 
        { CVIter<T> p2 = *this; p-=s; return p2; }

        inline CVIter<T>& operator+=(int n) { p += n*s; return *this; }
        inline CVIter<T>& operator-=(int n) { p -= n*s; return *this; }
        inline CVIter<T> operator+(int n) const 
        { return CVIter<T>(p+n*s,s,NonConj); }
        inline CVIter<T> operator-(int n) const 
        { return CVIter<T>(p-n*s,s,NonConj); }

        inline ptrdiff_t operator-(const CVIter<T>& rhs) const 
        { return (p-rhs.p)/s; }

        inline T operator[](int n) const { return *(p+n*s); }

        typedef std::random_access_iterator_tag iterator_category;
        typedef T                               value_type;
        typedef ptrdiff_t                       difference_type;
        typedef const T*                        pointer;
        typedef const T&                        reference;

    private :

        const T* p;
        const int s;

    };

    // These (C)VIt, rather than (C)VIter have the StepType and 
    // ConjType specified as templates.  When these can be known,
    // it is faster to let the compiler know which variety
    // of iterator you are dealing with.  
    // For example, if the step is 1, p++ is much faster than p+=s
    // (where s==1), so the S=Unit iterator does this.
    // Likewise when you know you do not have a Conjugate iterator, 
    // you can specify NonConj for the ConjType and it doesn't have
    // to deal with the complications of the VarConjRef type.

    // The general template specification is really that of
    // T = real, S = Unit, C = irrelevant
    // T = complex and S = Step are specialized below.
    template <class T, StepType S, ConjType C> class VIt 
    { 
    public :

        VIt() : p(0) TMV_DEFFIRSTLAST(0,0) {}

        VIt(T* inp, int TMV_DEBUGPARAM(step) TMV_PARAMFIRSTLAST(T) ) : 
            p(inp) TMV_DEFFIRSTLAST(_first,_last)
        { TMVAssert(step==1); }

        VIt(const VIt<T,S,C>& rhs) : 
            p(rhs.p)  TMV_DEFFIRSTLAST(rhs._first,rhs._last) {}

        VIt(const VIter<T>& rhs) : 
            p(rhs.getP())  TMV_DEFFIRSTLAST(rhs._first,rhs._last) 
        { TMVAssert(rhs.step()==1); }

        template <StepType S2> VIt(const VIt<T,S2,C>& rhs) : 
            p(rhs.getP())  TMV_DEFFIRSTLAST(rhs._first,rhs._last) 
        { TMVAssert(rhs.step()==1); }

        VIt<T,S,C>& operator=(const VIt<T,S,C>& rhs) 
        { 
            p=rhs.p; 
            TMV_SETFIRSTLAST(rhs._first,rhs._last);
            return *this; 
        }

        VIt<T,S,C>& operator=(const VIter<T>& rhs) 
        { 
            p=rhs.getP(); 
            TMVAssert(rhs.step()==1);
            TMV_SETFIRSTLAST(rhs._first,rhs._last);
            return *this; 
        }

        template <StepType S2> VIt<T,S,C>& operator=(const VIt<T,S2,C>& rhs) 
        { 
            p=rhs.getP(); 
            TMVAssert(rhs.step()==1);
            TMV_SETFIRSTLAST(rhs._first,rhs._last);
            return *this; 
        }

        ~VIt() {}

        T* getP() const { return p; }
        int step() const { return 1; }

        inline bool operator==(const VIt<T,S,C>& rhs) const 
        { return p == rhs.p; }
        inline bool operator!=(const VIt<T,S,C>& rhs) const 
        { return p != rhs.p; }
        inline bool operator<(const VIt<T,S,C>& rhs) const 
        { return (p < rhs.p); }

        inline T& operator*() const
        { 
#ifdef TMVFLDEBUG
            if (!(p>=_first && p<_last)) {
                std::cerr<<"p = "<<p<<std::endl;
                std::cerr<<"first,last = "<<_first<<"  "<<_last<<std::endl;
            }
            TMVAssert(p>=_first);
            TMVAssert(p<_last);
#endif
            return *p; 
        }

        inline VIt<T,S,C>& operator++() { ++p; return *this; }
        inline VIt<T,S,C>& operator--() { --p; return *this; }
        inline VIt<T,S,C> operator++(int) 
        { VIt<T,S,C> p2 = *this; ++p; return p2; }
        inline VIt<T,S,C> operator--(int) 
        { VIt<T,S,C> p2 = *this; --p; return p2; }

        inline VIt<T,S,C>& operator+=(int n) { p += n; return *this; }
        inline VIt<T,S,C>& operator-=(int n) { p -= n; return *this; }
        inline VIt<T,S,C> operator+(int n) const 
        { return VIt<T,S,C>(p+n,1 TMV_FIRSTLAST ); }
        inline VIt<T,S,C> operator-(int n) const 
        { return VIt<T,S,C>(p-n,1 TMV_FIRSTLAST ); }

        inline ptrdiff_t operator-(const VIt<T,S,C>& rhs) const 
        { return p-rhs.p; }

        inline T& operator[](int n) const
        {
#ifdef TMVFLDEBUG
            T* pn = p+n;
            if (!(p>=_first && p<_last)) {
                std::cerr<<"p = "<<p<<std::endl;
                std::cerr<<"first,last = "<<_first<<"  "<<_last<<std::endl;
            }
            TMVAssert(pn >= _first);
            TMVAssert(pn < _first);
            return *pn; 
#else
            return *(p+n);
#endif
        }

        typedef std::random_access_iterator_tag iterator_category;
        typedef T                               value_type;
        typedef ptrdiff_t                       difference_type;
        typedef T*                              pointer;
        typedef T&                              reference;

    private :

        T* p;

#ifdef TMVFLDEBUG
    public :
        const T* _first;
        const T* _last;
#endif

    };

    template <class T, StepType S, ConjType C> class CVIt 
    { 
    public :

        CVIt() : p(0) {}
        CVIt(const T* inp, int TMV_DEBUGPARAM(step)) : p(inp) 
        { TMVAssert(step==1); }
        CVIt(const CVIt<T,S,C>& rhs) : p(rhs.p) 
        { TMVAssert(rhs.step()==1); }
        CVIt(const VIt<T,S,C>& rhs) : p(rhs.getP()) 
        { TMVAssert(rhs.step()==1); }
        CVIt(const CVIter<T>& rhs) : p(rhs.getP()) 
        { TMVAssert(rhs.step()==1); }
        CVIt(const VIter<T>& rhs) : p(rhs.getP()) 
        { TMVAssert(rhs.step()==1); }
        template <StepType S2> CVIt(const CVIt<T,S2,C>& rhs) : p(rhs.getP())
        { TMVAssert(rhs.step()==1); }
        template <StepType S2> CVIt(const VIt<T,S2,C>& rhs) : p(rhs.getP())
        { TMVAssert(rhs.step()==1); }

        CVIt<T,S,C>& operator=(const CVIt<T,S,C>& rhs)
        { TMVAssert(rhs.step()==1); p = rhs.p; return *this; }
        CVIt<T,S,C>& operator=(const VIt<T,S,C>& rhs)
        { TMVAssert(rhs.step()==1); p = rhs.getP(); return *this; }
        CVIt<T,S,C>& operator=(const CVIter<T>& rhs)
        { TMVAssert(rhs.step()==1); p = rhs.getP(); return *this; }
        CVIt<T,S,C>& operator=(const VIter<T>& rhs)
        { TMVAssert(rhs.step()==1); p = rhs.getP(); return *this; }
        template <StepType S2> CVIt<T,S,C>& operator=(const CVIt<T,S2,C>& rhs)
        { TMVAssert(rhs.step()==1); p = rhs.getP(); return *this; }
        template <StepType S2> CVIt<T,S,C>& operator=(const VIt<T,S2,C>& rhs)
        { TMVAssert(rhs.step()==1); p = rhs.getP(); return *this; }
        ~CVIt() {}

        const T* getP() const { return p; }
        int step() const { return 1; }

        inline bool operator==(const CVIt<T,S,C>& rhs) const 
        { return p == rhs.p; }
        inline bool operator!=(const CVIt<T,S,C>& rhs) const 
        { return p != rhs.p; }
        inline bool operator<(const CVIt<T,S,C>& rhs) const 
        { return (p < rhs.p); }

        inline T operator*() const { return *p; }

        inline CVIt<T,S,C>& operator++() { ++p; return *this; }
        inline CVIt<T,S,C>& operator--() { --p; return *this; }
        inline CVIt<T,S,C> operator++(int) 
        { CVIt<T,S,C> p2 = *this; ++p; return p2; }
        inline CVIt<T,S,C> operator--(int) 
        { CVIt<T,S,C> p2 = *this; --p; return p2; }

        inline CVIt<T,S,C>& operator+=(int n) { p += n; return *this; }
        inline CVIt<T,S,C>& operator-=(int n) { p -= n; return *this; }
        inline CVIt<T,S,C> operator+(int n) const 
        { return CVIt<T,S,C>(p+n,1); }
        inline CVIt<T,S,C> operator-(int n) const 
        { return CVIt<T,S,C>(p-n,1); }

        inline ptrdiff_t operator-(const CVIt<T,S,C>& rhs) const 
        { return (p-rhs.p); }

        inline T operator[](int n) const { return *(p+n); }

        typedef std::random_access_iterator_tag iterator_category;
        typedef T                               value_type;
        typedef ptrdiff_t                       difference_type;
        typedef const T*                        pointer;
        typedef const T&                        reference;

    private :

        const T* p;
    };

    template <class T, ConjType C> class VIt<T,Step,C> 
    { 
    public :

        VIt() : p(0), s(0) TMV_DEFFIRSTLAST(0,0) {}
        VIt(T* inp, int instep TMV_PARAMFIRSTLAST(T) ) : 
            p(inp), s(instep) TMV_DEFFIRSTLAST(_first,_last) {}
        VIt(const VIt<T,Step,C>& rhs) : 
            p(rhs.p), s(rhs.s) TMV_DEFFIRSTLAST(rhs._first,rhs._last) {}
        VIt(const VIt<T,Unit,C>& rhs) : 
            p(rhs.getP()), s(1) TMV_DEFFIRSTLAST(rhs._first,rhs._last) {}
        VIt(const VIter<T>& rhs) : 
            p(rhs.getP()), s(rhs.step()) TMV_DEFFIRSTLAST(rhs._first,rhs._last)
        {}
        VIt<T,Step,C>& operator=(const VIt<T,Step,C>& rhs) 
        { 
            TMVAssert(s==rhs.s);
            p=rhs.p; 
            TMV_SETFIRSTLAST(rhs._first,rhs._last);
            return *this; 
        }
        VIt<T,Step,C>& operator=(const VIt<T,Unit,C>& rhs) 
        { 
            TMVAssert(s==1);
            p=rhs.getP(); 
            TMV_SETFIRSTLAST(rhs._first,rhs._last);
            return *this; 
        }
        VIt<T,Step,C>& operator=(const VIter<T>& rhs) 
        { 
            TMVAssert(s==rhs.step());
            p=rhs.getP(); 
            TMV_SETFIRSTLAST(rhs._first,rhs._last);
            return *this; 
        }
        ~VIt() {}

        T* getP() const { return p; }
        int step() const { return s; }

        inline bool operator==(const VIt<T,Step,C>& rhs) const 
        { return p == rhs.p; }
        inline bool operator!=(const VIt<T,Step,C>& rhs) const 
        { return p != rhs.p; }
        inline bool operator<(const VIt<T,Step,C>& rhs) const 
        { return (s > 0 ? p < rhs.p : p > rhs.p); }

        inline T& operator*() const
        { 
#ifdef TMVFLDEBUG
            if (!(p>=_first && p<_last)) {
                std::cerr<<"p = "<<p<<std::endl;
                std::cerr<<"first,last = "<<_first<<"  "<<_last<<std::endl;
            }
            TMVAssert(p>=_first);
            TMVAssert(p<_last);
#endif
            return *p; 
        }

        inline VIt<T,Step,C>& operator++() { p += s; return *this; }
        inline VIt<T,Step,C>& operator--() { p -= s; return *this; }
        inline VIt<T,Step,C> operator++(int) 
        { VIt<T,Step,C> p2 = *this; p+=s; return p2; }
        inline VIt<T,Step,C> operator--(int) 
        { VIt<T,Step,C> p2 = *this; p-=s; return p2; }

        inline VIt<T,Step,C>& operator+=(int n) { p += n*s; return *this; }
        inline VIt<T,Step,C>& operator-=(int n) { p -= n*s; return *this; }
        inline VIt<T,Step,C> operator+(int n) const 
        { return VIt<T,Step,C>(p+n*s,s TMV_FIRSTLAST ); }
        inline VIt<T,Step,C> operator-(int n) const 
        { return VIt<T,Step,C>(p-n*s,s TMV_FIRSTLAST ); }

        inline ptrdiff_t operator-(const VIt<T,Step,C>& rhs) const 
        { return (p-rhs.p)/s; }

        inline T& operator[](int n) const
        {
#ifdef TMVFLDEBUG
            T* pn = p+n*s;
            if (!(pn>=_first && pn<_last)) {
                std::cerr<<"pn = "<<pn<<std::endl;
                std::cerr<<"first,last = "<<_first<<"  "<<_last<<std::endl;
            }
            TMVAssert(pn >= _first);
            TMVAssert(pn < _last);
            return *pn; 
#else
            return *(p+n*s); 
#endif
        }

        typedef std::random_access_iterator_tag iterator_category;
        typedef T                               value_type;
        typedef ptrdiff_t                       difference_type;
        typedef T*                              pointer;
        typedef T&                              reference;

    private :

        T* p;
        const int s;

#ifdef TMVFLDEBUG
    public :
        const T* _first;
        const T* _last;
#endif

    };

    template <class T, ConjType C> class CVIt<T,Step,C> 
    { 
    public :

        CVIt() : p(0), s(0) {}
        CVIt(const T* inp, int instep) : p(inp), s(instep) {}
        CVIt(const CVIt<T,Step,C>& rhs) : p(rhs.p), s(rhs.s) {}
        CVIt(const VIt<T,Step,C>& rhs) : p(rhs.getP()), s(rhs.step()) {}
        CVIt(const CVIt<T,Unit,C>& rhs) : p(rhs.getP()), s(1) {}
        CVIt(const VIt<T,Unit,C>& rhs) : p(rhs.getP()), s(1) {}
        CVIt(const CVIter<T>& rhs) : p(rhs.getP()), s(rhs.step()) {}
        CVIt(const VIter<T>& rhs) : p(rhs.getP()), s(rhs.step()) {}
        CVIt<T,Step,C>& operator=(const CVIt<T,Step,C>& rhs)
        { TMVAssert(s==rhs.s); p = rhs.p; return *this; }
        CVIt<T,Step,C>& operator=(const VIt<T,Step,C>& rhs)
        { TMVAssert(s==rhs.step()); p = rhs.getP(); return *this; }
        CVIt<T,Step,C>& operator=(const CVIt<T,Unit,C>& rhs)
        { TMVAssert(s==1); p = rhs.getP(); return *this; }
        CVIt<T,Step,C>& operator=(const VIt<T,Unit,C>& rhs)
        { TMVAssert(s==1); p = rhs.getP(); return *this; }
        CVIt<T,Step,C>& operator=(const CVIter<T>& rhs)
        { TMVAssert(s==rhs.step()); p = rhs.getP(); return *this; }
        CVIt<T,Step,C>& operator=(const VIter<T>& rhs)
        { TMVAssert(s==rhs.step()); p = rhs.getP(); return *this; }
        ~CVIt() {}

        const T* getP() const { return p; }
        int step() const { return s; }

        inline bool operator==(const CVIt<T,Step,C>& rhs) const 
        { return p == rhs.p; }
        inline bool operator!=(const CVIt<T,Step,C>& rhs) const 
        { return p != rhs.p; }
        inline bool operator<(const CVIt<T,Step,C>& rhs) const 
        { return (s > 0 ? p < rhs.p : p > rhs.p); }

        inline T operator*() const { return *p; }

        inline CVIt<T,Step,C>& operator++() { p += s; return *this; }
        inline CVIt<T,Step,C>& operator--() { p -= s; return *this; }
        inline CVIt<T,Step,C> operator++(int) 
        { CVIt<T,Step,C> p2 = *this; p+=s; return p2; }
        inline CVIt<T,Step,C> operator--(int) 
        { CVIt<T,Step,C> p2 = *this; p-=s; return p2; }

        inline CVIt<T,Step,C>& operator+=(int n) 
        { p += n*s; return *this; }
        inline CVIt<T,Step,C>& operator-=(int n) 
        { p -= n*s; return *this; }
        inline CVIt<T,Step,C> operator+(int n) const 
        { return CVIt<T,Step,C>(p+n*s,s); }
        inline CVIt<T,Step,C> operator-(int n) const 
        { return CVIt<T,Step,C>(p-n*s,s); }

        inline ptrdiff_t operator-(const CVIt<T,Step,C>& rhs) const 
        { return (p-rhs.p)/s; }

        inline T operator[](int n) const { return *(p+n*s); }

        typedef std::random_access_iterator_tag iterator_category;
        typedef T                               value_type;
        typedef ptrdiff_t                       difference_type;
        typedef const T*                        pointer;
        typedef const T&                        reference;

    private :

        const T* p;
        const int s;
    };

    template <class T> 
    class ConjRef; // Undefined unless T is complex<T>
    
    template <class T> 
    class ConjRef<std::complex<T> >
    {
    public:

        typedef std::complex<T> CT;

        explicit ConjRef(CT& _val) : val(_val) {}
        ConjRef(const ConjRef<CT>& rhs) : val(rhs.val) {}
        ~ConjRef() {}

        inline operator CT() const { return std::conj(val); }
        inline CT& getRef() { return val; }
        inline CT conj() const { return val; }
        inline T real() const { return val.real(); }
        inline T imag() const { return -val.imag(); }
        inline CT operator-() const { return -std::conj(val); }

        inline ConjRef<CT>& operator=(const ConjRef<CT>& rhs)
        { val = rhs.val; return *this; }
        inline ConjRef<CT>& operator=(CT rhs)
        { val = std::conj(rhs); return *this; }
        inline ConjRef<CT>& operator=(T rhs)
        { val = rhs; return *this; }

        inline ConjRef<CT>& operator+=(const ConjRef<CT>& x2)
        { val += x2.val; return *this; }
        inline ConjRef<CT>& operator+=(CT x2)
        { val += std::conj(x2); return *this; }
        inline ConjRef<CT>& operator+=(T x2)
        { val += x2; return *this; }
        inline CT operator+(const ConjRef<CT>& x2)
        { return std::conj(val+x2.val); }
        inline friend CT operator+(const ConjRef<CT>& x1, CT x2)
        { return std::conj(x1.val)+x2; }
        inline friend CT operator+(ConjRef<CT>& x1, T x2)
        { return std::conj(x1.val)+x2; }
        inline friend CT operator+(CT x1, const ConjRef<CT>& x2)
        { return x1+std::conj(x2.val); }
        inline friend CT operator+(T x1, const ConjRef<CT>& x2)
        { return x1+std::conj(x2.val); }
        //inline friend CT& operator+=(CT& x1, const ConjRef<CT>& x2)
        //{ return x1+=std::conj(x2.val); }

        inline ConjRef<CT>& operator-=(const ConjRef<CT>& x2) 
        { val -= x2.val; return *this; }
        inline ConjRef<CT>& operator-=(CT x2) 
        { val -= std::conj(x2); return *this; }
        inline ConjRef<CT>& operator-=(T x2) 
        { val -= x2; return *this; }
        inline CT operator-(const ConjRef<CT>& x2)
        { return std::conj(val-x2.val); }
        inline friend CT operator-(const ConjRef<CT>& x1, CT x2)
        { return std::conj(x1.val)-x2; }
        inline friend CT operator-(const ConjRef<CT>& x1, T x2)
        { return std::conj(x1.val)-x2; }
        inline friend CT operator-(CT x1, const ConjRef<CT>& x2)
        { return x1-std::conj(x2.val); }
        inline friend CT operator-(T x1, const ConjRef<CT>& x2)
        { return x1-std::conj(x2.val); }
        //inline friend CT& operator-=(CT& x1, const ConjRef<CT>& x2)
        //{ return x1-=std::conj(x2.val); }

        inline ConjRef<CT>& operator*=(const ConjRef<CT>& x2) 
        { val *= x2.val; return *this; }
        inline ConjRef<CT>& operator*=(CT x2) 
        { val *= std::conj(x2); return *this; }
        inline ConjRef<CT>& operator*=(T x2) 
        { val *= x2; return *this; }
        inline CT operator*(const ConjRef<CT> x2)
        { return std::conj(val*x2.val); }
        inline friend CT operator*(const ConjRef<CT>& x1, CT x2)
        { return std::conj(x1.val)*x2; }
        inline friend CT operator*(const ConjRef<CT>& x1, T x2)
        { return std::conj(x1.val)*x2; }
        inline friend CT operator*(CT x1, const ConjRef<CT>& x2)
        { return x1*std::conj(x2.val); }
        inline friend CT operator*(T x1, const ConjRef<CT>& x2)
        { return x1*std::conj(x2.val); }
        //inline friend CT& operator*=(CT& x1, const ConjRef<CT>& x2)
        //{ return x1*=std::conj(x2.val); }

        inline ConjRef<CT>& operator/=(const ConjRef<CT>& x2) 
        { val /= x2.val; return *this; }
        inline ConjRef<CT>& operator/=(CT x2) 
        { val /= std::conj(x2); return *this; }
        inline ConjRef<CT>& operator/=(T x2) 
        { val /= x2; return *this; }
        inline CT operator/(const ConjRef<CT>& x2)
        { return std::conj(val/x2.val); }
        inline friend CT operator/(const ConjRef<CT>& x1, CT x2)
        { return std::conj(x1.val)/x2; }
        inline friend CT operator/(const ConjRef<CT>& x1, T x2)
        { return std::conj(x1.val)/x2; }
        inline friend CT operator/(CT x1, const ConjRef<CT>& x2)
        { return x1/std::conj(x2.val); }
        inline friend CT operator/(T x1, const ConjRef<CT>& x2)
        { return x1/std::conj(x2.val); }
        //inline friend CT& operator/=(CT& x1, const ConjRef<CT>& x2)
        //{ return x1/=std::conj(x2.val); }

        inline bool operator==(const ConjRef<CT>& x2) const
        { return val == x2.val; }
        inline bool operator==(CT x2) const 
        { return std::conj(val) == x2; }
        inline bool operator==(T x2) const 
        { return std::real(val) == x2 && std::imag(val) == T(0); }
        inline friend bool operator==(CT x1, const ConjRef<CT>& x2)
        { return x2==x1; }
        inline friend bool operator==(T x1, const ConjRef<CT>& x2)
        { return x2==x1; }
        inline bool operator!=(const ConjRef<CT>& x2) const
        { return !(operator==(x2)); }
        inline bool operator!=(CT x2) const 
        { return !(operator==(x2)); }
        inline bool operator!=(T x2) const 
        { return !(operator==(x2)); }
        inline friend bool operator!=(CT x1, const ConjRef<CT>& x2)
        { return !(x2==x1); }
        inline friend bool operator!=(T x1, const ConjRef<CT>& x2)
        { return !(x2==x1); }

        inline void swapWith(CT& x2)
        { 
            TMVAssert(&val != &x2);
            CT temp = x2; x2 = std::conj(val); val = std::conj(temp); 
        }
        inline void swapWith(ConjRef<CT> x2)
        { 
            TMVAssert(&val != &x2);
            CT temp = x2.val; x2.val = val; val = temp; 
        }

        inline friend std::ostream& operator<<(std::ostream& os, ConjRef<CT> x)
        { os << std::conj(x.val); return os; }
        inline friend std::istream& operator>>(std::istream& is, ConjRef<CT> x)
        { is >> x.val; x.val = std::conj(x.val); return is; }

    private:

        CT& val;
    };

    template <class T> inline T TMV_CONJ(const ConjRef<T>& x) 
    { return x.conj(); }
    template <class T> inline TMV_RealType(T) TMV_NORM(const ConjRef<T>& x) 
    { return norm(x.conj()); }
    template <class T> inline TMV_RealType(T) TMV_ABS(const ConjRef<T>& x) 
    { return std::abs(x.conj()); }
    template <class T> inline T TMV_SQR(const ConjRef<T>& x) 
    { return TMV_SQR(T(x)); }
    template <class T> inline T TMV_SQRT(const ConjRef<T>& x) 
    { return TMV_SQRT(T(x)); }
    template <class T> inline TMV_RealType(T) TMV_REAL(const ConjRef<T>& x) 
    { return x.real(); }
    template <class T> inline TMV_RealType(T) TMV_IMAG(const ConjRef<T>& x) 
    { return x.imag(); }
    template <class T> inline void TMV_SWAP(
        tmv::ConjRef<std::complex<T> > x1, tmv::ConjRef<std::complex<T> > x2)
    { return x1.swapWith(x2); }
    template <class T> inline void TMV_SWAP(
        std::complex<T>& x1, tmv::ConjRef<std::complex<T> > x2)
    { return x2.swapWith(x1); }
    template <class T> inline void TMV_SWAP(
        tmv::ConjRef<std::complex<T> > x1, std::complex<T>& x2)
    { return x1.swapWith(x2); }
}

namespace std {
    template <class T> inline T conj(const tmv::ConjRef<T>& x) 
    { return x.conj(); }
    template <class T> inline TMV_RealType(T) norm(const tmv::ConjRef<T>& x) 
    { return norm(x.conj()); }
    template <class T> inline TMV_RealType(T) real(const tmv::ConjRef<T>& x) 
    { return x.real(); }
    template <class T> inline TMV_RealType(T) imag(const tmv::ConjRef<T>& x) 
    { return x.imag(); }
    template <class T> inline void swap(
        tmv::ConjRef<std::complex<T> > x1, tmv::ConjRef<std::complex<T> > x2)
    { return x1.swapWith(x2); }
    template <class T> inline void swap(
        std::complex<T>& x1, tmv::ConjRef<std::complex<T> > x2)
    { return x2.swapWith(x1); }
    template <class T> inline void swap(
        tmv::ConjRef<std::complex<T> > x1, std::complex<T>& x2)
    { return x1.swapWith(x2); }
}

namespace tmv {

    template <class T> class VIt<std::complex<T>,Unit,Conj> 
    { 
        typedef std::complex<T> CT;
    public :

        VIt() : p(0) TMV_DEFFIRSTLAST(0,0) {}
        VIt(CT* inp, int TMV_DEBUGPARAM(step) 
            TMV_PARAMFIRSTLAST(CT)) :
            p(inp) TMV_DEFFIRSTLAST(_first,_last)
        { TMVAssert(step==1); }
        VIt(const VIt<CT,Unit,Conj>& rhs) : 
            p(rhs.p) TMV_DEFFIRSTLAST(rhs._first,rhs._last) {}
        VIt(const VIt<CT,Step,Conj>& rhs) : 
            p(rhs.getP())  TMV_DEFFIRSTLAST(rhs._first,rhs._last) 
        { TMVAssert(rhs.step() == 1); }
        VIt(const VIter<CT >& rhs) : 
            p(rhs.getP())  TMV_DEFFIRSTLAST(rhs._first,rhs._last) 
        { 
            TMVAssert(rhs.step() == 1); 
            TMVAssert(rhs.getC() == Conj);
        }
        VIt& operator=(const VIt<CT,Unit,Conj>& rhs) 
        { 
            p=rhs.p; 
            TMV_SETFIRSTLAST(rhs._first,rhs._last);
            return *this; 
        }
        VIt& operator=(const VIt<CT,Step,Conj>& rhs) 
        { 
            TMVAssert(rhs.step()==1);
            p=rhs.getP(); 
            TMV_SETFIRSTLAST(rhs._first,rhs._last);
            return *this; 
        }
        VIt& operator=(const VIter<CT >& rhs) 
        { 
            TMVAssert(rhs.step()==1);
            TMVAssert(rhs.getC()==Conj);
            p=rhs.getP(); 
            TMV_SETFIRSTLAST(rhs._first,rhs._last);
            return *this; 
        }
        ~VIt() {}

        CT* getP() const { return p; }
        int step() const { return 1; }

        inline bool operator==(const VIt<CT,Unit,Conj>& rhs) const 
        { return p == rhs.p; }
        inline bool operator!=(const VIt<CT,Unit,Conj>& rhs) const 
        { return p != rhs.p; }
        inline bool operator<(const VIt<CT,Unit,Conj>& rhs) const 
        { return (p < rhs.p); }

        inline ConjRef<CT > operator*() const
        { 
#ifdef TMVFLDEBUG
            if (!(p>=_first && p<_last)) {
                std::cerr<<"p = "<<p<<std::endl;
                std::cerr<<"first,last = "<<_first<<"  "<<_last<<std::endl;
            }
            TMVAssert(p>=_first);
            TMVAssert(p<_last);
#endif
            return ConjRef<CT >(*p);
        }

        inline VIt<CT,Unit,Conj>& operator++() 
        { ++p; return *this; }
        inline VIt<CT,Unit,Conj>& operator--() 
        { --p; return *this; }
        inline VIt<CT,Unit,Conj> operator++(int) 
        { VIt<CT,Unit,Conj> p2 = *this; ++p; return p2; }
        inline VIt<CT,Unit,Conj> operator--(int) 
        { VIt<CT,Unit,Conj> p2 = *this; --p; return p2; }

        inline VIt<CT,Unit,Conj>& operator+=(int n) 
        { p += n; return *this; }
        inline VIt<CT,Unit,Conj>& operator-=(int n) 
        { p -= n; return *this; }
        inline VIt<CT,Unit,Conj> operator+(int n) const 
        { return VIt<CT,Unit,Conj>(p+n,1 TMV_FIRSTLAST ); }
        inline VIt<CT,Unit,Conj> operator-(int n) const 
        { return VIt<CT,Unit,Conj>(p-n,1 TMV_FIRSTLAST ); }

        inline ptrdiff_t operator-(
            const VIt<CT,Unit,Conj>& rhs) const 
        { return (p-rhs.p); }

        inline CT& operator[](int n) const
        {
#ifdef TMVFLDEBUG
            CT* pn = p+n;
            if (!(pn>=_first && pn<_last)) {
                std::cerr<<"pn = "<<pn<<std::endl;
                std::cerr<<"first,last = "<<_first<<"  "<<_last<<std::endl;
            }
            TMVAssert(pn >= _first);
            TMVAssert(pn < _last);
            return ConjRef<CT >(*pn); 
#else
            return ConjRef<CT >(*(p+n)); 
#endif
        }

        typedef std::random_access_iterator_tag   iterator_category;
        typedef CT                                value_type;
        typedef ptrdiff_t                         difference_type;
        typedef VIt<CT,Unit,Conj>                 pointer;
        typedef ConjRef<CT >                      reference;

    private :

        CT* p;

#ifdef TMVFLDEBUG
    public :
        const CT* _first;
        const CT* _last;
#endif
    };

    template <class T> class CVIt<std::complex<T>,Unit,Conj> 
    { 
        typedef std::complex<T> CT;
    public :

        CVIt() : p(0) {}
        CVIt(const CT* inp, int TMV_DEBUGPARAM(step)) : p(inp) 
        { TMVAssert(step==1); }
        CVIt(const CVIt<CT,Unit,Conj>& rhs) : p(rhs.p) {}
        CVIt(const VIt<CT,Unit,Conj>& rhs) : p(rhs.getP()) {}
        CVIt(const CVIt<CT,Step,Conj>& rhs) : p(rhs.getP())
        { TMVAssert(rhs.step() == 1); }
        CVIt(const VIt<CT,Step,Conj>& rhs) : p(rhs.getP())
        { TMVAssert(rhs.step() == 1); }
        CVIt(const CVIter<CT >& rhs) : p(rhs.getP())
        { 
            TMVAssert(rhs.getC() == Conj);
            TMVAssert(rhs.step() == 1); 
        }
        CVIt(const VIter<CT >& rhs) : p(rhs.getP())
        { 
            TMVAssert(rhs.getC() == Conj);
            TMVAssert(rhs.step() == 1); 
        }
        CVIt<CT,Unit,Conj>& operator=(const CVIt<CT,Unit,Conj>& rhs)
        { p = rhs.p; return *this; }
        CVIt<CT,Unit,Conj>& operator=(const VIt<CT,Unit,Conj>& rhs)
        { p = rhs.getP(); return *this; }
        CVIt<CT,Unit,Conj>& operator=(const CVIt<CT,Step,Conj>& rhs)
        { TMVAssert(rhs.step()==1); p = rhs.getP(); return *this; }
        CVIt<CT,Unit,Conj>& operator=(const VIt<CT,Step,Conj>& rhs)
        { TMVAssert(rhs.step()==1); p = rhs.getP(); return *this; }
        CVIt<CT,Unit,Conj>& operator=(const CVIter<std::complex<T> >& rhs)
        { 
            TMVAssert(rhs.getC() == Conj);
            TMVAssert(rhs.step()==1); 
            p = rhs.getP(); 
            return *this; 
        }
        CVIt<CT,Unit,Conj>& operator=(const VIter<CT >& rhs)
        { 
            TMVAssert(rhs.getC() == Conj);
            TMVAssert(rhs.step()==1); 
            p = rhs.getP(); 
            return *this; 
        }
        ~CVIt() {}

        const CT* getP() const { return p; }
        int step() const { return 1; }

        inline bool operator==(const CVIt<CT,Unit,Conj>& rhs) const 
        { return p == rhs.p; }
        inline bool operator!=(const CVIt<CT,Unit,Conj>& rhs) const 
        { return p != rhs.p; }
        inline bool operator<(const CVIt<CT,Unit,Conj>& rhs) const 
        { return (p < rhs.p); }

        inline CT operator*() const { return std::conj(*p); }

        inline CVIt<CT,Unit,Conj>& operator++() { ++p; return *this; }
        inline CVIt<CT,Unit,Conj>& operator--() { --p; return *this; }
        inline CVIt<CT,Unit,Conj> operator++(int) 
        { CVIt<CT,Unit,Conj> p2 = *this; ++p; return p2; }
        inline CVIt<CT,Unit,Conj> operator--(int) 
        { CVIt<CT,Unit,Conj> p2 = *this; --p; return p2; }

        inline CVIt<CT,Unit,Conj>& operator+=(int n) { p += n; return *this; }
        inline CVIt<CT,Unit,Conj>& operator-=(int n) { p -= n; return *this; }
        inline CVIt<CT,Unit,Conj> operator+(int n) const 
        { return CVIt<CT,Unit,Conj>(p+n,1); }
        inline CVIt<CT,Unit,Conj> operator-(int n) const 
        { return CVIt<CT,Unit,Conj>(p-n,1); }

        inline ptrdiff_t operator-(const CVIt<CT,Unit,Conj>& rhs) const 
        { return (p-rhs.p); }

        inline CT operator[](int n) const 
        { return std::conj(*(p+n)); }

        typedef std::random_access_iterator_tag   iterator_category;
        typedef CT                                value_type;
        typedef ptrdiff_t                         difference_type;
        typedef CVIt<CT,Unit,Conj>                pointer;
        typedef const ConjRef<CT >                reference;

    private :

        const std::complex<T>* p;
    };

    template <class T> class VIt<std::complex<T>,Step,Conj> 
    { 
        typedef std::complex<T> CT;
    public :

        VIt() : p(0), s(0) TMV_DEFFIRSTLAST(0,0) {}
        VIt(CT* inp, int instep TMV_PARAMFIRSTLAST(CT)) :
            p(inp), s(instep) TMV_DEFFIRSTLAST(_first,_last) {}
        VIt(const VIt<CT,Step,Conj>& rhs) : 
            p(rhs.p), s(rhs.s) TMV_DEFFIRSTLAST(rhs._first,rhs._last) {}
        VIt(const VIt<CT,Unit,Conj>& rhs) : 
            p(rhs.getP()), s(1) TMV_DEFFIRSTLAST(rhs._first,rhs._last) {}
        VIt(const VIter<CT >& rhs) : 
            p(rhs.getP()), s(rhs.step()) TMV_DEFFIRSTLAST(rhs._first,rhs._last) 
        { TMVAssert(rhs.getC()==Conj); }
        VIt& operator=(const VIt<CT,Step,Conj>& rhs) 
        { 
            TMVAssert(s==rhs.s);
            p=rhs.p; 
            TMV_SETFIRSTLAST(rhs._first,rhs._last);
            return *this; 
        }
        VIt& operator=(const VIt<CT,Unit,Conj>& rhs) 
        { 
            TMVAssert(s==1);
            p=rhs.getP(); 
            TMV_SETFIRSTLAST(rhs._first,rhs._last);
            return *this; 
        }
        VIt& operator=(const VIter<CT >& rhs) 
        { 
            TMVAssert(s==rhs.step());
            TMVAssert(rhs.getC()==Conj);
            p=rhs.getP(); 
            TMV_SETFIRSTLAST(rhs._first,rhs._last);
            return *this; 
        }
        ~VIt() {}

        CT* getP() const { return p; }
        int step() const { return s; }

        inline bool operator==(const VIt<CT,Step,Conj>& rhs) const 
        { return p == rhs.p; }
        inline bool operator!=(const VIt<CT,Step,Conj>& rhs) const 
        { return p != rhs.p; }
        inline bool operator<(const VIt<CT,Step,Conj>& rhs) const 
        { return (s > 0 ? p < rhs.p : p > rhs.p); }

        inline ConjRef<CT > operator*() const
        { 
#ifdef TMVFLDEBUG
            if (!(p>=_first && p<_last)) {
                std::cerr<<"p = "<<p<<std::endl;
                std::cerr<<"first,last = "<<_first<<"  "<<_last<<std::endl;
            }
            TMVAssert(p>=_first);
            TMVAssert(p<_last);
#endif
            return ConjRef<CT >(*p);
        }

        inline VIt<CT,Step,Conj>& operator++() { p+=s; return *this; }
        inline VIt<CT,Step,Conj>& operator--() { p-=s; return *this; }
        inline VIt<CT,Step,Conj> operator++(int) 
        { VIt<CT,Step,Conj> p2 = *this; p+=s; return p2; }
        inline VIt<CT,Step,Conj> operator--(int) 
        { VIt<CT,Step,Conj> p2 = *this; p-=s; return p2; }

        inline VIt<CT,Step,Conj>& operator+=(int n) { p += n*s; return *this; }
        inline VIt<CT,Step,Conj>& operator-=(int n) { p -= n*s; return *this; }
        inline VIt<CT,Step,Conj> operator+(int n) const 
        { return VIt<CT,Step,Conj>(p+n*s,s TMV_FIRSTLAST ); }
        inline VIt<CT,Step,Conj> operator-(int n) const 
        { return VIt<CT,Step,Conj>(p-n*s,s TMV_FIRSTLAST ); }

        inline ptrdiff_t operator-(const VIt<CT,Step,Conj>& rhs) const 
        { return (p-rhs.p)/s; }

        inline CT& operator[](int n) const
        {
#ifdef TMVFLDEBUG
            CT* pn = p+n*s;
            if (!(pn>=_first && pn<_last)) {
                std::cerr<<"pn = "<<pn<<std::endl;
                std::cerr<<"first,last = "<<_first<<"  "<<_last<<std::endl;
            }
            TMVAssert(pn >= _first);
            TMVAssert(pn < _last);
            return ConjRef<CT >(*pn); 
#else
            return ConjRef<CT >(*(p+n*s)); 
#endif
        }

        typedef std::random_access_iterator_tag   iterator_category;
        typedef CT                                value_type;
        typedef ptrdiff_t                         difference_type;
        typedef VIt<CT,Step,Conj>                 pointer;
        typedef ConjRef<CT >                      reference;

    private :

        CT* p;
        const int s;

#ifdef TMVFLDEBUG
    public :
        const CT* _first;
        const CT* _last;
#endif
    };

    template <class T> class CVIt<std::complex<T>,Step,Conj> 
    { 
        typedef std::complex<T> CT;
    public :

        CVIt() : p(0), s(0) {}
        CVIt(const CT* inp, int instep) : p(inp), s(instep) {}
        CVIt(const CVIt<CT,Step,Conj>& rhs) : p(rhs.p), s(rhs.s) {}
        CVIt(const VIt<CT,Step,Conj>& rhs) : p(rhs.getP()), s(rhs.step()) {}
        CVIt(const CVIt<CT,Unit,Conj>& rhs) : p(rhs.getP()), s(1) {}
        CVIt(const VIt<CT,Unit,Conj>& rhs) : p(rhs.getP()), s(1) {}
        CVIt(const CVIter<CT >& rhs) : p(rhs.getP()), s(rhs.step()) 
        { TMVAssert(rhs.getC()==Conj); }
        CVIt(const VIter<CT >& rhs) : p(rhs.getP()), s(rhs.step()) 
        { TMVAssert(rhs.getC()==Conj); }
        CVIt<CT,Step,Conj>& operator=(const CVIt<CT,Step,Conj>& rhs)
        { TMVAssert(s==rhs.s); p = rhs.p; return *this; }
        CVIt<CT,Step,Conj>& operator=(const VIt<CT,Step,Conj>& rhs)
        { TMVAssert(s==rhs.s); p = rhs.getP(); return *this; }
        CVIt<CT,Step,Conj>& operator=(const CVIt<CT,Unit,Conj>& rhs)
        { TMVAssert(s==1); p = rhs.getP(); return *this; }
        CVIt<CT,Step,Conj>& operator=(const VIt<CT,Unit,Conj>& rhs)
        { TMVAssert(s==1); p = rhs.getP(); return *this; }
        CVIt<CT,Step,Conj>& operator=(const CVIter<CT >& rhs)
        { 
            TMVAssert(s==rhs.step()); 
            TMVAssert(rhs.getC()==Conj);
            p = rhs.getP(); return *this; 
        }
        CVIt<CT,Step,Conj>& operator=(const VIter<CT >& rhs)
        { 
            TMVAssert(s==rhs.step()); 
            TMVAssert(rhs.getC()==Conj);
            p = rhs.getP(); return *this; 
        }
        ~CVIt() {}

        const CT* getP() const { return p; }
        int step() const { return s; }

        inline bool operator==(const CVIt<CT,Step,Conj>& rhs) const 
        { return p == rhs.p; }
        inline bool operator!=(const CVIt<CT,Step,Conj>& rhs) const 
        { return p != rhs.p; }
        inline bool operator<(const CVIt<CT,Step,Conj>& rhs) const 
        { return (s > 0 ? p < rhs.p : p > rhs.p); }

        inline CT operator*() const { return std::conj(*p); }

        inline CVIt<CT,Step,Conj>& operator++() { p+=s; return *this; }
        inline CVIt<CT,Step,Conj>& operator--() { p-=s; return *this; }
        inline CVIt<CT,Step,Conj> operator++(int) 
        { CVIt<CT,Step,Conj> p2 = *this; p+=s; return p2; }
        inline CVIt<CT,Step,Conj> operator--(int) 
        { CVIt<CT,Step,Conj> p2 = *this; p-=s; return p2; }

        inline CVIt<CT,Step,Conj>& operator+=(int n) 
        { p += n*s; return *this; }
        inline CVIt<CT,Step,Conj>& operator-=(int n) 
        { p -= n*s; return *this; }
        inline CVIt<CT,Step,Conj> operator+(int n) const 
        { return CVIt<CT,Step,Conj>(p+n*s,s); }
        inline CVIt<CT,Step,Conj> operator-(int n) const 
        { return CVIt<CT,Step,Conj>(p-n*s,s); }

        inline ptrdiff_t operator-(const CVIt<CT,Step,Conj>& rhs) const 
        { return (p-rhs.p)/s; }

        inline CT operator[](int n) const 
        { return std::conj(*(p+n*s)); }

        typedef std::random_access_iterator_tag   iterator_category;
        typedef CT                                value_type;
        typedef ptrdiff_t                         difference_type;
        typedef CVIt<CT,Step,Conj>                pointer;
        typedef const ConjRef<CT >                reference;

    private :

        const CT* p;
        const int s;
    };

    template <class T> 
    class VarConjRef; // Undefined unless T is complex<T>

    template <class T> 
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

    template <class T> inline T TMV_CONJ(const VarConjRef<T>& x) 
    { return x.conj(); }
    template <class T> inline TMV_RealType(T) TMV_NORM(const VarConjRef<T>& x) 
    { return TMV_NORM(T(x)); }
    template <class T> inline TMV_RealType(T) TMV_ABS(const VarConjRef<T>& x) 
    { return TMV_ABS(T(x)); }
    template <class T> inline T TMV_SQR(const VarConjRef<T>& x) 
    { return TMV_SQR(T(x)); }
    template <class T> inline T TMV_SQRT(const VarConjRef<T>& x) 
    { return TMV_SQRT(T(x)); }
    template <class T> inline TMV_RealType(T) TMV_REAL(const VarConjRef<T>& x) 
    { return x.real(); }
    template <class T> inline TMV_RealType(T) TMV_IMAG(const VarConjRef<T>& x) 
    { return x.imag(); }
    template <class T> inline void TMV_SWAP(
        tmv::VarConjRef<std::complex<T> > x1,
        tmv::VarConjRef<std::complex<T> > x2)
    { return x1.swapWith(x2); }
    template <class T> inline void TMV_SWAP(
        std::complex<T>& x1, tmv::VarConjRef<std::complex<T> > x2)
    { return x2.swapWith(x1); }
    template <class T> inline void TMV_SWAP(
        tmv::VarConjRef<std::complex<T> > x1, std::complex<T>& x2)
    { return x1.swapWith(x2); }
}

namespace std {
    template <class T> inline T conj(const tmv::VarConjRef<T>& x) 
    { return x.conj(); }
    template <class T> inline TMV_RealType(T) norm(const tmv::VarConjRef<T>& x) 
    { return norm(T(x)); }
    template <class T> inline TMV_RealType(T) real(const tmv::VarConjRef<T>& x) 
    { return x.real(); }
    template <class T> inline TMV_RealType(T) imag(const tmv::VarConjRef<T>& x) 
    { return x.imag(); }
    template <class T> 
    inline void swap(
        tmv::VarConjRef<std::complex<T> > x1,
        tmv::VarConjRef<std::complex<T> > x2)
    { return x1.swapWith(x2); }
    template <class T> 
    inline void swap(
        std::complex<T>& x1, tmv::VarConjRef<std::complex<T> > x2)
    { return x2.swapWith(x1); }
    template <class T> 
    inline void swap(tmv::VarConjRef<std::complex<T> > x1, std::complex<T>& x2)
    { return x1.swapWith(x2); }
}

namespace tmv {

    template <class T> class VIter<std::complex<T> >
    { 
        typedef std::complex<T> CT;
    public :

        VIter() : p(0), s(0), c(NonConj) TMV_DEFFIRSTLAST(0,0) {}
        VIter(CT* inp, int instep, ConjType inc 
              TMV_PARAMFIRSTLAST(CT) ) : 
            p(inp), s(instep), c(inc) TMV_DEFFIRSTLAST(_first,_last) {}
        VIter(const VIter<CT >& rhs) : 
            p(rhs.p), s(rhs.s), c(rhs.c) TMV_DEFFIRSTLAST(rhs._first,rhs._last) 
        {}
        VIter<CT >& operator=(const VIter<CT >& rhs) 
        { 
            TMVAssert(s==rhs.s);
            TMVAssert(c==rhs.c);
            p=rhs.p; 
            TMV_SETFIRSTLAST(rhs._first,rhs._last);
            return *this; 
        }
        ~VIter() {}

        CT* getP() const { return p; }
        int step() const { return s; }
        ConjType getC() const { return c; }

        inline bool operator==(const VIter<CT >& rhs) const 
        { return p == rhs.p; }
        inline bool operator!=(const VIter<CT >& rhs) const 
        { return p != rhs.p; }
        inline bool operator<(const VIter<CT >& rhs) const 
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

        inline VIter<CT >& operator++() { p += s; return *this; }
        inline VIter<CT >& operator--() { p -= s; return *this; }
        inline VIter<CT > operator++(int) 
        { VIter<CT > p2 = *this; p+=s; return p2; }
        inline VIter<CT > operator--(int) 
        { VIter<CT > p2 = *this; p-=s; return p2; }

        inline VIter<CT >& operator+=(int n) 
        { if(s==1) ++p; else p += n*s; return *this; }
        inline VIter<CT >& operator-=(int n) 
        { if(s==1) --p; else p -= n*s; return *this; }
        inline VIter<CT > operator+(int n) const 
        { return VIter<CT >(s==1?p+n:p+n*s,s,c TMV_FIRSTLAST ); }
        inline VIter<CT > operator-(int n) const 
        { return VIter<CT >(s==1?p-n:p-n*s,s,c TMV_FIRSTLAST ); }

        inline ptrdiff_t operator-(const VIter<CT >& rhs) const 
        {
            TMVAssert(rhs.c==c);
            TMVAssert(rhs.s==s);
            return s==1 ? p-rhs.p : (p-rhs.p)/s; 
        }

        inline VarConjRef<CT > operator[](int n) const
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
        const int s;
        ConjType c;

#ifdef TMVFLDEBUG
    public :
        const CT* _first;
        const CT* _last;
#endif
    };

    template <class T> class CVIter<std::complex<T> >
    { 
        typedef std::complex<T> CT;
    public :

        CVIter() : p(0), s(0), c(NonConj) {}
        CVIter(const CT* inp, int instep, ConjType inc) :
            p(inp), s(instep), c(inc) {}
        CVIter(const CVIter<CT >& rhs) : p(rhs.p), s(rhs.s), c(rhs.c) {}
        CVIter(const VIter<CT >& rhs) : 
            p(rhs.getP()), s(rhs.step()), c(rhs.getC()) {}
        CVIter<CT >& operator=(const CVIter<CT >& rhs) 
        { 
            TMVAssert(s==rhs.s);
            TMVAssert(c==rhs.c);
            p=rhs.p; 
            return *this; 
        }
        CVIter<CT >& operator=(const VIter<CT >& rhs) 
        { 
            TMVAssert(s==rhs.step());
            TMVAssert(c==rhs.getC());
            p=rhs.getP(); 
            return *this; 
        }
        ~CVIter() {}

        const CT* getP() const { return p; }
        int step() const { return s; }
        ConjType getC() const { return c; }

        inline bool operator==(const CVIter<CT >& rhs) const 
        { return p == rhs.p; }
        inline bool operator!=(const CVIter<CT >& rhs) const 
        { return p != rhs.p; }
        inline bool operator<(const CVIter<CT >& rhs) const 
        { return (s > 0 ? p < rhs.p : p > rhs.p); }

        inline CT operator*() const 
        { return c==Conj ? std::conj(*p) : *p; }

        inline CVIter<CT >& operator++() { p += s; return *this; }
        inline CVIter<CT >& operator--() { p -= s; return *this; }
        inline CVIter<CT > operator++(int) 
        { CVIter<CT > p2 = *this; p+=s; return p2; }
        inline CVIter<CT > operator--(int) 
        { CVIter<CT > p2 = *this; p-=s; return p2; }

        inline CVIter<CT >& operator+=(int n) 
        { if(s==1) ++p; else p += n*s; return *this; }
        inline CVIter<CT >& operator-=(int n) 
        { if(s==1) --p; else p -= n*s; return *this; }
        inline CVIter<CT > operator+(int n) const 
        { return CVIter<CT >(s==1?p+n:p+n*s,s,c); }
        inline CVIter<CT > operator-(int n) const 
        { return CVIter<CT >(s==1?p-n:p-n*s,s,c); }

        inline ptrdiff_t operator-(const CVIter<CT >& rhs) const 
        {
            TMVAssert(rhs.c==c);
            TMVAssert(rhs.s==s);
            return s==1 ? p-rhs.p : (p-rhs.p)/s; 
        }

        inline CT operator[](int n) const 
        { return c==Conj ? std::conj(*(p+(s==1?n:n*s))) : *(p+(s==1?n:n*s)); }

        typedef std::random_access_iterator_tag   iterator_category;
        typedef CT                                value_type;
        typedef ptrdiff_t                         difference_type;
        typedef const CT*                         pointer;
        typedef const VarConjRef<CT >             reference;

    private :

        const CT* p;
        const int s;
        ConjType c;
    };

    template <class T> inline T& TMV_REF(T* vi, ConjType )
    { return *vi; }
    template <class T> inline VarConjRef<std::complex<T> > TMV_REF(
        std::complex<T>* vi, ConjType ct)
    { return VarConjRef<std::complex<T> >(*vi,ct); }

    template <class T> inline std::string TMV_Text(ConjRef<T>)
    { return std::string("ConjRef<") + TMV_Text(T()) + ">"; }

    template <class T> inline std::string TMV_Text(VarConjRef<T>)
    { return std::string("VarConjRef<") + TMV_Text(T()) + ">"; }

    template <class T, StepType S, ConjType C> 
    inline std::string TMV_Text(VIt<T,S,C>)
    { 
        return std::string("VIt<") + TMV_Text(T()) + "," +
            TMV_Text(S) + "," + TMV_Text(C) + ">"; 
    }

    template <class T, StepType S, ConjType C> 
    inline std::string TMV_Text(CVIt<T,S,C>)
    {
        return std::string("CVIt<") + TMV_Text(T()) + "," +
            TMV_Text(S) + "," + TMV_Text(C) + ">"; 
    }

    template <class T> inline std::string TMV_Text(VIter<T> it)
    { 
        return std::string("VIter<") + TMV_Text(T()) + "," +
            it.step()==1 ? TMV_Text(Unit) : TMV_Text(Step) + "," +
            TMV_Text(it.getC()) + ">"; 
    }

    template <class T> inline std::string TMV_Text(CVIter<T> it)
    { 
        return std::string("CVIter<") + TMV_Text(T()) + "," +
            it.step()==1 ? TMV_Text(Unit) : TMV_Text(Step) + "," +
            TMV_Text(it.getC()) + ">"; 
    }

} // namespace tmv

#endif
