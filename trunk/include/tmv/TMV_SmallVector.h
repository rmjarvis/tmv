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
// This file defines the SmallVector class.
//
// A SmallVector is a vector whose size is known at compile time.
// This allows for a lot of optimizations by the compiler in implementing
// the various operations on it. 
// 
// Since the instantiation needs to know the size of the vector, all 
// operations on the SmallVector are done inline, rather than precompiled.
// This gives another speed boost in most cases, since the compiler can
// avoid implementing a function call in most cases.
//
// Finally, another advantage is that we allocate the data on the stack
// rather than the heap, so we avoid new and delete calls as well.
// However, stack sizes are usually limited to hundreds of KB, 
// (mine is 8 MB), so we set a maximum size of 1KB for each SmallVector.
// (For double, this means up to N=128 will be allocated on the stack.)
// Any bigger than that, and the performance drop from using the
// heap is pretty irrelevant.  (See TMV_StackArray.h to change this.)
// 
// One drawback of using a SmallVector is that it does not do any
// alias checking in the aritmetic statements.  So a statement like
// v = m * v will not produce a correct answer. 
// Normally this is a feature, since the alias checks can be a 
// significant fraction of the calculation time for small vectors/matrices.
// 
// You can workaround this when necessary by explicitly making a copy.
// The easiest way is with the .copy() method.  e.g. v = m * v.copy().
//
//
// Constructors:
//
//    SmallVector<T,N,I>()  
//        Makes a Vector of size N with _uninitialized_ values
//
//    SmallVector<T,N,I>(T x)
//        Makes a Vector of size N with all values = x
//
//    SmallVector<T,N,I>(const T* v2)
//    SmallVector<T,N,I>(const vector<T>& v2)
//    SmallVector<T,N,I>(const BaseVector<V2>& v2)
//        Makes a SmallVector which copies the elements of v2.
//

#ifndef TMV_SmallVector_H
#define TMV_SmallVector_H


#include <vector>
#include <sstream>

#include "TMV_BaseVector.h"
#include "TMV_VIt.h"
#include "TMV_Array.h"

namespace tmv {

    //
    // SmallVector
    //

    template <class T, int N, IndexStyle I> 
    struct Traits<SmallVector<T,N,I> >
    {
        typedef T value_type;

        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;
        enum { visreal = Traits<T>::isreal };
        enum { viscomplex = Traits<T>::iscomplex };

        typedef SmallVector<T,N,I> type;
        typedef const type& calc_type; 
        typedef const type& eval_type; 
        typedef type copy_type;

        enum { vsize = N }; 
        enum { vfort = (I == FortranStyle) };
        enum { vcalc = true };
        enum { vstep = 1 }; 
        enum { vconj = false }; 

        typedef ConstVectorView<T,1,false,I> const_subvector_type;
        typedef ConstVectorView<T,UNKNOWN,false,I> const_subvector_step_type;
        typedef ConstSmallVectorView<T,N,1,false,I> const_view_type;
        typedef ConstSmallVectorView<T,N,1,false,CStyle> const_cview_type;
        typedef ConstSmallVectorView<T,N,1,false,FortranStyle> 
            const_fview_type;
        typedef ConstVectorView<T> const_xview_type;
        typedef ConstSmallVectorView<T,N,1,false,I> const_unitview_type;
        typedef ConstSmallVectorView<T,N,1,viscomplex,I> const_conjugate_type;
        typedef ConstSmallVectorView<T,N,-1,false,I> const_reverse_type;
        typedef ConstSmallVectorView<real_type,N,visreal?1:2,false,I> 
            const_realpart_type;
        typedef const_realpart_type const_imagpart_type;
        typedef ConstSmallVectorView<real_type,visreal?N:N*2,1,false,I> 
            const_flatten_type;
        typedef ConstSmallVectorView<T,N,1,false,I> const_nonconj_type;
        typedef SmallVectorView<T,N,1,false,I> nonconst_type;

        typedef CVIt<T,1,false> const_iterator;
        typedef CVIt<T,-1,false> const_reverse_iterator;

        typedef T& reference;

        typedef VectorView<T,1,false,I> subvector_type;
        typedef VectorView<T,UNKNOWN,false,I> subvector_step_type;
        typedef SmallVectorView<T,N,1,false,I> view_type;
        typedef SmallVectorView<T,N,1,false,CStyle> cview_type;
        typedef SmallVectorView<T,N,1,false,FortranStyle> fview_type;
        typedef VectorView<T> xview_type;
        typedef SmallVectorView<T,N,1,false,I> unitview_type;
        typedef SmallVectorView<T,N,1,viscomplex,I> conjugate_type;
        typedef SmallVectorView<T,N,-1,false,I> reverse_type;
        typedef SmallVectorView<real_type,N,visreal?1:2,false,I> realpart_type;
        typedef realpart_type imagpart_type;
        typedef SmallVectorView<real_type,visreal?N:N*2,1,false,I> 
            flatten_type;
        typedef SmallVectorView<T,N,1,false,I> nonconj_type;

        typedef VIt<T,1,false> iterator;
        typedef VIt<T,-1,false> reverse_iterator;
    };


    //#ifdef XTEST
#ifdef TMV_DEBUG
#ifndef XTEST_DEBUG
#define XTEST_DEBUG
#endif
#endif
    //#endif

    template <class T, int N, IndexStyle I> 
    class SmallVector : public BaseVector_Mutable<SmallVector<T,N,I> >
    {
    public:

        typedef SmallVector<T,N,I> type;
        typedef BaseVector_Mutable<type> base_mut;
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        typedef typename base_mut::iterator iterator;

        enum { vsize = Traits<type>::vsize };
        enum { vfort = Traits<type>::vfort };
        enum { vcalc = Traits<type>::vcalc };
        enum { visreal = Traits<type>::visreal };
        enum { viscomplex = Traits<type>::viscomplex };
        enum { vstep = Traits<type>::vstep };
        enum { vconj = Traits<type>::vconj };

        //
        // Constructors
        //

        inline SmallVector(size_t n=N) 
        {
            TMVStaticAssert(N > 0);
            TMVAssert(n == N);
#ifdef TMV_DEBUG
            this->setAllTo(T(888));
#endif
        }

        explicit inline SmallVector(size_t n, T x) 
        {
            TMVStaticAssert(N > 0);
            TMVAssert(n == N);
            this->setAllTo(x); 
        }

        explicit inline SmallVector(size_t n, const T* v) 
        {
            TMVStaticAssert(N > 0);
            TMVAssert(n == N);
#ifdef XTEST_DEBUG
            this->setAllTo(T(888));
#endif
            ConstSmallVectorView<T,N>(v).newAssignTo(*this);
        }

        explicit inline SmallVector(const std::vector<T>& v2) 
        {
            TMVStaticAssert(N > 0);
            TMVAssert(v2.size() == N);
#ifdef XTEST_DEBUG
            this->setAllTo(T(888));
#endif
            ConstSmallVectorView<T,N>(&v2[0]).newAssignTo(*this);
        }

        inline SmallVector(const type& v2) 
        {
            TMVStaticAssert(N > 0);
#ifdef XTEST_DEBUG
            this->setAllTo(T(888));
#endif
            v2.newAssignTo(*this);
        }

        template <class V2>
        inline SmallVector(const BaseVector<V2>& v2) 
        {
            TMVStaticAssert(N > 0);
            TMVAssert(v2.size() == N);
#ifdef XTEST_DEBUG
            this->setAllTo(T(888));
#endif
            v2.newAssignTo(*this);
        }

        inline ~SmallVector()
        {
#ifdef TMV_DEBUG
            this->setAllTo(T(999));
#endif
        }

        //
        // Op =
        //

        inline type& operator=(const type& v2)
        { 
            if (&v2 != this) base_mut::operator=(v2); 
            return *this; 
        }

        template <class V2>
        inline type& operator=(const BaseVector<V2>& v2) 
        { base_mut::operator=(v2); return *this; }


        //
        // Auxilliary Functions
        //

        inline const T* cptr() const { return itsv; }
        inline T* ptr() { return itsv; }
        inline T cref(int i) const  { return itsv[i]; }
        inline T& ref(int i) { return itsv[i]; }

        inline size_t size() const { return N; }
        inline int step() const { return 1; }
        inline bool isconj() const { return false; }


    protected :

        StackArray<T,N> itsv;

    }; // SmallVector

    template <class T, int N>
    class SmallVectorF : public SmallVector<T,N,FortranStyle>
    {
    public:

        typedef SmallVectorF<T,N> type;
        typedef SmallVector<T,N,FortranStyle> vtype;

        inline SmallVectorF() : type() {}
        explicit inline SmallVectorF(T x) : type(x) {}
        explicit inline SmallVectorF(const T* v2)  : type(v2) {}
        explicit inline SmallVectorF(const std::vector<T>& v2)  : type(v2) {}
        inline SmallVectorF(const type& v2)  : type(v2) {}
        template <class V2>
        inline SmallVectorF(const BaseVector<V2>& v2)  : type(v2) {}
        inline ~SmallVectorF() {}

        template <class V2>
        inline type& operator=(const BaseVector<V2>& v2) 
        { vtype::operator=(v2); return *this; }
        inline type& operator=(const type& v2)
        { vtype::operator=(v2); return *this; }
    }; // SmallVectorF


    //
    // ConstSmallVectorView
    //

    template <class T, int N, int S, bool C, IndexStyle I>
    struct Traits<ConstSmallVectorView<T,N,S,C,I> >
    {
        typedef T value_type;
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;
        enum { visreal = Traits<T>::isreal };
        enum { viscomplex = Traits<T>::iscomplex };

        typedef ConstSmallVectorView<T,N,S,C,I> type;
        typedef const type& calc_type;
        typedef const type& eval_type;

        enum { vsize = N }; 
        enum { vfort = (I==FortranStyle) }; 
        enum { vcalc = true };
        enum { vstep = S };
        enum { vconj = C };

        // In case N == UNKNOWN
        typedef typename VCopyHelper<T,N,vfort>::type copy_type;

        enum { negS = IntTraits<S>::negS };
        enum { twoS = visreal ? S : IntTraits<S>::twoS };
        enum { notC = !C && viscomplex };
        enum { twoN = visreal ? N : IntTraits<N>::twoS };

        typedef ConstVectorView<T,S,C,I> const_subvector_type;
        typedef ConstVectorView<T,UNKNOWN,C,I> const_subvector_step_type;
        typedef ConstSmallVectorView<T,N,S,C,I> const_view_type;
        typedef ConstSmallVectorView<T,N,S,C,CStyle> const_cview_type;
        typedef ConstSmallVectorView<T,N,S,C,FortranStyle> const_fview_type;
        typedef ConstVectorView<T,UNKNOWN,C> const_xview_type;
        typedef ConstSmallVectorView<T,N,1,C,I> const_unitview_type;
        typedef ConstSmallVectorView<T,N,S,notC,I> const_conjugate_type;
        typedef ConstSmallVectorView<T,N,negS,C,I> const_reverse_type;
        typedef ConstSmallVectorView<real_type,N,twoS,false,I> 
            const_realpart_type;
        typedef const_realpart_type const_imagpart_type;
        typedef ConstSmallVectorView<real_type,twoN,1,false,I> 
            const_flatten_type;
        typedef ConstSmallVectorView<T,N,S,false,I> const_nonconj_type;
        typedef SmallVectorView<T,N,S,C,I> nonconst_type;

        typedef CVIt<T,S,C> const_iterator;
        typedef CVIt<T,negS,C> const_reverse_iterator;
    };

    template <class T, int N, int S, bool C, IndexStyle I>
    class ConstSmallVectorView : 
        public BaseVector_Calc<ConstSmallVectorView<T,N,S,C,I> >
    {
    public:

        typedef ConstSmallVectorView<T,N,S,C,I> type;

        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        enum { vsize = Traits<type>::vsize };
        enum { vfort = Traits<type>::vfort };
        enum { vcalc = Traits<type>::vcalc };
        enum { visreal = Traits<type>::visreal };
        enum { viscomplex = Traits<type>::viscomplex };
        enum { vstep = Traits<type>::vstep };
        enum { vconj = Traits<type>::vconj };

        //
        // Constructors
        //

        inline ConstSmallVectorView(const T* v, size_t n, int s) : 
            itsv(v), itssize(n), itsstep(s) {}

        inline ConstSmallVectorView(const T* v, size_t n) : 
            itsv(v), itssize(n), itsstep(S) 
        { TMVStaticAssert(S != UNKNOWN); }

        inline ConstSmallVectorView(const T* v) :
            itsv(v), itssize(N), itsstep(S) 
        { TMVStaticAssert(N != UNKNOWN); TMVStaticAssert(S != UNKNOWN); }

        inline ConstSmallVectorView(const type& v2) : 
            itsv(v2.cptr()), itssize(v2.size()), itsstep(v2.step()) {}

        inline ConstSmallVectorView(const SmallVectorView<T,N,S,C,I>& v2) : 
            itsv(v2.cptr()), itssize(v2.size()), itsstep(v2.step()) {}

        template <int S2, IndexStyle I2>
        inline ConstSmallVectorView(const ConstVectorView<T,S2,C,I2>& v2) :
            itsv(v2.cptr()), itssize(v2.size()), itsstep(v2.step()) {}

        template <int S2, IndexStyle I2>
        inline ConstSmallVectorView(const VectorView<T,S2,C,I2>& v2) :
            itsv(v2.cptr()), itssize(v2.size()), itsstep(v2.step()) {}

        template <int N2, int S2, IndexStyle I2>
        inline ConstSmallVectorView(
            const ConstSmallVectorView<T,N2,S2,C,I2>& v2) :
            itsv(v2.cptr()), itssize(v2.size()), itsstep(v2.step()) {}

        template <int N2, int S2, IndexStyle I2>
        inline ConstSmallVectorView(const SmallVectorView<T,N2,S2,C,I2>& v2) :
            itsv(v2.cptr()), itssize(v2.size()), itsstep(v2.step()) {}

        inline ~ConstSmallVectorView() {
#ifdef TMV_DEBUG
            itsv = 0; 
#endif
        }

    private :
        inline void operator=(const type& v2);
    public :


        //
        // Auxilliary Functions
        //

        inline const T* cptr() const { return itsv; }
        inline T cref(int i) const  { return DoConj<C>(itsv[i*step()]); }

        inline size_t size() const { return itssize; }
        inline int step() const { return itsstep; }
        inline bool isconj() const { return C; }

    protected :

#ifdef TMV_DEBUG
        const T* itsv;
#else
        const T*const itsv;
#endif
        const CheckedInt<N> itssize;
        const CheckedInt<S> itsstep;

    }; // ConstSmallVectorView

    template <class T, int N, int S, bool C>
    class ConstSmallVectorViewF : 
        public ConstSmallVectorView<T,N,S,C,FortranStyle>
    {
    public:

        typedef ConstSmallVectorViewF<T,N,S,C> type;
        typedef ConstSmallVectorView<T,N,S,C,FortranStyle> vtype;

        inline ConstSmallVectorViewF(const T* v, size_t n, int s) :
            vtype(v,n,s) {}
        inline ConstSmallVectorViewF(const T* v, size_t n) : vtype(v,n) {}
        inline ConstSmallVectorViewF(const T* v) : vtype(v) {}
        inline ConstSmallVectorViewF(const type& v2) : vtype(v2) {}
        inline ConstSmallVectorViewF(
            const SmallVectorView<T,N,S,C,FortranStyle>& v2) : vtype(v2) {}
        template <int S2, IndexStyle I2>
        inline ConstSmallVectorViewF(const ConstVectorView<T,S2,C,I2>& v2) :
            vtype(v2) {}
        template <int S2, IndexStyle I2>
        inline ConstSmallVectorViewF(const VectorView<T,S2,C,I2>& v2) :
            vtype(v2) {}
        template <int N2, int S2, IndexStyle I2>
        inline ConstSmallVectorViewF(
            const ConstSmallVectorView<T,N2,S2,C,I2>& v2) :
            vtype(v2) {}
        template <int N2, int S2, IndexStyle I2>
        inline ConstSmallVectorViewF(const SmallVectorView<T,N2,S2,C,I2>& v2) :
            vtype(v2) {}
        inline ~ConstSmallVectorViewF() {}

    private :
        inline void operator=(const type& v2);
    }; // ConstSmallVectorView


    // 
    // SmallVectorView
    //

    template <class T, int N, int S, bool C, IndexStyle I>
    struct Traits<SmallVectorView<T,N,S,C,I> >
    {
        typedef T value_type;
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;
        enum { visreal = Traits<T>::isreal };
        enum { viscomplex = Traits<T>::iscomplex };

        typedef SmallVectorView<T,N,S,C,I> type;
        typedef const ConstSmallVectorView<T,N,S,C,I> calc_type;
        typedef calc_type eval_type;

        enum { vsize = N }; 
        enum { vfort = (I==FortranStyle) }; 
        enum { vcalc = true };
        enum { vstep = S };
        enum { vconj = C };

        // In case N == UNKNOWN
        typedef typename VCopyHelper<T,N,vfort>::type copy_type;

        enum { negS = IntTraits<S>::negS };
        enum { twoS = visreal ? S : IntTraits<S>::twoS };
        enum { notC = !C && viscomplex };
        enum { twoN = visreal ? N : IntTraits<N>::twoS };

        typedef ConstVectorView<T,S,C,I> const_subvector_type;
        typedef ConstVectorView<T,UNKNOWN,C,I> const_subvector_step_type;
        typedef ConstSmallVectorView<T,N,S,C,I> const_view_type;
        typedef ConstSmallVectorView<T,N,S,C,CStyle> const_cview_type;
        typedef ConstSmallVectorView<T,N,S,C,FortranStyle> const_fview_type;
        typedef ConstVectorView<T,UNKNOWN,C> const_xview_type;
        typedef ConstSmallVectorView<T,N,1,C,I> const_unitview_type;
        typedef ConstSmallVectorView<T,N,S,notC,I> const_conjugate_type;
        typedef ConstSmallVectorView<T,N,negS,C,I> const_reverse_type;
        typedef ConstSmallVectorView<real_type,N,twoS,false,I> 
            const_realpart_type;
        typedef const_realpart_type const_imagpart_type;
        typedef ConstSmallVectorView<real_type,twoN,1,false,I> 
            const_flatten_type;
        typedef ConstSmallVectorView<T,N,S,false,I> const_nonconj_type;
        typedef SmallVectorView<T,N,S,C,I> nonconst_type;

        typedef CVIt<T,S,C> const_iterator;
        typedef CVIt<T,negS,C> const_reverse_iterator;

        typedef typename AuxRef<T,C>::reference reference;

        typedef VectorView<T,S,C,I> subvector_type;
        typedef VectorView<T,UNKNOWN,C,I> subvector_step_type;
        typedef SmallVectorView<T,N,S,C,I> view_type;
        typedef SmallVectorView<T,N,S,C,CStyle> cview_type;
        typedef SmallVectorView<T,N,S,C,FortranStyle> fview_type;
        typedef VectorView<T,UNKNOWN,C> xview_type;
        typedef SmallVectorView<T,N,1,C,I> unitview_type;
        typedef SmallVectorView<T,N,S,notC,I> conjugate_type;
        typedef SmallVectorView<T,N,negS,C,I> reverse_type;
        typedef SmallVectorView<real_type,N,twoS,false,I> realpart_type;
        typedef realpart_type imagpart_type;
        typedef SmallVectorView<real_type,twoN,1,false,I> flatten_type;
        typedef SmallVectorView<T,N,S,false,I> nonconj_type;

        typedef VIt<T,S,C> iterator;
        typedef VIt<T,negS,C> reverse_iterator;
    };

    template <class T, int N, int S, bool C, IndexStyle I>
    class SmallVectorView : 
        public BaseVector_Mutable<SmallVectorView<T,N,S,C,I> >
    {
    public:

        typedef SmallVectorView<T,N,S,C,I> type;
        typedef BaseVector_Mutable<type> base_mut;
        typedef typename base_mut::reference reference;

        enum { vsize = Traits<type>::vsize };
        enum { vfort = Traits<type>::vfort };
        enum { vcalc = Traits<type>::vcalc };
        enum { visreal = Traits<type>::visreal };
        enum { viscomplex = Traits<type>::viscomplex };
        enum { vstep = Traits<type>::vstep };
        enum { vconj = Traits<type>::vconj };

        //
        // Constructors
        //

        inline SmallVectorView(T* v, size_t n, int s) :
            itsv(v), itssize(n), itsstep(s) {}

        inline SmallVectorView(T* v, size_t n) :
            itsv(v), itssize(n), itsstep(S) 
        { TMVAssert(S != UNKNOWN); }

        inline SmallVectorView(T* v) : itsv(v), itssize(N), itsstep(S) 
        { TMVAssert(N != UNKNOWN);  TMVAssert(S != UNKNOWN); }

        inline SmallVectorView(const type& v2) : 
            itsv(v2.itsv), itssize(v2.size()), itsstep(v2.step()) {}

        template <IndexStyle I2>
        inline SmallVectorView(SmallVectorView<T,N,S,C,I2> v2) :
            itsv(v2.ptr()), itssize(v2.size()), itsstep(v2.step()) {}

        template <int S2, IndexStyle I2>
        inline SmallVectorView(VectorView<T,S2,C,I2> v2) :
            itsv(v2.ptr()), itssize(v2.size()), itsstep(v2.step()) {}

        inline ~SmallVectorView() {
#ifdef TMV_DEBUG
            itsv = 0; 
#endif
        }


        //
        // Op =
        //

        template <class V2>
        inline type& operator=(const BaseVector<V2>& v2) 
        { base_mut::operator=(v2); return *this; }

        inline type& operator=(const type& v2)
        { base_mut::operator=(v2); return *this; }

        //
        // Auxilliary Functions
        //

        inline const T* cptr() const { return itsv; }
        inline T cref(int i) const  { return DoConj<C>(itsv[i*step()]); }
        inline T* ptr() { return itsv; }
        inline reference ref(int i) { return reference(itsv[i*step()]); }

        inline size_t size() const { return itssize; }
        inline int step() const { return itsstep; }
        inline bool isconj() const { return C; }

    protected :

#ifdef TMV_DEBUG
        T* itsv;
#else
        T*const itsv;
#endif
        const CheckedInt<N> itssize;
        const CheckedInt<S> itsstep;

    }; // SmallVectorView

    template <class T, int N, int S, bool C>
    class SmallVectorViewF : public SmallVectorView<T,N,S,C,FortranStyle>
    {
    public:
        typedef SmallVectorViewF<T,N,S,C> type;
        typedef SmallVectorView<T,N,S,C,FortranStyle> vtype;

        inline SmallVectorViewF(T* v, size_t n=N, int s=S) : vtype(v,n,s) {}
        inline SmallVectorViewF(T* v, size_t n=N) : vtype(v,n) {}
        inline SmallVectorViewF(T* v) : vtype(v) {}
        inline SmallVectorViewF(const type& v2) : vtype(v2) {}
        template <IndexStyle I2>
        inline SmallVectorViewF(SmallVectorView<T,N,S,C,I2> v2) : vtype(v2) {}
        template <int S2, IndexStyle I2>
        inline SmallVectorViewF(VectorView<T,S2,C,I2> v2) : vtype(v2) {}
        inline ~SmallVectorViewF() {}

        template <class V2>
        inline type& operator=(const BaseVector<V2>& v2) 
        { vtype::operator=(v2); return *this; }
        inline type& operator=(const type& v2)
        { vtype::operator=(v2); return *this; }

    }; // SmallVectorViewF


    // 
    // Swap
    //

    template <class V, class T, int N, int S, bool C, IndexStyle I>
    inline void Swap(BaseVector_Mutable<V>& v1, SmallVectorView<T,N,S,C,I> v2)
    { DoSwap(v1,v2); }
    template <class V, class T, int N, int S, bool C, IndexStyle I>
    inline void Swap(SmallVectorView<T,N,S,C,I> v1, BaseVector_Mutable<V>& v2)
    { DoSwap(v1,v2); }
    template <class T, int N, int S1, bool C1, IndexStyle I1, 
              int S2, bool C2, IndexStyle I2>
    inline void Swap(
        SmallVectorView<T,N,S1,C1,I1> v1, SmallVectorView<T,N,S2,C2,I2> v2)
    { DoSwap(v1,v2); }
    template <class T, int N, int S1, bool C1, IndexStyle I1, 
              int S2, bool C2, IndexStyle I2>
    inline void Swap(
        VectorView<T,S1,C1,I1> v1, SmallVectorView<T,N,S2,C2,I2> v2)
    { DoSwap(v1,v2); }
    template <class T, int N, int S1, bool C1, IndexStyle I1, 
              int S2, bool C2, IndexStyle I2>
    inline void Swap(
        SmallVectorView<T,N,S1,C1,I1> v1, VectorView<T,S2,C2,I2> v2)
    { DoSwap(v1,v2); }

    //
    // TMV_Text functions
    //

    template <class T, int N, IndexStyle I> 
    inline std::string TMV_Text(const SmallVector<T,N,I>& )
    {
        std::ostringstream s;
        s << "SmallVector<"<<TMV_Text(T())<<","<<N<<","<<TMV_Text(I)<<">";
        return s.str();
    }

    template <class T, int N, int S, bool C, IndexStyle I>
    inline std::string TMV_Text(const SmallVectorView<T,N,S,C,I>& v)
    {
        std::ostringstream s;
        s << "SmallVectorView<"<<TMV_Text(T());
        s << ","<<IntTraits<N>::text();
        if (N == UNKNOWN) s << "("<<v.size()<<")";
        s << ","<<IntTraits<S>::text();
        if (S == UNKNOWN) s << "("<<v.step()<<")";
        s << ","<<C<<","<<TMV_Text(I)<<">";
        return s.str();
    }

    template <class T, int N, int S, bool C, IndexStyle I>
    inline std::string TMV_Text(const ConstSmallVectorView<T,N,S,C,I>& v)
    {
        std::ostringstream s;
        s << "ConstSmallVectorView<"<<TMV_Text(T());
        s << ","<<IntTraits<N>::text();
        if (N == UNKNOWN) s << "("<<v.size()<<")";
        s << ","<<IntTraits<S>::text();
        if (S == UNKNOWN) s << "("<<v.step()<<")";
        s << ","<<C<<","<<TMV_Text(I)<<">";
        return s.str();
    }

} // namespace tmv

#endif
