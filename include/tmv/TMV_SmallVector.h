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
// heap is pretty irrelevant.  (See TMV_Array.h to change this.)
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
//    SmallVector<T,N,A>()  
//        Makes a Vector of size N with _uninitialized_ values
//
//    SmallVector<T,N,A>(T x)
//        Makes a Vector of size N with all values = x
//
//    SmallVector<T,N,A>(const T* v2)
//    SmallVector<T,N,A>(const vector<T>& v2)
//    SmallVector<T,N,A>(const BaseVector<V2>& v2)
//        Makes a SmallVector which copies the elements of v2.
//

#ifndef TMV_SmallVector_H
#define TMV_SmallVector_H


#include <vector>
#include "TMV_BaseVector.h"
#include "TMV_VIt.h"
#include "TMV_Array.h"

namespace tmv {

    //
    // SmallVector
    //

    template <class T, int N, int A0>
    struct Traits<SmallVector<T,N,A0> >
    {
        enum { A = (A0 & ~NoAlias) | Unit };
        enum { okA = (
                Attrib<A>::vectoronly &&
                !Attrib<A>::noalias &&
                !Attrib<A>::conj )};
        enum { _attrib = A };

        typedef T value_type;
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;
        enum { isreal = Traits<T>::isreal };
        enum { iscomplex = Traits<T>::iscomplex };

        typedef SmallVector<T,N,A0> type;
        typedef const type& calc_type; 
        typedef const type& eval_type; 
        typedef type copy_type;

        enum { _size = N }; 
        enum { _fort = Attrib<A>::fort };
        enum { _calc = true };
        enum { _step = 1 }; 
        enum { _conj = false }; 
        enum { _checkalias = Attrib<A>::checkalias };
        enum { _unit = true };
        enum { twoS = isreal ? 1 : 2 };

        enum { unitA = A };
        enum { nonunitA = A & ~Unit };
        enum { conjA = iscomplex ? (A ^ Conj) : int(A) };
        enum { nonconjA = A };
        enum { cstyleA = A & ~FortranStyle };
        enum { fstyleA = A | FortranStyle };
        enum { twosA = isreal ? int(A) : (A & ~Conj & ~Unit) };
        enum { An = _checkalias ? (A & ~CheckAlias) : (A | NoAlias) };
        enum { nonunitAn = An & ~Unit };

        typedef ConstVectorView<T,An> const_subvector_type;
        typedef ConstVectorView<T,nonunitAn> const_subvector_step_type;
        typedef ConstSmallVectorView<T,N,1,A> const_view_type;
        typedef ConstSmallVectorView<T,N,1,cstyleA> const_cview_type;
        typedef ConstSmallVectorView<T,N,1,fstyleA> const_fview_type;
        typedef ConstVectorView<T> const_xview_type;
        typedef ConstSmallVectorView<T,N,1,unitA> const_unitview_type;
        typedef ConstSmallVectorView<T,N,1,conjA> const_conjugate_type;
        typedef ConstSmallVectorView<T,N,-1,nonunitA> const_reverse_type;
        typedef ConstSmallVectorView<real_type,N,twoS,twosA>
            const_realpart_type;
        typedef const_realpart_type const_imagpart_type;
        typedef ConstSmallVectorView<real_type,isreal?N:N*2,1,unitA>
            const_flatten_type;
        typedef ConstSmallVectorView<T,N,1,nonconjA> const_nonconj_type;
        typedef SmallVectorView<T,N,1,A> nonconst_type;

        typedef CVIt<T,1,false> const_iterator;
        typedef CVIt<T,-1,false> const_reverse_iterator;

        typedef T& reference;

        typedef VectorView<T,An> subvector_type;
        typedef VectorView<T,nonunitAn> subvector_step_type;
        typedef SmallVectorView<T,N,1,A> view_type;
        typedef SmallVectorView<T,N,1,cstyleA> cview_type;
        typedef SmallVectorView<T,N,1,fstyleA> fview_type;
        typedef VectorView<T> xview_type;
        typedef SmallVectorView<T,N,1,unitA> unitview_type;
        typedef SmallVectorView<T,N,1,conjA> conjugate_type;
        typedef SmallVectorView<T,N,-1,nonunitA> reverse_type;
        typedef SmallVectorView<real_type,N,twoS,twosA> realpart_type;
        typedef realpart_type imagpart_type;
        typedef SmallVectorView<real_type,isreal?N:N*2,1,unitA>
            flatten_type;
        typedef SmallVectorView<T,N,1,A> nonconj_type;

        typedef VIt<T,1,false> iterator;
        typedef VIt<T,-1,false> reverse_iterator;
    };

    template <class T, int N, int A>
    class SmallVector : 
        public BaseVector_Mutable<SmallVector<T,N,A> >
    {
    public:

        typedef SmallVector<T,N,A> type;
        typedef BaseVector_Mutable<type> base_mut;
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        typedef typename base_mut::iterator iterator;

        enum { _size = Traits<type>::_size };
        enum { _fort = Traits<type>::_fort };
        enum { _calc = Traits<type>::_calc };
        enum { isreal = Traits<type>::isreal };
        enum { iscomplex = Traits<type>::iscomplex };
        enum { _step = Traits<type>::_step };
        enum { _conj = Traits<type>::_conj };
        enum { _checkalias = Traits<type>::_checkalias };
        enum { _unit = Traits<type>::_unit };
        enum { _attrib = Traits<type>::_attrib };

        //
        // Constructors
        //

        TMV_INLINE_ND SmallVector()
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(N >= 0);
#ifdef TMV_DEBUG
            this->flatten().setAllTo(Traits<real_type>::constr_value());
#endif
        }

        explicit SmallVector(T x) 
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(N >= 0);
            this->setAllTo(x); 
        }

        explicit SmallVector(const T* v) 
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(N >= 0);
            ConstSmallVectorView<T,N,1>(v).newAssignTo(*this);
        }

        explicit SmallVector(const std::vector<T>& v2) 
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(N >= 0);
            TMVAssert(v2.size() == N);
            ConstSmallVectorView<T,N,1>(&v2[0]).newAssignTo(*this);
        }

        SmallVector(const type& v2) 
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(N >= 0);
            v2.newAssignTo(*this);
        }

        template <class V2>
        SmallVector(const BaseVector<V2>& v2) 
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(N >= 0);
            TMVAssert(v2.size() == N);
            v2.newAssignTo(*this);
        }

        TMV_INLINE_ND ~SmallVector()
        {
#ifdef TMV_DEBUG
            this->flatten().setAllTo(Traits<real_type>::destr_value());
#endif
        }

        //
        // Op =
        //

        TMV_INLINE type& operator=(const type& v2)
        { if (this != &v2) base_mut::operator=(v2); return *this; }

        template <class V2>
        TMV_INLINE type& operator=(const BaseVector<V2>& v2) 
        { base_mut::operator=(v2); return *this; }


        //
        // Auxilliary Functions
        //

        TMV_INLINE const T* cptr() const { return itsv; }
        TMV_INLINE T* ptr() { return itsv; }
        T cref(int i) const  { return itsv[i]; }
        T& ref(int i) { return itsv[i]; }

        TMV_INLINE size_t size() const { return N; }
        TMV_INLINE int nElements() const { return N; }
        TMV_INLINE int step() const { return 1; }
        TMV_INLINE bool isconj() const { return false; }


    protected :

        StackArray<T,N> itsv;

    }; // SmallVector


    //
    // ConstSmallVectorView
    //

    template <class T, int N, int S, int A0>
    struct Traits<ConstSmallVectorView<T,N,S,A0> >
    {
        enum { A = (A0 & ~NoAlias) | (
                ( (S == 1) ? Unit : 0 ) )};
        enum { okA = (
                Attrib<A>::vectoronly &&
                !Attrib<A>::noalias &&
                ( Traits<T>::iscomplex || !Attrib<A>::conj ) &&
                ( Attrib<A>::unit == (S == 1) ) )};
        enum { _attrib = A };

        typedef T value_type;
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;
        enum { isreal = Traits<T>::isreal };
        enum { iscomplex = Traits<T>::iscomplex };

        typedef ConstSmallVectorView<T,N,S,A0> type;
        typedef const type& calc_type;
        typedef const type& eval_type;

        enum { _size = N }; 
        enum { _fort = Attrib<A>::fort };
        enum { _calc = true };
        enum { _step = S };
        enum { _conj = Attrib<A>::conj };
        enum { _checkalias = Attrib<A>::checkalias };
        enum { _unit = Attrib<A>::unit };
        enum { negS = IntTraits<S>::negS };
        enum { twoS = isreal ? S : IntTraits<S>::twoS };
        enum { twoN = isreal ? N : IntTraits<N>::twoS };

        // Use VCopyHelper for copy_type in case N == UNKNOWN
        typedef typename VCopyHelper<T,N,_fort>::type copy_type;

        enum { unitA = A | Unit };
        enum { nonunitA = A & ~Unit };
        enum { conjA = iscomplex ? (A ^ Conj) : int(A) };
        enum { nonconjA = A & ~Conj };
        enum { cstyleA = A & ~FortranStyle };
        enum { fstyleA = A | FortranStyle };
        enum { revA = (S == -1) ? int(unitA) : nonunitA };
        enum { twosA = isreal ? int(A) : (A & ~Conj & ~Unit) };
        enum { flatA = isreal ? int(A) : (unitA & ~Conj) };
        enum { An = _checkalias ? (A & ~CheckAlias) : (A | NoAlias) };
        enum { nonunitAn = An & ~Unit };

        typedef ConstVectorView<T,An> const_subvector_type;
        typedef ConstVectorView<T,nonunitAn> const_subvector_step_type;
        typedef ConstSmallVectorView<T,N,S,A> const_view_type;
        typedef ConstSmallVectorView<T,N,S,cstyleA> const_cview_type;
        typedef ConstSmallVectorView<T,N,S,fstyleA> const_fview_type;
        typedef ConstVectorView<T,_conj ? Conj : NonConj> const_xview_type;
        typedef ConstSmallVectorView<T,N,1,unitA> const_unitview_type;
        typedef ConstSmallVectorView<T,N,S,conjA> const_conjugate_type;
        typedef ConstSmallVectorView<T,N,negS,revA> const_reverse_type;
        typedef ConstSmallVectorView<real_type,N,twoS,twosA>
            const_realpart_type;
        typedef const_realpart_type const_imagpart_type;
        typedef ConstSmallVectorView<real_type,twoN,1,flatA>
            const_flatten_type;
        typedef ConstSmallVectorView<T,N,S,nonconjA> const_nonconj_type;
        typedef SmallVectorView<T,N,S,A> nonconst_type;

        typedef CVIt<T,S,_conj> const_iterator;
        typedef CVIt<T,negS,_conj> const_reverse_iterator;
    };

    template <class T, int N, int S, int A>
    class ConstSmallVectorView : 
        public BaseVector_Calc<ConstSmallVectorView<T,N,S,A> >
    {
    public:

        typedef ConstSmallVectorView<T,N,S,A> type;

        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        enum { _size = Traits<type>::_size };
        enum { _fort = Traits<type>::_fort };
        enum { _calc = Traits<type>::_calc };
        enum { isreal = Traits<type>::isreal };
        enum { iscomplex = Traits<type>::iscomplex };
        enum { _step = Traits<type>::_step };
        enum { _conj = Traits<type>::_conj };
        enum { _checkalias = Traits<type>::_checkalias };
        enum { _unit = Traits<type>::_unit };
        enum { _attrib = Traits<type>::_attrib };

        //
        // Constructors
        //

        ConstSmallVectorView(const T* v, size_t n, int s) : 
            itsv(v), itssize(n), itsstep(s)
        {
            TMVStaticAssert(Traits<type>::okA);
        }

        ConstSmallVectorView(const T* v, size_t n) : 
            itsv(v), itssize(n), itsstep(S) 
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(S != UNKNOWN); 
        }

        ConstSmallVectorView(const T* v) :
            itsv(v), itssize(N), itsstep(S) 
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(N != UNKNOWN); TMVStaticAssert(S != UNKNOWN); 
        }

        ConstSmallVectorView(const type& v2) : 
            itsv(v2.cptr()), itssize(v2.size()), itsstep(v2.step()) 
        {
            TMVStaticAssert(Traits<type>::okA);
        }

        ConstSmallVectorView(const SmallVectorView<T,N,S,A>& v2) : 
            itsv(v2.cptr()), itssize(v2.size()), itsstep(v2.step()) 
        {
            TMVStaticAssert(Traits<type>::okA);
        }

        template <int A2>
        ConstSmallVectorView(const ConstVectorView<T,A2>& v2) :
            itsv(v2.cptr()), itssize(v2.size()), itsstep(v2.step()) 
        {
            TMVStaticAssert(Traits<type>::okA);
        }

        template <int A2>
        ConstSmallVectorView(const VectorView<T,A2>& v2) :
            itsv(v2.cptr()), itssize(v2.size()), itsstep(v2.step()) 
        {
            TMVStaticAssert(Traits<type>::okA);
        }

        template <int N2, int S2, int A2>
        ConstSmallVectorView(
            const ConstSmallVectorView<T,N2,S2,A2>& v2) :
            itsv(v2.cptr()), itssize(v2.size()), itsstep(v2.step()) 
        {
            TMVStaticAssert(Traits<type>::okA);
        }

        template <int N2, int S2, int A2>
        ConstSmallVectorView(const SmallVectorView<T,N2,S2,A2>& v2) :
            itsv(v2.cptr()), itssize(v2.size()), itsstep(v2.step()) 
        {
            TMVStaticAssert(Traits<type>::okA);
        }

        ~ConstSmallVectorView() {
#ifdef TMV_DEBUG
            itsv = 0; 
#endif
        }

    private :
        void operator=(const type& v2);
    public :


        //
        // Auxilliary Functions
        //

        TMV_INLINE const T* cptr() const { return itsv; }
        T cref(int i) const  { return DoConj<_conj>(itsv[i*step()]); }

        TMV_INLINE size_t size() const { return itssize; }
        TMV_INLINE int nElements() const { return itssize; }
        TMV_INLINE int step() const { return itsstep; }
        TMV_INLINE bool isconj() const { return _conj; }

    protected :

        const T* itsv;
        const CheckedInt<N> itssize;
        const CheckedInt<S> itsstep;

    }; // ConstSmallVectorView


    // 
    // SmallVectorView
    //

    template <class T, int N, int S, int A0>
    struct Traits<SmallVectorView<T,N,S,A0> >
    {
        enum { A = (A0 & ~NoAlias) | (
                ( (S == 1) ? Unit : 0 ) )};
        enum { okA = (
                Attrib<A>::vectoronly &&
                !Attrib<A>::noalias &&
                ( Traits<T>::iscomplex || !Attrib<A>::conj ) &&
                ( Attrib<A>::unit == (S == 1) ) )};
        enum { _attrib = A };

        typedef T value_type;
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;
        enum { isreal = Traits<T>::isreal };
        enum { iscomplex = Traits<T>::iscomplex };

        typedef SmallVectorView<T,N,S,A0> type;
        typedef ConstSmallVectorView<T,N,S,A> calc_type;
        typedef calc_type eval_type;

        enum { _size = N }; 
        enum { _fort = Attrib<A>::fort };
        enum { _calc = true };
        enum { _step = S };
        enum { _conj = Attrib<A>::conj };
        enum { _checkalias = Attrib<A>::checkalias };
        enum { _unit = Attrib<A>::unit };
        enum { negS = IntTraits<S>::negS };
        enum { twoS = isreal ? S : IntTraits<S>::twoS };
        enum { twoN = isreal ? N : IntTraits<N>::twoS };

        typedef typename VCopyHelper<T,N,_fort>::type copy_type;

        enum { unitA = A | Unit };
        enum { nonunitA = A & ~Unit };
        enum { conjA = iscomplex ? (A ^ Conj) : int(A) };
        enum { nonconjA = A & ~Conj };
        enum { cstyleA = A & ~FortranStyle };
        enum { fstyleA = A | FortranStyle };
        enum { revA = (S == -1) ? int(unitA) : nonunitA };
        enum { twosA = isreal ? int(A) : (A & ~Conj & ~Unit) };
        enum { flatA = isreal ? int(A) : (unitA & ~Conj) };
        enum { An = _checkalias ? (A & ~CheckAlias) : (A | NoAlias) };
        enum { nonunitAn = An & ~Unit };

        typedef ConstVectorView<T,An> const_subvector_type;
        typedef ConstVectorView<T,nonunitAn> const_subvector_step_type;
        typedef ConstSmallVectorView<T,N,S,A> const_view_type;
        typedef ConstSmallVectorView<T,N,S,cstyleA> const_cview_type;
        typedef ConstSmallVectorView<T,N,S,fstyleA> const_fview_type;
        typedef ConstVectorView<T,_conj ? Conj : NonConj> const_xview_type;
        typedef ConstSmallVectorView<T,N,1,unitA> const_unitview_type;
        typedef ConstSmallVectorView<T,N,S,conjA> const_conjugate_type;
        typedef ConstSmallVectorView<T,N,negS,revA> const_reverse_type;
        typedef ConstSmallVectorView<real_type,N,twoS,twosA>
            const_realpart_type;
        typedef const_realpart_type const_imagpart_type;
        typedef ConstSmallVectorView<real_type,twoN,1,flatA>
            const_flatten_type;
        typedef ConstSmallVectorView<T,N,S,nonconjA> const_nonconj_type;
        typedef SmallVectorView<T,N,S,A> nonconst_type;

        typedef CVIt<T,S,_conj> const_iterator;
        typedef CVIt<T,negS,_conj> const_reverse_iterator;
    
        typedef typename AuxRef<T,_conj>::reference reference;

        typedef VectorView<T,An> subvector_type;
        typedef VectorView<T,nonunitAn> subvector_step_type;
        typedef SmallVectorView<T,N,S,A> view_type;
        typedef SmallVectorView<T,N,S,cstyleA> cview_type;
        typedef SmallVectorView<T,N,S,fstyleA> fview_type;
        typedef VectorView<T,_conj ? Conj : NonConj> xview_type;
        typedef SmallVectorView<T,N,1,unitA> unitview_type;
        typedef SmallVectorView<T,N,S,conjA> conjugate_type;
        typedef SmallVectorView<T,N,negS,revA> reverse_type;
        typedef SmallVectorView<real_type,N,twoS,twosA> realpart_type;
        typedef realpart_type imagpart_type;
        typedef SmallVectorView<real_type,twoN,1,flatA> flatten_type;
        typedef SmallVectorView<T,N,S,nonconjA> nonconj_type;

        typedef VIt<T,S,_conj> iterator;
        typedef VIt<T,negS,_conj> reverse_iterator;
    };

    template <class T, int N, int S, int A>
    class SmallVectorView : 
        public BaseVector_Mutable<SmallVectorView<T,N,S,A> >
    {
    public:

        typedef SmallVectorView<T,N,S,A> type;
        typedef BaseVector_Mutable<type> base_mut;
        typedef typename base_mut::reference reference;

        enum { _size = Traits<type>::_size };
        enum { _fort = Traits<type>::_fort };
        enum { _calc = Traits<type>::_calc };
        enum { isreal = Traits<type>::isreal };
        enum { iscomplex = Traits<type>::iscomplex };
        enum { _step = Traits<type>::_step };
        enum { _conj = Traits<type>::_conj };
        enum { _checkalias = Traits<type>::_checkalias };
        enum { _unit = Traits<type>::_unit };
        enum { _attrib = Traits<type>::_attrib };

        //
        // Constructors
        //

        SmallVectorView(T* v, size_t n, int s) :
            itsv(v), itssize(n), itsstep(s) 
        {
            TMVStaticAssert(Traits<type>::okA);
        }

        SmallVectorView(T* v, size_t n) :
            itsv(v), itssize(n), itsstep(S) 
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(S != UNKNOWN); 
        }

        SmallVectorView(T* v) : itsv(v), itssize(N), itsstep(S) 
        { 
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(N != UNKNOWN);  TMVStaticAssert(S != UNKNOWN); 
        }

        SmallVectorView(const type& v2) : 
            itsv(v2.itsv), itssize(v2.size()), itsstep(v2.step())
        {
            TMVStaticAssert(Traits<type>::okA);
        }

        template <int N2, int S2, int A2>
        SmallVectorView(SmallVectorView<T,N2,S2,A2> v2) :
            itsv(v2.ptr()), itssize(v2.size()), itsstep(v2.step())
        {
            TMVStaticAssert(Traits<type>::okA);
        }

        template <int A2>
        SmallVectorView(VectorView<T,A2> v2) :
            itsv(v2.ptr()), itssize(v2.size()), itsstep(v2.step())
        {
            TMVStaticAssert(Traits<type>::okA);
        }

        ~SmallVectorView() {
#ifdef TMV_DEBUG
            itsv = 0; 
#endif
        }


        //
        // Op =
        //

        TMV_INLINE type& operator=(const type& v2)
        { if (this != &v2) base_mut::operator=(v2); return *this; }

        template <class V2>
        TMV_INLINE type& operator=(const BaseVector<V2>& v2) 
        { base_mut::operator=(v2); return *this; }

        //
        // Auxilliary Functions
        //

        TMV_INLINE const T* cptr() const { return itsv; }
        T cref(int i) const  { return DoConj<_conj>(itsv[i*step()]); }
        TMV_INLINE T* ptr() { return itsv; }
        reference ref(int i) { return reference(itsv[i*step()]); }

        TMV_INLINE size_t size() const { return itssize; }
        TMV_INLINE int nElements() const { return itssize; }
        TMV_INLINE int step() const { return itsstep; }
        TMV_INLINE bool isconj() const { return _conj; }

    protected :

        T* itsv;
        const CheckedInt<N> itssize;
        const CheckedInt<S> itsstep;

    }; // SmallVectorView


    // 
    // Swap
    //

    template <class V, class T, int N, int S, int A>
    static TMV_INLINE void Swap(
        BaseVector_Mutable<V>& v1, SmallVectorView<T,N,S,A> v2)
    { DoSwap(v1,v2); }
    template <class V, class T, int N, int S, int A>
    static TMV_INLINE void Swap(
        SmallVectorView<T,N,S,A> v1, BaseVector_Mutable<V>& v2)
    { DoSwap(v1,v2); }
    template <class T, int N, int S1, int A1, int S2, int A2>
    static TMV_INLINE void Swap(
        SmallVectorView<T,N,S1,A1> v1, SmallVectorView<T,N,S2,A2> v2)
    { DoSwap(v1,v2); }
    template <class T, int N, int A1, int S2, int A2>
    static TMV_INLINE void Swap(
        VectorView<T,A1> v1, SmallVectorView<T,N,S2,A2> v2)
    { DoSwap(v1,v2); }
    template <class T, int N, int S1, int A1, int A2>
    static TMV_INLINE void Swap(
        SmallVectorView<T,N,S1,A1> v1, VectorView<T,A2> v2)
    { DoSwap(v1,v2); }


    //
    // Conjugate
    //
    
    template <class T, int N, int A>
    static TMV_INLINE typename SmallVector<T,N,A>::conjugate_type Conjugate(
        SmallVector<T,N,A>& v)
    { return v.conjugate(); }
    template <class T, int N, int S, int A>
    static TMV_INLINE typename SmallVectorView<T,N,S,A>::conjugate_type Conjugate(
        SmallVectorView<T,N,S,A> v)
    { return v.conjugate(); }


    //
    // TMV_Text functions
    //

#ifdef TMV_TEXT
    template <class T, int N, int A>
    static inline std::string TMV_Text(const SmallVector<T,N,A>& v)
    {
        std::ostringstream s;
        s << "SmallVector<"<<TMV_Text(T());
        s << ","<<N<<","<<Attrib<A>::vtext()<<">";
        s << "("<<v.size()<<","<<v.step()<<")";
        return s.str();
    }

    template <class T, int N, int S, int A>
    static inline std::string TMV_Text(const SmallVectorView<T,N,S,A>& v)
    {
        std::ostringstream s;
        s << "SmallVectorView<"<<TMV_Text(T());
        s << ","<<IntTraits<N>::text();
        s << ","<<IntTraits<S>::text();
        s << ","<<Attrib<A>::vtext()<<">";
        s << "("<<v.size()<<","<<v.step()<<")";
        return s.str();
    }

    template <class T, int N, int S, int A>
    static inline std::string TMV_Text(const ConstSmallVectorView<T,N,S,A>& v)
    {
        std::ostringstream s;
        s << "ConstSmallVectorView<"<<TMV_Text(T());
        s << ","<<IntTraits<N>::text();
        s << ","<<IntTraits<S>::text();
        s << ","<<Attrib<A>::vtext()<<">";
        s << "("<<v.size()<<","<<v.step()<<")";
        return s.str();
    }
#endif

} // namespace tmv

#endif
