

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
//    SmallVector<T,N,A>(const BaseVector<V2>& v2)
//        Makes a SmallVector which copies the elements of v2.
//

#ifndef TMV_SmallVector_H
#define TMV_SmallVector_H


#include "TMV_BaseVector.h"
#include "TMV_VIt.h"
#include "TMV_Array.h"

namespace tmv {

    //
    // SmallVector
    //

    template <class T, ptrdiff_t N, int A0>
    struct Traits<SmallVector<T,N,A0> >
    {
        enum { A = (A0 & ~NoAlias) | Unit };
        enum { okA = (
                Attrib<A>::vectoronly &&
                !Attrib<A>::conj &&
                !Attrib<A>::noalias )};
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
        enum { Ar = _checkalias ? (A & ~CheckAlias) : (A | NoAlias) };
        enum { nonunitAr = Ar & ~Unit };
        enum { An = (A & ~CheckAlias) };

        typedef ConstVectorView<T,Ar> const_subvector_type;
        typedef ConstVectorView<T,nonunitAr> const_subvector_step_type;
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

        typedef CVIt<T,1,NonConj> const_iterator;
        typedef CVIt<T,-1,NonConj> const_reverse_iterator;

        typedef T& reference;

        typedef VectorView<T,Ar> subvector_type;
        typedef VectorView<T,nonunitAr> subvector_step_type;
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
        typedef SmallVectorView<T,N,1,An> noalias_type;
        typedef SmallVectorView<T,N,1,An|CheckAlias> alias_type;

        typedef VIt<T,1,NonConj> iterator;
        typedef VIt<T,-1,NonConj> reverse_iterator;
    };

    template <class T, ptrdiff_t N, int A>
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
#ifdef TMV_EXTRA_DEBUG
            this->setAllTo(Traits<T>::constr_value());
#endif
        }

        explicit SmallVector(T x) 
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(N >= 0);
            this->setAllTo(x); 
        }

        SmallVector(const type& v2) 
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(N >= 0);
            typename Traits<type>::noalias_type na = this->noAlias();
            v2.assignTo(na);
        }

        template <class V2>
        SmallVector(const BaseVector<V2>& v2) 
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(N >= 0);
            TMVAssert(v2.size() == N);
            typename Traits<type>::noalias_type na = this->noAlias();
            v2.assignTo(na);
        }

        TMV_INLINE_ND ~SmallVector()
        {
#ifdef TMV_EXTRA_DEBUG
            this->setAllTo(Traits<T>::destr_value());
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
        T cref(ptrdiff_t i) const  { return itsv[i]; }
        T& ref(ptrdiff_t i) { return itsv[i]; }

        TMV_INLINE ptrdiff_t size() const { return N; }
        TMV_INLINE ptrdiff_t nElements() const { return N; }
        TMV_INLINE ptrdiff_t step() const { return 1; }
        TMV_INLINE bool isconj() const { return false; }


    protected :

        StackArray<T,N> itsv;

    }; // SmallVector


    //
    // ConstSmallVectorView
    //

    template <class T, ptrdiff_t N, ptrdiff_t S, int A0>
    struct Traits<ConstSmallVectorView<T,N,S,A0> >
    {
        enum { A = (A0 & ~NoAlias) | (
                ( (S == 1) ? Unit : 0 ) )};
        enum { okA = (
                Attrib<A>::vectoronly &&
                ( Traits<T>::iscomplex || !Attrib<A>::conj ) &&
                !Attrib<A>::noalias &&
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

        enum { known = N != Unknown };
        enum { copyA = _fort ? FortranStyle : CStyle };
        typedef typename TypeSelect<known, 
                SmallVector<T,N,copyA>,
                Vector<T,copyA|NoAlias> >::type copy_type;

        enum { unitA = A | Unit };
        enum { nonunitA = A & ~Unit };
        enum { conjA = iscomplex ? (A ^ Conj) : int(A) };
        enum { nonconjA = A & ~Conj };
        enum { cstyleA = A & ~FortranStyle };
        enum { fstyleA = A | FortranStyle };
        enum { revA = (S == -1) ? int(unitA) : nonunitA };
        enum { twosA = isreal ? int(A) : (A & ~Conj & ~Unit) };
        enum { flatA = isreal ? int(A) : (unitA & ~Conj) };
        enum { Ar = _checkalias ? (A & ~CheckAlias) : (A | NoAlias) };
        enum { nonunitAr = Ar & ~Unit };

        typedef ConstVectorView<T,Ar> const_subvector_type;
        typedef ConstVectorView<T,nonunitAr> const_subvector_step_type;
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

        typedef CVIt<T,S,_conj?Conj:NonConj> const_iterator;
        typedef CVIt<T,negS,_conj?Conj:NonConj> const_reverse_iterator;
    };

    template <class T, ptrdiff_t N, ptrdiff_t S, int A>
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

        ConstSmallVectorView(const T* v, ptrdiff_t n, ptrdiff_t s) : 
            itsv(v), itssize(n), itsstep(s)
        {
            TMVStaticAssert(Traits<type>::okA);
        }

        ConstSmallVectorView(const T* v, ptrdiff_t n) : 
            itsv(v), itssize(n), itsstep(S) 
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(S != Unknown); 
        }

        ConstSmallVectorView(const T* v) :
            itsv(v), itssize(N), itsstep(S) 
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(N != Unknown); TMVStaticAssert(S != Unknown); 
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
            TMVStaticAssert(Attrib<A>::conj == int(Attrib<A2>::conj)); 
        }

        template <int A2>
        ConstSmallVectorView(const VectorView<T,A2>& v2) :
            itsv(v2.cptr()), itssize(v2.size()), itsstep(v2.step()) 
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Attrib<A>::conj == int(Attrib<A2>::conj)); 
        }

        template <ptrdiff_t N2, ptrdiff_t S2, int A2>
        ConstSmallVectorView(
            const ConstSmallVectorView<T,N2,S2,A2>& v2) :
            itsv(v2.cptr()), itssize(v2.size()), itsstep(v2.step()) 
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Attrib<A>::conj == int(Attrib<A2>::conj)); 
        }

        template <ptrdiff_t N2, ptrdiff_t S2, int A2>
        ConstSmallVectorView(const SmallVectorView<T,N2,S2,A2>& v2) :
            itsv(v2.cptr()), itssize(v2.size()), itsstep(v2.step()) 
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Attrib<A>::conj == int(Attrib<A2>::conj)); 
        }

        ~ConstSmallVectorView() 
        {
#ifdef TMV_EXTRA_DEBUG
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
        T cref(ptrdiff_t i) const  { return DoConj<_conj>(itsv[i*step()]); }

        TMV_INLINE ptrdiff_t size() const { return itssize; }
        TMV_INLINE ptrdiff_t nElements() const { return itssize; }
        TMV_INLINE ptrdiff_t step() const { return itsstep; }
        TMV_INLINE bool isconj() const { return _conj; }

    protected :

        const T* itsv;
        const CheckedInt<N> itssize;
        const CheckedInt<S> itsstep;

    }; // ConstSmallVectorView


    // 
    // SmallVectorView
    //

    template <class T, ptrdiff_t N, ptrdiff_t S, int A0>
    struct Traits<SmallVectorView<T,N,S,A0> >
    {
        enum { A = (A0 & ~NoAlias) | (
                ( (S == 1) ? Unit : 0 ) )};
        enum { okA = (
                Attrib<A>::vectoronly &&
                ( Traits<T>::iscomplex || !Attrib<A>::conj ) &&
                !Attrib<A>::noalias &&
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

        enum { known = N != Unknown };
        enum { copyA = _fort ? FortranStyle : CStyle };
        typedef typename TypeSelect<known, 
                SmallVector<T,N,copyA>,
                Vector<T,copyA|NoAlias> >::type copy_type;

        enum { unitA = A | Unit };
        enum { nonunitA = A & ~Unit };
        enum { conjA = iscomplex ? (A ^ Conj) : int(A) };
        enum { nonconjA = A & ~Conj };
        enum { cstyleA = A & ~FortranStyle };
        enum { fstyleA = A | FortranStyle };
        enum { revA = (S == -1) ? int(unitA) : nonunitA };
        enum { twosA = isreal ? int(A) : (A & ~Conj & ~Unit) };
        enum { flatA = isreal ? int(A) : (unitA & ~Conj) };
        enum { Ar = _checkalias ? (A & ~CheckAlias) : (A | NoAlias) };
        enum { nonunitAr = Ar & ~Unit };
        enum { An = (A & ~CheckAlias) };

        typedef ConstVectorView<T,Ar> const_subvector_type;
        typedef ConstVectorView<T,nonunitAr> const_subvector_step_type;
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

        typedef CVIt<T,S,_conj?Conj:NonConj> const_iterator;
        typedef CVIt<T,negS,_conj?Conj:NonConj> const_reverse_iterator;
    
        typedef typename AuxRef<T,_conj>::reference reference;

        typedef VectorView<T,Ar> subvector_type;
        typedef VectorView<T,nonunitAr> subvector_step_type;
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
        typedef SmallVectorView<T,N,S,An> noalias_type;
        typedef SmallVectorView<T,N,S,An|CheckAlias> alias_type;

        typedef VIt<T,S,_conj?Conj:NonConj> iterator;
        typedef VIt<T,negS,_conj?Conj:NonConj> reverse_iterator;
    };

    template <class T, ptrdiff_t N, ptrdiff_t S, int A>
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

        SmallVectorView(T* v, ptrdiff_t n, ptrdiff_t s) :
            itsv(v), itssize(n), itsstep(s) 
        {
            TMVStaticAssert(Traits<type>::okA);
        }

        SmallVectorView(T* v, ptrdiff_t n) :
            itsv(v), itssize(n), itsstep(S) 
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(S != Unknown); 
        }

        SmallVectorView(T* v) : itsv(v), itssize(N), itsstep(S) 
        { 
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(N != Unknown);  TMVStaticAssert(S != Unknown); 
        }

        SmallVectorView(const type& v2) : 
            itsv(v2.itsv), itssize(v2.size()), itsstep(v2.step())
        {
            TMVStaticAssert(Traits<type>::okA);
        }

        template <ptrdiff_t N2, ptrdiff_t S2, int A2>
        SmallVectorView(SmallVectorView<T,N2,S2,A2> v2) :
            itsv(v2.ptr()), itssize(v2.size()), itsstep(v2.step())
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Attrib<A>::conj == int(Attrib<A2>::conj)); 
        }

        template <int A2>
        SmallVectorView(VectorView<T,A2> v2) :
            itsv(v2.ptr()), itssize(v2.size()), itsstep(v2.step())
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Attrib<A>::conj == int(Attrib<A2>::conj)); 
        }

        ~SmallVectorView() 
        {
#ifdef TMV_EXTRA_DEBUG
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
        T cref(ptrdiff_t i) const  { return DoConj<_conj>(itsv[i*step()]); }
        TMV_INLINE T* ptr() { return itsv; }
        reference ref(ptrdiff_t i) { return reference(itsv[i*step()]); }

        TMV_INLINE ptrdiff_t size() const { return itssize; }
        TMV_INLINE ptrdiff_t nElements() const { return itssize; }
        TMV_INLINE ptrdiff_t step() const { return itsstep; }
        TMV_INLINE bool isconj() const { return _conj; }

    protected :

        T* itsv;
        const CheckedInt<N> itssize;
        const CheckedInt<S> itsstep;

    }; // SmallVectorView


    // 
    // Swap
    //

    template <class V, class T, ptrdiff_t N, ptrdiff_t S, int A>
    TMV_INLINE void Swap(
        BaseVector_Mutable<V>& v1, SmallVectorView<T,N,S,A> v2)
    { DoSwap(v1,v2); }
    template <class V, class T, ptrdiff_t N, ptrdiff_t S, int A>
    TMV_INLINE void Swap(
        SmallVectorView<T,N,S,A> v1, BaseVector_Mutable<V>& v2)
    { DoSwap(v1,v2); }
    template <class T, ptrdiff_t N, ptrdiff_t S1, int A1, ptrdiff_t S2, int A2>
    TMV_INLINE void Swap(
        SmallVectorView<T,N,S1,A1> v1, SmallVectorView<T,N,S2,A2> v2)
    { DoSwap(v1,v2); }
    template <class T, ptrdiff_t N, int A1, ptrdiff_t S2, int A2>
    TMV_INLINE void Swap(
        VectorView<T,A1> v1, SmallVectorView<T,N,S2,A2> v2)
    { DoSwap(v1,v2); }
    template <class T, ptrdiff_t N, ptrdiff_t S1, int A1, int A2>
    TMV_INLINE void Swap(
        SmallVectorView<T,N,S1,A1> v1, VectorView<T,A2> v2)
    { DoSwap(v1,v2); }


    //
    // Conjugate
    //
    
    template <class T, ptrdiff_t N, int A>
    TMV_INLINE typename SmallVector<T,N,A>::conjugate_type Conjugate(
        SmallVector<T,N,A>& v)
    { return v.conjugate(); }
    template <class T, ptrdiff_t N, ptrdiff_t S, int A>
    TMV_INLINE typename SmallVectorView<T,N,S,A>::conjugate_type Conjugate(
        SmallVectorView<T,N,S,A> v)
    { return v.conjugate(); }


    //
    // TMV_Text functions
    //

    template <class T, ptrdiff_t N, int A>
    inline std::string TMV_Text(const SmallVector<T,N,A>& v)
    {
        std::ostringstream s;
        s << "SmallVector<"<<TMV_Text(T());
        s << ","<<N<<","<<Attrib<A>::vtext()<<">";
        s << "("<<v.size()<<","<<v.step()<<")";
        return s.str();
    }

    template <class T, ptrdiff_t N, ptrdiff_t S, int A>
    inline std::string TMV_Text(const SmallVectorView<T,N,S,A>& v)
    {
        std::ostringstream s;
        s << "SmallVectorView<"<<TMV_Text(T());
        s << ","<<IntTraits<N>::text();
        s << ","<<IntTraits<S>::text();
        s << ","<<Attrib<A>::vtext()<<">";
        s << "("<<v.size()<<","<<v.step()<<")";
        return s.str();
    }

    template <class T, ptrdiff_t N, ptrdiff_t S, int A>
    std::string TMV_Text(const ConstSmallVectorView<T,N,S,A>& v)
    {
        std::ostringstream s;
        s << "ConstSmallVectorView<"<<TMV_Text(T());
        s << ","<<IntTraits<N>::text();
        s << ","<<IntTraits<S>::text();
        s << ","<<Attrib<A>::vtext()<<">";
        s << "("<<v.size()<<","<<v.step()<<")";
        return s.str();
    }

} // namespace tmv

#endif
