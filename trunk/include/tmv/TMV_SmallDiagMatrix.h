

//---------------------------------------------------------------------------
//
// This file defines the TMV SmallDiagMatrix class.
//
// The SmallDiagMatrix class acts just like a DiagMatrix, except that
// the size of the matrix is provided as a template parameter.
// This makes some operations faster.
//
// Constructors:
//
//    SmallDiagMatrix<T,N,A>()
//        Makes a DiagMatrix with column size and row size = size
//        with _uninitialized_ values
//
//    SmallDiagMatrix<T,N,A>(T x)
//        Makes a DiagMatrix of size n with all values = x
//
//    SmallDiagMatrix<T,N,A>(const BaseVector<V>& vv)
//        Make a DiagMatrix which copies the elements of vv.
//
// All the other operations with a DiagMatrix work the same for 
// SmallDiagMatrix.
//

#ifndef TMV_SmallDiagMatrix_H
#define TMV_SmallDiagMatrix_H

#include "TMV_BaseMatrix_Diag.h"
#include "TMV_VIt.h"
#include "TMV_Array.h"
#include "TMV_DiagMatrix.h"

namespace tmv {

    //
    // SmallDiagMatrix
    //

    template <class T, int N, int A0>
    struct Traits<SmallDiagMatrix<T,N,A0> >
    {
        enum { A = (A0 & ~NoDivider & ~NoAlias) | Unit };
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

        typedef SmallDiagMatrix<T,N,A0> type;
        typedef const type& calc_type;
        typedef const type& eval_type;
        typedef type copy_type;

        enum { _colsize = N };
        enum { _rowsize = N };
        enum { _size = N };
        enum { _nlo = 0 };
        enum { _nhi = 0 };
        enum { _shape = Diag };
        enum { _fort = Attrib<A>::fort };
        enum { _calc = true };
        enum { _rowmajor = false };
        enum { _colmajor = false };
        enum { _step = 1 };
        enum { _diagstep = 1 };
        enum { _conj = false };
        enum { _checkalias = Attrib<A>::checkalias };
        enum { _unit = true };
        enum { twoS = isreal ? 1 : 2 };

        enum { _hasdivider = false };
        typedef QuotXM<1,real_type,type> inverse_type;

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

        typedef ConstSmallVectorView<T,N,_diagstep,A> const_diag_type;
        typedef ConstDiagMatrixView<T,Ar> const_subdiagmatrix_type;
        typedef ConstDiagMatrixView<T,nonunitAr> const_subdiagmatrix_step_type;
        typedef ConstSmallDiagMatrixView<T,N,_diagstep,A> const_view_type;
        typedef ConstSmallDiagMatrixView<T,N,_diagstep,cstyleA>
            const_cview_type;
        typedef ConstSmallDiagMatrixView<T,N,_diagstep,fstyleA>
            const_fview_type;
        typedef ConstDiagMatrixView<T> const_xview_type;
        typedef ConstSmallDiagMatrixView<T,N,1,unitA> const_unitview_type;
        typedef ConstSmallDiagMatrixView<T,N,_diagstep,conjA> 
            const_conjugate_type;
        typedef ConstSmallDiagMatrixView<T,N,_diagstep,A> const_transpose_type;
        typedef ConstSmallDiagMatrixView<T,N,_diagstep,conjA> 
            const_adjoint_type;
        typedef ConstSmallDiagMatrixView<real_type,N,twoS,twosA>
            const_realpart_type;
        typedef const_realpart_type const_imagpart_type;
        typedef ConstSmallDiagMatrixView<T,N,_diagstep,A> const_nonconj_type;
        typedef SmallDiagMatrixView<T,N,_diagstep,A> nonconst_type;

        typedef typename const_diag_type::const_iterator const_iterator;
        typedef const_iterator const_rowmajor_iterator;
        typedef const_iterator const_colmajor_iterator;

        typedef T& reference;

        typedef SmallVectorView<T,N,_diagstep,A> diag_type;
        typedef DiagMatrixView<T,Ar> subdiagmatrix_type;
        typedef DiagMatrixView<T,nonunitAr> subdiagmatrix_step_type;
        typedef SmallDiagMatrixView<T,N,_diagstep,A> view_type;
        typedef SmallDiagMatrixView<T,N,_diagstep,cstyleA> cview_type;
        typedef SmallDiagMatrixView<T,N,_diagstep,fstyleA> fview_type;
        typedef DiagMatrixView<T> xview_type;
        typedef SmallDiagMatrixView<T,N,1,unitA> unitview_type;
        typedef SmallDiagMatrixView<T,N,_diagstep,conjA> conjugate_type;
        typedef SmallDiagMatrixView<T,N,_diagstep,A> transpose_type;
        typedef SmallDiagMatrixView<T,N,_diagstep,conjA> adjoint_type;
        typedef SmallDiagMatrixView<real_type,N,twoS,twosA> realpart_type;
        typedef realpart_type imagpart_type;
        typedef SmallDiagMatrixView<T,N,_diagstep,A> nonconj_type;
        typedef SmallDiagMatrixView<T,N,_diagstep,An> noalias_type;
        typedef SmallDiagMatrixView<T,N,_diagstep,An|CheckAlias> alias_type;

        typedef typename diag_type::iterator iterator;
        typedef iterator rowmajor_iterator;
        typedef iterator colmajor_iterator;
    };

    template <class T, int N, int A>
    class SmallDiagMatrix : 
        public BaseMatrix_Diag_Mutable<SmallDiagMatrix<T,N,A> >
    {
    public:

        typedef SmallDiagMatrix<T,N,A> type;
        typedef BaseMatrix_Diag_Mutable<type> base_mut;

        typedef typename Traits<T>::real_type real_type;

        enum { _colsize = Traits<type>::_size };
        enum { _rowsize = Traits<type>::_size };
        enum { _size = Traits<type>::_size };
        enum { _fort = Traits<type>::_fort };
        enum { _shape = Traits<type>::_shape };
        enum { _rowmajor = Traits<type>::_rowmajor };
        enum { _colmajor = Traits<type>::_colmajor };
        enum { _calc = Traits<type>::_calc };
        enum { _step = Traits<type>::_step };
        enum { _diagstep = Traits<type>::_diagstep };
        enum { _conj = Traits<type>::_conj };
        enum { _checkalias = Traits<type>::_checkalias };
        enum { _unit = Traits<type>::_unit };
        enum { _attrib = Traits<type>::_attrib };

        //
        // Constructors
        //

        TMV_INLINE_ND SmallDiagMatrix()
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(N >= 0);
#ifdef TMV_DEBUG
            this->diag().setAllTo(Traits<T>::constr_value());
#endif
        }

        explicit SmallDiagMatrix(T x) 
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(N >= 0);
            this->setAllTo(x);
        }

        SmallDiagMatrix(const type& m2) 
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(N >= 0);
            this->noAlias() = m2;
        }

        template <class V2>
        SmallDiagMatrix(const BaseVector<V2>& v2) 
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(N >= 0);
            typename type::diag_type d = this->diag();
            this->diag().noAlias() = v2;
        }

        template <class M2>
        SmallDiagMatrix(const BaseMatrix<M2>& m2) 
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(N >= 0);
            TMVStaticAssert((Sizes<M2::_rowsize,M2::_colsize>::same));
            TMVAssert(m2.colsize() == m2.rowsize());
            DiagCopy<ShapeTraits2<M2::_shape,Diag>::assignable>::copy(
                m2,*this);
        }

        TMV_INLINE_ND ~SmallDiagMatrix() 
        {
#ifdef TMV_DEBUG
            this->diag().setAllTo(Traits<T>::destr_value());
#endif
        }

        //
        // Op=
        //

        TMV_INLINE type& operator=(const type& m2)
        { if (this != &m2) base_mut::operator=(m2); return *this; }

        template <class M2>
        TMV_INLINE type& operator=(const BaseMatrix<M2>& m2)
        { base_mut::operator=(m2); return *this; }

        TMV_INLINE type& operator=(T x)
        { base_mut::operator=(x); return *this; }


        //
        // Auxilliary Functions
        //

        TMV_INLINE const T* cptr() const { return itsm; }
        TMV_INLINE T* ptr() { return itsm; }

        T cref(int i, int j) const { return (i!=j ? T(0) : cref(i)); }
        T cref(int i) const { return itsm[i]; }
        T& ref(int i) { return itsm[i]; }

        TMV_INLINE size_t size() const { return N; }
        TMV_INLINE int nElements() const { return N; }
        TMV_INLINE int step() const { return 1; }
        TMV_INLINE bool isconj() const { return false; }
        TMV_INLINE bool isrm() const { return true; }
        TMV_INLINE bool iscm() const { return true; }

    private:

        StackArray<T,N> itsm;

    }; // SmallDiagMatrix

    template <class T, int N, int S, int A0>
    struct Traits<ConstSmallDiagMatrixView<T,N,S,A0> >
    {
        enum { A = (A0 & ~NoDivider & ~NoAlias) | (
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

        typedef ConstSmallDiagMatrixView<T,N,S,A0> type;
        typedef const type& calc_type;
        typedef const type& eval_type;

        enum { _colsize = N };
        enum { _rowsize = N };
        enum { _size = N };
        enum { _nlo = 0 };
        enum { _nhi = 0 };
        enum { _shape = Diag };
        enum { _fort = Attrib<A>::fort };
        enum { _calc = true };
        enum { _rowmajor = false };
        enum { _colmajor = false };
        enum { _step = S };
        enum { _diagstep = S };
        enum { _conj = Attrib<A>::conj };
        enum { _checkalias = Attrib<A>::checkalias };
        enum { _unit = Attrib<A>::unit };
        enum { twoS = isreal ? 1 : 2 };

        enum { known = (N != TMV_UNKNOWN) };
        enum { copyA = _fort ? FortranStyle : CStyle };
        typedef typename TypeSelect<known,
                SmallDiagMatrix<T,N,copyA|NoAlias>,
                DiagMatrix<T,copyA> >::type copy_type;

        enum { _hasdivider = false };
        typedef QuotXM<1,real_type,type> inverse_type;

        enum { unitA = A | Unit };
        enum { nonunitA = A & ~Unit };
        enum { conjA = iscomplex ? (A ^ Conj) : int(A) };
        enum { nonconjA = A & ~Conj };
        enum { cstyleA = A & ~FortranStyle };
        enum { fstyleA = A | FortranStyle };
        enum { twosA = isreal ? int(A) : (A & ~Conj & ~Unit) };
        enum { Ar = _checkalias ? (A & ~CheckAlias) : (A | NoAlias) };
        enum { nonunitAr = Ar & ~Unit };

        typedef ConstSmallVectorView<T,N,_diagstep,A> const_diag_type;
        typedef ConstDiagMatrixView<T,Ar> const_subdiagmatrix_type;
        typedef ConstDiagMatrixView<T,nonunitAr> const_subdiagmatrix_step_type;
        typedef ConstSmallDiagMatrixView<T,N,_diagstep,A> const_view_type;
        typedef ConstSmallDiagMatrixView<T,N,_diagstep,cstyleA>
            const_cview_type;
        typedef ConstSmallDiagMatrixView<T,N,_diagstep,fstyleA>
            const_fview_type;
        typedef ConstDiagMatrixView<T,_conj ? Conj : NonConj> const_xview_type;
        typedef ConstSmallDiagMatrixView<T,N,1,unitA> const_unitview_type;
        typedef ConstSmallDiagMatrixView<T,N,_diagstep,conjA> 
            const_conjugate_type;
        typedef ConstSmallDiagMatrixView<T,N,_diagstep,A> const_transpose_type;
        typedef ConstSmallDiagMatrixView<T,N,_diagstep,conjA> 
            const_adjoint_type;
        typedef ConstSmallDiagMatrixView<real_type,N,twoS,twosA>
            const_realpart_type;
        typedef const_realpart_type const_imagpart_type;
        typedef ConstSmallDiagMatrixView<T,N,_diagstep,A> const_nonconj_type;
        typedef SmallDiagMatrixView<T,N,_diagstep,A> nonconst_type;

        typedef typename const_diag_type::const_iterator const_iterator;
        typedef const_iterator const_rowmajor_iterator;
        typedef const_iterator const_colmajor_iterator;
    };

    template <class T, int N, int S, int A>
    class ConstSmallDiagMatrixView :
        public BaseMatrix_Diag<ConstSmallDiagMatrixView<T,N,S,A> >
    {
    public:

        typedef ConstSmallDiagMatrixView<T,N,S,A> type;

        enum { _colsize = Traits<type>::_size };
        enum { _rowsize = Traits<type>::_size };
        enum { _size = Traits<type>::_size };
        enum { _fort = Traits<type>::_fort };
        enum { _shape = Traits<type>::_shape };
        enum { _rowmajor = Traits<type>::_rowmajor };
        enum { _colmajor = Traits<type>::_colmajor };
        enum { _calc = Traits<type>::_calc };
        enum { _step = Traits<type>::_step };
        enum { _diagstep = Traits<type>::_diagstep };
        enum { _conj = Traits<type>::_conj };
        enum { _checkalias = Traits<type>::_checkalias };
        enum { _unit = Traits<type>::_unit };
        enum { _attrib = Traits<type>::_attrib };

        //
        // Constructors
        //

        ConstSmallDiagMatrixView(const T* m, size_t n, int s) :
            itsm(m), itssize(n), itsstep(s)
        {
            TMVStaticAssert(Traits<type>::okA);
        }

        ConstSmallDiagMatrixView(const T* m, size_t n) :
            itsm(m), itssize(n), itsstep(S)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(S != TMV_UNKNOWN); 
        }

        ConstSmallDiagMatrixView(const T* m) :
            itsm(m), itssize(N), itsstep(S)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(N != TMV_UNKNOWN); 
            TMVStaticAssert(S != TMV_UNKNOWN); 
        }

        ConstSmallDiagMatrixView(const type& m2) :
            itsm(m2.cptr()), itssize(m2.size()), itsstep(m2.step()) 
        {
            TMVStaticAssert(Traits<type>::okA);
        }

        template <int A2>
        ConstSmallDiagMatrixView(const ConstDiagMatrixView<T,A2>& m2) :
            itsm(m2.cptr()), itssize(m2.size()), itsstep(m2.step()) 
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Attrib<A>::conj == int(Attrib<A2>::conj)); 
        }

        template <int A2>
        ConstSmallDiagMatrixView(const DiagMatrixView<T,A2>& m2) :
            itsm(m2.cptr()), itssize(m2.size()), itsstep(m2.step()) 
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Attrib<A>::conj == int(Attrib<A2>::conj)); 
        }

        template <int N2, int S2, int A2>
        ConstSmallDiagMatrixView(
            const ConstSmallDiagMatrixView<T,N2,S2,A2>& m2) :
            itsm(m2.cptr()), itssize(m2.size()), itsstep(m2.step()) 
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Attrib<A>::conj == int(Attrib<A2>::conj)); 
        }

        template <int N2, int S2, int A2>
        ConstSmallDiagMatrixView(
            const SmallDiagMatrixView<T,N2,S2,A2>& m2) :
            itsm(m2.cptr()), itssize(m2.size()), itsstep(m2.step()) 
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Attrib<A>::conj == int(Attrib<A2>::conj)); 
        }

        ~ConstSmallDiagMatrixView() {
#ifdef TMV_DEBUG
            itsm = 0; 
#endif
        }

    private :
        void operator=(const type& m2);
    public :

        //
        // Auxilliary Functions
        //

        TMV_INLINE const T* cptr() const { return itsm; }

        T cref(int i, int j) const { return (i!=j ? T(0) : cref(i)); }
        T cref(int i) const
        { return DoConj<_conj>(itsm[i*step()]); }

        TMV_INLINE size_t size() const { return itssize; }
        TMV_INLINE int nElements() const { return itssize; }
        TMV_INLINE int step() const { return itsstep; }
        TMV_INLINE bool isconj() const { return _conj; }
        TMV_INLINE bool isrm() const { return step()==1; }
        TMV_INLINE bool iscm() const { return step()==1; }

    private :

        const T* itsm;
        const CheckedInt<N> itssize;
        const CheckedInt<S> itsstep;

    }; // ConstSmallDiagMatrixView

    template <class T, int N, int S, int A0>
    struct Traits<SmallDiagMatrixView<T,N,S,A0> >
    {
        enum { A = (A0 & ~NoDivider & ~NoAlias) | (
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

        typedef SmallDiagMatrixView<T,N,S,A0> type;
        typedef const type& calc_type;
        typedef const type& eval_type;

        enum { _colsize = N };
        enum { _rowsize = N };
        enum { _size = N };
        enum { _nlo = 0 };
        enum { _nhi = 0 };
        enum { _shape = Diag };
        enum { _fort = Attrib<A>::fort };
        enum { _calc = true };
        enum { _rowmajor = false };
        enum { _colmajor = false };
        enum { _step = S };
        enum { _diagstep = S };
        enum { _conj = Attrib<A>::conj };
        enum { _checkalias = Attrib<A>::checkalias };
        enum { _unit = Attrib<A>::unit };
        enum { twoS = isreal ? 1 : 2 };

        enum { known = (N != TMV_UNKNOWN) };
        enum { copyA = _fort ? FortranStyle : CStyle };
        typedef typename TypeSelect<known,
                SmallDiagMatrix<T,N,copyA|NoAlias>,
                DiagMatrix<T,copyA> >::type copy_type;

        enum { _hasdivider = false };
        typedef QuotXM<1,real_type,type> inverse_type;

        enum { unitA = A | Unit };
        enum { nonunitA = A & ~Unit };
        enum { conjA = iscomplex ? (A ^ Conj) : int(A) };
        enum { nonconjA = A & ~Conj };
        enum { cstyleA = A & ~FortranStyle };
        enum { fstyleA = A | FortranStyle };
        enum { twosA = isreal ? int(A) : (A & ~Conj & ~Unit) };
        enum { Ar = _checkalias ? (A & ~CheckAlias) : (A | NoAlias) };
        enum { nonunitAr = Ar & ~Unit };
        enum { An = (A & ~CheckAlias) };

        typedef ConstSmallVectorView<T,N,_diagstep,A> const_diag_type;
        typedef ConstDiagMatrixView<T,Ar> const_subdiagmatrix_type;
        typedef ConstDiagMatrixView<T,nonunitAr> const_subdiagmatrix_step_type;
        typedef ConstSmallDiagMatrixView<T,N,_diagstep,A> const_view_type;
        typedef ConstSmallDiagMatrixView<T,N,_diagstep,cstyleA>
            const_cview_type;
        typedef ConstSmallDiagMatrixView<T,N,_diagstep,fstyleA>
            const_fview_type;
        typedef ConstDiagMatrixView<T,_conj ? Conj : NonConj> const_xview_type;
        typedef ConstSmallDiagMatrixView<T,N,1,unitA> const_unitview_type;
        typedef ConstSmallDiagMatrixView<T,N,_diagstep,conjA> 
            const_conjugate_type;
        typedef ConstSmallDiagMatrixView<T,N,_diagstep,A> const_transpose_type;
        typedef ConstSmallDiagMatrixView<T,N,_diagstep,conjA> 
            const_adjoint_type;
        typedef ConstSmallDiagMatrixView<real_type,N,twoS,twosA>
            const_realpart_type;
        typedef const_realpart_type const_imagpart_type;
        typedef ConstSmallDiagMatrixView<T,N,_diagstep,A> const_nonconj_type;
        typedef SmallDiagMatrixView<T,N,_diagstep,A> nonconst_type;

        typedef typename const_diag_type::const_iterator const_iterator;
        typedef const_iterator const_rowmajor_iterator;
        typedef const_iterator const_colmajor_iterator;

        typedef typename AuxRef<T,_conj>::reference reference;

        typedef SmallVectorView<T,N,_diagstep,A> diag_type;
        typedef DiagMatrixView<T,Ar> subdiagmatrix_type;
        typedef DiagMatrixView<T,nonunitAr> subdiagmatrix_step_type;
        typedef SmallDiagMatrixView<T,N,_diagstep,A> view_type;
        typedef SmallDiagMatrixView<T,N,_diagstep,cstyleA> cview_type;
        typedef SmallDiagMatrixView<T,N,_diagstep,fstyleA> fview_type;
        typedef DiagMatrixView<T,_conj ? Conj : NonConj> xview_type;
        typedef SmallDiagMatrixView<T,N,1,unitA> unitview_type;
        typedef SmallDiagMatrixView<T,N,_diagstep,conjA> conjugate_type;
        typedef SmallDiagMatrixView<T,N,_diagstep,A> transpose_type;
        typedef SmallDiagMatrixView<T,N,_diagstep,conjA> adjoint_type;
        typedef SmallDiagMatrixView<real_type,N,twoS,twosA> realpart_type;
        typedef realpart_type imagpart_type;
        typedef SmallDiagMatrixView<T,N,_diagstep,A> nonconj_type;
        typedef SmallDiagMatrixView<T,N,_diagstep,An> noalias_type;
        typedef SmallDiagMatrixView<T,N,_diagstep,An|CheckAlias> alias_type;

        typedef typename diag_type::iterator iterator;
        typedef iterator rowmajor_iterator;
        typedef iterator colmajor_iterator;
    };

    template <class T, int N, int S, int A>
    class SmallDiagMatrixView :
        public BaseMatrix_Diag_Mutable<SmallDiagMatrixView<T,N,S,A> >
    {
    public:

        typedef SmallDiagMatrixView<T,N,S,A> type;
        typedef BaseMatrix_Diag_Mutable<type> base_mut;
        typedef typename base_mut::reference reference;

        enum { _colsize = Traits<type>::_size };
        enum { _rowsize = Traits<type>::_size };
        enum { _size = Traits<type>::_size };
        enum { _fort = Traits<type>::_fort };
        enum { _shape = Traits<type>::_shape };
        enum { _rowmajor = Traits<type>::_rowmajor };
        enum { _colmajor = Traits<type>::_colmajor };
        enum { _calc = Traits<type>::_calc };
        enum { _step = Traits<type>::_step };
        enum { _diagstep = Traits<type>::_diagstep };
        enum { _conj = Traits<type>::_conj };
        enum { _checkalias = Traits<type>::_checkalias };
        enum { _unit = Traits<type>::_unit };
        enum { _attrib = Traits<type>::_attrib };

        //
        // Constructors
        //

        SmallDiagMatrixView(T* m, size_t n, int s) :
            itsm(m), itssize(n), itsstep(s) 
        {
            TMVStaticAssert(Traits<type>::okA);
        }

        SmallDiagMatrixView(T* m, size_t n) :
            itsm(m), itssize(n), itsstep(S)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(S != TMV_UNKNOWN); 
        }

        SmallDiagMatrixView(T* m) : itsm(m), itssize(N), itsstep(S)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(N != TMV_UNKNOWN); TMVStaticAssert(S != TMV_UNKNOWN); 
        }

        SmallDiagMatrixView(const type& m2) :
            itsm(m2.itsm), itssize(m2.size()), itsstep(m2.step()) 
        {
            TMVStaticAssert(Traits<type>::okA);
        }

        template <int A2>
        SmallDiagMatrixView(DiagMatrixView<T,A2> m2) :
            itsm(m2.ptr()), itssize(m2.size()), itsstep(m2.step()) 
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Attrib<A>::conj == int(Attrib<A2>::conj)); 
        }

        template <int N2, int S2, int A2>
        SmallDiagMatrixView(SmallDiagMatrixView<T,N2,S2,A2> m2) :
            itsm(m2.ptr()), itssize(m2.size()), itsstep(m2.step()) 
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Attrib<A>::conj == int(Attrib<A2>::conj)); 
        }

        ~SmallDiagMatrixView() {
#ifdef TMV_DEBUG
            itsm = 0; 
#endif
        }


        //
        // Op = 
        //

        TMV_INLINE type& operator=(const type& m2)
        { if (this != &m2) base_mut::operator=(m2); return *this; }

        template <class M2>
        TMV_INLINE type& operator=(const BaseMatrix<M2>& m2)
        { base_mut::operator=(m2); return *this; }

        TMV_INLINE type& operator=(const T x)
        { base_mut::operator=(x); return *this; }


        //
        // Auxilliary Functions
        //

        TMV_INLINE const T* cptr() const { return itsm; }
        TMV_INLINE T* ptr() { return itsm; }

        T cref(int i, int j) const { return (i!=j ? T(0) : cref(i)); }
        T cref(int i) const
        { return DoConj<_conj>(itsm[i*step()]); }

        reference ref(int i) 
        { return reference(itsm[i*step()]); }

        TMV_INLINE size_t size() const { return itssize; }
        TMV_INLINE int nElements() const { return itssize; }
        TMV_INLINE int step() const { return itsstep; }
        TMV_INLINE bool isconj() const { return _conj; }
        TMV_INLINE bool isrm() const { return step()==1; }
        TMV_INLINE bool iscm() const { return step()==1; }

    private :

        T* itsm;
        const size_t itssize;
        const CheckedInt<S> itsstep;

    }; // SmallDiagMatrixView


    //
    // DiagMatrixViewOf
    //
    
    template <class V>
    struct DMVO
    {
        typedef typename V::value_type T;
        enum { N = V::_size };
        enum { S = V::_step };
        enum { A = (
                ( V::_conj ? Conj : NonConj ) |
                ( V::_fort ? FortranStyle : CStyle ) |
                ( V::_checkalias ? CheckAlias : NoAlias ) |
                ( S == 1 ? Unit : NonUnit ) )};

        typedef typename TypeSelect<
            ( N==TMV_UNKNOWN && (S==TMV_UNKNOWN || S==1) ),
            ConstDiagMatrixView<T,A> ,
            ConstSmallDiagMatrixView<T,N,S,A> >::type cv;
        typedef typename TypeSelect<
            ( N==TMV_UNKNOWN && (S==TMV_UNKNOWN || S==1) ),
            DiagMatrixView<T,A> ,
            SmallDiagMatrixView<T,N,S,A> >::type v;
    };

    template <class V>
    TMV_INLINE typename DMVO<V>::cv DiagMatrixViewOf(
        const BaseVector_Calc<V>& v)
    { return typename DMVO<V>::cv(v.cptr(),v.size(),v.step()); }

    template <class V>
    TMV_INLINE typename DMVO<V>::v DiagMatrixViewOf(
        BaseVector_Mutable<V>& v)
    { return typename DMVO<V>::v(v.ptr(),v.size(),v.step()); }

    template <class T, int A>
    TMV_INLINE typename DMVO<VectorView<T,A> >::v DiagMatrixViewOf(
        VectorView<T,A> v)
    { return typename DMVO<VectorView<T,A> >::v(v.ptr(),v.size(),v.step()); }

    template <class T, int N, int S, int A>
    TMV_INLINE typename DMVO<SmallVectorView<T,N,S,A> >::v DiagMatrixViewOf(
        SmallVectorView<T,N,S,A> v)
    {
        return typename DMVO<SmallVectorView<T,N,S,A> >::v(
            v.ptr(),v.size(),v.step()); 
    }

    
    //
    // Swap
    //

    template <class T, int N, int S, int A, class MM>
    TMV_INLINE void Swap(
        BaseMatrix_Diag_Mutable<MM>& m1, SmallDiagMatrixView<T,N,S,A> m2)
    { Swap(m1.diag(),m2.diag()); }
    template <class T, int N, int S, int A, class MM>
    TMV_INLINE void Swap(
        SmallDiagMatrixView<T,N,S,A> m1, BaseMatrix_Diag_Mutable<MM>& m2)
    { Swap(m1.diag(),m2.diag()); }
    template <class T, int N, int S1, int A1, int S2, int A2>
    TMV_INLINE void Swap(
        SmallDiagMatrixView<T,N,S1,A1> m1, SmallDiagMatrixView<T,N,S2,A2> m2)
    { Swap(m1.diag(),m2.diag()); }
    template <class T, int N, int S1, int A1, int A2>
    TMV_INLINE void Swap(
        SmallDiagMatrixView<T,N,S1,A1> m1, DiagMatrixView<T,A2> m2)
    { Swap(m1.diag(),m2.diag()); }
    template <class T, int N, int A1, int S2, int A2>
    TMV_INLINE void Swap(
        DiagMatrixView<T,A1> m1, SmallDiagMatrixView<T,N,S2,A2> m2)
    { Swap(m1.diag(),m2.diag()); }


    //
    // Conjugate, Transpose, Adjoint
    //
    
    template <class T, int N, int A>
    TMV_INLINE typename SmallDiagMatrix<T,N,A>::conjugate_type Conjugate(
        SmallDiagMatrix<T,N,A>& m)
    { return m.conjugate(); }
    template <class T, int N, int S, int A>
    TMV_INLINE typename SmallDiagMatrixView<T,N,S,A>::conjugate_type Conjugate(
        SmallDiagMatrixView<T,N,S,A> m)
    { return m.conjugate(); }

    template <class T, int N, int A>
    TMV_INLINE typename SmallDiagMatrix<T,N,A>::transpose_type Transpose(
        SmallDiagMatrix<T,N,A>& m)
    { return m.transpose(); }
    template <class T, int N, int S, int A>
    TMV_INLINE typename SmallDiagMatrixView<T,N,S,A>::transpose_type Transpose(
        SmallDiagMatrixView<T,N,S,A> m)
    { return m.transpose(); }

    template <class T, int N, int A>
    TMV_INLINE typename SmallDiagMatrix<T,N,A>::adjoint_type Adjoint(
        SmallDiagMatrix<T,N,A>& m)
    { return m.adjoint(); }
    template <class T, int N, int S, int A>
    TMV_INLINE typename SmallDiagMatrixView<T,N,S,A>::adjoint_type Adjoint(
        SmallDiagMatrixView<T,N,S,A> m)
    { return m.adjoint(); }


    //
    // TMV_Text 
    //

#ifdef TMV_TEXT
    template <class T, int N, int A>
    inline std::string TMV_Text(const SmallDiagMatrix<T,N,A>& m)
    {
        std::ostringstream s;
        s << "SmallDiagMatrix<"<<TMV_Text(T());
        s << ","<<N<<','<<Attrib<A>::vtext<<">";
        s << "("<<m.size()<<","<<m.step()<<")";
        return s.str();
    }

    template <class T, int N, int S, int A>
    inline std::string TMV_Text(
        const ConstSmallDiagMatrixView<T,N,S,A>& m)
    {
        std::ostringstream s;
        s << "ConstSmallDiagMatrixView<"<<TMV_Text(T());
        s << ","<<IntTraits<N>::text();
        s << ","<<IntTraits<S>::text();
        s << ","<<Attrib<A>::vtext<<">";
        s << "("<<m.size()<<","<<m.step()<<")";
        return s.str();
    }

    template <class T, int N, int S, int A>
    inline std::string TMV_Text(const SmallDiagMatrixView<T,N,S,A>& m)
    {
        std::ostringstream s;
        s << "SmallDiagMatrixView<"<<TMV_Text(T());
        s << ","<<IntTraits<N>::text();
        s << ","<<IntTraits<S>::text();
        s << ","<<Attrib<A>::vtext<<">";
        s << "("<<m.size()<<","<<m.step()<<")";
        return s.str();
    }
#endif

} // namespace tmv

#endif
