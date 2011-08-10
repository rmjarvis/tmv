

//---------------------------------------------------------------------------
//
// This file defines the TMV DiagMatrix class.
//
// The DiagMatrix class is provided for efficient storage of a diagonal
// matrix.  You can do most of the things that you can do with a 
// regular Matrix, but it will do them more efficiently.
//
// Constructors:
//
//    DiagMatrix<T>(size_t size)
//        Makes a DiagMatrix with column size and row size = size
//        with _uninitialized_ values
//
//    DiagMatrix<T>(size_t size, T x)
//        Makes a DiagMatrix of size n with all values = x
//
//    DiagMatrix<T>(size_t size, T* vv)
//    DiagMatrix<T>(size_t size, const std::vector<T>& vv)
//        Makes a DiagMatrix of size n which copies the values is vv
//
//    DiagMatrix<T>(const Vector<T>& vv)
//        Make a DiagMatrix which copies the elements of vv.
//
//    ConstDiagMatrixView<T>(const Vector<T>& v)
//        Make a constant DiagMatrix view with v as the diagonal.
//        While this view cannon be modified, changing the original v or m
//        will cause corresponding changes in this view.
//
//    DiagMatrixView<T>(Vector<T>& v)
//        Make a mutable DiagMatrix view with v as the diagonal.
//
//
// Access Functions
//
//    size_t colsize() const
//    size_t rowsize() const
//    size_t size() const
//        Return the dimensions of the DiagMatrix
//
//    T& operator()(int i)
//    T operator()(int i) const
//    T& operator()(int i, int j)
//    T operator()(int i, int j) const
//        Return the (i,j) element of the DiagMatrix
//        For the single paramter version, j=i
//
//    VectorView& diag()
//    ConstVectorView& diag() const
//        Return the diagonal of the DiagMatrix as a VectorView
//
//
// Modifying Functions - The same as the regular Matrix counterparts
//
//    DiagMatrix& setZero()
//    DiagMatrix& setAllTo(T x)
//    DiagMatrix& addToAll(T x)
//    DiagMatrix<T>& transposeSelf() 
//        (Does nothing.)
//    DiagMatrix& conjugateSelf()
//    DiagMatrix& setToIdentity(x = 1)
//    void Swap(DiagMatrix& m1, DiagMatrix& m2)
//
//
// SubDiagMatrix:
//
//    subDiagMatrix(int i1, int i2, int istep=1)
//        Returns a Sub-DiagMatrix which extends from i1 to i2 (step istep)
//        which refers to the same physical elements as the original.
//        As usual, i2 is the "one past the end" element.
//
//
// Functions of DiagMatrices - Same as for regular Matrices:
//
//    Det(m)
//    LogDet(m)
//    Trace(m)
//    Norm(m) or NormF(m)
//    NormSq(m)
//    Norm1(m) 
//    Norm2(m) 
//    NormInf(m) 
//    MaxAbsElement(m) 
//        (Note - for diagonal matrices, 
//        Norm1 = Norm2 = NormInf = MaxAbsElement.)
//    MaxAbs2Element(m) 
//    SumElements(m) 
//    SumAbsElements(m) 
//    SumAbs2Elements(m) 
//    Transpose(m)
//        (same as the original).
//    Conjugate(m)
//    Adjoint(m)
//
//    m.inverse()
//    Inverse(m)
//    m.invertSelf()
//    m.makeInverse(minv) (takes either Matrix or DiagMatrix argument)
//    m.makeInverseATA(invata) (takes either Matrix or DiagMatrix argument)
//
// I/O: 
//
//    os << d 
//        Writes d to ostream os as a full matrix
//
//    d.writeCompact(os)
//        Writes only the diagonal Vector to os
//
//    is >> d
//        Reads in d in the compact format
//
//

#ifndef TMV_DiagMatrix_H
#define TMV_DiagMatrix_H

#include <vector>
#include "TMV_BaseMatrix_Diag.h"
#include "TMV_Array.h"

namespace tmv {

    template <class T, int A0>
    struct Traits<DiagMatrix<T,A0> >
    {
        enum { A = (A0 & ~NoDivider) | Unit };
        enum { okA = (
                Attrib<A>::vectoronly && 
                !Attrib<A>::conj )};
        enum { _attrib = A };

        typedef T value_type;
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;
        enum { isreal = Traits<T>::isreal };
        enum { iscomplex = Traits<T>::iscomplex };

        typedef DiagMatrix<T,A0> type;
        typedef const type& calc_type;
        typedef const type& eval_type;
        typedef type copy_type;

        enum { _colsize = TMV_UNKNOWN };
        enum { _rowsize = TMV_UNKNOWN };
        enum { _size = TMV_UNKNOWN };
        enum { _fort = Attrib<A>::fort };
        enum { _shape = Diag };
        enum { _rowmajor = false };
        enum { _colmajor = false };
        enum { _calc = true };
        enum { _step = 1 };
        enum { _diagstep = 1 };
        enum { _conj = false };
        enum { _checkalias = !Attrib<A>::noalias };
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
        enum { Ac = _checkalias ? (A & CheckAlias) : (A & ~NoAlias) };
        enum { nonunitAc = Ac & ~Unit };
        enum { twosAc = isreal ? int(Ac) : (Ac & ~Conj & ~Unit) };

        typedef ConstVectorView<T,A> const_diag_type;
        typedef ConstDiagMatrixView<T,A> const_subdiagmatrix_type;
        typedef ConstDiagMatrixView<T,nonunitA> const_subdiagmatrix_step_type;
        typedef ConstDiagMatrixView<T,A> const_view_type;
        typedef ConstDiagMatrixView<T,cstyleA> const_cview_type;
        typedef ConstDiagMatrixView<T,fstyleA> const_fview_type;
        typedef ConstDiagMatrixView<T> const_xview_type;
        typedef ConstDiagMatrixView<T,unitA> const_unitview_type;
        typedef ConstDiagMatrixView<T,conjA> const_conjugate_type;
        typedef ConstDiagMatrixView<T,A> const_transpose_type;
        typedef ConstDiagMatrixView<T,conjA> const_adjoint_type;
        typedef typename TypeSelect< iscomplex ,
                ConstSmallDiagMatrixView<real_type,TMV_UNKNOWN,twoS,twosAc> ,
                ConstDiagMatrixView<real_type,twosA> >::type const_realpart_type;
        typedef const_realpart_type const_imagpart_type;
        typedef ConstDiagMatrixView<T,nonconjA> const_nonconj_type;
        typedef DiagMatrixView<T,A> nonconst_type;

        typedef T& reference;

        typedef VectorView<T,A> diag_type;
        typedef DiagMatrixView<T,A> subdiagmatrix_type;
        typedef DiagMatrixView<T,nonunitA> subdiagmatrix_step_type;
        typedef DiagMatrixView<T,A> view_type;
        typedef DiagMatrixView<T,cstyleA> cview_type;
        typedef DiagMatrixView<T,fstyleA> fview_type;
        typedef DiagMatrixView<T> xview_type;
        typedef DiagMatrixView<T,unitA> unitview_type;
        typedef DiagMatrixView<T,conjA> conjugate_type;
        typedef DiagMatrixView<T,A> transpose_type;
        typedef DiagMatrixView<T,conjA> adjoint_type;
        typedef typename TypeSelect< iscomplex ,
                SmallDiagMatrixView<real_type,TMV_UNKNOWN,twoS,twosAc> ,
                DiagMatrixView<real_type,twosA> >::type realpart_type;
        typedef realpart_type imagpart_type;
        typedef DiagMatrixView<T,nonconjA> nonconj_type;
    };

    // A Helper class to make a DiagMatrix from a BaseMatrix
    // If the BaseMatrix is assignable to Diag then we do the 
    // assign.  But if not, then we allow the construction - we
    // just make the DiagMatrix from the diagonal of the matrix.
    template <bool assignable_to_diag>
    struct DiagCopy // true
    {
        template <class M1, class M2>
        static TMV_INLINE void copy(
            const BaseMatrix<M1>& m1, BaseMatrix_Diag_Mutable<M2>& m2)
        { m1.newAssignTo(m2); }
    };
    template <>
    struct DiagCopy<false>
    {
        template <class M1, class M2>
        static TMV_INLINE void copy(
            const BaseMatrix<M1>& m1, BaseMatrix_Diag_Mutable<M2>& m2)
        {
            typename M2::diag_type m2d = m2.diag();
            m1.calc().diag().newAssignTo(m2d);
        }
    };

    template <class T, int A>
    class DiagMatrix : 
        public BaseMatrix_Diag_Mutable<DiagMatrix<T,A> >
    {
    public:
        typedef DiagMatrix<T,A> type;
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

#define NEW_SIZE(cs,rs) \

        DiagMatrix(size_t n=0) : itssize(n), itsm(n)
        {
            TMVStaticAssert(Traits<type>::okA);
#ifdef TMV_DEBUG
            this->diag().flatten().setAllTo(Traits<real_type>::constr_value());
#endif
        }

        DiagMatrix(size_t n, T x) : itssize(n), itsm(n)
        {
            TMVStaticAssert(Traits<type>::okA);
            this->setAllTo(x);
        }

        DiagMatrix(size_t n, const T* vv) : itssize(n), itsm(n)
        {
            TMVStaticAssert(Traits<type>::okA);
            ConstDiagMatrixView<T,Unit>(vv,n).newAssignTo(*this);
        }

        DiagMatrix(const std::vector<T>& vv) : 
            itssize(vv.size()), itsm(itssize)
        {
            TMVStaticAssert(Traits<type>::okA);
            ConstDiagMatrixView<T,Unit>(&vv[0],itssize).newAssignTo(*this);
        }

        DiagMatrix(const type& m2) : itssize(m2.itssize), itsm(itssize)
        {
            TMVStaticAssert(Traits<type>::okA);
            m2.newAssignTo(*this);
        }

        template <class V2>
        DiagMatrix(const BaseVector<V2>& v2) : 
            itssize(v2.size()), itsm(itssize)
        {
            TMVStaticAssert(Traits<type>::okA);
            typename type::diag_type d = this->diag();
            v2.newAssignTo(d);
        }

        template <class M2>
        DiagMatrix(const BaseMatrix<M2>& m2) :
            itssize(m2.colsize()), itsm(itssize)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert((Sizes<M2::_rowsize,M2::_colsize>::same));
            TMVAssert(m2.colsize() == m2.rowsize());
            DiagCopy<ShapeTraits2<M2::_shape,Diag>::assignable>::copy(
                m2,*this);
        }

        ~DiagMatrix() 
        {
#ifdef TMV_DEBUG
            this->diag().flatten().setAllTo(Traits<real_type>::destr_value());
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

        void swapWith(type& m2)
        {
            TMVAssert(m2.size() == size());
            if (itsm.get() == m2.itsm.get()) return;
            itsm.swapWith(m2.itsm);
        }

        void resize(const size_t n)
        {
#ifdef TMV_DEBUG
            this->diag().flatten().setAllTo(Traits<real_type>::destr_value());
#endif
            itssize = n;
            itsm.resize(n);
#ifdef TMV_DEBUG
            this->diag().flatten().setAllTo(Traits<real_type>::constr_value());
#endif
        }

        TMV_INLINE size_t size() const { return itssize; }
        TMV_INLINE int nElements() const { return itssize; }
        TMV_INLINE int step() const { return 1; }
        TMV_INLINE bool isconj() const { return false; }
        TMV_INLINE bool isrm() const { return true; }
        TMV_INLINE bool iscm() const { return true; }

    private:

        size_t itssize;
        AlignedArray<T> itsm;

    }; // DiagMatrix

    template <class T, int A0>
    struct Traits<ConstDiagMatrixView<T,A0> >
    {
        enum { A = (A0 & ~NoDivider) };
        enum { okA = (
                Attrib<A>::vectoronly &&
                ( Traits<T>::iscomplex || !Attrib<A>::conj ) )};
        enum { _attrib = A };

        typedef T value_type;
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;
        enum { isreal = Traits<T>::isreal };
        enum { iscomplex = Traits<T>::iscomplex };

        typedef ConstDiagMatrixView<T,A0> type;
        typedef const type& calc_type;
        typedef const type& eval_type;
        enum { copyA = Attrib<A>::fort ? FortranStyle : CStyle };
        typedef DiagMatrix<T,copyA> copy_type;

        enum { _colsize = TMV_UNKNOWN };
        enum { _rowsize = TMV_UNKNOWN };
        enum { _size = TMV_UNKNOWN };
        enum { _fort = Attrib<A>::fort };
        enum { _shape = Diag };
        enum { _rowmajor = false };
        enum { _colmajor = false };
        enum { _calc = true };
        enum { _step = Attrib<A>::unit ? 1 : TMV_UNKNOWN };
        enum { _diagstep = _step };
        enum { _conj = Attrib<A>::conj };
        enum { _checkalias = !Attrib<A>::noalias };
        enum { _unit = Attrib<A>::unit };
        enum { twoS = isreal ? int(_step) : IntTraits<_step>::twoS };

        enum { _hasdivider = false };
        typedef QuotXM<1,real_type,type> inverse_type;

        enum { unitA = A | Unit };
        enum { nonunitA = A & ~Unit };
        enum { conjA = iscomplex ? (A ^ Conj) : int(A) };
        enum { nonconjA = A & ~Conj };
        enum { cstyleA = A & ~FortranStyle };
        enum { fstyleA = A | FortranStyle };
        enum { twosA = isreal ? int(A) : (A & ~Conj & ~Unit) };
        enum { Ac = _checkalias ? (A & CheckAlias) : (A & ~NoAlias) };
        enum { nonunitAc = Ac & ~Unit };
        enum { twosAc = isreal ? int(Ac) : (Ac & ~Conj & ~Unit) };

        typedef ConstVectorView<T,A> const_diag_type;
        typedef ConstDiagMatrixView<T,A> const_subdiagmatrix_type;
        typedef ConstDiagMatrixView<T,nonunitA> const_subdiagmatrix_step_type;
        typedef ConstDiagMatrixView<T,A> const_view_type;
        typedef ConstDiagMatrixView<T,cstyleA> const_cview_type;
        typedef ConstDiagMatrixView<T,fstyleA> const_fview_type;
        typedef ConstDiagMatrixView<T,(_conj ? Conj : NonConj)> 
            const_xview_type;
        typedef ConstDiagMatrixView<T,unitA> const_unitview_type;
        typedef ConstDiagMatrixView<T,conjA> const_conjugate_type;
        typedef ConstDiagMatrixView<T,A> const_transpose_type;
        typedef ConstDiagMatrixView<T,conjA> const_adjoint_type;
        typedef typename TypeSelect< iscomplex ,
                ConstSmallDiagMatrixView<real_type,TMV_UNKNOWN,twoS,twosAc> ,
                ConstDiagMatrixView<real_type,twosA> >::type const_realpart_type;
        typedef const_realpart_type const_imagpart_type;
        typedef ConstDiagMatrixView<T,nonconjA> const_nonconj_type;
        typedef DiagMatrixView<T,A> nonconst_type;
    };

    template <class T, int A>
    class ConstDiagMatrixView :
        public BaseMatrix_Diag<ConstDiagMatrixView<T,A> >
    {
    public:
        typedef ConstDiagMatrixView<T,A> type;

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

        ConstDiagMatrixView(const T* m, size_t n, int s) :
            itsm(m), itssize(n), itsstep(s) 
        {
            TMVStaticAssert(Traits<type>::okA);
        }

        ConstDiagMatrixView(const T* m, size_t n) :
            itsm(m), itssize(n), itsstep(_step)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(_step != TMV_UNKNOWN); 
        }

        ConstDiagMatrixView(const type& m2) :
            itsm(m2.cptr()), itssize(m2.size()), itsstep(m2.step()) 
        {
            TMVStaticAssert(Traits<type>::okA);
        }

        template <int A2>
        ConstDiagMatrixView(const ConstDiagMatrixView<T,A2>& m2) :
            itsm(m2.cptr()), itssize(m2.size()), itsstep(m2.step()) 
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Attrib<A>::conj == int(Attrib<A2>::conj)); 
        }

        template <int A2>
        ConstDiagMatrixView(const DiagMatrixView<T,A2>& m2) :
            itsm(m2.cptr()), itssize(m2.size()), itsstep(m2.step()) 
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Attrib<A>::conj == int(Attrib<A2>::conj)); 
        }

        template <int N2, int S2, int A2>
        ConstDiagMatrixView(
            const ConstSmallDiagMatrixView<T,N2,S2,A2>& m2) :
            itsm(m2.cptr()), itssize(m2.size()), itsstep(m2.step()) 
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Attrib<A>::conj == int(Attrib<A2>::conj)); 
        }

        template <int N2, int S2, int A2>
        ConstDiagMatrixView(
            const SmallDiagMatrixView<T,N2,S2,A2>& m2) :
            itsm(m2.cptr()), itssize(m2.size()), itsstep(m2.step()) 
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Attrib<A>::conj == int(Attrib<A2>::conj)); 
        }

        ~ConstDiagMatrixView() {
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
        const size_t itssize;
        const CheckedInt<_step> itsstep;

    }; // ConstDiagMatrixView

    template <class T, int A0>
    struct Traits<DiagMatrixView<T,A0> >
    {
        enum { A = (A0 & ~NoDivider) };
        enum { okA = (
                Attrib<A>::vectoronly &&
                ( Traits<T>::iscomplex || !Attrib<A>::conj ) )};
        enum { _attrib = A };

        typedef T value_type;
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;
        enum { isreal = Traits<T>::isreal };
        enum { iscomplex = Traits<T>::iscomplex };

        typedef DiagMatrixView<T,A0> type;
        typedef ConstDiagMatrixView<T,A> calc_type;
        typedef const type& eval_type;
        enum { copyA = Attrib<A>::fort ? FortranStyle : CStyle };
        typedef DiagMatrix<T,copyA> copy_type;

        enum { _colsize = TMV_UNKNOWN };
        enum { _rowsize = TMV_UNKNOWN };
        enum { _size = TMV_UNKNOWN };
        enum { _fort = Attrib<A>::fort };
        enum { _shape = Diag };
        enum { _rowmajor = false };
        enum { _colmajor = false };
        enum { _calc = true };
        enum { _step = Attrib<A>::unit ? 1 : TMV_UNKNOWN };
        enum { _diagstep = _step };
        enum { _conj = Attrib<A>::conj };
        enum { _checkalias = !Attrib<A>::noalias };
        enum { _unit = Attrib<A>::unit };
        enum { twoS = isreal ? int(_step) : IntTraits<_step>::twoS };

        enum { _hasdivider = false };
        typedef QuotXM<1,real_type,type> inverse_type;

        enum { unitA = A | Unit };
        enum { nonunitA = A & ~Unit };
        enum { conjA = iscomplex ? (A ^ Conj) : int(A) };
        enum { nonconjA = A & ~Conj };
        enum { cstyleA = A & ~FortranStyle };
        enum { fstyleA = A | FortranStyle };
        enum { twosA = isreal ? int(A) : (A & ~Conj & ~Unit) };
        enum { Ac = _checkalias ? (A & CheckAlias) : (A & ~NoAlias) };
        enum { nonunitAc = Ac & ~Unit };
        enum { twosAc = isreal ? int(Ac) : (Ac & ~Conj & ~Unit) };

        typedef ConstVectorView<T,A> const_diag_type;
        typedef ConstDiagMatrixView<T,A> const_subdiagmatrix_type;
        typedef ConstDiagMatrixView<T,nonunitA> const_subdiagmatrix_step_type;
        typedef ConstDiagMatrixView<T,A> const_view_type;
        typedef ConstDiagMatrixView<T,cstyleA> const_cview_type;
        typedef ConstDiagMatrixView<T,fstyleA> const_fview_type;
        typedef ConstDiagMatrixView<T,(_conj ? Conj : NonConj)> 
            const_xview_type;
        typedef ConstDiagMatrixView<T,unitA> const_unitview_type;
        typedef ConstDiagMatrixView<T,conjA> const_conjugate_type;
        typedef ConstDiagMatrixView<T,A> const_transpose_type;
        typedef ConstDiagMatrixView<T,conjA> const_adjoint_type;
        typedef typename TypeSelect< iscomplex ,
                ConstSmallDiagMatrixView<real_type,TMV_UNKNOWN,twoS,twosAc> ,
                ConstDiagMatrixView<real_type,twosA> >::type const_realpart_type;
        typedef const_realpart_type const_imagpart_type;
        typedef ConstDiagMatrixView<T,nonconjA> const_nonconj_type;
        typedef DiagMatrixView<T,A> nonconst_type;

        typedef T& reference;

        typedef VectorView<T,A> diag_type;
        typedef DiagMatrixView<T,A> subdiagmatrix_type;
        typedef DiagMatrixView<T,nonunitA> subdiagmatrix_step_type;
        typedef DiagMatrixView<T,A> view_type;
        typedef DiagMatrixView<T,cstyleA> cview_type;
        typedef DiagMatrixView<T,fstyleA> fview_type;
        typedef DiagMatrixView<T,(_conj ? Conj : NonConj)> xview_type;
        typedef DiagMatrixView<T,unitA> unitview_type;
        typedef DiagMatrixView<T,conjA> conjugate_type;
        typedef DiagMatrixView<T,A> transpose_type;
        typedef DiagMatrixView<T,conjA> adjoint_type;
        typedef typename TypeSelect< iscomplex ,
                SmallDiagMatrixView<real_type,TMV_UNKNOWN,twoS,twosAc> ,
                DiagMatrixView<real_type,twosA> >::type realpart_type;
        typedef realpart_type imagpart_type;
        typedef DiagMatrixView<T,nonconjA> nonconj_type;
    };

    template <class T, int A>
    class DiagMatrixView :
        public BaseMatrix_Diag_Mutable<DiagMatrixView<T,A> >
    {
    public:
        typedef DiagMatrixView<T,A> type;
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

        DiagMatrixView(T* m, size_t n, int s) :
            itsm(m), itssize(n), itsstep(s) 
        {
            TMVStaticAssert(Traits<type>::okA);
        }

        DiagMatrixView(T* m, size_t n) :
            itsm(m), itssize(n), itsstep(_step)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(_step != TMV_UNKNOWN); 
        }

        DiagMatrixView(const type& m2) :
            itsm(m2.itsm), itssize(m2.size()), itsstep(m2.step()) 
        {
            TMVStaticAssert(Traits<type>::okA);
        }

        template <int A2>
        DiagMatrixView(DiagMatrixView<T,A2> m2) :
            itsm(m2.ptr()), itssize(m2.size()), itsstep(m2.step()) 
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Attrib<A>::conj == int(Attrib<A2>::conj)); 
        }

        template <int N2, int S2, int A2>
        DiagMatrixView(SmallDiagMatrixView<T,N2,S2,A2> m2) :
            itsm(m2.ptr()), itssize(m2.size()), itsstep(m2.step()) 
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Attrib<A>::conj == int(Attrib<A2>::conj)); 
        }

        ~DiagMatrixView() {
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
        const CheckedInt<_step> itsstep;

    }; // DiagMatrixView


    //
    // Special Constructors:
    //   DiagMatrixViewOf(v)
    //   DiagMatrixViewOf(T* v, n)
    //

    template <class T>
    static TMV_INLINE ConstDiagMatrixView<T,Unit> DiagMatrixViewOf(
        const T* v, size_t size)
    { return ConstDiagMatrixView<T,Unit>(v,size); }

    template <class T>
    static TMV_INLINE ConstDiagMatrixView<T> DiagMatrixViewOf(
        const T* v, size_t size, int step)
    { return ConstDiagMatrixView<T>(v,size,step); }

    template <class T>
    static TMV_INLINE DiagMatrixView<T,Unit> DiagMatrixViewOf(T* v, size_t size)
    { return DiagMatrixView<T,Unit>(v,size); }

    template <class T>
    static TMV_INLINE DiagMatrixView<T> DiagMatrixViewOf(
        T* v, size_t size, int step)
    { return DiagMatrixView<T>(v,size,step); }


    //
    // Swap
    //

    template <class T, int A>
    static TMV_INLINE void Swap(DiagMatrix<T,A>& m1, DiagMatrix<T,A>& m2)
    { m1.swapWith(m2); }
    template <class M, class T, int A>
    static TMV_INLINE void Swap(
        BaseMatrix_Diag_Mutable<M>& m1, DiagMatrixView<T,A> m2)
    { Swap(m1.diag(),m2.diag()); }
    template <class M, class T, int A>
    static TMV_INLINE void Swap(
        DiagMatrixView<T,A> m1, BaseMatrix_Diag_Mutable<M>& m2)
    { Swap(m1.diag(),m2.diag()); }
    template <class T, int A1, int A2>
    static TMV_INLINE void Swap(
        DiagMatrixView<T,A1> m1, DiagMatrixView<T,A2> m2)
    { Swap(m1.diag(),m2.diag()); }


    //
    // Conjugate, Transpose, Adjoint
    //

    template <class T, int A>
    static TMV_INLINE typename DiagMatrix<T,A>::conjugate_type Conjugate(
        DiagMatrix<T,A>& m)
    { return m.conjugate(); }
    template <class T, int A>
    static TMV_INLINE typename DiagMatrixView<T,A>::conjugate_type Conjugate(
        DiagMatrixView<T,A> m)
    { return m.conjugate(); }

    template <class T, int A>
    static TMV_INLINE typename DiagMatrix<T,A>::transpose_type Transpose(
        DiagMatrix<T,A>& m)
    { return m.transpose(); }
    template <class T, int A>
    static TMV_INLINE typename DiagMatrixView<T,A>::transpose_type Transpose(
        DiagMatrixView<T,A> m)
    { return m.transpose(); }

    template <class T, int A>
    static TMV_INLINE typename DiagMatrix<T,A>::adjoint_type Adjoint(
        DiagMatrix<T,A>& m)
    { return m.adjoint(); }
    template <class T, int A>
    static TMV_INLINE typename DiagMatrixView<T,A>::adjoint_type Adjoint(
        DiagMatrixView<T,A> m)
    { return m.adjoint(); }


    //
    // TMV_Text 
    //

#ifdef TMV_TEXT
    template <class T, int A>
    static inline std::string TMV_Text(const DiagMatrix<T,A>& m)
    {
        std::ostringstream s;
        s << "DiagMatrix<"<<TMV_Text(T());
        s << ","<<Attrib<A>::vtext()<<">";
        s << "("<<m.size()<<","<<m.step()<<")";
        return s.str();
    }

    template <class T, int A>
    static inline std::string TMV_Text(const ConstDiagMatrixView<T,A>& m)
    {
        std::ostringstream s;
        s << "ConstDiagMatrixView<"<<TMV_Text(T());
        s << ","<<Attrib<A>::vtext()<<">";
        s << "("<<m.size()<<","<<m.step()<<")";
        return s.str();
    }

    template <class T, int A>
    static inline std::string TMV_Text(const DiagMatrixView<T,A>& m)
    {
        std::ostringstream s;
        s << "DiagMatrixView<"<<TMV_Text(T());
        s << ","<<Attrib<A>::vtext()<<">";
        s << "("<<m.size()<<","<<m.step()<<")";
        return s.str();
    }
#endif

} // namespace tmv

#endif
