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
//    SmallDiagMatrix<T,N,I>()
//        Makes a DiagMatrix with column size and row size = size
//        with _uninitialized_ values
//
//    SmallDiagMatrix<T,N,I>(T x)
//        Makes a DiagMatrix of size n with all values = x
//
//    SmallDiagMatrix<T,N,I>(T* vv)
//    SmallDiagMatrix<T,N,I>(const std::vector<T>& vv)
//        Makes a DiagMatrix of size n which copies the values is vv
//
//    SmallDiagMatrix<T,N,I>(const BaseVector<V>& vv)
//        Make a DiagMatrix which copies the elements of vv.
//
// All the other operations with a DiagMatrix work the same for 
// SmallDiagMatrix.
//

#ifndef TMV_SmallDiagMatrix_H
#define TMV_SmallDiagMatrix_H

#include "TMV_BaseMatrix_Diag.h"
#include "TMV_VIt.h"
#include "TMV_DiagMatrix.h"
#include <vector>

namespace tmv {

    //
    // SmallDiagMatrix
    //

    template <class T, int N, IndexStyle I>
    struct Traits<SmallDiagMatrix<T,N,I> >
    {
        typedef T value_type;

        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;
        enum { isreal = Traits<T>::isreal };
        enum { iscomplex = Traits<T>::iscomplex };

        typedef SmallDiagMatrix<T,N,I> type;
        typedef const type& calc_type;
        typedef const type& eval_type;
        typedef type copy_type;

        enum { _colsize = N };
        enum { _rowsize = N };
        enum { _size = N };
        enum { _fort = (I == FortranStyle) };
        enum { _shape = Diag };
        enum { _rowmajor = false };
        enum { _colmajor = false };
        enum { _stor = ColMajor }; // arbitrary
        enum { _calc = true };
        enum { _diagstep = 1 };
        enum { _conj = false };
        enum { twoS = 2 };
        enum { notC = iscomplex };
        enum { _hasdivider = false };

        typedef ConstSmallVectorView<T,N,_diagstep,false,I> const_diag_type;

        typedef ConstDiagMatrixView<T,_diagstep,false,I> const_subdiagmatrix_type;
        typedef ConstDiagMatrixView<T,UNKNOWN,false,I>
            const_subdiagmatrix_step_type;
        typedef ConstSmallDiagMatrixView<T,N,_diagstep,false,I> const_view_type;
        typedef ConstSmallDiagMatrixView<T,N,_diagstep,false,CStyle>
            const_cview_type;
        typedef ConstSmallDiagMatrixView<T,N,_diagstep,false,FortranStyle>
            const_fview_type;
        typedef ConstDiagMatrixView<T> const_xview_type;
        typedef ConstSmallDiagMatrixView<T,N,1,false> const_unitview_type;
        typedef ConstSmallDiagMatrixView<T,N,_diagstep,notC,I>
            const_conjugate_type;
        typedef ConstSmallDiagMatrixView<T,N,_diagstep,false,I>
            const_transpose_type;
        typedef ConstSmallDiagMatrixView<T,N,_diagstep,notC,I> const_adjoint_type;
        typedef ConstSmallDiagMatrixView<real_type,N,twoS,false,I>
            const_realpart_type;
        typedef const_realpart_type const_imagpart_type;
        typedef ConstSmallDiagMatrixView<T,N,_diagstep,false,I> const_nonconj_type;
        typedef SmallDiagMatrixView<T,N,_diagstep,false,I> nonconst_type;

        typedef QuotXM<1,real_type,type> inverse_type;

        typedef T& reference;

        typedef SmallVectorView<T,N,_diagstep,false,I> diag_type;

        typedef DiagMatrixView<T,_diagstep,false,I> subdiagmatrix_type;
        typedef DiagMatrixView<T,UNKNOWN,false,I> subdiagmatrix_step_type;
        typedef SmallDiagMatrixView<T,N,_diagstep,false,I> view_type;
        typedef SmallDiagMatrixView<T,N,_diagstep,false,CStyle> cview_type;
        typedef SmallDiagMatrixView<T,N,_diagstep,false,FortranStyle> fview_type;
        typedef DiagMatrixView<T> xview_type;
        typedef SmallDiagMatrixView<T,N,1,false> unitview_type;
        typedef SmallDiagMatrixView<T,N,_diagstep,notC,I> conjugate_type;
        typedef SmallDiagMatrixView<T,N,_diagstep,false,I> transpose_type;
        typedef SmallDiagMatrixView<T,N,_diagstep,notC,I> adjoint_type;
        typedef SmallDiagMatrixView<real_type,N,twoS,false,I> realpart_type;
        typedef realpart_type imagpart_type;
        typedef SmallDiagMatrixView<T,N,_diagstep,false,I> nonconj_type;
    };

    template <class T, int N, IndexStyle I>
    class SmallDiagMatrix : 
        public BaseMatrix_Diag_Mutable<SmallDiagMatrix<T,N,I> >
    {
    public:

        typedef SmallDiagMatrix<T,N,I> type;
        typedef BaseMatrix_Diag_Mutable<type> base_mut;

        enum { _colsize = Traits<type>::_size };
        enum { _rowsize = Traits<type>::_size };
        enum { _size = Traits<type>::_size };
        enum { _fort = Traits<type>::_fort };
        enum { _shape = Traits<type>::_shape };
        enum { _rowmajor = Traits<type>::_rowmajor };
        enum { _colmajor = Traits<type>::_colmajor };
        enum { _calc = Traits<type>::_calc };
        enum { _diagstep = Traits<type>::_diagstep };
        enum { _conj = Traits<type>::_conj };

        //
        // Constructors
        //

        SmallDiagMatrix(size_t n=N)
        {
            TMVStaticAssert(N >= 0);
            TMVAssert(n==N);
#ifdef TMV_DEBUG
            this->setAllTo(T(888));
#endif
        }

        explicit SmallDiagMatrix(size_t n, T x) 
        {
            TMVStaticAssert(N >= 0);
            TMVAssert(n==N);
            this->setAllTo(x);
        }

        explicit SmallDiagMatrix(size_t n, const T* vv) 
        {
            TMVStaticAssert(N >= 0);
            TMVAssert(n==N);
            ConstSmallDiagMatrixView<T,N,1>(vv).newAssignTo(*this);
        }

        SmallDiagMatrix(const std::vector<T>& vv)
        {
            TMVStaticAssert(N >= 0);
            TMVAssert(vv.size() == N);
            ConstSmallDiagMatrixView<T,N,1>(&vv[0]).newAssignTo(*this);
        }

        SmallDiagMatrix(const type& m2) 
        {
            TMVStaticAssert(N >= 0);
            m2.newAssignTo(*this);
        }

        template <class V2>
        SmallDiagMatrix(const BaseVector<V2>& v2) 
        {
            TMVStaticAssert(N >= 0);
            typename type::diag_type d = this->diag();
            v2.newAssignTo(d);
        }

        template <class M2>
        SmallDiagMatrix(const BaseMatrix<M2>& m2) 
        {
            TMVStaticAssert(N >= 0);
            TMVStaticAssert((Sizes<M2::_rowsize,M2::_colsize>::same));
            TMVAssert(m2.colsize() == m2.rowsize());
            DiagCopy<ShapeTraits2<M2::_shape,Diag>::assignable>::copy(
                m2,*this);
        }

        ~SmallDiagMatrix() 
        {
#ifdef TMV_DEBUG
            this->setAllTo(T(999));
#endif
        }

        //
        // Op=
        //

        type& operator=(const type& m2)
        {
            if (&m2 != this) base_mut::operator=(m2);
            return *this; 
        }

        template <class M2>
        type& operator=(const BaseMatrix<M2>& m2)
        {
            base_mut::operator=(m2);
            return *this; 
        }

        type& operator=(T x)
        {
            base_mut::operator=(x);
            return *this; 
        }


        //
        // Auxilliary Functions
        //

        const T* cptr() const { return itsm; }
        T* ptr() { return itsm; }

        T cref(int i, int j) const { return (i!=j ? T(0) : cref(i)); }
        T cref(int i) const { return itsm[i]; }
        T& ref(int i) { return itsm[i]; }

        size_t size() const { return N; }
        int nElements() const { return N; }
        int step() const { return 1; }
        bool isconj() const { return false; }
        bool isrm() const { return true; }
        bool iscm() const { return true; }

    private:

        StackArray<T,N> itsm;

    }; // SmallDiagMatrix

    template <class T, int N>
    class SmallDiagMatrixF : 
        public SmallDiagMatrix<T,N,FortranStyle>
    {
    public:
        typedef SmallDiagMatrixF<T,N> type;
        typedef SmallDiagMatrix<T,N,FortranStyle> mtype;

        SmallDiagMatrixF() {}
        explicit SmallDiagMatrixF(T x) : mtype(x) {}
        explicit SmallDiagMatrixF(const T* vv) : mtype(vv) {}
        explicit SmallDiagMatrixF(const std::vector<T>& vv) :
            mtype(vv) {}
        SmallDiagMatrixF(const type& m2) : mtype(m2) {}
        template <class M2>
        SmallDiagMatrixF(const BaseMatrix_Diag<M2>& m2) : mtype(m2) {}
        ~SmallDiagMatrixF() {}

        type& operator=(const type& m2)
        { mtype::operator=(m2); return *this; }
        template <class M2>
        type& operator=(const BaseMatrix<M2>& m2)
        { mtype::operator=(m2); return *this; }
        type& operator=(T x)
        { mtype::operator=(x); return *this; }

    }; // SmallDiagMatrixF

    template <class T, int N, int S, bool C, IndexStyle I>
    struct Traits<ConstSmallDiagMatrixView<T,N,S,C,I> >
    {
        typedef T value_type;

        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;
        enum { isreal = Traits<T>::isreal };
        enum { iscomplex = Traits<T>::iscomplex };

        typedef ConstSmallDiagMatrixView<T,N,S,C,I> type;
        typedef const type& calc_type;
        typedef const type& eval_type;
        typedef SmallDiagMatrix<T,I> copy_type;

        enum { _colsize = N };
        enum { _rowsize = N };
        enum { _size = N };
        enum { _fort = (I == FortranStyle) };
        enum { _shape = Diag };
        enum { _rowmajor = (S==1) };
        enum { _colmajor = (S==1) };
        enum { _stor = ColMajor }; // arbitrary
        enum { _calc = true };
        enum { _diagstep = S };
        enum { _conj = C };
        enum { twoS = isreal ? S : IntTraits<S>::twoS };
        enum { notC = !C && iscomplex };
        enum { _hasdivider = false };

        typedef ConstSmallVectorView<T,N,S,C,I> const_diag_type;

        typedef ConstDiagMatrixView<T,S,C,I> const_subdiagmatrix_type;
        typedef ConstDiagMatrixView<T,UNKNOWN,C,I>
            const_subdiagmatrix_step_type;
        typedef ConstSmallDiagMatrixView<T,N,S,C,I> const_view_type;
        typedef ConstSmallDiagMatrixView<T,N,S,C,CStyle> const_cview_type;
        typedef ConstSmallDiagMatrixView<T,N,S,C,FortranStyle>
            const_fview_type;
        typedef ConstDiagMatrixView<T,UNKNOWN,C> const_xview_type;
        typedef ConstSmallDiagMatrixView<T,N,1,C> const_unitview_type;
        typedef ConstSmallDiagMatrixView<T,N,S,notC,I> const_conjugate_type;
        typedef ConstSmallDiagMatrixView<T,N,S,C,I> const_transpose_type;
        typedef ConstSmallDiagMatrixView<T,N,S,notC,I> const_adjoint_type;
        typedef ConstSmallDiagMatrixView<real_type,N,twoS,C,I>
            const_realpart_type;
        typedef const_realpart_type const_imagpart_type;
        typedef ConstSmallDiagMatrixView<T,N,S,false,I> const_nonconj_type;
        typedef SmallDiagMatrixView<T,N,S,C,I> nonconst_type;

        typedef QuotXM<1,real_type,type> inverse_type;
    };

    template <class T, int N, int S, bool C, IndexStyle I>
    class ConstSmallDiagMatrixView :
        public BaseMatrix_Diag<ConstSmallDiagMatrixView<T,N,S,C,I> >
    {
    public:

        typedef ConstSmallDiagMatrixView<T,N,S,C,I> type;

        enum { _colsize = Traits<type>::_size };
        enum { _rowsize = Traits<type>::_size };
        enum { _size = Traits<type>::_size };
        enum { _fort = Traits<type>::_fort };
        enum { _shape = Traits<type>::_shape };
        enum { _rowmajor = Traits<type>::_rowmajor };
        enum { _colmajor = Traits<type>::_colmajor };
        enum { _calc = Traits<type>::_calc };
        enum { _diagstep = Traits<type>::_diagstep };
        enum { _conj = Traits<type>::_conj };

        //
        // Constructors
        //

        ConstSmallDiagMatrixView(const T* m, size_t n, int s) :
            itsm(m), itssize(n), itsstep(s) {}

        ConstSmallDiagMatrixView(const T* m, size_t n) :
            itsm(m), itssize(n), itsstep(S)
        { TMVStaticAssert(S != UNKNOWN); }

        ConstSmallDiagMatrixView(const T* m) :
            itsm(m), itssize(N), itsstep(S)
        {
            TMVStaticAssert(N != UNKNOWN); 
            TMVStaticAssert(S != UNKNOWN); 
        }

        ConstSmallDiagMatrixView(const type& m2) :
            itsm(m2.cptr()), itssize(m2.size()), itsstep(m2.step()) {}

        template <int S2, IndexStyle I2>
        ConstSmallDiagMatrixView(
            const ConstDiagMatrixView<T,S2,C,I2>& m2
        ) :
            itsm(m2.cptr()), itssize(m2.size()), itsstep(m2.step()) {}

        template <int S2, IndexStyle I2>
        ConstSmallDiagMatrixView(const DiagMatrixView<T,S2,C,I2>& m2) :
            itsm(m2.cptr()), itssize(m2.size()), itsstep(m2.step()) {}

        template <int N2, int S2, IndexStyle I2>
        ConstSmallDiagMatrixView(
            const ConstSmallDiagMatrixView<T,N2,S2,C,I2>& m2) :
            itsm(m2.cptr()), itssize(m2.size()), itsstep(m2.step()) {}

        template <int N2, int S2, IndexStyle I2>
        ConstSmallDiagMatrixView(
            const SmallDiagMatrixView<T,N2,S2,C,I2>& m2) :
            itsm(m2.cptr()), itssize(m2.size()), itsstep(m2.step()) {}

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

        const T* cptr() const { return itsm; }

        T cref(int i, int j) const { return (i!=j ? T(0) : cref(i)); }
        T cref(int i) const
        { return DoConj<C>(itsm[i*step()]); }

        size_t size() const { return itssize; }
        int nElements() const { return itssize; }
        int step() const { return itsstep; }
        bool isconj() const { return C; }
        bool isrm() const { return step()==1; }
        bool iscm() const { return step()==1; }

    private :

        const T* itsm;
        const CheckedInt<N> itssize;
        const CheckedInt<S> itsstep;

    }; // ConstSmallDiagMatrixView

    template <class T, int N, int S=1, bool C=false>
    class ConstSmallDiagMatrixViewF :
        public ConstSmallDiagMatrixView<T,N,S,C,FortranStyle>
    {
    public:

        typedef ConstSmallDiagMatrixViewF<T,N,S,C> type;
        typedef ConstSmallDiagMatrixView<T,N,S,C,FortranStyle> mtype;

        ConstSmallDiagMatrixViewF(const T* m, size_t n, int s) : 
            mtype(m,n,s) {}
        ConstSmallDiagMatrixViewF(const T* m, size_t n) : mtype(m,n) {}
        ConstSmallDiagMatrixViewF(const T* m) : mtype(m) {}
        ConstSmallDiagMatrixViewF(const type& m2) : mtype(m2) {}
        template <int S2, IndexStyle I2>
        ConstSmallDiagMatrixViewF(
            const ConstSmallDiagMatrixView<T,S2,C,I2>& m2) : mtype(m2) {}
        template <int S2, IndexStyle I2>
        ConstSmallDiagMatrixViewF(
            const SmallDiagMatrixView<T,S2,C,I2>& m2
        ) :
            mtype(m2) {}
        template <int N2, int S2, IndexStyle I2>
        ConstSmallDiagMatrixViewF(
            const ConstSmallDiagMatrixView<T,N2,S2,C,I2>& m2) : mtype(m2) {}
        template <int N2, int S2, IndexStyle I2>
        ConstSmallDiagMatrixViewF(
            const SmallDiagMatrixView<T,N2,S2,C,I2>& m2) :
            mtype(m2) {}
        ~ConstSmallDiagMatrixViewF() {}

    private :
        void operator=(const type& m2);

    }; // ConstSmallDiagMatrixViewF


    template <class T, int N, int S, bool C, IndexStyle I>
    struct Traits<SmallDiagMatrixView<T,N,S,C,I> >
    {
        typedef T value_type;

        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;
        enum { isreal = Traits<T>::isreal };
        enum { iscomplex = Traits<T>::iscomplex };

        typedef SmallDiagMatrixView<T,N,S,C,I> type;
        typedef const ConstSmallDiagMatrixView<T,N,S,C,I> calc_type;
        typedef calc_type eval_type;
        typedef SmallDiagMatrix<T,N,I> copy_type;

        enum { _colsize = N };
        enum { _rowsize = N };
        enum { _size = N };
        enum { _fort = (I == FortranStyle) };
        enum { _shape = Diag };
        enum { _rowmajor = (S==1) }; 
        enum { _colmajor = (S==1) }; 
        enum { _stor = ColMajor }; // arbitrary
        enum { _calc = true };
        enum { _diagstep = S };
        enum { _conj = C };
        enum { twoS = isreal ? S : IntTraits<S>::twoS };
        enum { notC = !C && iscomplex };
        enum { _hasdivider = false };

        typedef ConstSmallVectorView<T,N,S,C,I> const_diag_type;

        typedef ConstDiagMatrixView<T,S,C,I> const_subdiagmatrix_type;
        typedef ConstDiagMatrixView<T,UNKNOWN,C,I>
            const_subdiagmatrix_step_type;
        typedef ConstSmallDiagMatrixView<T,N,S,C,I> const_view_type;
        typedef ConstSmallDiagMatrixView<T,N,S,C,CStyle> const_cview_type;
        typedef ConstSmallDiagMatrixView<T,N,S,C,FortranStyle>
            const_fview_type;
        typedef ConstDiagMatrixView<T,UNKNOWN,C> const_xview_type;
        typedef ConstSmallDiagMatrixView<T,N,1,C> const_unitview_type;
        typedef ConstSmallDiagMatrixView<T,N,S,notC,I> const_conjugate_type;
        typedef ConstSmallDiagMatrixView<T,N,S,C,I> const_transpose_type;
        typedef ConstSmallDiagMatrixView<T,N,S,notC,I> const_adjoint_type;
        typedef ConstSmallDiagMatrixView<real_type,N,twoS,C,I>
            const_realpart_type;
        typedef const_realpart_type const_imagpart_type;
        typedef ConstSmallDiagMatrixView<T,N,S,false,I> const_nonconj_type;
        typedef SmallDiagMatrixView<T,N,S,C,I> nonconst_type;

        typedef QuotXM<1,real_type,type> inverse_type;

        typedef typename AuxRef<T,C>::reference reference;

        typedef SmallVectorView<T,N,S,C,I> diag_type;

        typedef DiagMatrixView<T,S,C,I> subdiagmatrix_type;
        typedef DiagMatrixView<T,UNKNOWN,C,I> subdiagmatrix_step_type;
        typedef SmallDiagMatrixView<T,N,S,C,I> view_type;
        typedef SmallDiagMatrixView<T,N,S,C,CStyle> cview_type;
        typedef SmallDiagMatrixView<T,N,S,C,FortranStyle> fview_type;
        typedef DiagMatrixView<T,UNKNOWN,C> xview_type;
        typedef SmallDiagMatrixView<T,N,1,C> unitview_type;
        typedef SmallDiagMatrixView<T,N,S,notC,I> conjugate_type;
        typedef SmallDiagMatrixView<T,N,S,C,I> transpose_type;
        typedef SmallDiagMatrixView<T,N,S,notC,I> adjoint_type;
        typedef SmallDiagMatrixView<real_type,N,twoS,C,I> realpart_type;
        typedef realpart_type imagpart_type;
        typedef SmallDiagMatrixView<T,N,S,false,I> nonconj_type;
    };

    template <class T, int N, int S, bool C, IndexStyle I>
    class SmallDiagMatrixView :
        public BaseMatrix_Diag_Mutable<SmallDiagMatrixView<T,N,S,C,I> >
    {
    public:

        typedef SmallDiagMatrixView<T,N,S,C,I> type;
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
        enum { _diagstep = Traits<type>::_diagstep };
        enum { _conj = Traits<type>::_conj };

        //
        // Constructors
        //

        SmallDiagMatrixView(T* m, size_t n, int s) :
            itsm(m), itssize(n), itsstep(s) {}

        SmallDiagMatrixView(T* m, size_t n) :
            itsm(m), itssize(n), itsstep(S)
        { TMVStaticAssert(S != UNKNOWN); }

        SmallDiagMatrixView(T* m) : itsm(m), itssize(N), itsstep(S)
        { TMVStaticAssert(N != UNKNOWN); TMVStaticAssert(S != UNKNOWN); }

        SmallDiagMatrixView(const type& m2) :
            itsm(m2.itsm), itssize(m2.size()), itsstep(m2.step()) {}

        template <int S2, IndexStyle I2>
        SmallDiagMatrixView(DiagMatrixView<T,S2,C,I2> m2) :
            itsm(m2.ptr()), itssize(m2.size()), itsstep(m2.step()) {}

        template <int N2, int S2, IndexStyle I2>
        SmallDiagMatrixView(SmallDiagMatrixView<T,N2,S2,C,I2> m2) :
            itsm(m2.ptr()), itssize(m2.size()), itsstep(m2.step()) {}

        ~SmallDiagMatrixView() {
#ifdef TMV_DEBUG
            itsm = 0; 
#endif
        }


        //
        // Op = 
        //

        template <class M2>
        type& operator=(const BaseMatrix<M2>& m2)
        {
            base_mut::operator=(m2); 
            return *this;
        }

        type& operator=(const type& m2)
        {
            base_mut::operator=(m2); 
            return *this;
        }

        type& operator=(const T x)
        {
            base_mut::operator=(x); 
            return *this;
        }


        //
        // Auxilliary Functions
        //

        const T* cptr() const { return itsm; }
        T* ptr() { return itsm; }

        T cref(int i, int j) const { return (i!=j ? T(0) : cref(i)); }
        T cref(int i) const
        { return DoConj<C>(itsm[i*step()]); }

        reference ref(int i) 
        { return reference(itsm[i*step()]); }

        size_t size() const { return itssize; }
        int nElements() const { return itssize; }
        int step() const { return itsstep; }
        bool isconj() const { return C; }
        bool isrm() const { return step()==1; }
        bool iscm() const { return step()==1; }

    private :

        T* itsm;
        const size_t itssize;
        const CheckedInt<S> itsstep;

    }; // SmallDiagMatrixView

    template <class T, int N, int S=1, bool C=false>
    class SmallDiagMatrixViewF :
        public SmallDiagMatrixView<T,N,S,C,FortranStyle>
    {
    public:

        typedef SmallDiagMatrixViewF<T,N,S,C> type;
        typedef SmallDiagMatrixView<T,N,S,C,FortranStyle> mtype;

        SmallDiagMatrixViewF(T* m, size_t n, int s) : mtype(m,n,s) {}
        SmallDiagMatrixViewF(T* m, size_t n) : mtype(m,n) {}
        SmallDiagMatrixViewF(T* m) : mtype(m) {}

        SmallDiagMatrixViewF(const type& m2) : mtype(m2) {}
        template <int S2, IndexStyle I2>
        SmallDiagMatrixViewF(DiagMatrixView<T,S2,C,I2> m2) :
            mtype(m2) {}
        template <int N2, int S2, IndexStyle I2>
        SmallDiagMatrixViewF(SmallMatrixView<T,N2,S2,C,I2> m2) :
            mtype(m2) {}
        ~SmallDiagMatrixViewF() {}

        template <class M2>
        type& operator=(const BaseMatrix<M2>& m2)
        { mtype::operator=(m2); return *this; }
        type& operator=(const type& m2)
        { mtype::operator=(m2); return *this; }
        type& operator=(const T x)
        { mtype::operator=(x); return *this; }

    }; // SmallDiagMatrixViewF


    //
    // DiagMatrixViewOf
    //

    template <class T, int N, int S, bool C, IndexStyle I>
    static inline ConstSmallDiagMatrixView<T,N,S,C,I> DiagMatrixViewOf(
        const ConstSmallVectorView<T,N,S,C,I>& v)
    { return ConstSmallDiagMatrixView<T,N,S,C,I>(v.cptr(),v.size(),v.step()); }

    template <class T, int N, int S, bool C, IndexStyle I>
    static inline SmallDiagMatrixView<T,N,S,C,I> DiagMatrixViewOf(
        SmallVectorView<T,N,S,C,I> v)
    { return SmallDiagMatrixView<T,N,S,C,I>(v.ptr(),v.size(),v.step()); }

    template <class T, int N, IndexStyle I>
    static inline SmallDiagMatrixView<T,N,1,false,I> DiagMatrixViewOf(
        SmallVector<T,N,I>& v)
    { return SmallDiagMatrixView<T,N,1,false,I>(v.ptr()); }
    
    
    //
    // Swap
    //

    template <class T, int N, int S, bool C, IndexStyle I, class MM>
    static inline void Swap(
        BaseMatrix_Diag_Mutable<MM>& m1, SmallDiagMatrixView<T,N,S,C,I> m2)
    { Swap(m1.diag(),m2.diag()); }
    template <class T, int N, int S, bool C, IndexStyle I, class MM>
    static inline void Swap(
        SmallDiagMatrixView<T,N,S,C,I> m1, BaseMatrix_Diag_Mutable<MM>& m2)
    { Swap(m1.diag(),m2.diag()); }
    template <class T, int N, int S1, bool C1, IndexStyle I1, int S2, bool C2, IndexStyle I2>
    static inline void Swap(
        SmallDiagMatrixView<T,N,S1,C1,I1> m1,
        SmallDiagMatrixView<T,N,S2,C2,I2> m2)
    { Swap(m1.diag(),m2.diag()); }
    template <class T, int N, int S1, bool C1, IndexStyle I1, int S2, bool C2, IndexStyle I2>
    static inline void Swap(
        SmallDiagMatrixView<T,N,S1,C1,I1> m1, 
        DiagMatrixView<T,S2,C2,I2> m2)
    { Swap(m1.diag(),m2.diag()); }
    template <class T, int N, int S1, bool C1, IndexStyle I1, int S2, bool C2, IndexStyle I2>
    static inline void Swap(
        DiagMatrixView<T,S1,C1,I1> m1,
        SmallDiagMatrixView<T,N,S2,C2,I2> m2)
    { Swap(m1.diag(),m2.diag()); }


    //
    // Conjugate, Transpose, Adjoint
    //
    
    template <class T, int N, IndexStyle I>
    static inline typename SmallDiagMatrix<T,N,I>::conjugate_type Conjugate(
        SmallDiagMatrix<T,N,I>& m)
    { return m.conjugate(); }
    template <class T, int N, int S, bool C, IndexStyle I>
    static inline typename SmallDiagMatrixView<T,N,S,C,I>::conjugate_type Conjugate(
        SmallDiagMatrixView<T,N,S,C,I> m)
    { return m.conjugate(); }

    template <class T, int N, IndexStyle I>
    static inline typename SmallDiagMatrix<T,N,I>::transpose_type Transpose(
        SmallDiagMatrix<T,N,I>& m)
    { return m.transpose(); }
    template <class T, int N, int S, bool C, IndexStyle I>
    static inline typename SmallDiagMatrixView<T,N,S,C,I>::transpose_type Transpose(
        SmallDiagMatrixView<T,N,S,C,I> m)
    { return m.transpose(); }

    template <class T, int N, IndexStyle I>
    static inline typename SmallDiagMatrix<T,N,I>::adjoint_type Adjoint(
        SmallDiagMatrix<T,N,I>& m)
    { return m.adjoint(); }
    template <class T, int N, int S, bool C, IndexStyle I>
    static inline typename SmallDiagMatrixView<T,N,S,C,I>::adjoint_type Adjoint(
        SmallDiagMatrixView<T,N,S,C,I> m)
    { return m.adjoint(); }


    //
    // TMV_Text 
    //

    template <class T, int N, IndexStyle I>
    static inline std::string TMV_Text(const SmallDiagMatrix<T,N,I>& )
    {
        std::ostringstream s;
        s << "SmallDiagMatrix<"<<TMV_Text(T())<<","<<N<<','<<TMV_Text(I)<<">";
        return s.str();
    }

    template <class T, int N, int S, bool C, IndexStyle I>
    static inline std::string TMV_Text(
        const ConstSmallDiagMatrixView<T,N,S,C,I>& m)
    {
        std::ostringstream s;
        s << "ConstSmallDiagMatrixView<"<<TMV_Text(T());
        s << ","<<IntTraits<N>::text();
        if (N == UNKNOWN) s << "("<<m.size()<<")";
        s << ","<<IntTraits<S>::text();
        if (S == UNKNOWN) s << "("<<m.step()<<")";
        s << ","<<C<<","<<TMV_Text(I)<<">";
        return s.str();
    }

    template <class T, int N, int S, bool C, IndexStyle I>
    static inline std::string TMV_Text(const SmallDiagMatrixView<T,N,S,C,I>& m)
    {
        std::ostringstream s;
        s << "SmallDiagMatrixView<"<<TMV_Text(T());
        s << ","<<IntTraits<N>::text();
        if (N == UNKNOWN) s << "("<<m.size()<<")";
        s << ","<<IntTraits<S>::text();
        if (S == UNKNOWN) s << "("<<m.step()<<")";
        s << ","<<C<<","<<TMV_Text(I)<<">";
        return s.str();
    }

} // namespace tmv

#endif
