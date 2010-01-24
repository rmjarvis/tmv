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
#include "TMV_StackArray.h"
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
        enum { misreal = Traits<T>::isreal };
        enum { miscomplex = Traits<T>::iscomplex };

        typedef SmallDiagMatrix<T,N,I> type;
        typedef const type& calc_type;
        typedef const type& eval_type;
        typedef type copy_type;

        enum { mcolsize = N };
        enum { mrowsize = N };
        enum { msize = N };
        enum { mfort = (I == FortranStyle) };
        enum { mshape = Diag };
        enum { mrowmajor = false };
        enum { mcolmajor = false };
        enum { mstor = ColMajor }; // arbitrary
        enum { mcalc = true };
        enum { mstep = 1 };
        enum { mconj = false };
        enum { twoS = 2 };
        enum { notC = miscomplex };

        typedef ConstSmallVectorView<T,N,mstep,false,I> const_diag_type;

        typedef ConstDiagMatrixView<T,mstep,false,I> const_subdiagmatrix_type;
        typedef ConstDiagMatrixView<T,UNKNOWN,false,I> 
            const_subdiagmatrix_step_type;
        typedef ConstSmallDiagMatrixView<T,N,mstep,false,I> const_view_type;
        typedef ConstSmallDiagMatrixView<T,N,mstep,false,CStyle> 
            const_cview_type;
        typedef ConstSmallDiagMatrixView<T,N,mstep,false,FortranStyle> 
            const_fview_type;
        typedef ConstDiagMatrixView<T> const_xview_type;
        typedef ConstSmallDiagMatrixView<T,N,1,false> const_rmview_type;
        typedef ConstSmallDiagMatrixView<T,N,1,false> const_cmview_type;
        typedef ConstSmallDiagMatrixView<T,N,mstep,notC,I> 
            const_conjugate_type;
        typedef ConstSmallDiagMatrixView<T,N,mstep,false,I> 
            const_transpose_type;
        typedef ConstSmallDiagMatrixView<T,N,mstep,notC,I> const_adjoint_type;
        typedef ConstSmallDiagMatrixView<real_type,N,twoS,false,I> 
            const_realpart_type;
        typedef const_realpart_type const_imagpart_type;
        typedef ConstSmallDiagMatrixView<T,N,mstep,false,I> const_nonconj_type;
        typedef SmallDiagMatrixView<T,N,mstep,false,I> nonconst_type;

        typedef QuotXM<1,real_type,type> inverse_type;

        typedef T& reference;

        typedef SmallVectorView<T,N,mstep,false,I> diag_type;

        typedef DiagMatrixView<T,mstep,false,I> subdiagmatrix_type;
        typedef DiagMatrixView<T,UNKNOWN,false,I> subdiagmatrix_step_type;
        typedef SmallDiagMatrixView<T,N,mstep,false,I> view_type;
        typedef SmallDiagMatrixView<T,N,mstep,false,CStyle> cview_type;
        typedef SmallDiagMatrixView<T,N,mstep,false,FortranStyle> fview_type;
        typedef DiagMatrixView<T> xview_type;
        typedef SmallDiagMatrixView<T,N,1,false> rmview_type;
        typedef SmallDiagMatrixView<T,N,1,false> cmview_type;
        typedef SmallDiagMatrixView<T,N,mstep,notC,I> conjugate_type;
        typedef SmallDiagMatrixView<T,N,mstep,false,I> transpose_type;
        typedef SmallDiagMatrixView<T,N,mstep,notC,I> adjoint_type;
        typedef SmallDiagMatrixView<real_type,N,twoS,false,I> realpart_type;
        typedef realpart_type imagpart_type;
        typedef SmallDiagMatrixView<T,N,mstep,false,I> nonconj_type;
    };

#ifdef XTEST
#ifdef TMV_DEBUG
#define XTEST_DEBUG
#endif
#endif

    template <class T, int N, IndexStyle I>
    class SmallDiagMatrix : 
        public BaseMatrix_Diag_Mutable<SmallDiagMatrix<T,N,I> >
    {
    public:

        typedef SmallDiagMatrix<T,N,I> type;
        typedef BaseMatrix_Diag_Mutable<type> base_mut;

        enum { mcolsize = Traits<type>::msize };
        enum { mrowsize = Traits<type>::msize };
        enum { msize = Traits<type>::msize };
        enum { mfort = Traits<type>::mfort };
        enum { mshape = Traits<type>::mshape };
        enum { mrowmajor = Traits<type>::mrowmajor };
        enum { mcolmajor = Traits<type>::mcolmajor };
        enum { mcalc = Traits<type>::mcalc };
        enum { mstep = Traits<type>::mstep };
        enum { mconj = Traits<type>::mconj };

        //
        // Constructors
        //

        inline SmallDiagMatrix()
        {
            TMVStaticAssert(N > 0);
#ifdef TMV_DEBUG
            this->setAllTo(T(888));
#endif
        }

        explicit inline SmallDiagMatrix(T x) 
        {
            TMVStaticAssert(N > 0);
            this->setAllTo(x);
        }

        explicit inline SmallDiagMatrix(const T* vv) 
        {
            TMVStaticAssert(N > 0);
#ifdef XTEST_DEBUG
            this->setAllTo(T(888));
#endif
            ConstSmallDiagMatrixView<T,N,1>(vv).newAssignTo(*this);
        }

        inline SmallDiagMatrix(const std::vector<T>& vv)
        {
            TMVStaticAssert(N > 0);
#ifdef XTEST_DEBUG
            this->setAllTo(T(888));
#endif
            TMVAssert(vv.size() == N);
            ConstSmallDiagMatrixView<T,N,1>(&vv[0]).newAssignTo(*this);
        }

        inline SmallDiagMatrix(const type& m2) 
        {
            TMVStaticAssert(N > 0);
#ifdef XTEST_DEBUG
            this->setAllTo(T(888));
#endif
            m2.newAssignTo(*this);
        }

        template <class V2>
        inline SmallDiagMatrix(const BaseVector<V2>& v2) 
        {
            TMVStaticAssert(N > 0);
#ifdef XTEST_DEBUG
            this->setAllTo(T(888));
#endif
            typename type::diag_type d = this->diag();
            v2.newAssignTo(d);
        }

        template <class M2>
        inline SmallDiagMatrix(const BaseMatrix<M2>& m2) 
        {
            TMVStaticAssert(N > 0);
            TMVStaticAssert((Sizes<M2::mrowsize,M2::mcolsize>::same));
            TMVAssert(m2.colsize() == m2.rowsize());
#ifdef XTEST_DEBUG
            this->setAllTo(T(888));
#endif
            DiagCopy<ShapeTraits2<M2::mshape,Diag>::assignable>::copy(
                m2,*this);
        }

        inline ~SmallDiagMatrix() 
        {
#ifdef TMV_DEBUG
            this->setAllTo(T(999));
#endif
        }

        //
        // Op=
        //

        inline type& operator=(const type& m2)
        {
            if (&m2 != this) base_mut::operator=(m2);
            return *this; 
        }

        template <class M2>
        inline type& operator=(const BaseMatrix<M2>& m2)
        {
            base_mut::operator=(m2);
            return *this; 
        }

        inline type& operator=(T x)
        {
            base_mut::operator=(x);
            return *this; 
        }


        //
        // Auxilliary Functions
        //

        inline const T* cptr() const { return itsm; }
        inline T* ptr() { return itsm; }

        inline T cref(int i, int j) const { return (i!=j ? T(0) : cref(i)); }
        inline T cref(int i) const { return itsm[i]; }
        inline T& ref(int i) { return itsm[i]; }

        inline size_t size() const { return N; }
        inline int step() const { return 1; }
        inline bool isconj() const { return false; }
        inline bool isrm() const { return true; }
        inline bool iscm() const { return true; }

    private:

        StackArray<T,N> itsm;

    }; // SmallDiagMatrix

    template <class T, int N>
    class SmallDiagMatrixF : public SmallDiagMatrix<T,N,FortranStyle>
    {
    public:
        typedef SmallDiagMatrixF<T,N> type;
        typedef SmallDiagMatrix<T,N,FortranStyle> mtype;

        inline SmallDiagMatrixF() {}
        explicit inline SmallDiagMatrixF(T x) : mtype(x) {}
        explicit inline SmallDiagMatrixF(const T* vv) : mtype(vv) {}
        explicit inline SmallDiagMatrixF(const std::vector<T>& vv) :
            mtype(vv) {}
        inline SmallDiagMatrixF(const type& m2) : mtype(m2) {}
        template <class M2>
        inline SmallDiagMatrixF(const BaseMatrix_Diag<M2>& m2) : mtype(m2) {}
        inline ~SmallDiagMatrixF() {}

        inline type& operator=(const type& m2)
        { mtype::operator=(m2); return *this; }
        template <class M2>
        inline type& operator=(const BaseMatrix<M2>& m2)
        { mtype::operator=(m2); return *this; }
        inline type& operator=(T x)
        { mtype::operator=(x); return *this; }

    }; // SmallDiagMatrixF

    template <class T, int N, int S, bool C, IndexStyle I>
    struct Traits<ConstSmallDiagMatrixView<T,N,S,C,I> >
    {
        typedef T value_type;

        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;
        enum { misreal = Traits<T>::isreal };
        enum { miscomplex = Traits<T>::iscomplex };

        typedef ConstSmallDiagMatrixView<T,N,S,C,I> type;
        typedef const type& calc_type;
        typedef const type& eval_type;
        typedef SmallDiagMatrix<T,I> copy_type;

        enum { mcolsize = N };
        enum { mrowsize = N };
        enum { msize = N };
        enum { mfort = (I == FortranStyle) };
        enum { mshape = Diag };
        enum { mrowmajor = (S==1) };
        enum { mcolmajor = (S==1) };
        enum { mstor = ColMajor }; // arbitrary
        enum { mcalc = true };
        enum { mstep = S };
        enum { mconj = C };
        enum { twoS = misreal ? S : IntTraits<S>::twoS };
        enum { notC = !C && miscomplex };

        typedef ConstSmallVectorView<T,N,S,C,I> const_diag_type;

        typedef ConstDiagMatrixView<T,S,C,I> const_subdiagmatrix_type;
        typedef ConstDiagMatrixView<T,UNKNOWN,C,I> 
            const_subdiagmatrix_step_type;
        typedef ConstSmallDiagMatrixView<T,N,S,C,I> const_view_type;
        typedef ConstSmallDiagMatrixView<T,N,S,C,CStyle> const_cview_type;
        typedef ConstSmallDiagMatrixView<T,N,S,C,FortranStyle> 
            const_fview_type;
        typedef ConstDiagMatrixView<T,UNKNOWN,C> const_xview_type;
        typedef ConstSmallDiagMatrixView<T,N,1,C> const_rmview_type;
        typedef ConstSmallDiagMatrixView<T,N,1,C> const_cmview_type;
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

        enum { mcolsize = Traits<type>::msize };
        enum { mrowsize = Traits<type>::msize };
        enum { msize = Traits<type>::msize };
        enum { mfort = Traits<type>::mfort };
        enum { mshape = Traits<type>::mshape };
        enum { mrowmajor = Traits<type>::mrowmajor };
        enum { mcolmajor = Traits<type>::mcolmajor };
        enum { mcalc = Traits<type>::mcalc };
        enum { mstep = Traits<type>::mstep };
        enum { mconj = Traits<type>::mconj };

        //
        // Constructors
        //

        inline ConstSmallDiagMatrixView(const T* m, size_t n, int s) :
            itsm(m), itssize(n), itsstep(s) {}

        inline ConstSmallDiagMatrixView(const T* m, size_t n) :
            itsm(m), itssize(n), itsstep(S)
        { TMVStaticAssert(S != UNKNOWN); }

        inline ConstSmallDiagMatrixView(const T* m) :
            itsm(m), itssize(N), itsstep(S)
        { 
            TMVStaticAssert(N != UNKNOWN); 
            TMVStaticAssert(S != UNKNOWN); 
        }

        inline ConstSmallDiagMatrixView(const type& m2) :
            itsm(m2.cptr()), itssize(m2.size()), itsstep(m2.step()) {}

        template <int S2, IndexStyle I2>
        inline ConstSmallDiagMatrixView(
            const ConstDiagMatrixView<T,S2,C,I2>& m2
        ) :
            itsm(m2.cptr()), itssize(m2.size()), itsstep(m2.step()) {}

        template <int S2, IndexStyle I2>
        inline ConstSmallDiagMatrixView(const DiagMatrixView<T,S2,C,I2>& m2) :
            itsm(m2.cptr()), itssize(m2.size()), itsstep(m2.step()) {}

        template <int N2, int S2, IndexStyle I2>
        inline ConstSmallDiagMatrixView(
            const ConstSmallDiagMatrixView<T,N2,S2,C,I2>& m2) :
            itsm(m2.cptr()), itssize(m2.size()), itsstep(m2.step()) {}

        template <int N2, int S2, IndexStyle I2>
        inline ConstSmallDiagMatrixView(
            const SmallDiagMatrixView<T,N2,S2,C,I2>& m2) :
            itsm(m2.cptr()), itssize(m2.size()), itsstep(m2.step()) {}

        inline ~ConstSmallDiagMatrixView() { 
#ifdef TMV_DEBUG
            itsm = 0; 
#endif
        }

    private :
        inline void operator=(const type& m2);
    public :

        //
        // Auxilliary Functions
        //

        inline const T* cptr() const { return itsm; }

        inline T cref(int i, int j) const { return (i!=j ? T(0) : cref(i)); }
        inline T cref(int i) const
        { return DoConj<C>(itsm[i*step()]); }

        inline size_t size() const { return itssize; }
        inline int step() const { return itsstep; }
        inline bool isconj() const { return C; }
        inline bool isrm() const { return step()==1; }
        inline bool iscm() const { return step()==1; }

    private :

#ifdef TMV_DEBUG
        const T* itsm;
#else
        const T*const itsm;
#endif
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

        inline ConstSmallDiagMatrixViewF(const T* m, size_t n, int s) : 
            mtype(m,n,s) {}
        inline ConstSmallDiagMatrixViewF(const T* m, size_t n) : mtype(m,n) {}
        inline ConstSmallDiagMatrixViewF(const T* m) : mtype(m) {}
        inline ConstSmallDiagMatrixViewF(const type& m2) : mtype(m2) {}
        template <int S2, IndexStyle I2>
        inline ConstSmallDiagMatrixViewF(
            const ConstSmallDiagMatrixView<T,S2,C,I2>& m2) : mtype(m2) {}
        template <int S2, IndexStyle I2>
        inline ConstSmallDiagMatrixViewF(
            const SmallDiagMatrixView<T,S2,C,I2>& m2
        ) :
            mtype(m2) {}
        template <int N2, int S2, IndexStyle I2>
        inline ConstSmallDiagMatrixViewF(
            const ConstSmallDiagMatrixView<T,N2,S2,C,I2>& m2) : mtype(m2) {}
        template <int N2, int S2, IndexStyle I2>
        inline ConstSmallDiagMatrixViewF(
            const SmallDiagMatrixView<T,N2,S2,C,I2>& m2) :
            mtype(m2) {}
        inline ~ConstSmallDiagMatrixViewF() {}

    private :
        inline void operator=(const type& m2);

    }; // ConstSmallDiagMatrixViewF


    template <class T, int N, int S, bool C, IndexStyle I>
    struct Traits<SmallDiagMatrixView<T,N,S,C,I> >
    {
        typedef T value_type;

        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;
        enum { misreal = Traits<T>::isreal };
        enum { miscomplex = Traits<T>::iscomplex };

        typedef SmallDiagMatrixView<T,N,S,C,I> type;
        typedef const ConstSmallDiagMatrixView<T,N,S,C,I> calc_type;
        typedef calc_type eval_type;
        typedef SmallDiagMatrix<T,N,I> copy_type;

        enum { mcolsize = N };
        enum { mrowsize = N };
        enum { msize = N };
        enum { mfort = (I == FortranStyle) };
        enum { mshape = Diag };
        enum { mrowmajor = (S==1) }; 
        enum { mcolmajor = (S==1) }; 
        enum { mstor = ColMajor }; // arbitrary
        enum { mcalc = true };
        enum { mstep = S };
        enum { mconj = C };
        enum { twoS = misreal ? S : IntTraits<S>::twoS };
        enum { notC = !C && miscomplex };

        typedef ConstSmallVectorView<T,N,S,C,I> const_diag_type;

        typedef ConstDiagMatrixView<T,S,C,I> const_subdiagmatrix_type;
        typedef ConstDiagMatrixView<T,UNKNOWN,C,I> 
            const_subdiagmatrix_step_type;
        typedef ConstSmallDiagMatrixView<T,N,S,C,I> const_view_type;
        typedef ConstSmallDiagMatrixView<T,N,S,C,CStyle> const_cview_type;
        typedef ConstSmallDiagMatrixView<T,N,S,C,FortranStyle> 
            const_fview_type;
        typedef ConstDiagMatrixView<T,UNKNOWN,C> const_xview_type;
        typedef ConstSmallDiagMatrixView<T,N,1,C> const_rmview_type;
        typedef ConstSmallDiagMatrixView<T,N,1,C> const_cmview_type;
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
        typedef SmallDiagMatrixView<T,N,1,C> rmview_type;
        typedef SmallDiagMatrixView<T,N,1,C> cmview_type;
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

        enum { mcolsize = Traits<type>::msize };
        enum { mrowsize = Traits<type>::msize };
        enum { msize = Traits<type>::msize };
        enum { mfort = Traits<type>::mfort };
        enum { mshape = Traits<type>::mshape };
        enum { mrowmajor = Traits<type>::mrowmajor };
        enum { mcolmajor = Traits<type>::mcolmajor };
        enum { mcalc = Traits<type>::mcalc };
        enum { mstep = Traits<type>::mstep };
        enum { mconj = Traits<type>::mconj };

        //
        // Constructors
        //

        inline SmallDiagMatrixView(T* m, size_t n, int s) :
            itsm(m), itssize(n), itsstep(s) {}

        inline SmallDiagMatrixView(T* m, size_t n) :
            itsm(m), itssize(n), itsstep(S)
        { TMVStaticAssert(S != UNKNOWN); }

        inline SmallDiagMatrixView(T* m) : itsm(m), itssize(N), itsstep(S)
        { TMVStaticAssert(N != UNKNOWN); TMVStaticAssert(S != UNKNOWN); }

        inline SmallDiagMatrixView(const type& m2) :
            itsm(m2.itsm), itssize(m2.size()), itsstep(m2.step()) {}

        template <int S2, IndexStyle I2>
        inline SmallDiagMatrixView(DiagMatrixView<T,S2,C,I2> m2) :
            itsm(m2.ptr()), itssize(m2.size()), itsstep(m2.step()) {}

        template <int N2, int S2, IndexStyle I2>
        inline SmallDiagMatrixView(SmallDiagMatrixView<T,N2,S2,C,I2> m2) :
            itsm(m2.ptr()), itssize(m2.size()), itsstep(m2.step()) {}

        inline ~SmallDiagMatrixView() { 
#ifdef TMV_DEBUG
            itsm = 0; 
#endif
        }


        //
        // Op = 
        //

        template <class M2>
        inline type& operator=(const BaseMatrix<M2>& m2)
        {
            base_mut::operator=(m2); 
            return *this;
        }

        inline type& operator=(const type& m2)
        {
            base_mut::operator=(m2); 
            return *this;
        }

        inline type& operator=(const T x)
        {
            base_mut::operator=(x); 
            return *this;
        }


        //
        // Auxilliary Functions
        //

        inline const T* cptr() const { return itsm; }
        inline T* ptr() { return itsm; }

        inline T cref(int i, int j) const { return (i!=j ? T(0) : cref(i)); }
        inline T cref(int i) const
        { return DoConj<C>(itsm[i*step()]); }

        inline reference ref(int i) 
        { return reference(itsm[i*step()]); }

        inline size_t size() const { return itssize; }
        inline int step() const { return itsstep; }
        inline bool isconj() const { return C; }
        inline bool isrm() const { return step()==1; }
        inline bool iscm() const { return step()==1; }

    private :

#ifdef TMV_DEBUG
        T* itsm;
#else
        T*const itsm;
#endif
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

        inline SmallDiagMatrixViewF(T* m, size_t n, int s) : mtype(m,n,s) {}
        inline SmallDiagMatrixViewF(T* m, size_t n) : mtype(m,n) {}
        inline SmallDiagMatrixViewF(T* m) : mtype(m) {}

        inline SmallDiagMatrixViewF(const type& m2) : mtype(m2) {}
        template <int S2, IndexStyle I2>
        inline SmallDiagMatrixViewF(DiagMatrixView<T,S2,C,I2> m2) :
            mtype(m2) {}
        template <int N2, int S2, IndexStyle I2>
        inline SmallDiagMatrixViewF(SmallMatrixView<T,N2,S2,C,I2> m2) :
            mtype(m2) {}
        inline ~SmallDiagMatrixViewF() {}

        template <class M2>
        inline type& operator=(const BaseMatrix<M2>& m2)
        { mtype::operator=(m2); return *this; }
        inline type& operator=(const type& m2)
        { mtype::operator=(m2); return *this; }
        inline type& operator=(const T x)
        { mtype::operator=(x); return *this; }

    }; // SmallDiagMatrixViewF


    //
    // DiagMatrixViewOf
    //

    template <class T, int N, int S, bool C, IndexStyle I>
    inline ConstSmallDiagMatrixView<T,N,S,C,I> DiagMatrixViewOf(
        const ConstSmallVectorView<T,N,S,C,I>& v)
    { return ConstSmallDiagMatrixView<T,N,S,C,I>(v.cptr(),v.size(),v.step()); }

    template <class T, int N, int S, bool C, IndexStyle I>
    inline SmallDiagMatrixView<T,N,S,C,I> DiagMatrixViewOf(
        SmallVectorView<T,N,S,C,I> v)
    { return SmallDiagMatrixView<T,N,S,C,I>(v.ptr(),v.size(),v.step()); }

    template <class T, int N, IndexStyle I>
    inline SmallDiagMatrixView<T,N,1,false,I> DiagMatrixViewOf(
        SmallVector<T,N,I>& v)
    { return SmallDiagMatrixView<T,N,1,false,I>(v.ptr()); }


    //
    // TMV_Text 
    //

    template <class T, int N, IndexStyle I>
    inline std::string TMV_Text(const SmallDiagMatrix<T,N,I>& )
    {
        std::ostringstream s;
        s << "SmallDiagMatrix<"<<TMV_Text(T())<<","<<N<<','<<TMV_Text(I)<<">";
        return s.str();
    }

    template <class T, int N, int S, bool C, IndexStyle I>
    inline std::string TMV_Text(const ConstSmallDiagMatrixView<T,N,S,C,I>& m)
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
    inline std::string TMV_Text(const SmallDiagMatrixView<T,N,S,C,I>& m)
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
