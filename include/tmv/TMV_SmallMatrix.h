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
// This file defines the TMV SmallMatrix class.
//
// A SmallMatrix is a matrix whose size is known at compile time.
// This allows for a lot of optimizations by the compiler in implementing
// the various operations on it. 
// 
// Since the instantiation needs to know the size of the matrix, all 
// operations on the SmallMatrix are done inline, rather than precompiled.
// This gives another speed boost in most cases, since the compiler can
// avoid implementing a function call in most cases.
//
// Finally, another advantage is that we allocate the data on the stack
// rather than the heap, so we avoid new and delete calls as well.
// However, stack sizes are usually limited to hundreds of KB, 
// (mine is 8 MB), so we set a maximum size of 1KB for each SmallMatrix.
// (For double, this means up to 11x11 will be allocated on the stack.)
// Any bigger than that, and the performance drop from using the
// heap is pretty irrelevant.  (See TMV_Array.h to change this.)
// 
// One drawback of using a SmallMatrix is that it does not do any
// alias checking in the aritmetic statements.  So a statement like
// m = m2 * m will not produce a correct answer. 
// Normally this is a feature, since the alias checks can be a 
// significant fraction of the calculation time for small matrices.
// 
// You can workaround this when necessary by explicitly making a copy.
// The easiest way is with the .copy() method.  e.g. m = m2 * m.copy().
// 
//
// Constructors:
//
//    SmallMatrix<T,M,N,stor,I>()
//        Makes a SmallMatrix with column size = M and row size = N
//        with _uninitialized_ values
//
//    SmallMatrix<T,M,N,stor,I>(T x)
//        Makes a SmallMatrix of size n with all values = x
//
//    SmallMatrix<T,M,N,stor,I>(const vector<vector<T> >& m)
//        Makes a SmallMatrix with a_ij = m[i][j]
//
//    SmallMatrix<T,M,N,stor,I>(const T* m)
//    SmallMatrix<T,M,N,stor,I>(const vector<T>& m)
//    SmallMatrix<T,M,N,stor,I>(const GenMatrix<T>& m)
//        Make a SmallMatrix which copies the elements of m.
//
// Special Creators:
//
//    SmallMatrixView RowVectorViewOf(Vector& v)
//    SmallMatrixView RowVectorViewOf(VectorView v)
//    ConstSmallMatrixView RowVectorViewOf(const Vector& v)
//    ConstSmallMatrixView RowVectorViewOf(const ConstVectorView v)
//        Returns a 1xn MatrixView with v in the (only) row.
//        Unlike creating a Matrix with RowVector(v), this uses the same
//        data storage as the actual Vector v.
//        For a const Vector or a ConstVectorView, this returns a 
//        ConstMatrixView.
//
//    SmallMatrixView ColVectorViewOf(Vector& v)
//    SmallMatrixView ColVectorViewOf(VectorView v)
//    ConstSmallMatrixView ColVectorViewOf(const Vector& v)
//    ConstSmallMatrixView ColVectorViewOf(const ConstVectorView& v)
//        Returns an nx1 MatrixView with v in the (only) column.


#ifndef TMV_SmallMatrix_H
#define TMV_SmallMatrix_H

#include <vector>
#include "TMV_BaseMatrix_Rec.h"
#include "TMV_BaseMatrix_Tri.h"
#include "TMV_VIt.h"
#include "TMV_Array.h"
#include "TMV_Matrix.h"

namespace tmv {

    template <class M> class LUD;
    template <class T> class InstLUD;
    template <class M> class QRD;
    template <class T> class InstQRD;
    template <class M> class QRPD;
    template <class T> class InstQRPD;
    template <class M> class SVD;
    template <class T> class InstSVD;

    //
    // SmallMatrix
    //

    template <class T, int M, int N, int A0, int A1>
    struct Traits<SmallMatrix<T,M,N,A0,A1> >
    {
        enum { A01 = A0 | A1 };
        enum { A = (A01 & ~NoDivider) | (
                (Attrib<A01>::rowmajor ? 0 : ColMajor) )};
        enum { okA = (
                !Attrib<A>::conj &&
                (Attrib<A>::rowmajor || Attrib<A>::colmajor) &&
                (Attrib<A>::rowmajor != int(Attrib<A>::colmajor)) &&
                !Attrib<A>::diagmajor &&
                !Attrib<A>::nonunitdiag &&
                !Attrib<A>::unitdiag &&
                !Attrib<A>::zerodiag &&
                !Attrib<A>::packed &&
                !Attrib<A>::lower &&
                !Attrib<A>::withdivider )};
        enum { _attrib = A };

        typedef T value_type;
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;
        enum { isreal = Traits<T>::isreal };
        enum { iscomplex = Traits<T>::iscomplex };

        typedef SmallMatrix<T,M,N,A0,A1> type;
        typedef const type& calc_type;
        typedef const type& eval_type;
        typedef SmallMatrix<T,M,N,A01> copy_type;

        enum { _colsize = M };
        enum { _rowsize = N };
        enum { _shape = Rec };
        enum { _fort = Attrib<A>::fort };
        enum { _calc = true };
        enum { _rowmajor = Attrib<A>::rowmajor };
        enum { _colmajor = Attrib<A>::colmajor };
        enum { _stepi = _rowmajor ? N : 1 };
        enum { _stepj = _rowmajor ? 1 : M };
        enum { _diagstep = _stepi + _stepj };
        enum { _conj = false };
        enum { _checkalias = Attrib<A>::checkalias };
        enum { _canlin = true };
        enum { _linsize = IntTraits2<M,N>::prod };
        enum { twoSi = isreal ? int(_stepi) : int(IntTraits<_stepi>::twoS) };
        enum { twoSj = isreal ? int(_stepj) : int(IntTraits<_stepj>::twoS) };
        enum { minMN = IntTraits2<M,N>::min };

        enum { _hasdivider = false };
        typedef QuotXM<1,real_type,type> inverse_type;

        typedef typename TypeSelect<_hasdivider ,
                const InstLUD<value_type>& , LUD<copy_type> >::type lud_type;
        typedef typename TypeSelect<_hasdivider ,
                const InstQRD<value_type>& , QRD<copy_type> >::type qrd_type;
        typedef typename TypeSelect<_hasdivider ,
                const InstQRPD<value_type>& , QRPD<copy_type> >::type qrpd_type;
        typedef InvalidType svd_type;

        enum { vecA = (
                (_fort ? FortranStyle : 0) |
                (_checkalias ? CheckAlias : 0) )};
        enum { colA = vecA | (_colmajor ? Unit : 0) };
        enum { rowA = vecA | (_rowmajor ? Unit : 0) };
        enum { diagA = vecA | (_diagstep == 1 ? Unit : 0) };
        enum { cstyleA = A & ~FortranStyle };
        enum { fstyleA = A | FortranStyle };
        enum { nmA = (A & ~AllStorageType) };
        enum { cmA = nmA | ColMajor };
        enum { rmA = nmA | RowMajor };
        enum { colpairA = A & ~RowMajor };
        enum { rowpairA = A & ~ColMajor };
        enum { conjA = iscomplex ? (A ^ Conj) : int(A) };
        enum { trA = _colmajor ? rmA : _rowmajor ? cmA : int(A) };
        enum { adjA = iscomplex ? (trA ^ Conj) : int(trA) };
        enum { nonconjA = A };
        enum { twosA = isreal ? nonconjA : int(nonconjA & ~AllStorageType) };
        enum { An = _checkalias ? (A & ~CheckAlias) : (A | NoAlias) };
        enum { vecAn = _checkalias ? (vecA & ~CheckAlias) : (vecA | NoAlias) };
        enum { colAn = vecAn | (_colmajor ? Unit : 0) };
        enum { rowAn = vecAn | (_rowmajor ? Unit : 0) };
        enum { diagAn = vecAn | (_diagstep == 1 ? Unit : 0) };
        enum { nmAn = (An & ~AllStorageType) };

        typedef ConstSmallVectorView<T,M,_stepi,colA> const_col_type;
        typedef ConstVectorView<T,colAn> const_col_sub_type;
        typedef ConstSmallVectorView<T,N,_stepj,rowA> const_row_type;
        typedef ConstVectorView<T,rowAn> const_row_sub_type;
        typedef ConstSmallVectorView<T,minMN,_diagstep,diagA> const_diag_type;
        typedef ConstVectorView<T,diagAn> const_diag_sub_type;

        typedef ConstMatrixView<T,An> const_submatrix_type;
        typedef ConstMatrixView<T,nmAn> const_submatrix_step_type;
        typedef ConstVectorView<T,vecAn> const_subvector_type;
        typedef ConstSmallMatrixView<T,M,2,_stepi,UNKNOWN,colpairA>
            const_colpair_type;
        typedef ConstSmallMatrixView<T,2,N,UNKNOWN,_stepj,rowpairA>
            const_rowpair_type;
        typedef ConstSmallMatrixView<T,M,UNKNOWN,_stepi,_stepj,A>
            const_colrange_type;
        typedef ConstSmallMatrixView<T,UNKNOWN,N,_stepi,_stepj,A>
            const_rowrange_type;
        typedef ConstSmallMatrixView<T,M,N,_stepi,_stepj,A>
            const_view_type;
        typedef ConstSmallMatrixView<T,M,N,_stepi,_stepj,cstyleA>
            const_cview_type;
        typedef ConstSmallMatrixView<T,M,N,_stepi,_stepj,fstyleA>
            const_fview_type;
        typedef ConstMatrixView<T> const_xview_type;
        typedef ConstSmallMatrixView<T,M,N,1,_stepj,cmA> const_cmview_type;
        typedef ConstSmallMatrixView<T,M,N,_stepi,1,rmA> const_rmview_type;
        typedef ConstSmallMatrixView<T,M,N,_stepi,_stepj,conjA>
            const_conjugate_type;
        typedef ConstSmallMatrixView<T,N,M,_stepj,_stepi,trA>
            const_transpose_type;
        typedef ConstSmallMatrixView<T,N,M,_stepj,_stepi,adjA>
            const_adjoint_type;
        typedef ConstSmallUpperTriMatrixView<T,N,_stepi,_stepj,(A|NonUnitDiag)> 
            const_uppertri_type;
        typedef ConstSmallUpperTriMatrixView<T,N,_stepi,_stepj,(A|UnitDiag)>
            const_unit_uppertri_type;
        typedef ConstSmallUpperTriMatrixView<T,N,_stepi,_stepj,(A|UnknownDiag)>
            const_unknown_uppertri_type;
        typedef ConstSmallLowerTriMatrixView<T,M,_stepi,_stepj,(A|NonUnitDiag)>
            const_lowertri_type;
        typedef ConstSmallLowerTriMatrixView<T,M,_stepi,_stepj,(A|UnitDiag)>
            const_unit_lowertri_type;
        typedef ConstSmallLowerTriMatrixView<T,N,_stepi,_stepj,(A|UnknownDiag)>
            const_unknown_lowertri_type;
        typedef ConstSmallMatrixView<real_type,M,N,twoSi,twoSj,twosA>
            const_realpart_type;
        typedef const_realpart_type const_imagpart_type;
        typedef ConstSmallVectorView<T,_linsize,1,(vecA|Unit)> 
            const_linearview_type;
        typedef ConstSmallMatrixView<T,M,N,_stepi,_stepj,nonconjA>
            const_nonconj_type;
        typedef SmallMatrixView<T,M,N,_stepi,_stepj,A> nonconst_type;

        typedef T& reference;
        typedef VIt<T,1,false> linear_iterator;

        typedef SmallVectorView<T,M,_stepi,colA> col_type;
        typedef VectorView<T,colAn> col_sub_type;
        typedef SmallVectorView<T,N,_stepj,rowA> row_type;
        typedef VectorView<T,rowAn> row_sub_type;
        typedef SmallVectorView<T,minMN,_diagstep,diagA> diag_type;
        typedef VectorView<T,diagAn> diag_sub_type;

        typedef MatrixView<T,An> submatrix_type;
        typedef MatrixView<T,nmAn> submatrix_step_type;
        typedef VectorView<T,vecAn> subvector_type;
        typedef SmallMatrixView<T,M,2,_stepi,UNKNOWN,colpairA> colpair_type;
        typedef SmallMatrixView<T,2,N,UNKNOWN,_stepj,rowpairA> rowpair_type;
        typedef SmallMatrixView<T,M,UNKNOWN,_stepi,_stepj,A> colrange_type;
        typedef SmallMatrixView<T,UNKNOWN,N,_stepi,_stepj,A> rowrange_type;
        typedef SmallMatrixView<T,M,N,_stepi,_stepj,A> view_type;
        typedef SmallMatrixView<T,M,N,_stepi,_stepj,cstyleA> cview_type;
        typedef SmallMatrixView<T,M,N,_stepi,_stepj,fstyleA> fview_type;
        typedef MatrixView<T> xview_type;
        typedef SmallMatrixView<T,M,N,1,_stepj,cmA> cmview_type;
        typedef SmallMatrixView<T,M,N,_stepi,1,rmA> rmview_type;
        typedef SmallMatrixView<T,M,N,_stepi,_stepj,conjA> conjugate_type;
        typedef SmallMatrixView<T,N,M,_stepj,_stepi,trA> transpose_type;
        typedef SmallMatrixView<T,N,M,_stepj,_stepi,adjA> adjoint_type;
        typedef SmallUpperTriMatrixView<T,N,_stepi,_stepj,(A|NonUnitDiag)> 
            uppertri_type;
        typedef SmallUpperTriMatrixView<T,N,_stepi,_stepj,(A|UnitDiag)>
            unit_uppertri_type;
        typedef SmallUpperTriMatrixView<T,N,_stepi,_stepj,(A|UnknownDiag)>
            unknown_uppertri_type;
        typedef SmallLowerTriMatrixView<T,M,_stepi,_stepj,(A|NonUnitDiag)>
            lowertri_type;
        typedef SmallLowerTriMatrixView<T,M,_stepi,_stepj,(A|UnitDiag)>
            unit_lowertri_type;
        typedef SmallLowerTriMatrixView<T,N,_stepi,_stepj,(A|UnknownDiag)>
            unknown_lowertri_type;
        typedef SmallMatrixView<real_type,M,N,twoSi,twoSj,twosA> realpart_type;
        typedef realpart_type imagpart_type;
        typedef SmallVectorView<T,_linsize,1,(vecA|Unit)> linearview_type;
        typedef SmallMatrixView<T,M,N,_stepi,_stepj,nonconjA> nonconj_type;
    };

    template <class T, int M, int N, int A0, int A1>
    class SmallMatrix : 
        public BaseMatrix_Rec_Mutable<SmallMatrix<T,M,N,A0,A1> >,
        public MatrixDivHelper<SmallMatrix<T,M,N,A0,A1> >
    {
    public:

        typedef SmallMatrix<T,M,N,A0,A1> type;
        typedef BaseMatrix_Rec_Mutable<type> base_mut;
        typedef typename Traits<type>::lud_type lud_type;
        typedef typename Traits<type>::qrd_type qrd_type;
        typedef typename Traits<type>::qrpd_type qrpd_type;
        typedef typename Traits<type>::svd_type svd_type;

        typedef typename Traits<T>::real_type real_type;

        enum { _colsize = Traits<type>::_colsize };
        enum { _rowsize = Traits<type>::_rowsize };
        enum { _shape = Traits<type>::_shape };
        enum { _fort = Traits<type>::_fort };
        enum { _rowmajor = Traits<type>::_rowmajor };
        enum { _colmajor = Traits<type>::_colmajor };
        enum { _calc = Traits<type>::_calc };
        enum { _stepi = Traits<type>::_stepi };
        enum { _stepj = Traits<type>::_stepj };
        enum { _diagstep = Traits<type>::_diagstep };
        enum { _conj = Traits<type>::_conj };
        enum { _checkalias = Traits<type>::_checkalias };
        enum { _canlin = Traits<type>::_canlin };
        enum { _linsize = Traits<type>::_linsize };
        enum { _attrib = Traits<type>::_attrib };

        //
        // Constructors
        //

        SmallMatrix() 
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(M>=0);
            TMVStaticAssert(N>=0);
#ifdef TMV_DEBUG
            this->linearView().flatten().setAllTo(
                Traits<real_type>::constr_value());
#endif
        }

        explicit SmallMatrix(T x) 
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(M>=0);
            TMVStaticAssert(N>=0);
            this->setAllTo(x);
        }

        explicit SmallMatrix(const T* vv) 
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(M>=0);
            TMVStaticAssert(N>=0);
            typename type::linearview_type lv = this->linearView();
            ConstSmallVectorView<T,_linsize,1>(vv).newAssignTo(lv);
        }

        explicit SmallMatrix(const std::vector<T>& vv) 
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(M>=0);
            TMVStaticAssert(N>=0);
            TMVAssert(vv.size() == _linsize);
            typename type::linearview_type lv = this->linearView();
            ConstSmallVectorView<T,_linsize,1>(&vv[0]).newAssignTo(lv);
        }

        explicit SmallMatrix(const std::vector<std::vector<T> >& vv) 
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(M>=0);
            TMVStaticAssert(N>=0);
            TMVAssert(vv.size() == M);
            for(int i=0;i<M;++i) {
                TMVAssert(vv[i].size() == N);
                typename type::row_type mi = this->row(i);
                ConstSmallVectorView<T,N,1>(&vv[i][0]).newAssignTo(mi);
            }
        }

        SmallMatrix(const type& m2) 
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(M>=0);
            TMVStaticAssert(N>=0);
            m2.newAssignTo(*this);
        }

        template <class M2>
        SmallMatrix(const BaseMatrix<M2>& m2) 
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(M>=0);
            TMVStaticAssert(N>=0);
            TMVStaticAssert((ShapeTraits2<M2::_shape,_shape>::assignable));
            TMVAssert(m2.colsize() == M);
            TMVAssert(m2.rowsize() == N);
            m2.newAssignTo(*this);
        }

        ~SmallMatrix()
        {
#ifdef TMV_DEBUG
            this->linearView().flatten().setAllTo(
                Traits<real_type>::destr_value());
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

        T cref(int i, int j) const
        { return itsm[i*stepi()+j*stepj()]; }

        T& ref(int i, int j)
        { return itsm[i*stepi()+j*stepj()]; }

        TMV_INLINE size_t ls() const { return _linsize; }
        TMV_INLINE size_t colsize() const { return M; }
        TMV_INLINE size_t rowsize() const { return N; }
        TMV_INLINE int nElements() const { return _linsize; }
        TMV_INLINE int stepi() const { return _stepi; }
        TMV_INLINE int stepj() const { return _stepj; }
        TMV_INLINE bool isconj() const { return false; }
        TMV_INLINE bool isrm() const { return _rowmajor; }
        TMV_INLINE bool iscm() const { return _colmajor; }

    private:

        StackArray<T,_linsize> itsm;

    }; // SmallMatrix

    template <class T, int M, int N, int Si, int Sj, int A0>
    struct Traits<ConstSmallMatrixView<T,M,N,Si,Sj,A0> >
    {
        enum { A = (A0 & ~NoDivider) | (
                ( Attrib<A0>::colmajor ? 0 :
                  Attrib<A0>::rowmajor ? 0 :
                  Si == 1 ? ColMajor :
                  Sj == 1 ? RowMajor : 0 ) )};
        enum { okA = (
                !Attrib<A>::diagmajor &&
                !Attrib<A>::nonunitdiag &&
                !Attrib<A>::unitdiag &&
                !Attrib<A>::zerodiag &&
                !Attrib<A>::packed &&
                !Attrib<A>::lower &&
                !Attrib<A>::withdivider &&
                ( Traits<T>::iscomplex || !Attrib<A>::conj ) &&
                !( Si != UNKNOWN && Si != 1 && Attrib<A>::colmajor ) &&
                !( Sj != UNKNOWN && Sj != 1 && Attrib<A>::rowmajor ) )};
        enum { _attrib = A };

        typedef T value_type;
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;
        enum { isreal = Traits<T>::isreal };
        enum { iscomplex = Traits<T>::iscomplex };

        typedef ConstSmallMatrixView<T,M,N,Si,Sj,A0> type;
        typedef const type& calc_type;
        typedef const type& eval_type;

        enum { _colsize = M };
        enum { _rowsize = N };
        enum { _shape = Rec };
        enum { _fort = Attrib<A>::fort };
        enum { _calc = true };
        enum { _rowmajor = Attrib<A>::rowmajor };
        enum { _colmajor = Attrib<A>::colmajor };
        enum { _stepi = Si };
        enum { _stepj = Sj };
        enum { _diagstep = IntTraits2<Si,Sj>::sum };
        enum { _conj = Attrib<A>::conj };
        enum { _checkalias = Attrib<A>::checkalias };
        enum { _canlin = (
                (Si == 1 && M != UNKNOWN && Sj == M) || 
                (Sj == 1 && N != UNKNOWN && Si == N) )};
        enum { _linsize = IntTraits2<M,N>::prod };
        enum { twoSi = isreal ? int(_stepi) : int(IntTraits<_stepi>::twoS) };
        enum { twoSj = isreal ? int(_stepj) : int(IntTraits<_stepj>::twoS) };
        enum { minMN = IntTraits2<M,N>::min };

        // Use MCopyHelper for copy_type in case M or N == UNKNOWN
        typedef typename MCopyHelper<T,Rec,M,N,_rowmajor,_fort>::type copy_type;

        enum { _hasdivider = false };
        typedef QuotXM<1,real_type,type> inverse_type;

        typedef typename TypeSelect<_hasdivider ,
                const InstLUD<value_type>& , LUD<copy_type> >::type lud_type;
        typedef typename TypeSelect<_hasdivider ,
                const InstQRD<value_type>& , QRD<copy_type> >::type qrd_type;
        typedef typename TypeSelect<_hasdivider ,
                const InstQRPD<value_type>& , QRPD<copy_type> >::type qrpd_type;
        typedef InvalidType svd_type;

        enum { vecA = (
                (_fort ? FortranStyle : 0) |
                (_conj ? Conj : 0 ) |
                (_checkalias ? CheckAlias : 0) )};
        enum { colA = vecA | (_colmajor ? Unit : 0) };
        enum { rowA = vecA | (_rowmajor ? Unit : 0) };
        enum { diagA = vecA | (_diagstep == 1 ? Unit : 0) };
        enum { cstyleA = A & ~FortranStyle };
        enum { fstyleA = A | FortranStyle };
        enum { nmA = (A & ~AllStorageType) };
        enum { cmA = nmA | ColMajor };
        enum { rmA = nmA | RowMajor };
        enum { colpairA = A & ~RowMajor };
        enum { rowpairA = A & ~ColMajor };
        enum { conjA = iscomplex ? (A ^ Conj) : int(A) };
        enum { trA = _colmajor ? rmA : _rowmajor ? cmA : int(A) };
        enum { adjA = iscomplex ? (trA ^ Conj) : int(trA) };
        enum { nonconjA = A & ~Conj };
        enum { twosA = isreal ? nonconjA : int(nonconjA & ~AllStorageType) };
        enum { An = _checkalias ? (A & ~CheckAlias) : (A | NoAlias) };
        enum { vecAn = _checkalias ? (vecA & ~CheckAlias) : (vecA | NoAlias) };
        enum { colAn = vecAn | (_colmajor ? Unit : 0) };
        enum { rowAn = vecAn | (_rowmajor ? Unit : 0) };
        enum { diagAn = vecAn | (_diagstep == 1 ? Unit : 0) };
        enum { nmAn = (An & ~AllStorageType) };

        typedef ConstSmallVectorView<T,M,_stepi,colA> const_col_type;
        typedef ConstVectorView<T,colAn> const_col_sub_type;
        typedef ConstSmallVectorView<T,N,_stepj,rowA> const_row_type;
        typedef ConstVectorView<T,rowAn> const_row_sub_type;
        typedef ConstSmallVectorView<T,minMN,_diagstep,diagA> const_diag_type;
        typedef ConstVectorView<T,diagAn> const_diag_sub_type;

        typedef ConstMatrixView<T,An> const_submatrix_type;
        typedef ConstMatrixView<T,nmAn> const_submatrix_step_type;
        typedef ConstVectorView<T,vecAn> const_subvector_type;
        typedef ConstSmallMatrixView<T,M,2,_stepi,UNKNOWN,colpairA>
            const_colpair_type;
        typedef ConstSmallMatrixView<T,2,N,UNKNOWN,_stepj,rowpairA>
            const_rowpair_type;
        typedef ConstSmallMatrixView<T,M,UNKNOWN,_stepi,_stepj,A>
            const_colrange_type;
        typedef ConstSmallMatrixView<T,UNKNOWN,N,_stepi,_stepj,A>
            const_rowrange_type;
        typedef ConstSmallMatrixView<T,M,N,_stepi,_stepj,A>
            const_view_type;
        typedef ConstSmallMatrixView<T,M,N,_stepi,_stepj,cstyleA>
            const_cview_type;
        typedef ConstSmallMatrixView<T,M,N,_stepi,_stepj,fstyleA>
            const_fview_type;
        typedef ConstMatrixView<T,(_conj ? Conj : NonConj)> const_xview_type;
        typedef ConstSmallMatrixView<T,M,N,1,_stepj,cmA> const_cmview_type;
        typedef ConstSmallMatrixView<T,M,N,_stepi,1,rmA> const_rmview_type;
        typedef ConstSmallMatrixView<T,M,N,_stepi,_stepj,conjA>
            const_conjugate_type;
        typedef ConstSmallMatrixView<T,N,M,_stepj,_stepi,trA>
            const_transpose_type;
        typedef ConstSmallMatrixView<T,N,M,_stepj,_stepi,adjA>
            const_adjoint_type;
        typedef ConstSmallUpperTriMatrixView<T,N,_stepi,_stepj,(A|NonUnitDiag)> 
            const_uppertri_type;
        typedef ConstSmallUpperTriMatrixView<T,N,_stepi,_stepj,(A|UnitDiag)>
            const_unit_uppertri_type;
        typedef ConstSmallUpperTriMatrixView<T,N,_stepi,_stepj,(A|UnknownDiag)>
            const_unknown_uppertri_type;
        typedef ConstSmallLowerTriMatrixView<T,M,_stepi,_stepj,(A|NonUnitDiag)>
            const_lowertri_type;
        typedef ConstSmallLowerTriMatrixView<T,M,_stepi,_stepj,(A|UnitDiag)>
            const_unit_lowertri_type;
        typedef ConstSmallLowerTriMatrixView<T,N,_stepi,_stepj,(A|UnknownDiag)>
            const_unknown_lowertri_type;
        typedef ConstSmallMatrixView<real_type,M,N,twoSi,twoSj,twosA>
            const_realpart_type;
        typedef const_realpart_type const_imagpart_type;
        typedef ConstSmallVectorView<T,_linsize,1,(vecA|Unit)> 
            const_linearview_type;
        typedef ConstSmallMatrixView<T,M,N,_stepi,_stepj,nonconjA>
            const_nonconj_type;
        typedef SmallMatrixView<T,M,N,_stepi,_stepj,A> nonconst_type;
    };

    template <class T, int M, int N, int Si, int Sj, int A>
    class ConstSmallMatrixView :
        public BaseMatrix_Rec<ConstSmallMatrixView<T,M,N,Si,Sj,A> >,
        public MatrixDivHelper<ConstSmallMatrixView<T,M,N,Si,Sj,A> >
    {
    public:

        typedef ConstSmallMatrixView<T,M,N,Si,Sj,A> type;
        typedef typename Traits<type>::lud_type lud_type;
        typedef typename Traits<type>::qrd_type qrd_type;
        typedef typename Traits<type>::qrpd_type qrpd_type;
        typedef typename Traits<type>::svd_type svd_type;

        enum { _colsize = Traits<type>::_colsize };
        enum { _rowsize = Traits<type>::_rowsize };
        enum { _shape = Traits<type>::_shape };
        enum { _fort = Traits<type>::_fort };
        enum { _rowmajor = Traits<type>::_rowmajor };
        enum { _colmajor = Traits<type>::_colmajor };
        enum { _calc = Traits<type>::_calc };
        enum { _stepi = Traits<type>::_stepi };
        enum { _stepj = Traits<type>::_stepj };
        enum { _diagstep = Traits<type>::_diagstep };
        enum { _conj = Traits<type>::_conj };
        enum { _checkalias = Traits<type>::_checkalias };
        enum { _canlin = Traits<type>::_canlin };
        enum { _linsize = Traits<type>::_linsize };
        enum { _attrib = Traits<type>::_attrib };

        //
        // Constructors
        //

        TMV_INLINE ConstSmallMatrixView(
            const T* m, size_t cs, size_t rs, int si, int sj) :
            itsm(m), itscs(cs), itsrs(rs), itssi(si), itssj(sj) 
        {
            TMVStaticAssert(Traits<type>::okA);
        }

        TMV_INLINE ConstSmallMatrixView(
            const T* m, size_t cs, size_t rs, int si) :
            itsm(m), itscs(cs), itsrs(rs), itssi(si), itssj(Sj)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Sj != UNKNOWN); 
        }

        TMV_INLINE ConstSmallMatrixView(const T* m, size_t cs, size_t rs) :
            itsm(m), itscs(cs), itsrs(rs), itssi(Si), itssj(Sj)
        { TMVStaticAssert(Si != UNKNOWN); TMVStaticAssert(Sj != UNKNOWN); }

        TMV_INLINE ConstSmallMatrixView(const T* m, size_t cs) :
            itsm(m), itscs(cs), itsrs(N), itssi(Si), itssj(Sj)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(N != UNKNOWN);
            TMVStaticAssert(Si != UNKNOWN); TMVStaticAssert(Sj != UNKNOWN); 
        }

        TMV_INLINE ConstSmallMatrixView(const T* m) :
            itsm(m), itscs(M), itsrs(N), itssi(Si), itssj(Sj)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(M != UNKNOWN); TMVStaticAssert(N != UNKNOWN);
            TMVStaticAssert(Si != UNKNOWN); TMVStaticAssert(Sj != UNKNOWN); 
        }

        TMV_INLINE ConstSmallMatrixView(const type& m2) :
            itsm(m2.cptr()), itscs(m2.colsize()), itsrs(m2.rowsize()),
            itssi(m2.stepi()), itssj(m2.stepj()) 
        {
            TMVStaticAssert(Traits<type>::okA);
        }

        template <int A2>
        TMV_INLINE ConstSmallMatrixView(const ConstMatrixView<T,A2>& m2) :
            itsm(m2.cptr()), itscs(m2.colsize()), itsrs(m2.rowsize()),
            itssi(m2.stepi()), itssj(m2.stepj()) 
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Attrib<A>::conj == int(Attrib<A2>::conj)); 
        }

        template <int A2>
        TMV_INLINE ConstSmallMatrixView(const MatrixView<T,A2>& m2) :
            itsm(m2.cptr()), itscs(m2.colsize()), itsrs(m2.rowsize()),
            itssi(m2.stepi()), itssj(m2.stepj()) 
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Attrib<A>::conj == int(Attrib<A2>::conj)); 
        }

        template <int M2, int N2, int Si2, int Sj2, int A2>
        TMV_INLINE ConstSmallMatrixView(
            const ConstSmallMatrixView<T,M2,N2,Si2,Sj2,A2>& m2) :
            itsm(m2.cptr()), itscs(m2.colsize()), itsrs(m2.rowsize()),
            itssi(m2.stepi()), itssj(m2.stepj()) 
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Attrib<A>::conj == int(Attrib<A2>::conj)); 
        }

        template <int M2, int N2, int Si2, int Sj2, int A2>
        TMV_INLINE ConstSmallMatrixView(
            const SmallMatrixView<T,M2,N2,Si2,Sj2,A2>& m2) :
            itsm(m2.cptr()), itscs(m2.colsize()), itsrs(m2.rowsize()),
            itssi(m2.stepi()), itssj(m2.stepj()) 
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Attrib<A>::conj == int(Attrib<A2>::conj)); 
        }

        TMV_INLINE ~ConstSmallMatrixView() {
#ifdef TMV_DEBUG
            itsm = 0; 
#endif
        }

    private:
        void operator=(const type& v2);
    public:


        // 
        // Auxilliary Functions
        //

        TMV_INLINE const T* cptr() const { return itsm; }

        T cref(int i, int j) const
        { return DoConj<_conj>(itsm[i*stepi()+j*stepj()]); }

        TMV_INLINE size_t ls() const { return int(itscs)*int(itsrs); }
        TMV_INLINE size_t colsize() const { return itscs; }
        TMV_INLINE size_t rowsize() const { return itsrs; }
        TMV_INLINE int nElements() const { return int(itscs)*int(itsrs); }
        TMV_INLINE int stepi() const { return itssi; }
        TMV_INLINE int stepj() const { return itssj; }
        TMV_INLINE bool isconj() const { return _conj; }
        TMV_INLINE bool isrm() const 
        { return _rowmajor || (!_colmajor &&  stepj() == 1); }
        TMV_INLINE bool iscm() const 
        { return _colmajor || (!_rowmajor &&  stepi() == 1); }

    private :

        const T* itsm;
        const CheckedInt<M> itscs;
        const CheckedInt<N> itsrs;
        const CheckedInt<Si> itssi;
        const CheckedInt<Sj> itssj;

    }; // ConstSmallMatrixView

    template <class T, int M, int N, int Si, int Sj, int A0>
    struct Traits<SmallMatrixView<T,M,N,Si,Sj,A0> >
    {
        enum { A = (A0 & ~NoDivider) | (
                ( Attrib<A0>::colmajor ? 0 :
                  Attrib<A0>::rowmajor ? 0 :
                  Si == 1 ? ColMajor :
                  Sj == 1 ? RowMajor : 0 ) )};
        enum { okA = (
                !Attrib<A>::diagmajor &&
                !Attrib<A>::nonunitdiag &&
                !Attrib<A>::unitdiag &&
                !Attrib<A>::zerodiag &&
                !Attrib<A>::packed &&
                !Attrib<A>::lower &&
                !Attrib<A>::withdivider &&
                ( Traits<T>::iscomplex || !Attrib<A>::conj ) &&
                !( Si != UNKNOWN && Si != 1 && Attrib<A>::colmajor ) &&
                !( Sj != UNKNOWN && Sj != 1 && Attrib<A>::rowmajor ) )};
        enum { _attrib = A };

        typedef T value_type;
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;
        enum { isreal = Traits<T>::isreal };
        enum { iscomplex = Traits<T>::iscomplex };

        typedef SmallMatrixView<T,M,N,Si,Sj,A0> type;
        typedef ConstSmallMatrixView<T,M,N,Si,Sj,A> calc_type;
        typedef const type& eval_type;

        enum { _colsize = M };
        enum { _rowsize = N };
        enum { _shape = Rec };
        enum { _fort = Attrib<A>::fort };
        enum { _calc = true };
        enum { _rowmajor = Attrib<A>::rowmajor };
        enum { _colmajor = Attrib<A>::colmajor };
        enum { _stepi = Si };
        enum { _stepj = Sj };
        enum { _diagstep = IntTraits2<Si,Sj>::sum };
        enum { _conj = Attrib<A>::conj };
        enum { _checkalias = Attrib<A>::checkalias };
        enum { _canlin = (
                (Si == 1 && M != UNKNOWN && Sj == M) || 
                (Sj == 1 && N != UNKNOWN && Si == N) )};
        enum { _linsize = IntTraits2<M,N>::prod };
        enum { twoSi = isreal ? int(_stepi) : int(IntTraits<_stepi>::twoS) };
        enum { twoSj = isreal ? int(_stepj) : int(IntTraits<_stepj>::twoS) };
        enum { minMN = IntTraits2<M,N>::min };

        typedef typename MCopyHelper<T,Rec,M,N,_rowmajor,_fort>::type copy_type;

        enum { _hasdivider = false };
        typedef QuotXM<1,real_type,type> inverse_type;

        typedef typename TypeSelect<_hasdivider ,
                const InstLUD<value_type>& , LUD<copy_type> >::type lud_type;
        typedef typename TypeSelect<_hasdivider ,
                const InstQRD<value_type>& , QRD<copy_type> >::type qrd_type;
        typedef typename TypeSelect<_hasdivider ,
                const InstQRPD<value_type>& , QRPD<copy_type> >::type qrpd_type;
        typedef InvalidType svd_type;

        enum { vecA = (
                (_fort ? FortranStyle : 0) |
                (_conj ? Conj : 0 ) |
                (_checkalias ? CheckAlias : 0) )};
        enum { colA = vecA | (_colmajor ? Unit : 0) };
        enum { rowA = vecA | (_rowmajor ? Unit : 0) };
        enum { diagA = vecA | (_diagstep == 1 ? Unit : 0) };
        enum { cstyleA = A & ~FortranStyle };
        enum { fstyleA = A | FortranStyle };
        enum { nmA = (A & ~AllStorageType) };
        enum { cmA = nmA | ColMajor };
        enum { rmA = nmA | RowMajor };
        enum { colpairA = A & ~RowMajor };
        enum { rowpairA = A & ~ColMajor };
        enum { conjA = iscomplex ? (A ^ Conj) : int(A) };
        enum { trA = _colmajor ? rmA : _rowmajor ? cmA : int(A) };
        enum { adjA = iscomplex ? (trA ^ Conj) : int(trA) };
        enum { nonconjA = A & ~Conj };
        enum { twosA = isreal ? nonconjA : int(nonconjA & ~AllStorageType) };
        enum { An = _checkalias ? (A & ~CheckAlias) : (A | NoAlias) };
        enum { vecAn = _checkalias ? (vecA & ~CheckAlias) : (vecA | NoAlias) };
        enum { colAn = vecAn | (_colmajor ? Unit : 0) };
        enum { rowAn = vecAn | (_rowmajor ? Unit : 0) };
        enum { diagAn = vecAn | (_diagstep == 1 ? Unit : 0) };
        enum { nmAn = (An & ~AllStorageType) };

        typedef ConstSmallVectorView<T,M,_stepi,colA> const_col_type;
        typedef ConstVectorView<T,colAn> const_col_sub_type;
        typedef ConstSmallVectorView<T,N,_stepj,rowA> const_row_type;
        typedef ConstVectorView<T,rowAn> const_row_sub_type;
        typedef ConstSmallVectorView<T,minMN,_diagstep,diagA> const_diag_type;
        typedef ConstVectorView<T,diagAn> const_diag_sub_type;

        typedef ConstMatrixView<T,An> const_submatrix_type;
        typedef ConstMatrixView<T,nmAn> const_submatrix_step_type;
        typedef ConstVectorView<T,vecAn> const_subvector_type;
        typedef ConstSmallMatrixView<T,M,2,_stepi,UNKNOWN,colpairA>
            const_colpair_type;
        typedef ConstSmallMatrixView<T,2,N,UNKNOWN,_stepj,rowpairA>
            const_rowpair_type;
        typedef ConstSmallMatrixView<T,M,UNKNOWN,_stepi,_stepj,A>
            const_colrange_type;
        typedef ConstSmallMatrixView<T,UNKNOWN,N,_stepi,_stepj,A>
            const_rowrange_type;
        typedef ConstSmallMatrixView<T,M,N,_stepi,_stepj,A>
            const_view_type;
        typedef ConstSmallMatrixView<T,M,N,_stepi,_stepj,cstyleA>
            const_cview_type;
        typedef ConstSmallMatrixView<T,M,N,_stepi,_stepj,fstyleA>
            const_fview_type;
        typedef ConstMatrixView<T,(_conj ? Conj : NonConj)> const_xview_type;
        typedef ConstSmallMatrixView<T,M,N,1,_stepj,cmA> const_cmview_type;
        typedef ConstSmallMatrixView<T,M,N,_stepi,1,rmA> const_rmview_type;
        typedef ConstSmallMatrixView<T,M,N,_stepi,_stepj,conjA>
            const_conjugate_type;
        typedef ConstSmallMatrixView<T,N,M,_stepj,_stepi,trA>
            const_transpose_type;
        typedef ConstSmallMatrixView<T,N,M,_stepj,_stepi,adjA>
            const_adjoint_type;
        typedef ConstSmallUpperTriMatrixView<T,N,_stepi,_stepj,(A|NonUnitDiag)> 
            const_uppertri_type;
        typedef ConstSmallUpperTriMatrixView<T,N,_stepi,_stepj,(A|UnitDiag)>
            const_unit_uppertri_type;
        typedef ConstSmallUpperTriMatrixView<T,N,_stepi,_stepj,(A|UnknownDiag)>
            const_unknown_uppertri_type;
        typedef ConstSmallLowerTriMatrixView<T,M,_stepi,_stepj,(A|NonUnitDiag)>
            const_lowertri_type;
        typedef ConstSmallLowerTriMatrixView<T,M,_stepi,_stepj,(A|UnitDiag)>
            const_unit_lowertri_type;
        typedef ConstSmallLowerTriMatrixView<T,N,_stepi,_stepj,(A|UnknownDiag)>
            const_unknown_lowertri_type;
        typedef ConstSmallMatrixView<real_type,M,N,twoSi,twoSj,twosA>
            const_realpart_type;
        typedef const_realpart_type const_imagpart_type;
        typedef ConstSmallVectorView<T,_linsize,1,(vecA|Unit)> 
            const_linearview_type;
        typedef ConstSmallMatrixView<T,M,N,_stepi,_stepj,nonconjA>
            const_nonconj_type;
        typedef SmallMatrixView<T,M,N,_stepi,_stepj,A> nonconst_type;

        typedef typename AuxRef<T,_conj>::reference reference;
        typedef VIt<T,1,_conj> linear_iterator;

        typedef SmallVectorView<T,M,_stepi,colA> col_type;
        typedef VectorView<T,colAn> col_sub_type;
        typedef SmallVectorView<T,N,_stepj,rowA> row_type;
        typedef VectorView<T,rowAn> row_sub_type;
        typedef SmallVectorView<T,minMN,_diagstep,diagA> diag_type;
        typedef VectorView<T,diagAn> diag_sub_type;

        typedef MatrixView<T,An> submatrix_type;
        typedef MatrixView<T,nmAn> submatrix_step_type;
        typedef VectorView<T,vecAn> subvector_type;
        typedef SmallMatrixView<T,M,2,_stepi,UNKNOWN,colpairA> colpair_type;
        typedef SmallMatrixView<T,2,N,UNKNOWN,_stepj,rowpairA> rowpair_type;
        typedef SmallMatrixView<T,M,UNKNOWN,_stepi,_stepj,A> colrange_type;
        typedef SmallMatrixView<T,UNKNOWN,N,_stepi,_stepj,A> rowrange_type;
        typedef SmallMatrixView<T,M,N,_stepi,_stepj,A> view_type;
        typedef SmallMatrixView<T,M,N,_stepi,_stepj,cstyleA> cview_type;
        typedef SmallMatrixView<T,M,N,_stepi,_stepj,fstyleA> fview_type;
        typedef MatrixView<T,(_conj ? Conj : NonConj)> xview_type;
        typedef SmallMatrixView<T,M,N,1,_stepj,cmA> cmview_type;
        typedef SmallMatrixView<T,M,N,_stepi,1,rmA> rmview_type;
        typedef SmallMatrixView<T,M,N,_stepi,_stepj,conjA> conjugate_type;
        typedef SmallMatrixView<T,N,M,_stepj,_stepi,trA> transpose_type;
        typedef SmallMatrixView<T,N,M,_stepj,_stepi,adjA> adjoint_type;
        typedef SmallUpperTriMatrixView<T,N,_stepi,_stepj,(A|NonUnitDiag)> 
            uppertri_type;
        typedef SmallUpperTriMatrixView<T,N,_stepi,_stepj,(A|UnitDiag)>
            unit_uppertri_type;
        typedef SmallUpperTriMatrixView<T,N,_stepi,_stepj,(A|UnknownDiag)>
            unknown_uppertri_type;
        typedef SmallLowerTriMatrixView<T,M,_stepi,_stepj,(A|NonUnitDiag)>
            lowertri_type;
        typedef SmallLowerTriMatrixView<T,M,_stepi,_stepj,(A|UnitDiag)>
            unit_lowertri_type;
        typedef SmallLowerTriMatrixView<T,N,_stepi,_stepj,(A|UnknownDiag)>
            unknown_lowertri_type;
        typedef SmallMatrixView<real_type,M,N,twoSi,twoSj,twosA> realpart_type;
        typedef realpart_type imagpart_type;
        typedef SmallVectorView<T,_linsize,1,(vecA|Unit)> linearview_type;
        typedef SmallMatrixView<T,M,N,_stepi,_stepj,nonconjA> nonconj_type;
    };

    template <class T, int M, int N, int Si, int Sj, int A>
    class SmallMatrixView :
        public BaseMatrix_Rec_Mutable<SmallMatrixView<T,M,N,Si,Sj,A> >,
        public MatrixDivHelper<SmallMatrixView<T,M,N,Si,Sj,A> >
    {
    public:

        typedef SmallMatrixView<T,M,N,Si,Sj,A> type;
        typedef BaseMatrix_Rec_Mutable<type> base_mut;
        typedef typename base_mut::reference reference;
        typedef typename Traits<type>::lud_type lud_type;
        typedef typename Traits<type>::qrd_type qrd_type;
        typedef typename Traits<type>::qrpd_type qrpd_type;
        typedef typename Traits<type>::svd_type svd_type;

        enum { _colsize = Traits<type>::_colsize };
        enum { _rowsize = Traits<type>::_rowsize };
        enum { _shape = Traits<type>::_shape };
        enum { _fort = Traits<type>::_fort };
        enum { _rowmajor = Traits<type>::_rowmajor };
        enum { _colmajor = Traits<type>::_colmajor };
        enum { _calc = Traits<type>::_calc };
        enum { _stepi = Traits<type>::_stepi };
        enum { _stepj = Traits<type>::_stepj };
        enum { _diagstep = Traits<type>::_diagstep };
        enum { _conj = Traits<type>::_conj };
        enum { _checkalias = Traits<type>::_checkalias };
        enum { _canlin = Traits<type>::_canlin };
        enum { _linsize = Traits<type>::_linsize };
        enum { _attrib = Traits<type>::_attrib };

        //
        // Constructors
        //

        TMV_INLINE SmallMatrixView(T* m, size_t cs, size_t rs, int si, int sj) :
            itsm(m), itscs(cs), itsrs(rs), itssi(si), itssj(sj) 
        {
            TMVStaticAssert(Traits<type>::okA);
        }

        TMV_INLINE SmallMatrixView(T* m, size_t cs, size_t rs, int si) :
            itsm(m), itscs(cs), itsrs(rs), itssi(si), itssj(Sj)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Sj != UNKNOWN); 
        }

        TMV_INLINE SmallMatrixView(T* m, size_t cs, size_t rs) :
            itsm(m), itscs(cs), itsrs(rs), itssi(Si), itssj(Sj)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Si != UNKNOWN); TMVStaticAssert(Sj != UNKNOWN); 
        }

        TMV_INLINE SmallMatrixView(T* m, size_t cs) :
            itsm(m), itscs(cs), itsrs(N), itssi(Si), itssj(Sj)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(N != UNKNOWN);
            TMVStaticAssert(Si != UNKNOWN); TMVStaticAssert(Sj != UNKNOWN); 
        }

        TMV_INLINE SmallMatrixView(T* m) :
            itsm(m), itscs(M), itsrs(N), itssi(Si), itssj(Sj)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(M != UNKNOWN); TMVStaticAssert(N != UNKNOWN);
            TMVStaticAssert(Si != UNKNOWN); TMVStaticAssert(Sj != UNKNOWN); 
        }

        TMV_INLINE SmallMatrixView(const type& m2) :
            itsm(m2.itsm), itscs(m2.colsize()), itsrs(m2.rowsize()),
            itssi(m2.stepi()), itssj(m2.stepj()) 
        {
            TMVStaticAssert(Traits<type>::okA);
        }

        template <int A2>
        TMV_INLINE SmallMatrixView(MatrixView<T,A2> m2) :
            itsm(m2.ptr()), itscs(m2.colsize()), itsrs(m2.rowsize()),
            itssi(m2.stepi()), itssj(m2.stepj()) 
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Attrib<A>::conj == int(Attrib<A2>::conj)); 
        }

        template <int M2, int N2, int Si2, int Sj2, int A2>
        TMV_INLINE SmallMatrixView(SmallMatrixView<T,M2,N2,Si2,Sj2,A2> m2) :
            itsm(m2.ptr()), itscs(m2.colsize()), itsrs(m2.rowsize()),
            itssi(m2.stepi()), itssj(m2.stepj()) 
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Attrib<A>::conj == int(Attrib<A2>::conj)); 
        }

        TMV_INLINE ~SmallMatrixView() {
#ifdef TMV_DEBUG
            itsm = 0; 
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

        TMV_INLINE type& operator=(const T x)
        { base_mut::operator=(x); return *this; }


        //
        // Auxilliary Functions
        //

        TMV_INLINE const T* cptr() const { return itsm; }
        TMV_INLINE T* ptr() { return itsm; }

        T cref(int i, int j) const
        { return DoConj<_conj>(itsm[i*stepi()+j*stepj()]); }

        reference ref(int i, int j)
        { return reference(itsm[i*stepi()+j*stepj()]); }

        TMV_INLINE size_t ls() const { return int(itscs)*int(itsrs); }
        TMV_INLINE size_t colsize() const { return itscs; }
        TMV_INLINE size_t rowsize() const { return itsrs; }
        TMV_INLINE int nElements() const { return int(itscs)*int(itsrs); }
        TMV_INLINE int stepi() const { return itssi; }
        TMV_INLINE int stepj() const { return itssj; }
        TMV_INLINE bool isconj() const { return _conj; }
        TMV_INLINE bool isrm() const 
        { return _rowmajor || (!_colmajor && stepj() == 1); }
        TMV_INLINE bool iscm() const 
        { return _colmajor || (!_rowmajor && stepi() == 1); }

    private :

        T* itsm;
        const CheckedInt<M> itscs;
        const CheckedInt<N> itsrs;
        const CheckedInt<Si> itssi;
        const CheckedInt<Sj> itssj;

    }; // SmallMatrixView


    // Special Creators: 
    //   RowVectorViewOf(v) = 1xn Matrix with v in only row - Same Storage
    //   ColVectorViewOf(v) = nx1 Matrix with v in only col - Same Storage

    // A quick helper class to get the return types correct for
    // RowVectorView and ColVectorView.
    template <class V>
    struct VVO
    {
        typedef typename V::value_type T;
        enum { N = V::_size };
        enum { S = V::_step };
        enum { vecA = (
                ( V::_conj ? Conj : NonConj ) |
                ( V::_fort ? FortranStyle : CStyle ) |
                ( Traits<V>::A & AllAliasStatus ) )};
        enum { colA = vecA | (S == 1 ? ColMajor : 0 ) };
        enum { rowA = vecA | (S == 1 ? RowMajor : 0 ) };
        typedef ConstSmallMatrixView<T,1,N,N,S,rowA> crv;
        typedef SmallMatrixView<T,1,N,N,S,rowA> rv;
        typedef ConstSmallMatrixView<T,N,1,S,N,colA> ccv;
        typedef SmallMatrixView<T,N,1,S,N,colA> cv;
    };

    template <class V>
    static TMV_INLINE typename VVO<V>::crv RowVectorViewOf(
        const BaseVector_Calc<V>& v)
    {
        return typename VVO<V>::crv(
            v.cptr(),1,v.size(),v.size(),v.step());
    }

    template <class V>
    static TMV_INLINE typename VVO<V>::rv RowVectorViewOf(
        BaseVector_Mutable<V>& v)
    {
        return typename VVO<V>::rv(
            v.ptr(),1,v.size(),v.size(),v.step());
    }

    template <class V>
    static TMV_INLINE typename VVO<V>::ccv ColVectorViewOf(
        const BaseVector_Calc<V>& v)
    {
        return typename VVO<V>::ccv(
            v.cptr(),v.size(),1,v.step(),v.size());
    }

    template <class V>
    static TMV_INLINE typename VVO<V>::cv ColVectorViewOf(
        BaseVector_Mutable<V>& v)
    {
        return typename VVO<V>::cv(
            v.ptr(),v.size(),1,v.step(),v.size());
    }

    // Repeat for VectorView and SmallVectorView so we can have them 
    // work without requiring the non-const reference.
    template <class T, int A>
    static TMV_INLINE typename VVO<VectorView<T,A> >::rv RowVectorViewOf(
        VectorView<T,A> v)
    {
        return typename VVO<VectorView<T,A> >::rv(
            v.ptr(),1,v.size(),v.size(),v.step());
    }

    template <class T, int N, int S, int A>
    static TMV_INLINE typename VVO<SmallVectorView<T,N,S,A> >::rv RowVectorViewOf(
        SmallVectorView<T,N,S,A> v)
    {
        return typename VVO<SmallVectorView<T,N,S,A> >::rv(
            v.ptr(),1,v.size(),v.size(),v.step());
    }

    template <class T, int S, int A>
    static TMV_INLINE typename VVO<VectorView<T,A> >::cv ColVectorViewOf(
        VectorView<T,A> v)
    {
        return typename VVO<VectorView<T,A> >::cv(
            v.ptr(),v.size(),1,v.step(),v.size());
    }

    template <class T, int N, int S, int A>
    static TMV_INLINE typename VVO<SmallVectorView<T,N,S,A> >::cv ColVectorViewOf(
        SmallVectorView<T,N,S,A> v)
    {
        return typename VVO<SmallVectorView<T,N,S,A> >::cv(
            v.ptr(),v.size(),1,v.step(),v.size());
    }


    //
    // Swap
    //

    template <class T, int M, int N, int Si, int Sj, int A, class MM>
    static TMV_INLINE void Swap(
        BaseMatrix_Rec_Mutable<MM>& m1, SmallMatrixView<T,M,N,Si,Sj,A> m2)
    { DoSwap(m1,m2); }
    template <class T, int M, int N, int Si, int Sj, int A, class MM>
    static TMV_INLINE void Swap(
        SmallMatrixView<T,M,N,Si,Sj,A> m1, BaseMatrix_Rec_Mutable<MM>& m2)
    { DoSwap(m1,m2); }
    template <class T, int M, int N, int Si1, int Sj1, int A1, int Si2, int Sj2, int A2>
    static TMV_INLINE void Swap(
        SmallMatrixView<T,M,N,Si1,Sj1,A1> m1,
        SmallMatrixView<T,M,N,Si2,Sj2,A2> m2)
    { DoSwap(m1,m2); }
    template <class T, int M, int N, int Si1, int Sj1, int A1, int A2>
    static TMV_INLINE void Swap(
        SmallMatrixView<T,M,N,Si1,Sj1,A1> m1, 
        MatrixView<T,A2> m2)
    { DoSwap(m1,m2); }
    template <class T, int M, int N, int A1, int Si2, int Sj2, int A2>
    static TMV_INLINE void Swap(
        MatrixView<T,A1> m1,
        SmallMatrixView<T,M,N,Si2,Sj2,A2> m2)
    { DoSwap(m1,m2); }


    //
    // Conjugate, Transpose, Adjoint
    //

    template <class T, int M, int N, int A0, int A1>
    static TMV_INLINE typename SmallMatrix<T,M,N,A0,A1>::conjugate_type Conjugate(
        SmallMatrix<T,M,N,A0,A1>& m)
    { return m.conjugate(); }
    template <class T, int M, int N, int Si, int Sj, int A>
    static TMV_INLINE typename SmallMatrixView<T,M,N,Si,Sj,A>::conjugate_type Conjugate(
        SmallMatrixView<T,M,N,Si,Sj,A> m)
    { return m.conjugate(); }

    template <class T, int M, int N, int A0, int A1>
    static TMV_INLINE typename SmallMatrix<T,M,N,A0,A1>::transpose_type Transpose(
        SmallMatrix<T,M,N,A0,A1>& m)
    { return m.transpose(); }
    template <class T, int M, int N, int Si, int Sj, int A>
    static TMV_INLINE typename SmallMatrixView<T,M,N,Si,Sj,A>::transpose_type Transpose(
        SmallMatrixView<T,M,N,Si,Sj,A> m)
    { return m.transpose(); }

    template <class T, int M, int N, int A0, int A1>
    static TMV_INLINE typename SmallMatrix<T,M,N,A0,A1>::adjoint_type Adjoint(
        SmallMatrix<T,M,N,A0,A1>& m)
    { return m.adjoint(); }
    template <class T, int M, int N, int Si, int Sj, int A>
    static TMV_INLINE typename SmallMatrixView<T,M,N,Si,Sj,A>::adjoint_type Adjoint(
        SmallMatrixView<T,M,N,Si,Sj,A> m)
    { return m.adjoint(); }


    // 
    // TMV_Text
    //

#ifdef TMV_TEXT
    template <class T, int M, int N, int A0, int A1>
    static inline std::string TMV_Text(const SmallMatrix<T,M,N,A0,A1>& m)
    {
        const int A = A0 | A1;
        std::ostringstream s;
        s << "SmallMatrix<"<<TMV_Text(T());
        s << ','<<M<<','<<N;
        s << ","<<Attrib<A>::text()<<">";
        s << "("<<m.colsize()<<","<<m.rowsize()<<",";
        s << m.stepi()<<","<<m.stepj()<<")";
        return s.str();
    }

    template <class T, int M, int N, int Si, int Sj, int A>
    static inline std::string TMV_Text(
        const ConstSmallMatrixView<T,M,N,Si,Sj,A>& m)
    {
        std::ostringstream s;
        s << "ConstSmallMatrixView<"<<TMV_Text(T());
        s << ','<<IntTraits<M>::text()<<','<<IntTraits<N>::text();
        s << ','<<IntTraits<Si>::text()<<','<<IntTraits<Sj>::text();
        s << ","<<Attrib<A>::text()<<">";
        s << "("<<m.colsize()<<","<<m.rowsize()<<",";
        s << m.stepi()<<","<<m.stepj()<<")";
        return s.str();
    }

    template <class T, int M, int N, int Si, int Sj, int A>
    static inline std::string TMV_Text(
        const SmallMatrixView<T,M,N,Si,Sj,A>& m)
    {
        std::ostringstream s;
        s << "SmallMatrixView<"<<TMV_Text(T());
        s << ','<<IntTraits<M>::text()<<','<<IntTraits<N>::text();
        s << ','<<IntTraits<Si>::text()<<','<<IntTraits<Sj>::text();
        s << ","<<Attrib<A>::text()<<">";
        s << "("<<m.colsize()<<","<<m.rowsize()<<",";
        s << m.stepi()<<","<<m.stepj()<<")";
        return s.str();
    }
#endif

} // namespace tmv

#endif
