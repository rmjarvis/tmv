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
// This file defines the BaseMatrix_Diag and BaseMatrix_Diag_Mutable classes.
//
// See TMV_DiagMatrix.h for the functions that are defined for these objects.
//

#ifndef TMV_BaseMatrix_Diag_H
#define TMV_BaseMatrix_Diag_H

#include "TMV_BaseMatrix.h"
#include "TMV_BaseMatrix_Rec.h"

namespace tmv {

    template <class M>
    class BaseMatrix_Diag;
    template <class M>
    class BaseMatrix_Diag_Mutable;

    // Defined in TMV_DiagMatrixIO.h
    template <class M>
    static void Read(std::istream& is, BaseMatrix_Diag_Mutable<M>& m);

    // Defined in InvertD.h
    template <int ix, class T, class M1, class M2>
    static void Invert(
        const Scaling<ix,T>& x,
        const BaseMatrix_Diag<M1>& m1, BaseMatrix_Mutable<M2>& m2);
    template <int ix, class T, class M1, class M2>
    static void Invert(
        const Scaling<ix,T>& x,
        const BaseMatrix_Diag<M1>& m1, BaseMatrix_Diag_Mutable<M2>& m2);

    // Defined in ElemMultVV.h
    template <bool add, int ix, class T, class V1, class V2, class V3>
    static void NoAliasElemMultVV(
        const Scaling<ix,T>& x1, const BaseVector_Calc<V1>& v1,
        const BaseVector_Calc<V2>& v2, BaseVector_Mutable<V3>& v3);

    // Defined below:
    template <class M1, class M2> 
    static void Copy(
        const BaseMatrix_Diag<M1>& m1, BaseMatrix_Rec_Mutable<M2>& m2);
    template <class M1, class M2> 
    static void Copy(
        const BaseMatrix_Diag<M1>& m1, BaseMatrix_Diag_Mutable<M2>& m2);
    template <class M1, class M2> 
    static void NoAliasCopy(
        const BaseMatrix_Diag<M1>& m1, BaseMatrix_Rec_Mutable<M2>& m2);
    template <class M1, class M2> 
    static void NoAliasCopy(
        const BaseMatrix_Diag<M1>& m1, BaseMatrix_Diag_Mutable<M2>& m2);

    template <class M>
    class BaseMatrix_Diag : 
        public BaseMatrix_Calc<M>
    {
    public:
        enum { _colsize = Traits<M>::_size };
        enum { _rowsize = Traits<M>::_size };
        enum { _size = Traits<M>::_size };
        enum { _fort = Traits<M>::_fort };
        enum { _shape = Traits<M>::_shape };
        enum { _rowmajor = Traits<M>::_rowmajor };
        enum { _colmajor = Traits<M>::_colmajor };
        enum { _calc = Traits<M>::_calc };
        enum { _diagstep = Traits<M>::_diagstep };
        enum { _conj = Traits<M>::_conj };

        typedef M type;

        typedef typename Traits<M>::value_type value_type;

        typedef typename Traits<M>::calc_type calc_type;
        typedef typename Traits<M>::eval_type eval_type;
        typedef typename Traits<M>::copy_type copy_type;

        typedef typename Traits<M>::const_diag_type const_diag_type;

        typedef typename Traits<M>::const_subdiagmatrix_type 
            const_subdiagmatrix_type;
        typedef typename Traits<M>::const_subdiagmatrix_step_type 
            const_subdiagmatrix_step_type;
        typedef typename Traits<M>::const_view_type const_view_type;
        typedef typename Traits<M>::const_cview_type const_cview_type;
        typedef typename Traits<M>::const_fview_type const_fview_type;
        typedef typename Traits<M>::const_xview_type const_xview_type;
        typedef typename Traits<M>::const_cmview_type const_cmview_type;
        typedef typename Traits<M>::const_rmview_type const_rmview_type;
        typedef typename Traits<M>::const_transpose_type 
            const_transpose_type;
        typedef typename Traits<M>::const_conjugate_type 
            const_conjugate_type;
        typedef typename Traits<M>::const_adjoint_type const_adjoint_type;
        typedef typename Traits<M>::const_realpart_type const_realpart_type;
        typedef typename Traits<M>::const_imagpart_type const_imagpart_type;
        typedef typename Traits<M>::const_nonconj_type const_nonconj_type;
        typedef typename Traits<M>::nonconst_type nonconst_type;

        // Derived values:
        typedef typename Traits<value_type>::real_type real_type;
        typedef typename Traits<real_type>::float_type float_type;
        typedef typename Traits<value_type>::float_type zfloat_type;
        typedef typename Traits<value_type>::complex_type complex_type;


        //
        // Constructor
        //

        BaseMatrix_Diag() {}
        BaseMatrix_Diag(const BaseMatrix_Diag<M>&) {}
        ~BaseMatrix_Diag() {}

    private:
        void operator=(const BaseMatrix_Diag<M>&);
    public:


        //
        // Access 
        //

        value_type operator()(int i, int j) const
        {
            CheckRowIndex<_fort>(i,size());
            CheckColIndex<_fort>(j,size());
            TMVAssert(i == j);
            return cref(i);
        }

        value_type operator()(int i) const
        {
            CheckIndex<_fort>(i,size());
            return cref(i);
        }

        const_diag_type diag() const
        { return const_diag_type(cptr(),size(),step()); }


        //
        // Functions
        //

        value_type sumElements() const
        { return diag().sumElements(); }

        float_type sumAbsElements() const
        { return diag().sumAbsElements(); }

        real_type sumAbs2Elements() const
        { return diag().sumAbs2Elements(); }

        float_type maxAbsElement() const
        { return diag().maxAbsElement(); }

        real_type maxAbs2Element() const
        { return diag().maxAbs2Element(); }

        real_type normSq() const
        { return diag().normSq(); }

        float_type normSq(const float_type scale) const
        { return diag().normSq(scale); }

        float_type normF() const 
        { return diag().norm2(); }

        float_type norm() const
        { return normF(); }

        float_type norm1() const
        { return diag().maxAbsElement(); }

        float_type normInf() const
        { return diag().maxAbsElement(); }

        template <class ret_type, class F>
        ret_type sumElements(const F& f) const
        { return diag().sumElements(f); }



        //
        // Division Functions
        //

        float_type norm2() const
        { return diag().maxAbsElement(); }

        float_type condition() const
        { return diag().maxAbsElement() / diag().minAbsElement(); }


        //
        // subDiagMatrix, etc.
        //

        // These versions always uses CStyle
        const_subdiagmatrix_type cSubDiagMatrix(int i1, int i2) const
        {
            return const_subdiagmatrix_type(
                cptr()+i1*step(), i2-i1, step());
        }

        const_subdiagmatrix_step_type cSubDiagMatrix(
            int i1, int i2, int istep) const
        {
            return const_subdiagmatrix_step_type(
                cptr()+i1*step(), (i2-i1)/istep, istep*step());
        }


        // These check the indices according the the indexing style being
        // used, and then calls the above CStyle versions.
        const_subdiagmatrix_type subDiagMatrix(int i1, int i2) const
        {
            CheckRange<_fort>(i1,i2,size());
            return cSubDiagMatrix(i1,i2);
        }

        const_subdiagmatrix_step_type subDiagMatrix(
            int i1, int i2, int istep) const
        {
            CheckRange<_fort>(i1,i2,istep,size());
            return cSubDiagMatrix(i1,i2,istep);
        }


        //
        // Views
        //

        const_view_type view() const
        { return const_view_type(cptr(),size(),step()); }

        const_cview_type cView() const
        { return view(); }

        const_fview_type fView() const
        { return view(); }

        const_xview_type xView() const
        { return view(); }

        const_cmview_type cmView() const
        {
            TMVAssert(step()==1 && "Called cmView on DiagMatrix with step!=1");
            return view(); 
        }

        const_rmview_type rmView() const
        {
            TMVAssert(step()==1 && "Called rmView on DiagMatrix with step!=1");
            return view(); 
        }

        const_view_type constView() const
        { return view(); }

        const_transpose_type transpose() const
        { return view(); }

        const_conjugate_type conjugate() const
        { return const_conjugate_type(cptr(),size(),step()); }

        const_adjoint_type adjoint() const
        { return conjugate(); }

        const_realpart_type realPart() const
        {
            const bool isreal = Traits<value_type>::isreal;
            return const_realpart_type(
                reinterpret_cast<const real_type*>(cptr()), size(),
                isreal ? step() : 2*step());
        }

        const_imagpart_type imagPart() const
        {
            const bool isreal = Traits<value_type>::isreal;
            TMVStaticAssert(Traits<value_type>::iscomplex);
            return const_imagpart_type(
                reinterpret_cast<const real_type*>(cptr())+1, size(),
                isreal ? step() : 2*step());
        }

        const_nonconj_type nonConj() const
        { return const_nonconj_type(cptr(),size(),step()); }

        nonconst_type nonConst() const
        { return nonconst_type(const_cast<value_type*>(cptr()),size(),step()); }


        //
        // I/O
        //

        void writeCompact(std::ostream& os) const
        { os << "D " << diag() << std::endl; }
        void writeCompact(std::ostream& os, float_type thresh) const
        { os << "D "; diag().write(os,thresh); os << std::endl; }


        //
        // Auxilliary routines
        //

        template <class M2>
        void assignTo(BaseMatrix_Mutable<M2>& m2) const
        { 
            TMVStaticAssert((ShapeTraits2<_shape,M2::_shape>::assignable));
            tmv::Copy(mat(),m2.mat()); 
        }

        template <class M2>
        void newAssignTo(BaseMatrix_Mutable<M2>& m2) const
        { 
            TMVStaticAssert((ShapeTraits2<_shape,M2::_shape>::assignable));
            tmv::NoAliasCopy(mat(),m2.mat()); 
        }

        const type& mat() const
        { return *static_cast<const type*>(this); }

        bool isconj() const { return _conj; }

        // Note that these last functions need to be defined in a more derived
        // class than this, or an infinite loop will result.
        // Also, cref, get_row and get_col from BaseMatrix.

        size_t colsize() const { return mat().size(); }
        size_t rowsize() const { return mat().size(); }
        size_t size() const { return mat().size(); }
        int step() const { return mat().step(); }

        const value_type* cptr() const { return mat().cptr(); }
        value_type cref(int i, int j) const 
        { return (i!=j ? value_type(0) : cref(i)); }
        value_type cref(int i) const 
        { return mat().cref(i); }

    };

    // Specify ExactSameStorage for diagonal matrices:
    template <class M1, class M2>
    bool ExactSameStorage(
        const BaseMatrix_Diag<M1>& m1, const BaseMatrix_Diag<M2>& m2)
    {
        typedef typename M1::value_type T1;
        typedef typename M2::value_type T2;
        return Traits2<T1,T2>::sametype && (m1.step() == m2.step()); 
    }

    template <class M>
    class BaseMatrix_Diag_Mutable : 
        public BaseMatrix_Diag<M>,
        public BaseMatrix_Mutable<M>
    {
    public:
        enum { _colsize = Traits<M>::_size };
        enum { _rowsize = Traits<M>::_size };
        enum { _size = Traits<M>::_size };
        enum { _fort = Traits<M>::_fort };
        enum { _shape = Traits<M>::_shape };
        enum { _rowmajor = Traits<M>::_rowmajor };
        enum { _colmajor = Traits<M>::_colmajor };
        enum { _calc = Traits<M>::_calc };
        enum { _diagstep = Traits<M>::_diagstep };
        enum { _conj = Traits<M>::_conj };

        typedef M type;
        typedef BaseMatrix_Diag<M> base_diag;

        typedef typename Traits<M>::value_type value_type;

        typedef typename Traits<M>::calc_type calc_type;
        typedef typename Traits<M>::eval_type eval_type;
        typedef typename Traits<M>::copy_type copy_type;

        typedef typename Traits<M>::const_diag_type const_diag_type;

        typedef typename Traits<M>::const_subdiagmatrix_type 
            const_subdiagmatrix_type;
        typedef typename Traits<M>::const_subdiagmatrix_step_type 
            const_subdiagmatrix_step_type;
        typedef typename Traits<M>::const_view_type const_view_type;
        typedef typename Traits<M>::const_cview_type const_cview_type;
        typedef typename Traits<M>::const_fview_type const_fview_type;
        typedef typename Traits<M>::const_xview_type const_xview_type;
        typedef typename Traits<M>::const_cmview_type const_cmview_type;
        typedef typename Traits<M>::const_rmview_type const_rmview_type;
        typedef typename Traits<M>::const_transpose_type 
            const_transpose_type;
        typedef typename Traits<M>::const_conjugate_type 
            const_conjugate_type;
        typedef typename Traits<M>::const_adjoint_type const_adjoint_type;
        typedef typename Traits<M>::const_realpart_type const_realpart_type;
        typedef typename Traits<M>::const_imagpart_type const_imagpart_type;
        typedef typename Traits<M>::const_nonconj_type const_nonconj_type;

        typedef typename Traits<M>::diag_type diag_type;
        typedef typename Traits<M>::subdiagmatrix_type subdiagmatrix_type;
        typedef typename Traits<M>::subdiagmatrix_step_type 
            subdiagmatrix_step_type;
        typedef typename Traits<M>::view_type view_type;
        typedef typename Traits<M>::cview_type cview_type;
        typedef typename Traits<M>::fview_type fview_type;
        typedef typename Traits<M>::xview_type xview_type;
        typedef typename Traits<M>::cmview_type cmview_type;
        typedef typename Traits<M>::rmview_type rmview_type;
        typedef typename Traits<M>::transpose_type transpose_type;
        typedef typename Traits<M>::conjugate_type conjugate_type;
        typedef typename Traits<M>::adjoint_type adjoint_type;
        typedef typename Traits<M>::realpart_type realpart_type;
        typedef typename Traits<M>::imagpart_type imagpart_type;
        typedef typename Traits<M>::nonconj_type nonconj_type;

        typedef typename Traits<M>::reference reference;

        // Derived values:
        typedef typename Traits<value_type>::real_type real_type;
        typedef typename Traits<real_type>::float_type float_type;
        typedef typename Traits<value_type>::float_type zfloat_type;
        typedef typename Traits<value_type>::complex_type complex_type;

        //
        // Constructor
        //

        BaseMatrix_Diag_Mutable() {}
        BaseMatrix_Diag_Mutable(const BaseMatrix_Diag_Mutable<M>&) {}
        ~BaseMatrix_Diag_Mutable() {}


        //
        // Access 
        //

        reference operator()(int i, int j)
        {
            TMVAssert(i == j);
            CheckIndex<_fort>(i,size());
            return ref(i);
        }

        reference operator()(int i)
        {
            CheckIndex<_fort>(i,size());
            return ref(i);
        }

        diag_type diag() 
        { return diag_type(ptr(),size(),step()); }


        // We need to repeat the const versions so the non-const ones
        // don't clobber them.
        value_type operator()(int i, int j) const
        { return base_diag::operator()(i,j); }
        const_diag_type diag() const
        { return base_diag::diag(); }


        //
        // Op =
        //

        type& operator=(BaseMatrix_Diag_Mutable<M>& m2) 
        {
            TMVAssert(colsize() == m2.colsize());
            TMVAssert(rowsize() == m2.rowsize());
            m2.assignTo(mat());
            return mat(); 
        }

        template <class M2>
        type& operator=(const BaseMatrix<M2>& m2) 
        {
            TMVStaticAssert((Sizes<_colsize,M2::_colsize>::same));
            TMVStaticAssert((Sizes<_rowsize,M2::_rowsize>::same));
            TMVAssert(colsize() == m2.colsize());
            TMVAssert(rowsize() == m2.rowsize());
            m2.assignTo(mat());
            return mat(); 
        }

        type& operator=(const value_type x)
        {
            TMVStaticAssert((Sizes<_rowsize,_colsize>::same));
            TMVAssert(colsize() == rowsize());
            setToIdentity(x);
            return mat();
        }


        //
        // Modifying Functions
        //

        type& setZero() 
        { return setAllTo(value_type(0)); }

        type& setAllTo(value_type x) 
        { diag().setAllTo(x); return mat(); }

        type& addToAll(value_type x) 
        { diag().addToAll(x); return mat(); }

        type& clip(float_type thresh) 
        { diag().clip(thresh); return mat(); }

        template <class F>
        type& applyToAll(const F& f)
        { diag().applyToAll(f); return mat(); }

        type& conjugateSelf() 
        { diag().conjugateSelf(); return mat(); }

        type& transposeSelf() 
        { return mat(); }

        type& invertSelf() 
        {
            const int n=size();
            for(int i=0;i<n;++i) ref(i) = real_type(1)/cref(i);
            return mat(); 
        }

        type& setToIdentity(const value_type x=value_type(1))
        { diag().setAllTo(x); return mat(); }

        type& swap(int i1, int i2) 
        {
            CheckIndex<_fort>(i1,size());
            CheckIndex<_fort>(i2,size());
            diag().cSwap(i1,i2);
            return mat();
        }

        type& cPermute(const int* p, int i1, int i2) 
        { diag().permute(p,i1,i2); return mat(); }
        type& permute(const int* p, int i1, int i2) 
        {
            CheckRange<_fort>(i1,i2,size());
            return cPermute(p,i1,i2);
        }
        type& permute(const int* p) 
        { cPermute(p,0,size()); return mat(); }

        type& cReversePermute(const int* p, int i1, int i2) 
        { diag().reversePermute(p,i1,i2); }
        type& reversePermute(const int* p, int i1, int i2) 
        {
            CheckRange<_fort>(i1,i2,size());
            return cReversePermute(p,i1,i2);
        }
        type& reversePermute(const int* p) 
        { cReversePermute(p,0,size()); return mat(); }


        //
        // subDiagMatrix, etc.
        //

        // These versions always uses CStyle
        subdiagmatrix_type cSubDiagMatrix(int i1, int i2) 
        { return subdiagmatrix_type(ptr()+i1*step(), i2-i1, step()); }

        subdiagmatrix_step_type cSubDiagMatrix(
            int i1, int i2, int istep) 
        {
            return subdiagmatrix_step_type(
                ptr()+i1*step(), (i2-i1)/istep, istep*step());
        }

        // These check the indices according the the indexing style being
        // used, and then calls the above CStyle versions.
        subdiagmatrix_type subDiagMatrix(int i1, int i2) 
        {
            CheckRange<_fort>(i1,i2,size());
            return cSubDiagMatrix(i1,i2);
        }

        subdiagmatrix_step_type subMatrix(
            int i1, int i2, int istep) 
        {
            CheckRange<_fort>(i1,i2,istep,size());
            return cSubDiagMatrix(i1,i2,istep);
        }


        // Repeat the const versions:
        const_subdiagmatrix_type cSubDiagMatrix(int i1, int i2) const
        { return base_diag::cSubDiagMatrix(i1,i2); }
        const_subdiagmatrix_step_type cSubDiagMatrix(
            int i1, int i2, int istep) const
        { return base_diag::cSubDiagMatrix(i1,i2,istep); }

        const_subdiagmatrix_type subDiagMatrix(int i1, int i2) const
        { return base_diag::subDiagMatrix(i1,i2); }
        const_subdiagmatrix_step_type subDiagMatrix(
            int i1, int i2, int istep) const
        { return base_diag::subDiagMatrix(i1,i2,istep); }


        //
        // Views
        //

        view_type view() 
        { return view_type(ptr(),size(),step()); }

        cview_type cView() 
        { return view(); }

        fview_type fView() 
        { return view(); }

        xview_type xView() 
        { return view(); }

        cmview_type cmView() 
        {
            TMVAssert(step()==1 && "Called cmView on DiagMatrix with step!=1");
            return view(); 
        }

        rmview_type rmView() 
        {
            TMVAssert(step()==1 && "Called rmView on DiagMatrix with step!=1");
            return view(); 
        }

        transpose_type transpose() 
        { return view(); }

        conjugate_type conjugate() 
        { return conjugate_type(ptr(),size(),step()); }

        adjoint_type adjoint() 
        { return conjugate(); }

        realpart_type realPart() 
        {
            const bool isreal = Traits<value_type>::isreal;
            return realpart_type(
                reinterpret_cast<real_type*>(ptr()), size(),
                isreal ? step() : 2*step());
        }

        imagpart_type imagPart() 
        {
            const bool isreal = Traits<value_type>::isreal;
            TMVStaticAssert(Traits<value_type>::iscomplex);
            return imagpart_type(
                reinterpret_cast<real_type*>(ptr())+1, size(),
                isreal ? step() : 2*step());
        }

        nonconj_type nonConj()
        { return nonconj_type(ptr(),size(),step()); }


        // Repeat the const versions:
        const_view_type view() const
        { return base_diag::view(); }
        const_cview_type cView() const
        { return base_diag::cView(); }
        const_fview_type fView() const
        { return base_diag::fView(); }
        const_xview_type xView() const
        { return base_diag::xView(); }
        const_cmview_type cmView() const
        { return base_diag::cmView(); }
        const_rmview_type rmView() const
        { return base_diag::rmView(); }
        const_transpose_type transpose() const
        { return base_diag::transpose(); }
        const_conjugate_type conjugate() const
        { return base_diag::conjugate(); }
        const_adjoint_type adjoint() const
        { return base_diag::adjoint(); }
        const_realpart_type realPart() const
        { return base_diag::realPart(); }
        const_imagpart_type imagPart() const
        { return base_diag::imagPart(); }
        const_nonconj_type nonConj() const
        { return base_diag::nonConj(); }

        //
        // I/O
        //

        void read(std::istream& is)
        { tmv::Read(is,mat()); }

        //
        // Auxilliary routines
        //

        const type& mat() const
        { return *static_cast<const type*>(this); }
        type& mat()
        { return *static_cast<type*>(this); }


        // Note that these last functions need to be defined in a more derived
        // class than this, or an infinite loop will result when compiling.
        // Also, cref and cptr from above.

        size_t colsize() const { return mat().size(); }
        size_t rowsize() const { return mat().size(); }
        size_t size() const { return mat().size(); }
        int step() const { return mat().step(); }

        value_type* ptr() { return mat().ptr(); }
        reference ref(int i) { return mat().ref(i); }
        value_type cref(int i) { return mat().cref(i); }

    }; // BaseMatrix_Diag_Mutable

    // The BaseMatrix Trace call is efficient for composite types, since
    // it avoid calculating all the elements to do the sum.
    // But if we do have the elements calculated, this overloaded 
    // version will be faster:
    template <class M>
    static typename M::value_type DoTrace(const BaseMatrix_Diag<M>& m)
    { return m.diag().sumElements(); }


    template <class T, IndexStyle I=CStyle>
    class DiagMatrix;
    template <class T, int S=UNKNOWN, bool C=false, IndexStyle I=CStyle>
    class ConstDiagMatrixView;
    template <class T, int S=UNKNOWN, bool C=false, IndexStyle I=CStyle>
    class DiagMatrixView;
    template <class T, int N, IndexStyle I=CStyle>
    class SmallDiagMatrix;
    template <class T, int N, int S, bool C=false, IndexStyle I=CStyle>
    class ConstSmallDiagMatrixView;
    template <class T, int N, int S, bool C=false, IndexStyle I=CStyle>
    class SmallDiagMatrixView;

    // This helper class helps decide calc_type for composite classes:
    template <class T, int cs, int rs, bool rm, bool fort>
    struct MCopyHelper<T,Diag,cs,rs,rm,fort>
    { typedef SmallDiagMatrix<T,cs,fort?FortranStyle:CStyle> type; };
    template <class T, int rs, bool rm, bool fort>
    struct MCopyHelper<T,Diag,UNKNOWN,rs,rm,fort>
    { typedef SmallDiagMatrix<T,rs,fort?FortranStyle:CStyle> type; };
    template <class T, int cs, bool rm, bool fort>
    struct MCopyHelper<T,Diag,cs,UNKNOWN,rm,fort>
    { typedef SmallDiagMatrix<T,cs,fort?FortranStyle:CStyle> type; };
    template <class T, bool rm, bool fort>
    struct MCopyHelper<T,Diag,UNKNOWN,UNKNOWN,rm,fort>
    { typedef DiagMatrix<T,fort?FortranStyle:CStyle> type; };

    //
    // Copy Matrices
    //

    template <class M1, class M2> 
    static void Copy(
        const BaseMatrix_Diag<M1>& m1, BaseMatrix_Diag_Mutable<M2>& m2)
    { 
        typename M1::const_diag_type m1d = m1.diag();
        typename M2::diag_type m2d = m2.diag();
        Copy(m1d,m2d);
    }

    template <class M1, class M2> 
    static void NoAliasCopy(
        const BaseMatrix_Diag<M1>& m1, BaseMatrix_Diag_Mutable<M2>& m2)
    { 
        typename M1::const_diag_type m1d = m1.diag();
        typename M2::diag_type m2d = m2.diag();
        NoAliasCopy(m1d,m2d);
    }

    //
    // M = D
    //

    template <class M1, class M2> 
    static void Copy(
        const BaseMatrix_Diag<M1>& m1, BaseMatrix_Rec_Mutable<M2>& m2)
    {
        typename M1::const_diag_type m1d = m1.diag();
        typename M2::diag_type m2d = m2.diag();
        if (SameStorage(m1,m2)) {
            AliasCopy(m1d,m2d);
            m2.upperTri().offDiag().setZero();
            m2.lowerTri().offDiag().setZero();
        } else {
            m2.setZero();
            NoAliasCopy(m1d,m2d);
        }
    }

    template <class M1, class M2> 
    static void NoAliasCopy(
        const BaseMatrix_Diag<M1>& m1, BaseMatrix_Rec_Mutable<M2>& m2)
    {
        typename M1::const_diag_type m1d = m1.diag();
        typename M2::diag_type m2d = m2.diag();
        m2.setZero();
        NoAliasCopy(m1d,m2d);
    }

    template <class M1, class M2>
    static void AliasCopy(
        const BaseMatrix_Diag<M1>& m1, BaseMatrix_Rec_Mutable<M2>& m2)
    { Copy(m1,m2); }


    //
    // Swap Matrices
    //

    template <class M1, class M2> 
    static void Swap(
        BaseMatrix_Diag_Mutable<M1>& m1, BaseMatrix_Diag_Mutable<M2>& m2)
    { Swap(m1.diag(),m2.diag()); }


    //
    // TMV_Text 
    //

    template <class M>
    static std::string TMV_Text(const BaseMatrix_Diag<M>& m)
    {
        std::ostringstream s;
        s << "BaseMatrix_Diag< "<<TMV_Text(m.mat())<<" >";
        return s.str();
    }

    template <class M>
    static std::string TMV_Text(const BaseMatrix_Diag_Mutable<M>& m)
    {
        std::ostringstream s;
        s << "BaseMatrix_Diag_Mutable< "<<TMV_Text(m.mat())<<" >";
        return s.str();
    }

} // namespace tmv

#endif
