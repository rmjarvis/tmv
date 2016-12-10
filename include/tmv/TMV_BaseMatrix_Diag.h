

//---------------------------------------------------------------------------
//
// This file defines the BaseMatrix_Diag and BaseMatrix_Diag_Mutable classes.
//
// See TMV_DiagMatrix.h for the functions that are defined for these objects.
//

#ifndef TMV_BaseMatrix_Diag_H
#define TMV_BaseMatrix_Diag_H

#include "TMV_BaseMatrix.h"

namespace tmv {

    template <class M>
    class BaseMatrix_Diag;
    template <class M>
    class BaseMatrix_Diag_Mutable;

    // Defined in InvertD.h
    template <class M1>
    inline void InvertSelf(BaseMatrix_Diag_Mutable<M1>& m1);

    // Defined in ElemMultVV.h
    template <bool add, int ix, class T, class V1, class V2, class V3>
    inline void ElemMultVV(
        const Scaling<ix,T>& x1, const BaseVector_Calc<V1>& v1,
        const BaseVector_Calc<V2>& v2, BaseVector_Mutable<V3>& v3);

    // Defined below:
    template <class M1, class M2>
    static void Copy(
        const BaseMatrix_Diag<M1>& m1, BaseMatrix_Rec_Mutable<M2>& m2);
    template <class M1, class M2>
    inline void Copy(
        const BaseMatrix_Diag<M1>& m1, BaseMatrix_Diag_Mutable<M2>& m2);

    // A helper class for returning views without necessarily
    // making a new object.
    template <bool ref, class type, class view_type>
    struct MakeDiagView_Helper;

    template <class type, class view_type>
    struct MakeDiagView_Helper<true,type,view_type>
    {
        typedef type& ret_type;
        typedef const type& const_ret_type;
        static TMV_INLINE ret_type call(type& m) { return m; }
        static TMV_INLINE const_ret_type call(const type& m) { return m; }
    };

    template <class type, class view_type>
    struct MakeDiagView_Helper<false,type,view_type>
    {
        typedef view_type ret_type;
        typedef view_type const_ret_type;
        static TMV_INLINE ret_type call(type& m) 
        { return view_type(m.ptr(),m.size(),m.step()); }
        static TMV_INLINE const_ret_type call(const type& m) 
        { return view_type(m.cptr(),m.size(),m.step()); }
    };

    template <class type, class view_type>
    struct MakeDiagView
    {
        enum { ref = Traits2<type,view_type>::sametype };
        typedef MakeDiagView_Helper<ref,type,view_type> helper;

        static TMV_INLINE typename helper::ret_type call(type& m)
        { return helper::call(m); }
        static TMV_INLINE typename helper::const_ret_type call(const type& m)
        { return helper::call(m); }
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

 
    //
    // BaseVector
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
        typedef BaseMatrix_Calc<M> base;

        typedef typename base::calc_type calc_type;
        typedef typename base::eval_type eval_type;
        typedef typename base::copy_type copy_type;
        typedef typename base::inverse_type inverse_type;
        typedef typename base::value_type value_type;
        typedef typename base::real_type real_type;
        typedef typename base::complex_type complex_type;
        typedef typename base::float_type float_type;
        typedef typename base::zfloat_type zfloat_type;

        typedef typename base::const_view_type const_view_type;
        typedef typename base::const_cview_type const_cview_type;
        typedef typename base::const_fview_type const_fview_type;
        typedef typename base::const_xview_type const_xview_type;
        typedef typename base::const_transpose_type const_transpose_type;
        typedef typename base::const_conjugate_type const_conjugate_type;
        typedef typename base::const_adjoint_type const_adjoint_type;
        typedef typename base::const_realpart_type const_realpart_type;
        typedef typename base::const_imagpart_type const_imagpart_type;
        typedef typename base::const_nonconj_type const_nonconj_type;
        typedef typename base::nonconst_type nonconst_type;

        typedef typename Traits<M>::const_diag_type const_diag_type;

        typedef typename Traits<M>::const_unitview_type const_unitview_type;

        typedef typename Traits<M>::const_subdiagmatrix_type 
            const_subdiagmatrix_type;
        typedef typename Traits<M>::const_subdiagmatrix_step_type 
            const_subdiagmatrix_step_type;

        typedef typename Traits<M>::const_iterator const_iterator;
        typedef typename Traits<M>::const_rowmajor_iterator 
            const_rowmajor_iterator;
        typedef typename Traits<M>::const_colmajor_iterator 
            const_colmajor_iterator;

        //
        // Constructor
        //

    protected:
        TMV_INLINE BaseMatrix_Diag() {}
        TMV_INLINE BaseMatrix_Diag(const BaseMatrix_Diag<M>&) {}
        TMV_INLINE ~BaseMatrix_Diag() {}

    private:
        void operator=(const BaseMatrix_Diag<M>&);
    public:


        //
        // Access 
        //

        TMV_INLINE_ND value_type operator()(ptrdiff_t i, ptrdiff_t j) const
        {
            CheckRowIndex<_fort>(i,size());
            CheckColIndex<_fort>(j,size());
            return cref(i,j);
        }

        TMV_INLINE_ND value_type operator()(ptrdiff_t i) const
        {
            CheckIndex<_fort>(i,size());
            return cref(i);
        }

        TMV_INLINE const_diag_type diag() const
        { return const_diag_type(cptr(),size(),step()); }


        //
        // Functions
        //

        TMV_INLINE value_type sumElements() const
        { return diag().sumElements(); }

        TMV_INLINE float_type sumAbsElements() const
        { return diag().sumAbsElements(); }

        TMV_INLINE real_type sumAbs2Elements() const
        { return diag().sumAbs2Elements(); }

        TMV_INLINE float_type maxAbsElement() const
        { return diag().maxAbsElement(); }

        TMV_INLINE real_type maxAbs2Element() const
        { return diag().maxAbs2Element(); }

        TMV_INLINE real_type normSq() const
        { return diag().normSq(); }

        TMV_INLINE float_type normSq(const float_type scale) const
        { return diag().normSq(scale); }

        TMV_INLINE float_type normF() const 
        { return diag().norm2(); }

        TMV_INLINE float_type norm() const
        { return normF(); }

        TMV_INLINE float_type norm1() const
        { return diag().maxAbsElement(); }

        TMV_INLINE float_type norm2() const
        { return diag().maxAbsElement(); }

        TMV_INLINE float_type normInf() const
        { return diag().maxAbsElement(); }

        TMV_INLINE float_type condition() const
        { return diag().maxAbsElement() / diag().minAbsElement(); }

        template <class ret_type, class F>
        TMV_INLINE ret_type sumElements(const F& f) const
        { return diag().sumElements(f); }



        //
        // subDiagMatrix, etc.
        //

        // These versions always uses CStyle
        TMV_INLINE const_subdiagmatrix_type cSubDiagMatrix(ptrdiff_t i1, ptrdiff_t i2) const
        {
            return const_subdiagmatrix_type(
                cptr()+i1*step(), i2-i1, step());
        }

        TMV_INLINE const_subdiagmatrix_step_type cSubDiagMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t istep) const
        {
            return const_subdiagmatrix_step_type(
                cptr()+i1*step(), (i2-i1)/istep, istep*step());
        }


        // These check the indices according the the indexing style being
        // used, and then calls the above CStyle versions.
        TMV_INLINE_ND const_subdiagmatrix_type subDiagMatrix(
            ptrdiff_t i1, ptrdiff_t i2) const
        {
            CheckRange<_fort>(i1,i2,size());
            return cSubDiagMatrix(i1,i2);
        }

        TMV_INLINE_ND const_subdiagmatrix_step_type subDiagMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t istep) const
        {
            CheckRange<_fort>(i1,i2,istep,size());
            return cSubDiagMatrix(i1,i2,istep);
        }


        //
        // Views
        //

        TMV_INLINE TMV_MAYBE_CREF(type,const_view_type) view() const
        { return MakeDiagView<type,const_view_type>::call(mat()); }

        TMV_INLINE TMV_MAYBE_CREF(type,const_cview_type) cView() const
        { return MakeDiagView<type,const_cview_type>::call(mat()); }

        TMV_INLINE TMV_MAYBE_CREF(type,const_fview_type) fView() const
        { return MakeDiagView<type,const_fview_type>::call(mat()); }

        TMV_INLINE TMV_MAYBE_CREF(type,const_xview_type) xView() const
        { return MakeDiagView<type,const_xview_type>::call(mat()); }

        TMV_INLINE_ND TMV_MAYBE_CREF(type,const_unitview_type) unitView() const
        {
            TMVAssert(step()==1&&"Called unitView on DiagMatrix with step!=1");
            return MakeDiagView<type,const_unitview_type>::call(mat()); 
        }

        TMV_INLINE TMV_MAYBE_CREF(type,const_view_type) constView() const
        { return MakeDiagView<type,const_view_type>::call(mat()); }

        TMV_INLINE TMV_MAYBE_CREF(type,const_transpose_type) transpose() const
        { return MakeDiagView<type,const_transpose_type>::call(mat()); }

        TMV_INLINE TMV_MAYBE_CREF(type,const_conjugate_type) conjugate() const
        { return MakeDiagView<type,const_conjugate_type>::call(mat()); }

        TMV_INLINE TMV_MAYBE_CREF(type,const_adjoint_type) adjoint() const
        { return MakeDiagView<type,const_adjoint_type>::call(mat()); }

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

        TMV_INLINE TMV_MAYBE_CREF(type,const_nonconj_type) nonConj() const
        { return MakeDiagView<type,const_nonconj_type>::call(mat()); }

        TMV_INLINE nonconst_type nonConst() const
        {
            return nonconst_type(
                const_cast<value_type*>(cptr()),size(),step()); 
        }

        //
        // Auxilliary routines
        //

        template <class M2>
        TMV_INLINE void assignTo(BaseMatrix_Mutable<M2>& m2) const
        {
            TMVStaticAssert((ShapeTraits2<_shape,M2::_shape>::assignable));
            tmv::Copy(mat(),m2.mat()); 
        }

        TMV_INLINE const type& mat() const
        { return static_cast<const type&>(*this); }

        TMV_INLINE bool isconj() const { return _conj; }

        // Note that these last functions need to be defined in a more derived
        // class than this, or an infinite loop will result.
        // Also, cref, get_row and get_col from BaseMatrix.

        TMV_INLINE ptrdiff_t colsize() const { return mat().size(); }
        TMV_INLINE ptrdiff_t rowsize() const { return mat().size(); }
        TMV_INLINE ptrdiff_t size() const { return mat().size(); }
        TMV_INLINE ptrdiff_t nlo() const { return 0; }
        TMV_INLINE ptrdiff_t nhi() const { return 0; }
        TMV_INLINE ptrdiff_t step() const { return mat().step(); }

        TMV_INLINE ptrdiff_t rowstart(ptrdiff_t i) const { return i; }
        TMV_INLINE ptrdiff_t rowend(ptrdiff_t i) const { return i+1; }
        TMV_INLINE ptrdiff_t colstart(ptrdiff_t j) const { return j; }
        TMV_INLINE ptrdiff_t colend(ptrdiff_t j) const { return j+1; }

        TMV_INLINE const value_type* cptr() const { return mat().cptr(); }
        value_type cref(ptrdiff_t i, ptrdiff_t j) const 
        { return (i!=j ? value_type(0) : cref(i)); }
        TMV_INLINE value_type cref(ptrdiff_t i) const 
        { return mat().cref(i); }

        TMV_INLINE const_iterator begin() const
        { return diag().begin(); }
        TMV_INLINE const_iterator end() const
        { return diag().end(); }

        TMV_INLINE const_rowmajor_iterator rowmajor_begin() const
        { return begin(); }
        TMV_INLINE const_rowmajor_iterator rowmajor_end() const
        { return end(); }

        TMV_INLINE const_colmajor_iterator colmajor_begin() const
        { return begin(); }
        TMV_INLINE const_colmajor_iterator colmajor_end() const
        { return end(); }

    };

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
        typedef BaseMatrix_Diag<M> base;
        typedef BaseMatrix_Mutable<M> base_mut;

        typedef typename base::calc_type calc_type;
        typedef typename base::eval_type eval_type;
        typedef typename base::copy_type copy_type;
        typedef typename base::inverse_type inverse_type;
        typedef typename base::value_type value_type;
        typedef typename base::real_type real_type;
        typedef typename base::complex_type complex_type;
        typedef typename base::float_type float_type;
        typedef typename base::zfloat_type zfloat_type;

        typedef typename base::const_view_type const_view_type;
        typedef typename base::const_cview_type const_cview_type;
        typedef typename base::const_fview_type const_fview_type;
        typedef typename base::const_xview_type const_xview_type;
        typedef typename base::const_transpose_type const_transpose_type;
        typedef typename base::const_conjugate_type const_conjugate_type;
        typedef typename base::const_adjoint_type const_adjoint_type;
        typedef typename base::const_realpart_type const_realpart_type;
        typedef typename base::const_imagpart_type const_imagpart_type;
        typedef typename base::const_nonconj_type const_nonconj_type;
        typedef typename base::nonconst_type nonconst_type;

        typedef typename base_mut::view_type view_type;
        typedef typename base_mut::cview_type cview_type;
        typedef typename base_mut::fview_type fview_type;
        typedef typename base_mut::xview_type xview_type;
        typedef typename base_mut::transpose_type transpose_type;
        typedef typename base_mut::conjugate_type conjugate_type;
        typedef typename base_mut::adjoint_type adjoint_type;
        typedef typename base_mut::realpart_type realpart_type;
        typedef typename base_mut::imagpart_type imagpart_type;
        typedef typename base_mut::nonconj_type nonconj_type;
        typedef typename base_mut::noalias_type noalias_type;
        typedef typename base_mut::alias_type alias_type;
        typedef typename base_mut::reference reference;

        typedef typename base::const_diag_type const_diag_type;
        typedef typename base::const_unitview_type const_unitview_type;
        typedef typename base::const_subdiagmatrix_type 
            const_subdiagmatrix_type;
        typedef typename base::const_subdiagmatrix_step_type 
            const_subdiagmatrix_step_type;

        typedef typename base::const_iterator const_iterator;
        typedef typename base::const_rowmajor_iterator const_rowmajor_iterator;
        typedef typename base::const_colmajor_iterator const_colmajor_iterator;

        typedef typename Traits<M>::diag_type diag_type;

        typedef typename Traits<M>::unitview_type unitview_type;

        typedef typename Traits<M>::subdiagmatrix_type subdiagmatrix_type;
        typedef typename Traits<M>::subdiagmatrix_step_type 
            subdiagmatrix_step_type;

        typedef typename Traits<M>::iterator iterator;
        typedef typename Traits<M>::rowmajor_iterator rowmajor_iterator;
        typedef typename Traits<M>::colmajor_iterator colmajor_iterator;


        //
        //
        // Constructor
        //

    protected:
        TMV_INLINE BaseMatrix_Diag_Mutable() {}
        TMV_INLINE BaseMatrix_Diag_Mutable(const BaseMatrix_Diag_Mutable<M>&) {}
        TMV_INLINE ~BaseMatrix_Diag_Mutable() {}
    public:


        //
        // Access 
        //

        TMV_INLINE_ND reference operator()(ptrdiff_t i, ptrdiff_t j)
        {
            TMVAssert(i == j);
            CheckIndex<_fort>(i,size());
            return ref(i);
        }

        TMV_INLINE_ND reference operator()(ptrdiff_t i)
        {
            CheckIndex<_fort>(i,size());
            return ref(i);
        }

        TMV_INLINE diag_type diag() 
        { return diag_type(ptr(),size(),step()); }


        // We need to repeat the const versions so the non-const ones
        // don't clobber them.
        TMV_INLINE value_type operator()(ptrdiff_t i, ptrdiff_t j) const
        { return base::operator()(i,j); }
        TMV_INLINE const_diag_type diag() const
        { return base::diag(); }


        //
        // Op =
        //

        TMV_INLINE_ND type& operator=(const BaseMatrix_Diag_Mutable<M>& m2) 
        {
            TMVAssert(colsize() == m2.colsize());
            TMVAssert(rowsize() == m2.rowsize());
            m2.assignTo(mat());
            return mat(); 
        }

        template <class M2>
        TMV_INLINE_ND type& operator=(const BaseMatrix<M2>& m2) 
        {
            TMVStaticAssert((Sizes<_colsize,M2::_colsize>::same));
            TMVStaticAssert((Sizes<_rowsize,M2::_rowsize>::same));
            TMVAssert(colsize() == m2.colsize());
            TMVAssert(rowsize() == m2.rowsize());
            m2.assignTo(mat());
            return mat(); 
        }

        TMV_INLINE_ND type& operator=(const value_type x)
        {
            TMVStaticAssert((Sizes<_rowsize,_colsize>::same));
            TMVAssert(colsize() == rowsize());
            setToIdentity(x);
            return mat();
        }


        //
        // Modifying Functions
        //

        TMV_INLINE type& setZero() 
        { return setAllTo(value_type(0)); }

        TMV_INLINE type& setAllTo(value_type x) 
        { diag().setAllTo(x); return mat(); }

        TMV_INLINE type& addToAll(value_type x) 
        { diag().addToAll(x); return mat(); }

        TMV_INLINE type& clip(float_type thresh) 
        { diag().clip(thresh); return mat(); }

        template <class F>
        TMV_INLINE type& applyToAll(const F& f)
        { diag().applyToAll(f); return mat(); }

        TMV_INLINE type& conjugateSelf() 
        { diag().conjugateSelf(); return mat(); }

        TMV_INLINE type& transposeSelf() 
        { return mat(); }

        TMV_INLINE type& invertSelf() 
        { tmv::InvertSelf(mat()); return mat(); }

        TMV_INLINE type& setToIdentity(const value_type x=value_type(1))
        { diag().setAllTo(x); return mat(); }

        TMV_INLINE_ND type& swap(ptrdiff_t i1, ptrdiff_t i2) 
        {
            CheckIndex<_fort>(i1,size());
            CheckIndex<_fort>(i2,size());
            diag().cSwap(i1,i2);
            return mat();
        }

        TMV_INLINE type& cPermute(const ptrdiff_t* p, ptrdiff_t i1, ptrdiff_t i2) 
        { diag().permute(p,i1,i2); return mat(); }
        TMV_INLINE_ND type& permute(const ptrdiff_t* p, ptrdiff_t i1, ptrdiff_t i2) 
        {
            CheckRange<_fort>(i1,i2,size());
            return cPermute(p,i1,i2);
        }
        TMV_INLINE type& permute(const ptrdiff_t* p) 
        { cPermute(p,0,size()); return mat(); }

        TMV_INLINE type& cReversePermute(const ptrdiff_t* p, ptrdiff_t i1, ptrdiff_t i2) 
        { diag().reversePermute(p,i1,i2); }
        TMV_INLINE_ND type& reversePermute(const ptrdiff_t* p, ptrdiff_t i1, ptrdiff_t i2) 
        {
            CheckRange<_fort>(i1,i2,size());
            return cReversePermute(p,i1,i2);
        }
        TMV_INLINE type& reversePermute(const ptrdiff_t* p) 
        { cReversePermute(p,0,size()); return mat(); }


        //
        // subDiagMatrix, etc.
        //

        // These versions always uses CStyle
        TMV_INLINE subdiagmatrix_type cSubDiagMatrix(ptrdiff_t i1, ptrdiff_t i2) 
        { return subdiagmatrix_type(ptr()+i1*step(), i2-i1, step()); }

        TMV_INLINE subdiagmatrix_step_type cSubDiagMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t istep) 
        {
            return subdiagmatrix_step_type(
                ptr()+i1*step(), (i2-i1)/istep, istep*step());
        }

        // These check the indices according the the indexing style being
        // used, and then calls the above CStyle versions.
        TMV_INLINE_ND subdiagmatrix_type subDiagMatrix(ptrdiff_t i1, ptrdiff_t i2) 
        {
            CheckRange<_fort>(i1,i2,size());
            return cSubDiagMatrix(i1,i2);
        }

        TMV_INLINE_ND subdiagmatrix_step_type subDiagMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t istep) 
        {
            CheckRange<_fort>(i1,i2,istep,size());
            return cSubDiagMatrix(i1,i2,istep);
        }


        // Repeat the const versions:
        TMV_INLINE const_subdiagmatrix_type cSubDiagMatrix(ptrdiff_t i1, ptrdiff_t i2) const
        { return base::cSubDiagMatrix(i1,i2); }
        TMV_INLINE const_subdiagmatrix_step_type cSubDiagMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t istep) const
        { return base::cSubDiagMatrix(i1,i2,istep); }

        TMV_INLINE const_subdiagmatrix_type subDiagMatrix(ptrdiff_t i1, ptrdiff_t i2) const
        { return base::subDiagMatrix(i1,i2); }
        TMV_INLINE const_subdiagmatrix_step_type subDiagMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t istep) const
        { return base::subDiagMatrix(i1,i2,istep); }


        //
        // Views
        //

        TMV_INLINE TMV_MAYBE_REF(type,view_type) view() 
        { return MakeDiagView<type,view_type>::call(mat()); }

        TMV_INLINE TMV_MAYBE_REF(type,cview_type) cView() 
        { return MakeDiagView<type,cview_type>::call(mat()); }

        TMV_INLINE TMV_MAYBE_REF(type,fview_type) fView() 
        { return MakeDiagView<type,fview_type>::call(mat()); }

        TMV_INLINE TMV_MAYBE_REF(type,xview_type) xView() 
        { return MakeDiagView<type,xview_type>::call(mat()); }

        TMV_INLINE_ND TMV_MAYBE_REF(type,unitview_type) unitView() 
        {
            TMVAssert(step()==1&&"Called unitView on DiagMatrix with step!=1");
            return MakeDiagView<type,unitview_type>::call(mat()); 
        }

        TMV_INLINE TMV_MAYBE_REF(type,transpose_type) transpose()
        { return MakeDiagView<type,transpose_type>::call(mat()); }

        TMV_INLINE TMV_MAYBE_REF(type,conjugate_type) conjugate() 
        { return MakeDiagView<type,conjugate_type>::call(mat()); }

        TMV_INLINE TMV_MAYBE_REF(type,adjoint_type) adjoint()
        { return MakeDiagView<type,adjoint_type>::call(mat()); }

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

        TMV_INLINE TMV_MAYBE_REF(type,nonconj_type) nonConj()
        { return MakeDiagView<type,nonconj_type>::call(mat()); }

        TMV_INLINE TMV_MAYBE_REF(type,noalias_type) noAlias()
        { return MakeDiagView<type,noalias_type>::call(mat()); }

        TMV_INLINE TMV_MAYBE_REF(type,alias_type) alias()
        { return MakeDiagView<type,alias_type>::call(mat()); }


        // Repeat the const versions:
        TMV_INLINE TMV_MAYBE_CREF(type,const_view_type) view() const
        { return base::view(); }
        TMV_INLINE TMV_MAYBE_CREF(type,const_cview_type) cView() const
        { return base::cView(); }
        TMV_INLINE TMV_MAYBE_CREF(type,const_fview_type) fView() const
        { return base::fView(); }
        TMV_INLINE TMV_MAYBE_CREF(type,const_xview_type) xView() const
        { return base::xView(); }
        TMV_INLINE TMV_MAYBE_CREF(type,const_unitview_type) unitView() const
        { return base::unitView(); }
        TMV_INLINE TMV_MAYBE_CREF(type,const_transpose_type) transpose() const
        { return base::transpose(); }
        TMV_INLINE TMV_MAYBE_CREF(type,const_conjugate_type) conjugate() const
        { return base::conjugate(); }
        TMV_INLINE TMV_MAYBE_CREF(type,const_adjoint_type) adjoint() const
        { return base::adjoint(); }
        TMV_INLINE const_realpart_type realPart() const
        { return base::realPart(); }
        TMV_INLINE const_imagpart_type imagPart() const
        { return base::imagPart(); }
        TMV_INLINE TMV_MAYBE_CREF(type,const_nonconj_type) nonConj() const
        { return base::nonConj(); }

        //
        // Auxilliary routines
        //

        TMV_INLINE const type& mat() const
        { return static_cast<const type&>(*this); }
        TMV_INLINE type& mat()
        { return static_cast<type&>(*this); }


        // Note that these last functions need to be defined in a more derived
        // class than this, or an infinite loop will result when compiling.
        // Also, cref and cptr from above.

        TMV_INLINE ptrdiff_t colsize() const { return mat().size(); }
        TMV_INLINE ptrdiff_t rowsize() const { return mat().size(); }
        TMV_INLINE ptrdiff_t size() const { return mat().size(); }
        TMV_INLINE ptrdiff_t step() const { return mat().step(); }

        TMV_INLINE value_type* ptr() { return mat().ptr(); }
        TMV_INLINE reference ref(ptrdiff_t i, ptrdiff_t j) { return mat().ref(i); }
        TMV_INLINE reference ref(ptrdiff_t i) { return mat().ref(i); }
        TMV_INLINE value_type cref(ptrdiff_t i) { return mat().cref(i); }

        TMV_INLINE iterator begin() 
        { return diag().begin(); }
        TMV_INLINE iterator end() 
        { return diag().end(); }

        TMV_INLINE const_iterator begin() const
        { return diag().begin(); }
        TMV_INLINE const_iterator end() const
        { return diag().end(); }

        TMV_INLINE rowmajor_iterator rowmajor_begin() 
        { return begin(); }
        TMV_INLINE rowmajor_iterator rowmajor_end() 
        { return end(); }

        TMV_INLINE colmajor_iterator colmajor_begin() 
        { return begin(); }
        TMV_INLINE colmajor_iterator colmajor_end() 
        { return end(); }

        TMV_INLINE const_rowmajor_iterator rowmajor_begin() const
        { return begin(); }
        TMV_INLINE const_rowmajor_iterator rowmajor_end() const
        { return end(); }

        TMV_INLINE const_colmajor_iterator colmajor_begin() const
        { return begin(); }
        TMV_INLINE const_colmajor_iterator colmajor_end() const
        { return end(); }

    }; // BaseMatrix_Diag_Mutable


    template <class T, int A=0>
    class DiagMatrix;
    template <class T, int A=0>
    class ConstDiagMatrixView;
    template <class T, int A=0>
    class DiagMatrixView;
    template <class T, ptrdiff_t N, int A=0>
    class SmallDiagMatrix;
    template <class T, ptrdiff_t N, ptrdiff_t S=Unknown, int A=0>
    class ConstSmallDiagMatrixView;
    template <class T, ptrdiff_t N, ptrdiff_t S=Unknown, int A=0>
    class SmallDiagMatrixView;

    // This helper class helps decide calc_type for composite classes:
    template <class T, ptrdiff_t cs, ptrdiff_t rs, int A>
    struct MCopyHelper<T,Diag,cs,rs,A>
    { typedef SmallDiagMatrix<T,cs,A> type; };
    template <class T, ptrdiff_t rs, int A>
    struct MCopyHelper<T,Diag,Unknown,rs,A>
    { typedef SmallDiagMatrix<T,rs,A> type; };
    template <class T, ptrdiff_t cs, int A>
    struct MCopyHelper<T,Diag,cs,Unknown,A>
    { typedef SmallDiagMatrix<T,cs,A> type; };
    template <class T, int A>
    struct MCopyHelper<T,Diag,Unknown,Unknown,A>
    { typedef DiagMatrix<T,A|NoAlias> type; };

    template <class T, ptrdiff_t cs, ptrdiff_t rs, ptrdiff_t si, ptrdiff_t sj, int A>
    struct MViewHelper<T,Diag,cs,rs,si,sj,A>
    { 
        typedef SmallDiagMatrixView<T,cs,si,A> type; 
        typedef ConstSmallDiagMatrixView<T,cs,si,A> ctype; 
    };
    template <class T, ptrdiff_t si, ptrdiff_t sj, int A>
    struct MViewHelper<T,Diag,Unknown,Unknown,si,sj,A>
    {
        enum { A2 = A | (si == 1 ? Unit : NonUnit) | NoAlias };
        typedef DiagMatrixView<T,A2> type; 
        typedef ConstDiagMatrixView<T,A2> ctype; 
    };


    //
    // Copy Matrices
    //

    template <class M1, class M2>
    inline void Copy(
        const BaseMatrix_Diag<M1>& m1, BaseMatrix_Diag_Mutable<M2>& m2)
    {
        typename M1::const_diag_type m1d = m1.diag();
        typename M2::diag_type m2d = m2.diag();
        Copy(m1d,m2d);
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
            Copy(m1d,m2d);
            if (m1.size() > 1) {
                m2.upperTri().offDiag().setZero();
                m2.lowerTri().offDiag().setZero();
            }
        } else {
            m2.setZero();
            m2d.noAlias() = m1d;
        }
    }


    //
    // Swap Matrices
    //

    template <class M1, class M2>
    TMV_INLINE void Swap(
        BaseMatrix_Diag_Mutable<M1>& m1, BaseMatrix_Diag_Mutable<M2>& m2)
    { Swap(m1.diag(),m2.diag()); }


    //
    // TMV_Text 
    //

    template <class M>
    inline std::string TMV_Text(const BaseMatrix_Diag<M>& m)
    {
        std::ostringstream s;
        s << "BaseMatrix_Diag< "<<TMV_Text(m.mat())<<" >";
        return s.str();
    }

    template <class M>
    inline std::string TMV_Text(const BaseMatrix_Diag_Mutable<M>& m)
    {
        std::ostringstream s;
        s << "BaseMatrix_Diag_Mutable< "<<TMV_Text(m.mat())<<" >";
        return s.str();
    }

} // namespace tmv

#endif
