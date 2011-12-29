
//-----------------------------------------------------------------------------
//
// This file defines the BaseMatrix, BaseMatrix_Calc, and BaseMatrix_Mutable
// classes.  These define the functionality for matrices of all the
// different shapes, including dense rectangular, upper/lower triangular,
// diagonal, banded, symmetric, hermitian, band-symmetric and band-hermitian.
//
// BaseMatrix is the base class for all of the various types of matrices.
// It defines all of the functions that you can use on any matrix.
//
// There is a single template argument, which is the type of the most
// derived class for the object.
// So, for example, Matrix<double> inherits from BaseMatrix<Matrix<double> >.
//
// The methods that are defined for BaseMatrix are described in 
// TMV_Matrix.h and include all of the constant access methods,
// the non-modifying functions, write, and arithmetic operators.
//
//
// BaseMatrix_Calc is the base class for all matrices that have their
// data calculated in memory, possibly being the conjugate of the 
// underlying data.
//
// The methods that are defined for BaseMatrix_Calc are described in
// TMV_Matrix.h and include row, col, diag, transpose, conjugate, 
// and adjoint (among others).
//
//
// BaseMatrix_Mutable is the base class for all matrices with calculated
// data (like BaseMatrix_Calc) whose values are allowed to be modified.
// 
// The methods that are defined for BaseMatrix_Mutable are described in
// TMV_Matrix.h and include the non-const versions of transpose, conjugate,
// etc. and some functions such as setZero(), setAllTo(), and addToAll()
// that work for any shape matrix.  This class is also the type of the
// argument of assignTo(). 
//
//
#ifndef TMV_BaseMatrix_H
#define TMV_BaseMatrix_H

#include "TMV_BaseVector.h"

namespace tmv {

    // BaseMatrix is the base class for all matrix classes.
    // All non-modifying functions are defined for BaseMatrix, such as
    // norm, trace, det, etc.
    template <class M>
    class BaseMatrix;

    // BaseMatrix gets a lot of the relevant information about M 
    // from the Traits class: Traits<M>
    // The types are defined as typedef statements.
    // The integer (or boolean) values are defined as enum statements.
    // See TMV_Matrix.h or TMV_SmallMatrix.h for some concrete examples.
    //
    //
    //  value_type = The type of the individual elements
    //
    //  type = shorthand for the derived type 
    //
    //  calc_type = The type of the calculated version of the matrix.
    //  i.e. where all of the values are stored in memory somewhere.
    //  This is the return type of calc()
    //  Use this when you will access elements of a composite matrix 
    //  multiple times or need all of the elements.
    //  Also use this when you will use cptr().
    //
    //  eval_type = The type of the evaluated version of the matrix
    //  This is the return type of eval()
    //  Use this when you will only access a few elements of a 
    //  composite matrix one time each.
    //
    //  copy_type = The type of a new copy of the matrix
    //  This is the return type of copy()
    //
    //  inverse_type = The type of the inverse of the matrix
    //  This is the return type of inverse()
    //
    //  _colsize = column size of matrix (aka number of rows)
    //  _rowsize = row size of matrix (aka number of columns)
    //  (Use TMV_UNKNOWN if unknown at compile time)
    //
    //  _nlo = number of non-zero subdiagonals
    //  _nhi = number of non-zero superdiagonals
    //
    //  _shape = The shape of the non-zero elements of the matrix
    //
    //  _fort = does the indexing use fortran style?
    //
    //  _calc = are the element values already calculated in memory?


    // BaseMatrix_Calc is derived from BaseMatrix, and is used
    // for matrices that have their values already calculated somewhere.
    // So composite classes inherit directly from BaseMatrix, rather
    // than BaseMatrix_Calc.
    template <class M>
    class BaseMatrix_Calc;

    // BaseMatrix_Calc adds some more requirements to the Traits<M> class:
    //
    //  _conj = is the matrix the conjugate of the underlying data?
    //  _checkalias = do we need to check this matrix for aliases?
    //  _rowmajor = is the matrix RowMajor?
    //  _colmajor = is the matrix ColMajor?
    //  _hasdivider = does this matrix have a divider object?
    //
    //  const_view_type = return type from view() const
    //  const_cview_type = return type from cView() const
    //  const_fview_type = return type from fView() const
    //  const_xview_type = return type from xView() const
    //
    //  const_conjugate_type = return type from conjugate() const
    //  const_transpose_type = return type from transpose() const
    //  const_adjoint_type = return type from adjoint() const
    //  const_realpart_type = return type from realPart() const
    //  const_imagpart_type = return type from imagPart() const
    //  const_nonconj_type = return type from nonConj() const
    //  nonconst_type = return type from nonConst() const
    //
    //  const_rowmajor_iterator = return type from rowmajor_begin() const,
    //                            rowmajor_end() const
    //  const_colmajor_iterator = return type from colmajor_begin() const,
    //                            colmajor_end() const


    // BaseMatrix_Mutable is used for matrices that are allowed to have 
    // their data modified.
    // It would naturally derive from BaseMatrix_Calc, except that we want
    // to avoid having to deal with diamond inheritance patterns.
    // For example BaseMatrix_Diag derives from BaseMatrix_Calc, but not
    // from BaseMatrix_Mutable, since it represents an immutable DiagMatrix.
    // But then BaseMatrix_Diag_Mutable derives from both BaseMatrix_Diag
    // and BaseMatrix_Mutable to represent a mutable DiagMatrix.
    // So, to avoid dealing with the sticky issues of diamond inheritance,
    // we do not derive BaseMatrix_Mutable from BaseMatrix_Calc.
    template <class M>
    class BaseMatrix_Mutable;

    // BaseMatrix_Mutable adds:
    //
    //  reference = return type of m(i,j)
    //
    //  view_type = return type from view() 
    //  cview_type = return type from cView()
    //  fview_type = return type from fView()
    //  xview_type = return type from xView()
    //
    //  conjugate_type = return type from conjugate() 
    //  transpose_type = return type from transpose()
    //  adjoint_type = return type from adjoint()
    //  realpart_type = return type from realPart()
    //  imagpart_type = return type from imagPart()
    //  nonconj_type = return type from nonConj()
    //  noalias_type = return type from noAlias()
    //  alias_type = return type from alias()
    //
    //  rowmajor_iterator = return type from rowmajor_begin(), rowmajor_end()
    //  colmajor_iterator = return type from colmajor_begin(), colmajor_end()


    // Other BaseMatrxi varieties declared here without definition.
    template <class M> class BaseMatrix_Rec;
    template <class M> class BaseMatrix_Rec_Mutable;
    template <class M> class BaseMatrix_Diag;
    template <class M> class BaseMatrix_Diag_Mutable;
    template <class M> class BaseMatrix_Tri;
    template <class M> class BaseMatrix_Tri_Mutable;
    template <class M> class BaseMatrix_Band;
    template <class M> class BaseMatrix_Band_Mutable;
    template <class M> class BaseMatrix_Sym;
    template <class M> class BaseMatrix_Sym_Mutable;
    template <class M> class BaseMatrix_SymBand;
    template <class M> class BaseMatrix_SymBand_Mutable;

    //
    // Helper functions and values:
    //

    // These helper functions check the validity of indices according
    // to whether the matrix uses CStyle or FortranStyle indexing.
    // They also update the indices to be consistent with CStyle.
    template <bool _fort>
    TMV_INLINE_ND void CheckRowIndex(int& i, int m)
    { // CStyle
        TMVAssert(i >= 0 && "row index must be in matrix");
        TMVAssert(i < m && "row index must be in matrix");
    }
    template <bool _fort>
    TMV_INLINE_ND void CheckColIndex(int& j, int n)
    { // CStyle
        TMVAssert(j >= 0 && "column index must be in matrix");
        TMVAssert(j < n && "column index must be in matrix");
    }
    template <>
    TMV_INLINE_ND void CheckRowIndex<true>(int& i, int m)
    { // FortranStyle
        TMVAssert(i >= 1 && "row index must be in matrix");
        TMVAssert(i <= m && "row index must be in matrix");
        --i;
    }
    template <>
    TMV_INLINE_ND void CheckColIndex<true>(int& j, int n)
    { // FortranStyle
        TMVAssert(j >= 1 && "column index must be in matrix");
        TMVAssert(j <= n && "column index must be in matrix");
        --j;
    }

    // Override SameStorage for Matrix objects:
    template <class M1, class M2>
    TMV_INLINE bool SameStorage(
        const BaseMatrix<M1>& v1, const BaseMatrix<M2>& m2)
    { return false; }
    template <class V1, class M2>
    TMV_INLINE bool SameStorage(
        const BaseVector<V1>& v1, const BaseMatrix<M2>& m2)
    { return false; }
    template <class M1, class V2>
    TMV_INLINE bool SameStorage(
        const BaseMatrix<M1>& v1, const BaseVector<V2>& m2)
    { return false; }
#ifndef TMV_NO_ALIAS_CHECK
    template <class V1, class M2>
    TMV_INLINE bool SameStorage(
        const BaseVector_Calc<V1>& v1, const BaseMatrix_Calc<M2>& m2)
    {
        return 
            static_cast<const void*>(v1.vec().cptr()) == 
            static_cast<const void*>(m2.mat().cptr()); 
    }
    template <class M1, class V2>
    TMV_INLINE bool SameStorage(
        const BaseMatrix_Calc<M1>& m1, const BaseVector_Calc<V2>& v2)
    {
        return 
            static_cast<const void*>(m1.mat().cptr()) == 
            static_cast<const void*>(v2.vec().cptr()); 
    }
    template <class M1, class M2>
    TMV_INLINE bool SameStorage(
        const BaseMatrix_Calc<M1>& m1, const BaseMatrix_Calc<M2>& m2)
    {
        return 
            static_cast<const void*>(m1.mat().cptr()) == 
            static_cast<const void*>(m2.mat().cptr()); 
    }
#endif

    template <class M1, class M2>
    TMV_INLINE bool ExactSameStorage(
        const BaseMatrix<M1>& m1, const BaseMatrix<M2>& m2)
    { return false; }
    template <class M1, class M2>
    TMV_INLINE bool OppositeStorage(
        const BaseMatrix<M1>& m1, const BaseMatrix<M2>& m2)
    { return false; }
    // Step checks are done for specific shapes in the various files
    // for each shaped matrix.

    // This helper class determines if there is ExactSameStorage
    // or OppositeStorage for two matrices at compile time:
    template <class M1, class M2>
    struct MStepHelper
    {
        typedef typename M1::value_type T1;
        typedef typename M2::value_type T2;
        enum { known = (
                M1::_stepi != TMV_UNKNOWN &&
                M1::_stepj != TMV_UNKNOWN &&
                M2::_stepi != TMV_UNKNOWN &&
                M2::_stepj != TMV_UNKNOWN ) };
        enum { same = (
                Traits2<T1,T2>::sametype &&
                known &&
                M1::_stepi == int(M2::_stepi) &&
                M1::_stepj == int(M2::_stepj) ) };
        enum { opp = (
                Traits2<T1,T2>::sametype &&
                known &&
                M1::_stepi == int(M2::_stepj) &&
                M1::_stepj == int(M2::_stepi) ) };
    };

    // This helper class helps decide calc_type for composite classes:
    // We don't define anything, since it needs to be specialized
    // differently for each shape.
    template <class T, int shape, int cs, int rs, int A=0>
    struct MCopyHelper;

    // This is similar - it defines the right view type when some
    // sizes or steps might be known.
    template <class T, int shape, int cs, int rs, int si, int sj, int A=0>
    struct MViewHelper;

    template <class M>
    inline typename M::value_type DoTrace(const BaseMatrix<M>& m);

    // Defined in TMV_QuotXM.h
    template <int ix, class T, class M>
    class QuotXM;

    // Defined in TMV_Det.h
    template <class M>
    inline typename M::value_type Det(const BaseMatrix_Calc<M>& m);
    template <class M>
    inline typename M::float_type LogDet(
        const BaseMatrix_Calc<M>& m, typename M::zfloat_type* sign=0);
    template <class M>
    inline bool IsSingular(const BaseMatrix_Calc<M>& m);

    // Defined in TMV_InvertM.h
    template <int ix, class T, class M1, class M2>
    inline void MakeInverse(
        const Scaling<ix,T>& x, const BaseMatrix_Calc<M1>& m1,
        BaseMatrix_Mutable<M2>& minv);

    template <class M1, class M2>
    inline void MakeInverseATA(
        const BaseMatrix_Calc<M1>& m1, BaseMatrix_Mutable<M2>& mata);


    template <class M>
    class BaseMatrix
    {
    public:
        enum { _colsize = Traits<M>::_colsize };
        enum { _rowsize = Traits<M>::_rowsize };
        enum { _nlo = Traits<M>::_nlo };
        enum { _nhi = Traits<M>::_nhi };
        enum { _shape = Traits<M>::_shape };
        enum { _fort = Traits<M>::_fort };
        enum { _calc = Traits<M>::_calc };

        typedef M type;
        typedef typename Traits<M>::calc_type calc_type;
        typedef typename Traits<M>::eval_type eval_type;
        typedef typename Traits<M>::copy_type copy_type;
        typedef typename Traits<M>::inverse_type inverse_type;

        typedef typename Traits<M>::value_type value_type;
        typedef typename Traits<value_type>::real_type real_type;
        typedef typename Traits<value_type>::complex_type complex_type;
        typedef typename Traits<real_type>::float_type float_type;
        typedef typename Traits<value_type>::float_type zfloat_type;

        enum { isreal = Traits<value_type>::isreal };
        enum { iscomplex = Traits<value_type>::iscomplex };

        //
        // Constructor
        //

    protected:
        TMV_INLINE BaseMatrix() {}
        TMV_INLINE BaseMatrix(const BaseMatrix<M>&) {}
        TMV_INLINE ~BaseMatrix() {}

    private :
        void operator=(const BaseMatrix<M>& m2);
    public :


        //
        // Access
        //

        TMV_INLINE value_type operator()(int i, int j) const 
        {
            CheckRowIndex<_fort>(i,colsize());
            CheckColIndex<_fort>(j,rowsize());
            return cref(i,j);
        }


        //
        // Functions
        //

        TMV_INLINE value_type trace() const
        { return tmv::DoTrace(eval()); }

        TMV_INLINE value_type sumElements() const
        { return calc().sumElements(); }

        TMV_INLINE float_type sumAbsElements() const
        { return calc().sumAbsElements(); }

        TMV_INLINE real_type sumAbs2Elements() const
        { return calc().sumAbs2Elements(); }

        TMV_INLINE float_type maxAbsElement() const 
        { return calc().maxAbsElement(); }

        TMV_INLINE real_type maxAbs2Element() const 
        { return calc().maxAbs2Element(); }

        TMV_INLINE real_type normSq() const
        { return calc().normSq(); }

        TMV_INLINE float_type normSq(const float_type scale) const
        { return calc().normSq(scale); }

        TMV_INLINE float_type normF() const
        { return calc().normF(); }

        TMV_INLINE float_type norm() const
        { return normF(); }

        TMV_INLINE float_type norm1() const 
        { return calc().norm1(); }

        TMV_INLINE float_type normInf() const 
        { return calc().normInf(); }

        TMV_INLINE real_type norm2() const
        { return calc().norm2(); }

        TMV_INLINE real_type condition() const
        { return calc().condition(); }

        template <class ret_type, class F>
        TMV_INLINE ret_type sumElements(const F& f) const
        { return calc().sumElements(f); }


        //
        // Auxilliary routines
        //

        TMV_INLINE const type& mat() const 
        { return static_cast<const type&>(*this); }

        TMV_INLINE calc_type calc() const 
        { return static_cast<calc_type>(mat()); }

        TMV_INLINE eval_type eval() const 
        { return static_cast<eval_type>(mat()); }

        TMV_INLINE copy_type copy() const 
        { return static_cast<copy_type>(mat()); }

        TMV_INLINE int nrows() const { return colsize(); }
        TMV_INLINE int ncols() const { return rowsize(); }
        TMV_INLINE bool isSquare() const { return colsize() == rowsize(); }

        // Note that these last functions need to be defined in a more derived
        // class than this, or an infinite loop will result when compiling.

        TMV_INLINE int colsize() const { return mat().colsize(); }
        TMV_INLINE int rowsize() const { return mat().rowsize(); }
        TMV_INLINE int nlo() const { return mat().nlo(); }
        TMV_INLINE int nhi() const { return mat().nhi(); }
        TMV_INLINE int nElements() const { return mat().nElements(); }

        TMV_INLINE value_type cref(int i, int j) const  
        { return eval().cref(i,j); }

        template <class M2>
        TMV_INLINE void assignTo(BaseMatrix_Mutable<M2>& m2) const
        { mat().assignTo(m2.mat()); }

    }; // BaseMatrix

    template <class M>
    class BaseMatrix_Calc : 
        public BaseMatrix<M>
    {
    public:
        enum { _colsize = Traits<M>::_colsize };
        enum { _rowsize = Traits<M>::_rowsize };
        enum { _shape = Traits<M>::_shape };
        enum { _fort = Traits<M>::_fort };
        enum { _calc = Traits<M>::_calc };
        enum { _conj = Traits<M>::_conj };
        enum { _checkalias = Traits<M>::_checkalias };
        enum { _rowmajor = Traits<M>::_rowmajor }; 
        enum { _colmajor = Traits<M>::_colmajor }; 
        enum { _hasdivider = Traits<M>::_hasdivider };

        typedef M type;
        typedef BaseMatrix<M> base;

        typedef typename base::calc_type calc_type;
        typedef typename base::eval_type eval_type;
        typedef typename base::copy_type copy_type;
        typedef typename base::inverse_type inverse_type;

        typedef typename base::value_type value_type;
        typedef typename base::real_type real_type;
        typedef typename base::complex_type complex_type;
        typedef typename base::float_type float_type;
        typedef typename base::zfloat_type zfloat_type;

        typedef typename Traits<M>::const_view_type const_view_type;
        typedef typename Traits<M>::const_cview_type const_cview_type;
        typedef typename Traits<M>::const_fview_type const_fview_type;
        typedef typename Traits<M>::const_xview_type const_xview_type;

        typedef typename Traits<M>::const_transpose_type const_transpose_type;
        typedef typename Traits<M>::const_conjugate_type const_conjugate_type;
        typedef typename Traits<M>::const_adjoint_type const_adjoint_type;
        typedef typename Traits<M>::const_realpart_type const_realpart_type;
        typedef typename Traits<M>::const_imagpart_type const_imagpart_type;
        typedef typename Traits<M>::const_nonconj_type const_nonconj_type;
        typedef typename Traits<M>::nonconst_type nonconst_type;

        typedef typename Traits<M>::const_rowmajor_iterator 
            const_rowmajor_iterator;
        typedef typename Traits<M>::const_colmajor_iterator 
            const_colmajor_iterator;

        //
        // Constructor
        //

    protected:
        TMV_INLINE BaseMatrix_Calc() {}
        TMV_INLINE BaseMatrix_Calc(const BaseMatrix_Calc<M>&) {}
        TMV_INLINE ~BaseMatrix_Calc() {}

    private:
        void operator=(const BaseMatrix_Calc<M>&);
    public:


        // All of these functions are implemented in a derived class,
        // but these are the things that all of the differently shaped
        // matrices should be able to implement in their base class.
        // e.g. BaseMatrix_Rec, BaseMatrix_Diag, etc.

        //
        // Views
        //

        TMV_INLINE TMV_MAYBE_CREF(type,const_view_type) view() const
        { return mat().view(); }

        TMV_INLINE TMV_MAYBE_CREF(type,const_cview_type) cView() const
        { return mat().view(); }

        TMV_INLINE TMV_MAYBE_CREF(type,const_fview_type) fView() const
        { return mat().view(); }

        TMV_INLINE TMV_MAYBE_CREF(type,const_xview_type) xView() const
        { return mat().view(); }

        TMV_INLINE const_transpose_type transpose() const
        { return mat().transpose(); }

        TMV_INLINE TMV_MAYBE_CREF(type,const_conjugate_type) conjugate() const
        { return mat().conjugate(); }

        TMV_INLINE const_adjoint_type adjoint() const
        { return mat().adjoint(); }

        TMV_INLINE const_realpart_type realPart() const
        { return mat().realPart(); }

        TMV_INLINE const_imagpart_type imagPart() const
        { TMVStaticAssert(type::iscomplex); return mat().imagPart(); }

        TMV_INLINE TMV_MAYBE_CREF(type,const_nonconj_type) nonConj() const
        { return mat().nonConj(); }

        TMV_INLINE nonconst_type nonConst() const
        { return mat().nonConst(); } 


        //
        // Division Routines
        //

        TMV_INLINE inverse_type inverse() const
        { return inverse_type(real_type(1),mat()); }

        TMV_INLINE_ND value_type det() const
        {
            TMVStaticAssert((Sizes<_rowsize,_colsize>::same));
            TMVAssert(colsize() == rowsize());
            return tmv::Det(*this);
        }

        TMV_INLINE_ND float_type logDet(zfloat_type* sign=0) const
        {
            TMVStaticAssert((Sizes<_rowsize,_colsize>::same));
            TMVAssert(colsize() == rowsize());
            return tmv::LogDet(*this,sign);
        }

        TMV_INLINE_ND bool isSingular() const
        {
            TMVStaticAssert((Sizes<_rowsize,_colsize>::same));
            TMVAssert(colsize() == rowsize());
            return tmv::IsSingular(*this);
        }

        template <class M2>
        TMV_INLINE void makeInverse(BaseMatrix_Mutable<M2>& minv) const
        { tmv::MakeInverse(Scaling<1,real_type>(),*this,minv); }

        template <class M2>
        TMV_INLINE void makeInverseATA(BaseMatrix_Mutable<M2>& mata) const
        { tmv::MakeInverseATA(*this,mata); }

        TMV_INLINE real_type norm2() const
        { return mat().norm2(); }

        TMV_INLINE real_type condition() const
        { return mat().condition(); }

        //
        // Auxilliary routines
        //

        TMV_INLINE bool isconj() const { return _conj; }

        TMV_INLINE bool isrm() const { return mat().isrm(); }
        TMV_INLINE bool iscm() const { return mat().iscm(); }

        TMV_INLINE const type& mat() const
        { return static_cast<const type&>(*this); }

        TMV_INLINE int colsize() const { return mat().colsize(); }
        TMV_INLINE int rowsize() const { return mat().rowsize(); }

        TMV_INLINE int rowstart(int i) const { return mat().rowstart(i); }
        TMV_INLINE int rowend(int i) const { return mat().rowend(i); }
        TMV_INLINE int colstart(int j) const { return mat().colstart(j); }
        TMV_INLINE int colend(int j) const { return mat().colend(j); }

        TMV_INLINE const_rowmajor_iterator rowmajor_begin() const
        { return mat().rowmajor_begin(); }
        TMV_INLINE const_colmajor_iterator colmajor_begin() const
        { return mat().colmajor_begin(); }

    }; // BaseMatrix_Calc

    template <class M>
    class BaseMatrix_Mutable 
    {
    public:

        enum { _colsize = Traits<M>::_colsize };
        enum { _rowsize = Traits<M>::_rowsize };
        enum { _shape = Traits<M>::_shape };
        enum { _fort = Traits<M>::_fort };
        enum { _calc = Traits<M>::_calc };
        enum { _conj = Traits<M>::_conj };
        enum { _checkalias = Traits<M>::_checkalias };
        enum { _rowmajor = Traits<M>::_rowmajor }; 
        enum { _colmajor = Traits<M>::_colmajor }; 

        typedef M type;

        typedef typename Traits<M>::value_type value_type;
        typedef typename Traits<M>::real_type real_type;
        typedef typename Traits<real_type>::float_type float_type;

        typedef typename Traits<M>::view_type view_type;
        typedef typename Traits<M>::cview_type cview_type;
        typedef typename Traits<M>::fview_type fview_type;
        typedef typename Traits<M>::xview_type xview_type;

        typedef typename Traits<M>::transpose_type transpose_type;
        typedef typename Traits<M>::conjugate_type conjugate_type;
        typedef typename Traits<M>::adjoint_type adjoint_type;
        typedef typename Traits<M>::realpart_type realpart_type;
        typedef typename Traits<M>::imagpart_type imagpart_type;
        typedef typename Traits<M>::nonconj_type nonconj_type;
        typedef typename Traits<M>::noalias_type noalias_type;
        typedef typename Traits<M>::alias_type alias_type;
 
        typedef typename Traits<M>::reference reference;

        typedef typename Traits<M>::rowmajor_iterator rowmajor_iterator;
        typedef typename Traits<M>::colmajor_iterator colmajor_iterator;
        typedef ListAssigner<value_type,rowmajor_iterator> list_assigner;

        //
        // Constructor
        //

    protected:
        TMV_INLINE BaseMatrix_Mutable() {}
        TMV_INLINE BaseMatrix_Mutable(const BaseMatrix_Mutable<M>&) {}
        TMV_INLINE ~BaseMatrix_Mutable() {}
    public :


        //
        // Access 
        //

        TMV_INLINE_ND reference operator()(int i, int j)
        {
            CheckIndex<_fort>(i,colsize());
            CheckIndex<_fort>(j,rowsize());
            return ref(i,j);
        }


        //
        // Op =
        //

        TMV_INLINE_ND type& operator=(const BaseMatrix_Mutable<M>& m2) 
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

        TMV_INLINE list_assigner operator<<(const value_type& x) 
        { return list_assigner(rowmajor_begin(),mat().nElements(),x); }


        //
        // Modifying Functions
        //

        TMV_INLINE type& setZero() 
        { return mat().setZero(); }

        TMV_INLINE type& setAllTo(value_type val) 
        { return mat().setAllTo(val); }

        TMV_INLINE type& addToAll(value_type val) 
        { return mat().addToAll(val); }

        TMV_INLINE type& clip(float_type thresh) 
        { return mat().clip(thresh); }

        template <class F>
        TMV_INLINE type& applyToAll(const F& f)
        { return mat().applyToAll(f); }

        TMV_INLINE type& conjugateSelf() 
        { return mat().conjugateSelf(); }


        //
        // Views
        //

        TMV_INLINE TMV_MAYBE_REF(type,view_type) view() 
        { return mat().view(); }

        TMV_INLINE TMV_MAYBE_REF(type,cview_type) cView() 
        { return mat().cView(); }

        TMV_INLINE TMV_MAYBE_REF(type,fview_type) fView() 
        { return mat().fView(); }

        TMV_INLINE TMV_MAYBE_REF(type,xview_type) xView() 
        { return mat().xView(); }

        TMV_INLINE transpose_type transpose() 
        { return mat().transpose(); }

        TMV_INLINE TMV_MAYBE_REF(type,conjugate_type) conjugate() 
        { return mat().conjugate(); }

        TMV_INLINE adjoint_type adjoint() 
        { return mat().adjiont(); }

        TMV_INLINE realpart_type realPart()
        { return mat().realPart(); }

        TMV_INLINE imagpart_type imagPart()
        { TMVStaticAssert(type::iscomplex); return mat().imagPart(); }

        TMV_INLINE TMV_MAYBE_REF(type,nonconj_type) nonConj()
        { return mat().nonConj(); }

        TMV_INLINE TMV_MAYBE_REF(type,nonconj_type) noAlias()
        { return mat().noAlias(); }

        TMV_INLINE TMV_MAYBE_REF(type,nonconj_type) alias()
        { return mat().alias(); }


        //
        // Arithmetic
        //

        // These operators need to be here, rather than just defining
        // non-member operator*=, etc., since the argument to a non-member
        // function would have to be BaseMatrix_Mutable& (ie. a reference).
        // But then you couldn't write something like:
        // m.transpose() += m2;
        // since the m.transpose() function returns a view by value, which is
        // not castable to the non-const reference argument of operator+=.
        // So we define all these here with unspecified right hand sides.
        // They just have to be valid arguments to a MultEq, AddEq, etc.
        // function defined as non-member functions for various objects.

        template <class X2>
        TMV_INLINE type& operator+=(const X2& x2)
        { AddEq(mat(),x2); return mat(); }

        template <class X2>
        TMV_INLINE type& operator-=(const X2& x2)
        { SubtractEq(mat(),x2); return mat(); }

        template <class X2>
        TMV_INLINE type& operator*=(const X2& x2)
        { MultEq(mat(),x2); return mat(); }

        template <class X2>
        TMV_INLINE type& operator/=(const X2& x2)
        { LDivEq(mat(),x2); return mat(); }

        template <class X2>
        TMV_INLINE type& operator%=(const X2& x2)
        { RDivEq(mat(),x2); return mat(); }


        //
        // Auxilliary routines
        //

        TMV_INLINE const type& mat() const
        { return static_cast<const type&>(*this); }
        TMV_INLINE type& mat()
        { return static_cast<type&>(*this); }

        TMV_INLINE int colsize() const { return mat().colsize(); }
        TMV_INLINE int rowsize() const { return mat().rowsize(); }
        TMV_INLINE reference ref(int i, int j) { return mat().ref(i,j); }
        TMV_INLINE value_type cref(int i, int j) const  
        { return mat().cref(i,j); }

        TMV_INLINE rowmajor_iterator rowmajor_begin()
        { return mat().rowmajor_begin(); }
        TMV_INLINE colmajor_iterator colmajor_begin()
        { return mat().colmajor_begin(); }

    }; // BaseMatrix_Mutable

    template <class T>
    class MatrixSizer;

    template <class T>
    struct Traits<MatrixSizer<T> >
    {
        typedef T value_type;
        typedef MatrixSizer<T> type;
        typedef InvalidType calc_type;
        typedef InvalidType eval_type;
        typedef InvalidType copy_type;
        typedef InvalidType inverse_type;

        enum { _colsize = TMV_UNKNOWN };
        enum { _rowsize = TMV_UNKNOWN };
        enum { _nlo = TMV_UNKNOWN };
        enum { _nhi = TMV_UNKNOWN };
        enum { _shape = Null };
        enum { _fort = false };
        enum { _calc = false };
    };

    template <class T>
    class MatrixSizer : 
        public BaseMatrix<MatrixSizer<T> >
    {
    public:
        TMV_INLINE MatrixSizer(const int _cs, const int _rs) :
            cs(_cs), rs(_rs) {}
        TMV_INLINE int colsize() const { return cs; }
        TMV_INLINE int rowsize() const { return rs; }
        TMV_INLINE int nlo() const { return TMV_MAX(cs-1,0); }
        TMV_INLINE int nhi() const { return TMV_MAX(rs-1,0); }
        TMV_INLINE int nElements() const { return cs * rs; }

        TMV_INLINE T cref(int , int ) const  { return T(0); }

        template <class M2>
        TMV_INLINE void assignTo(BaseMatrix_Mutable<M2>& ) const {}

    private :
        const int cs, rs;
    }; // MatrixSizer

    //
    // Trace
    //

    // This one is BaseMatrix, not BaseMatrix_Calc, 
    // since it should really be called with an eval() object, not calc(), 
    // since you don't need to calculate most of the elements.
    // This is also why we need DoTrace, rather than simply Trace, since
    // we have to make sure eval() is called.

    template <int algo, int size, class M>
    struct Trace_Helper;

    // algo 0: Zero size, so trace is defined as 0.
    template <int size, class M>
    struct Trace_Helper<0,size,M>
    {
        static inline typename M::value_type call(const M& m) 
        { return typename M::value_type(0); }
    };

    // algo 1: Values are calculated.  Call diag().sumElements().
    template <int size, class M>
    struct Trace_Helper<1,size,M>
    {
        static inline typename M::value_type call(const M& m) 
        { return m.diag().sumElements(); }
    };

    // algo 2: Values are not calculated.  Do sum using cref(i,i).
    template <int size, class M>
    struct Trace_Helper<2,size,M>
    {
        static inline typename M::value_type call(const M& m) 
        {
            const int n = size == TMV_UNKNOWN ? m.colsize() : size;
            typename M::value_type sum(0);
            for (int i=0;i<n;++i) sum += m.cref(i,i);
            return sum;
        }
    };

    // algo 3: Special for triangular with unit diagonal
    template <int size, class M>
    struct Trace_Helper<3,size,M>
    {
        static inline typename M::value_type call(const M& m) 
        {
            const int n = size == TMV_UNKNOWN ? m.colsize() : size;
            return typename M::value_type(n); 
        }
    };

    // algo 4: Check for UnknownDiag
    template <int size, class M>
    struct Trace_Helper<4,size,M>
    {
        static inline typename M::value_type call(const M& m) 
        {
            if (m.isunit()) return Trace_Helper<3,size,M>::call(m);
            else return Trace_Helper<1,size,M>::call(m);
        }
    };

    // algo 11: Triangular matrix.  Check _unit and _unknowndiag.
    // (Can't do this below, since _unit and _unknowndiag aren't 
    //  valid members of most matrix classes.)
    template <int size, class M>
    struct Trace_Helper<11,size,M>
    {
        static inline typename M::value_type call(const M& m) 
        {
            const int algo =
                M::_unit ? 3 :
                M::_unknowndiag ? 4 :
                1;
            return Trace_Helper<algo,size,M>::call(m);
        }
    };


    template <class M>
    inline typename M::value_type DoTrace(const BaseMatrix<M>& m)
    {
        TMVStaticAssert((Sizes<M::_rowsize,M::_colsize>::same));
        TMVAssert(m.colsize() == m.rowsize());
        const int size = Sizes<M::_rowsize,M::_colsize>::size;
        const int algo = 
            size == 0 ? 0 :
            ShapeTraits<M::_shape>::unit ? 3 :
            !M::_calc ? 2 : 
            // Check for possible UnknownDiag
            ShapeTraits<M::_shape>::unit_shape != int(InvalidShape) ? 11 :
            1;
        return Trace_Helper<algo,size,M>::call(m.mat());
    }


    //
    // Matrix ==, != Matrix
    //

    // I don't make any effort to optimize this, since it's 
    // probably only going to be used in assert statements and 
    // the like.  If anyone really wants to test equality, the 
    // better way to do it is something like 
    // if (Norm(a-b)<=1.e-10) {...}
    //
    // Also, using crefs as I do here means that I don't have to 
    // write different versions for each kind of matrix.
    // This will work for every possible pairing.
    template <bool rm, int cs, int rs, class M1, class M2>
    struct EqMM_Helper;

    template <int cs, int rs, class M1, class M2>
    struct EqMM_Helper<true,cs,rs,M1,M2> // rm = true
    {
        static bool eq(const M1& m1, const M2& m2)
        {
            const int M = cs == TMV_UNKNOWN ? m1.colsize() : cs;
            const int N = rs == TMV_UNKNOWN ? m1.rowsize() : rs;
            for(int i=0;i<M;++i) {
                for(int j=0;j<N;++j) {
                    if (m1.cref(i,j) != m2.cref(i,j)) return false;
                }
            }
            return true;
        }
    };

    template <int cs, int rs, class M1, class M2>
    struct EqMM_Helper<false,cs,rs,M1,M2> // rm = false
    {
        static bool eq(const M1& m1, const M2& m2)
        {
            const int M = cs == TMV_UNKNOWN ? m1.colsize() : cs;
            const int N = rs == TMV_UNKNOWN ? m1.rowsize() : rs;
            for(int j=0;j<N;++j) {
                for(int i=0;i<M;++i) {
                    if (m1.cref(i,j) != m2.cref(i,j)) return false;
                }
            }
            return true;
        }
    };

    template <class M1, class M2>
    inline bool CallEq(
        const BaseMatrix_Calc<M1>& m1, const BaseMatrix_Calc<M2>& m2)
    {
        TMVStaticAssert((Sizes<M1::_colsize,M2::_colsize>::same)); 
        TMVStaticAssert((Sizes<M1::_rowsize,M2::_rowsize>::same)); 
        TMVAssert(m1.colsize() == m2.colsize());
        TMVAssert(m1.rowsize() == m2.rowsize());
        const int cs = Sizes<M1::_colsize,M2::_colsize>::size;
        const int rs = Sizes<M1::_rowsize,M2::_rowsize>::size;
        const bool rm = M2::_rowmajor || (M1::_rowmajor && !M2::_colmajor);
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::const_cview_type M2v;
        TMV_MAYBE_CREF(M1,M1v) m1v = m1.cView();
        TMV_MAYBE_CREF(M2,M2v) m2v = m2.cView();
        return EqMM_Helper<rm,cs,rs,M1v,M2v>::eq(m1v,m2v);
    }

    template <class M1, class M2>
    TMV_INLINE bool operator==(
        const BaseMatrix<M1>& m1, const BaseMatrix<M2>& m2)
    { return CallEq(m1.calc(),m2.calc()); }

    template <class M1, class M2>
    TMV_INLINE bool operator!=(
        const BaseMatrix<M1>& m1, const BaseMatrix<M2>& m2)
    { return !(m1 == m2); }


    //
    // Other Functions of Matrices
    // (These just call the method mersion.)
    //

    template <class M>
    TMV_INLINE typename M::value_type Trace(const BaseMatrix<M>& m)
    { return m.trace(); }

    template <class M>
    TMV_INLINE typename M::float_type Norm(const BaseMatrix<M>& m)
    { return m.norm(); }
    template <class M>
    TMV_INLINE typename M::float_type NormF(const BaseMatrix<M>& m)
    { return m.normF(); }
    template <class M>
    TMV_INLINE typename M::real_type NormSq(const BaseMatrix<M>& m)
    { return m.normSq(); }
    template <class M>
    TMV_INLINE typename M::float_type Norm1(const BaseMatrix<M>& m)
    { return m.norm1(); }
    template <class M>
    TMV_INLINE typename M::float_type NormInf(const BaseMatrix<M>& m)
    { return m.normInf(); }
    template <class M>
    TMV_INLINE typename M::float_type Norm2(const BaseMatrix<M>& m)
    { return m.norm2(); }
    template <class M>
    TMV_INLINE typename M::float_type Condition(const BaseMatrix<M>& m)
    { return m.condition(); }

    template <class M>
    TMV_INLINE typename M::float_type MaxAbsElement(const BaseMatrix<M>& m)
    { return m.maxAbsElement(); }
    template <class M>
    TMV_INLINE typename M::real_type MaxAbs2Element(const BaseMatrix<M>& m)
    { return m.maxAbs2Element(); }

    template <class M>
    TMV_INLINE typename M::value_type SumElements(const BaseMatrix<M>& m)
    { return m.sumElements(); }
    template <class M>
    TMV_INLINE typename M::float_type SumAbsElements(const BaseMatrix<M>& m)
    { return m.sumAbsElements(); }
    template <class M>
    TMV_INLINE typename M::real_type SumAbs2Elements(const BaseMatrix<M>& m)
    { return m.sumAbs2Elements(); }

    template <class M>
    TMV_INLINE typename M::const_conjugate_type Conjugate(
        const BaseMatrix_Calc<M>& m)
    { return m.conjugate(); }
    template <class M>
    TMV_INLINE typename M::const_transpose_type Transpose(
        const BaseMatrix_Calc<M>& m)
    { return m.transpose(); }
    template <class M>
    TMV_INLINE typename M::const_adjoint_type Adjoint(
        const BaseMatrix_Calc<M>& m)
    { return m.adjoint(); }

    template <class M>
    TMV_INLINE typename M::inverse_type Inverse(const BaseMatrix<M>& m)
    { return m.inverse(); }

    // These next few allow expressions like c = Transpose(a*b)
    template <class M>
    TMV_INLINE typename M::calc_type::const_conjugate_type::copy_type Conjugate(
        const BaseMatrix<M>& m)
    { return m.calc().conjugate(); }
    template <class M>
    TMV_INLINE typename M::calc_type::const_transpose_type::copy_type Transpose(
        const BaseMatrix<M>& m)
    { return m.calc().transpose(); }
    template <class M>
    TMV_INLINE typename M::calc_type::const_adjoint_type::copy_type Adjoint(
        const BaseMatrix<M>& m)
    { return m.calc().adjoint(); }


    //
    // TMV_Text 
    //

#ifdef TMV_TEXT
    template <class M>
    inline std::string TMV_Text(const BaseMatrix<M>& m)
    {
        std::ostringstream s;
        s << "BaseMatrix< "<<TMV_Text(m.mat())<<" >";
        return s.str();
    }

    template <class M>
    inline std::string TMV_Text(const BaseMatrix_Calc<M>& m)
    {
        std::ostringstream s;
        s << "BaseMatrix_Calc< "<<TMV_Text(m.mat())<<" >";
        return s.str();
    }

    template <class M>
    inline std::string TMV_Text(const BaseMatrix_Mutable<M>& m)
    {
        std::ostringstream s;
        s << "BaseMatrix_Mutable< "<<TMV_Text(m.mat())<<" >";
        return s.str();
    }
#endif

} // namespace tmv

#endif
