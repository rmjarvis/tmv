
//-----------------------------------------------------------------------------
//
// This file defines the BaseMatrix_Band and BaseMatrix_Band_Mutable
// classes.  These are the base classes for the banded
// matrices, BandMatrix, ThinBandMatrix, and SmallBandMatrix.

#ifndef TMV_BaseMatrix_Band_H
#define TMV_BaseMatrix_Band_H

#include "TMV_BaseMatrix.h"
#include "TMV_BaseMatrix_Rec.h"
#include "TMV_BaseMatrix_Tri.h"
#include "TMV_BaseMatrix_Diag.h"

namespace tmv {

    // BaseMatrix_Band is derived from BaseMatrix_Calc, and is used
    // for banded matrices.
    template <class M>
    class BaseMatrix_Band;

    // BaseMatrix_Band adds the following requirements to Traits<M>:
    //
    //  _nlo = the number of (non-zero) bands below the main diagonal if known
    //  _nhi = the number of (non-zero) bands above the main diagonal if known
    //  _stepi = the step size along column if known (else Unknown)
    //  _stepj = the step size along row if known (else Unknown)
    //  _diagstep = the step size along row if known (else Unknown)
    //
    //  const_diag_type = return type from diag() const
    //  const_row_sub_type = return type from row(i,j1,j2) const
    //  const_col_sub_type = return type from col(j,i1,i2) const
    //  const_diag_sub_type = return type from diag(i) and diag(i,j1,j2) const
    //
    //  const_submatrix_type = return type from subMatrix(i1,i2,j1,j2) const
    //  const_submatrix_step_type = return type from 
    //      subMatrix(i1,i2,j1,j2,istep,jstep) const
    //  const_subvector_type = return type from 
    //      subVector(i,j,istep,jstep,size) const
    //  const_subbandmatrix_type = return type from 
    //      subBandMatrix(i1,i2,j1,j2) const and
    //      subBandMatrix(i1,i2,j1,j2,newnlo,newnhi) const
    //  const_subbandmatrix_step_type = return type from 
    //      subBandMatrix(i1,i2,j1,j2,newnlo,newnhi,istep,jstep) const
    //  const_colrange_type = return type from colRange(j1,j2) const
    //  const_rowrange_type = return type from rowRange(i1,i2) const
    //  const_diagrange_type = return type from diagRange(k1,k2) const
    //
    //  const_upperband_type = return type from upperBand() const
    //  const_lowerband_type = return type from lowerBand() const
    //  const_upperbandoff_type = return type from upperBandOff() const
    //  const_lowerbandoff_type = return type from lowerBandOff() const


    // BaseMatrix_Band_Mutable is derived from both BaseMatrix_Band and
    // BaseMatrix_Mutable, and is used for mutable, banded matrices.
    template <class M>
    class BaseMatrix_Band_Mutable;

    // BaseMatrix_Band_Mutable adds:
    //
    //  diag_type = return type from diag() 
    //  row_sub_type = return type from row(i,j1,j2) 
    //  col_sub_type = return type from col(j,i1,i2) 
    //  diag_sub_type = return type from diag(i) and diag(i,j1,j2) 
    //
    //  submatrix_type = return type from subMatrix(i1,i2,j1,j2) 
    //  submatrix_step_type = return type from 
    //      subMatrix(i1,i2,j1,j2,istep,jstep) 
    //  subvector_type = return type from subVector(i,j,istep,jstep,size)
    //  subbandmatrix_type = return type from 
    //      subBandMatrix(i1,i2,j1,j2) and
    //      subBandMatrix(i1,i2,j1,j2,newnlo,newnhi) 
    //  subbandmatrix_step_type = return type from 
    //      subBandMatrix(i1,i2,j1,j2,istep,jstep) 
    //  colrange_type = return type from colRange(j1,j2) 
    //  rowrange_type = return type from rowRange(i1,i2) 
    //
    //  upperband_type = return type from upperBand() 
    //  lowerband_type = return type from lowerBand() 
    //  upperbandoff_type = return type from upperBandOff() 
    //  lowerbandoff_type = return type from lowerBandOff() 

    // The following all derive from BaseMatrix_Band or BaseMatrix_Band_Mutable.
    // See TMV_BandMatrix.h and TMV_SmallBandMatrix.h for their definitions:
    template <class T, int A=0>
    class BandMatrix;
    template <class T, int A=0>
    class ConstBandMatrixView;
    template <class T, int A=0>
    class BandMatrixView;
    template <class T, ptrdiff_t LO, ptrdiff_t HI, int A=0>
    class ThinBandMatrix;
    template <class T, ptrdiff_t M, ptrdiff_t N, ptrdiff_t LO, ptrdiff_t HI, int A=0>
    class SmallBandMatrix;
    template <class T, ptrdiff_t M, ptrdiff_t N, ptrdiff_t LO, ptrdiff_t HI, ptrdiff_t Si=Unknown, ptrdiff_t Sj=Unknown, int A=0>
    class ConstSmallBandMatrixView;
    template <class T, ptrdiff_t M, ptrdiff_t N, ptrdiff_t LO, ptrdiff_t HI, ptrdiff_t Si=Unknown, ptrdiff_t Sj=Unknown, int A=0>
    class SmallBandMatrixView;

    // In TMV_Norm.h
    template <class M>
    inline typename M::float_type DoNormF(const BaseMatrix_Band<M>& m);
    template <class M>
    inline typename M::float_type DoNorm2(const BaseMatrix<M>& m);
    template <class M>
    inline typename M::float_type DoCondition(const BaseMatrix<M>& m);

    //
    // Helper functions and values:
    //

    // Specify ExactSameStorage for banded matrices:
    template <class M1, class M2>
    TMV_INLINE bool ExactSameStorage(
        const BaseMatrix_Band<M1>& m1, const BaseMatrix_Band<M2>& m2)
    {
        typedef typename M1::value_type T1;
        typedef typename M2::value_type T2;
        return Traits2<T1,T2>::sametype && 
            (m1.stepi() == m2.stepi() && m1.stepj() == m2.stepj()) &&
            (m1.nhi() == m2.nhi() && m1.nlo() == m2.nlo()); 
    }

    inline bool InBand(ptrdiff_t i, ptrdiff_t j, ptrdiff_t nlo, ptrdiff_t nhi)
    { return (j+nlo >= i) && (i+nhi >= j); }

    // These helper functions check the validity of indices according
    // to whether the matrix uses CStyle or FortranStyle indexing.
    // They also update the indices to be consistent with CStyle.
    template <bool _fort>
    TMV_INLINE_ND void CheckRowRange_Band(
        ptrdiff_t& j, ptrdiff_t& i1, ptrdiff_t i2, ptrdiff_t m, ptrdiff_t nlo, ptrdiff_t nhi)
    { // CStyle or FortranStyle
        CheckRowRange<_fort>(i1,i2,m);
        TMVAssert(InBand(i1,j,nlo,nhi) && "first element must be in band");
        TMVAssert(InBand(i2-1,j,nlo,nhi) && "last element must be in band");
    }
    template <bool _fort>
    TMV_INLINE_ND void CheckColRange_Band(
        ptrdiff_t& i, ptrdiff_t& j1, ptrdiff_t j2, ptrdiff_t n, ptrdiff_t nlo, ptrdiff_t nhi)
    { // CStyle or FortranStyle
        CheckColRange<_fort>(j1,j2,n);
        TMVAssert(InBand(i,j1,nlo,nhi) && "first element must be in band");
        TMVAssert(InBand(i,j2-1,nlo,nhi) && "last element must be in band");
    }
    template <bool _fort>
    TMV_INLINE_ND void CheckDiagIndex_Band(
        ptrdiff_t& i, ptrdiff_t m, ptrdiff_t n, ptrdiff_t nlo, ptrdiff_t nhi)
    { // CStyle
        TMVAssert(i >= -nlo && "negative diag index must be <= nlo");
        TMVAssert(i <= nhi+1 && "positive diag index must be <= nhi+1");
        CheckDiagIndex<false>(i,m,n);
    }
    template <bool _fort>
    TMV_INLINE_ND void CheckDiagIndex_Band(
        ptrdiff_t& i, ptrdiff_t& j1, ptrdiff_t& j2, ptrdiff_t m, ptrdiff_t n, ptrdiff_t nlo, ptrdiff_t nhi)
    { // CStyle
        TMVAssert(i >= -nlo && "negative diag index must be <= nlo");
        TMVAssert(i <= nhi+1 && "positive diag index must be <= nhi+1");
        CheckDiagIndex<false>(i,j1,j2,m,n);
    }
    template <bool _fort>
    TMV_INLINE_ND void CheckMatSubVector_Band(
        ptrdiff_t& i, ptrdiff_t& j, ptrdiff_t istep, ptrdiff_t jstep, ptrdiff_t size, ptrdiff_t m, ptrdiff_t n,
        ptrdiff_t nlo, ptrdiff_t nhi)
    { // CStyle or FortranStyle
        TMVAssert(InBand(i,j,nlo,nhi) && 
                  "first element must be in band");
        TMVAssert(InBand(i+istep*(size-1),j+jstep*(size-1),nlo,nhi) && 
                  "last element must be in band");
        CheckMatSubVector<_fort>(i,j,istep,jstep,size,m,n);
    }
    template <bool _fort>
    TMV_INLINE_ND void CheckDiagRange_Band(
        ptrdiff_t k1, ptrdiff_t& k2, ptrdiff_t nlo, ptrdiff_t nhi)
    { // CStyle
        TMVAssert(k1 <= k2 && 
                  "range must have a non-negative number of diagonals");
        TMVAssert(k1 >= -nlo && "first diagonal in range must be >= -nlo");
        TMVAssert(k2 <= nhi+1 && "end of diagonal range must be <= nhi+1");
    }

    template <>
    TMV_INLINE_ND void CheckDiagIndex_Band<true>(
        ptrdiff_t& i, ptrdiff_t m, ptrdiff_t n, ptrdiff_t nlo, ptrdiff_t nhi)
    { // FortranStyle
        TMVAssert(i >= -nlo && "negative diag index must be <= nlo");
        TMVAssert(i <= nhi && "positive diag index must be <= nhi");
        CheckDiagIndex<true>(i,m,n);
    }
    template <>
    TMV_INLINE_ND void CheckDiagIndex_Band<true>(
        ptrdiff_t& i, ptrdiff_t& j1, ptrdiff_t& j2, ptrdiff_t m, ptrdiff_t n, ptrdiff_t nlo, ptrdiff_t nhi)
    { // FortranStyle
        TMVAssert(i >= -nlo && "negative diag index must be <= nlo");
        TMVAssert(i <= nhi && "positive diag index must be <= nhi");
        CheckDiagIndex<true>(i,j1,j2,m,n);
    }
    template <>
    TMV_INLINE_ND void CheckDiagRange_Band<true>(
        ptrdiff_t k1, ptrdiff_t& k2, ptrdiff_t nlo, ptrdiff_t nhi)
    { // FortranStyle
        TMVAssert(k1 <= k2 && 
                  "range must have a positive number of diagonals");
        TMVAssert(k1 >= -nlo && "first diagonal in range must be >= -nlo");
        TMVAssert(k2 <= nhi && "last diagonal in range must be <= nhi");
        ++k2;
    }

    TMV_INLINE_ND void CheckInBand(
        ptrdiff_t i, ptrdiff_t j, ptrdiff_t nlo, ptrdiff_t nhi)
    { TMVAssert(InBand(i,j,nlo,nhi) && "Reference value must be in band"); }
    TMV_INLINE_ND void CheckSubMatrix_Band(
        ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2, ptrdiff_t nlo, ptrdiff_t nhi)
    {
        TMVAssert(InBand(i1,j1,nlo,nhi) && 
                  "Upper left corner must be in band");
        TMVAssert(InBand(i1,j2-1,nlo,nhi) && 
                  "Upper right corner must be in band");
        TMVAssert(InBand(i2-1,j1,nlo,nhi) && 
                  "Lower left corner must be in band");
        TMVAssert(InBand(i2-1,j2-1,nlo,nhi) && 
                  "Lower right corner must be in band");
    }
    TMV_INLINE_ND void CheckSubMatrix_Band(
        ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2, ptrdiff_t istep, ptrdiff_t jstep, ptrdiff_t nlo, ptrdiff_t nhi)
    {
        TMVAssert(InBand(i1,j1,nlo,nhi) && 
                  "Upper left corner must be in band");
        TMVAssert(InBand(i1,j2-jstep,nlo,nhi) && 
                  "Upper right corner must be in band");
        TMVAssert(InBand(i2-istep,j1,nlo,nhi) && 
                  "Lower left corner must be in band");
        TMVAssert(InBand(i2-istep,j2-jstep,nlo,nhi) && 
                  "Lower right corner must be in band");
    }

    TMV_INLINE_ND void CheckSubBandMatrix(
        ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2, ptrdiff_t newnlo, ptrdiff_t newnhi,
        ptrdiff_t nlo, ptrdiff_t nhi)
    {
        TMVAssert(InBand(i1,j1,nlo,nhi) && 
                  "Upper left corner must be in band");
        TMVAssert(InBand(i1+newnlo,j1,nlo,nhi) &&
                  "New subdiagonals must be in band");
        TMVAssert(InBand(i1,j1+newnhi,nlo,nhi) &&
                  "New superdiagonals must be in band");
        TMVAssert(newnhi < j2-j1 &&
                  "New nhi must be less than new rowsize");
        TMVAssert(newnlo < i2-i1 &&
                  "New nlo must be less than new colsize");
    }
    TMV_INLINE_ND void CheckSubBandMatrix(
        ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2, ptrdiff_t newnlo, ptrdiff_t newnhi,
        ptrdiff_t istep, ptrdiff_t jstep, ptrdiff_t nlo, ptrdiff_t nhi)
    {
        TMVAssert(InBand(i1,j1,nlo,nhi) && 
                  "Upper left corner must be in band");
        TMVAssert(InBand(i1+newnlo*istep,j1,nlo,nhi) &&
                  "New subdiagonals must be in band");
        TMVAssert(InBand(i1,j1+newnhi*jstep,nlo,nhi) &&
                  "New superdiagonals must be in band");
        TMVAssert(newnhi < j2-j1 &&
                  "New nhi must be less than new rowsize");
        TMVAssert(newnlo < i2-i1 &&
                  "New nlo must be less than new colsize");
    }

    TMV_INLINE_ND void CheckUpperBandOff(ptrdiff_t nhi, ptrdiff_t n)
    {
        TMVAssert(n>0 && "upperBandOff not possible for zero-sized matrix");
        TMVAssert(nhi>0 && 
                  "upperBandOff not possible for band matrix with nhi=0");
    }
    TMV_INLINE_ND void CheckLowerBandOff(ptrdiff_t nlo, ptrdiff_t m)
    {
        TMVAssert(m>0 && "lowerBandOff not possible for zero-sized matrix");
        TMVAssert(nlo>0 && 
                  "lowerBandOff not possible for band matrix with nlo=0");
    }
    TMV_INLINE_ND void CheckUpperBandOff(ptrdiff_t i, ptrdiff_t nhi, ptrdiff_t n)
    {
        CheckUpperBandOff(nhi,n);
        TMVAssert(i>0 && i<=nhi && "upperBandOff index is not valid");
    }
    TMV_INLINE_ND void CheckLowerBandOff(ptrdiff_t i, ptrdiff_t nlo, ptrdiff_t m)
    {
        CheckLowerBandOff(nlo,m);
        TMVAssert(i>0 && i<=nlo && "lowerBandOff index is not valid");
    }


    // Since SmallBandMatrix needs to know nlo and nhi at compile time,
    // we always use BandMatrix here. 
    // BCopyHelper takes lo,hi params.
    template <class T, ptrdiff_t cs, ptrdiff_t rs, int A>
    struct MCopyHelper<T,Band,cs,rs,A>
    {
        enum { A2 = A | NoDivider | NoAlias };
        typedef BandMatrix<T,A2> type;
    };

    template <class T, ptrdiff_t cs, ptrdiff_t rs, ptrdiff_t si, ptrdiff_t sj, int A>
    struct MViewHelper<T,Band,cs,rs,si,sj,A>
    { 
        enum { xx = Unknown };
        typedef SmallBandMatrixView<T,cs,rs,xx,xx,si,sj,A> type; 
        typedef ConstSmallBandMatrixView<T,cs,rs,xx,xx,si,sj,A> ctype; 
    };
    template <class T, ptrdiff_t si, ptrdiff_t sj, int A>
    struct MViewHelper<T,Band,Unknown,Unknown,si,sj,A>
    {
        enum { A2 = A |
            (si == 1 ? ColMajor : sj == 1 ? RowMajor : NonMajor) | NoAlias };
        typedef BandMatrixView<T,A2> type; 
        typedef ConstBandMatrixView<T,A2> ctype; 
    };


    template <class T, int shape, ptrdiff_t cs, ptrdiff_t rs, ptrdiff_t lo, ptrdiff_t hi, int A>
    struct BCopyHelper;
    template <class T, ptrdiff_t cs, ptrdiff_t rs, ptrdiff_t lo, ptrdiff_t hi, int A>
    struct BCopyHelper<T,Band,cs,rs,lo,hi,A>
    {
        typedef SmallBandMatrix<T,cs,rs,lo,hi,A> type;
    };
    template <class T, ptrdiff_t cs, ptrdiff_t rs, ptrdiff_t lo, int A>
    struct BCopyHelper<T,Band,cs,rs,lo,Unknown,A>
    {
        typedef typename MCopyHelper<T,Band,cs,rs,A>::type type;
    };
    template <class T, ptrdiff_t cs, ptrdiff_t rs, ptrdiff_t hi, int A>
    struct BCopyHelper<T,Band,cs,rs,Unknown,hi,A>
    {
        typedef typename MCopyHelper<T,Band,cs,rs,A>::type type;
    };
    template <class T, ptrdiff_t cs, ptrdiff_t rs, int A>
    struct BCopyHelper<T,Band,cs,rs,Unknown,Unknown,A>
    {
        typedef typename MCopyHelper<T,Band,cs,rs,A>::type type;
    };

    template <class T, int shape,  ptrdiff_t cs, ptrdiff_t rs, ptrdiff_t lo, ptrdiff_t hi, ptrdiff_t si, ptrdiff_t sj, int A>
    struct BViewHelper;
    template <class T, ptrdiff_t cs, ptrdiff_t rs, ptrdiff_t lo, ptrdiff_t hi, ptrdiff_t si, ptrdiff_t sj, int A>
    struct BViewHelper<T,Band,cs,rs,lo,hi,si,sj,A>
    { 
        enum { xx = Unknown };
        typedef SmallBandMatrixView<T,cs,rs,lo,hi,si,sj,A> type; 
        typedef ConstSmallBandMatrixView<T,cs,rs,lo,hi,si,sj,A> ctype; 
    };
    template <class T, ptrdiff_t si, ptrdiff_t sj, int A>
    struct BViewHelper<T,Band,Unknown,Unknown,Unknown,Unknown,si,sj,A>
    {
        enum { A2 = A |
            (si == 1 ? ColMajor : sj == 1 ? RowMajor : NonMajor) | NoAlias };
        typedef BandMatrixView<T,A2> type; 
        typedef ConstBandMatrixView<T,A2> ctype; 
    };


    // Defined in TMV_CopyB.h
    template <class M1, class M2>
    inline void Copy(
        const BaseMatrix_Band<M1>& m1, BaseMatrix_Band_Mutable<M2>& m2);

    template <class M1, class M2>
    inline void Copy(
        const BaseMatrix_Band<M1>& m1, BaseMatrix_Rec_Mutable<M2>& m2);

    template <class M1, class M2>
    inline void Copy(
        const BaseMatrix_Band<M1>& m1, BaseMatrix_Tri_Mutable<M2>& m2);

    template <class M1, class M2>
    inline void Copy(
        const BaseMatrix_Tri<M1>& m1, BaseMatrix_Band_Mutable<M2>& m2);

    template <class M1, class M2>
    inline void Copy(
        const BaseMatrix_Diag<M1>& m1, BaseMatrix_Band_Mutable<M2>& m2);

    // Defined in TMV_SwapB.h
    template <class M1, class M2>
    inline void Swap(
        BaseMatrix_Band_Mutable<M1>& m1, BaseMatrix_Band_Mutable<M2>& m2);
    template <class M>
    inline void TransposeSelf(BaseMatrix_Band_Mutable<M>& m);

    // Defined in TMV_NormB.h
    template <class M>
    inline typename M::real_type DoNormSq(const BaseMatrix_Band<M>& m);
    template <class M>
    inline typename M::float_type DoNormSq(
        const BaseMatrix_Band<M>& m, const typename M::float_type scale);
    template <class M>
    inline typename M::float_type DoNorm1(const BaseMatrix_Band<M>& m);
    template <class M>
    inline typename M::float_type DoNormInf(const BaseMatrix_Band<M>& m);
    template <class M>
    inline typename M::value_type DoSumElements(const BaseMatrix_Band<M>& m);
    template <class M>
    inline typename M::float_type DoSumAbsElements(const BaseMatrix_Band<M>& m);
    template <class M>
    inline typename M::real_type DoSumAbs2Elements(const BaseMatrix_Band<M>& m);
    template <class M>
    inline typename M::float_type DoMaxAbsElement(const BaseMatrix_Band<M>& m);
    template <class M>
    inline typename M::real_type DoMaxAbs2Element(const BaseMatrix_Band<M>& m);

    // Defined below:
    template <class M>
    inline void SetZero(BaseMatrix_Band_Mutable<M>& m);
    template <class M, class RT>
    inline void Clip(BaseMatrix_Band_Mutable<M>& m, const RT& thresh);
    template <class M, class T>
    inline void SetAllTo(BaseMatrix_Band_Mutable<M>& m, const T& val);
    template <class M, class T>
    inline void AddToAll(BaseMatrix_Band_Mutable<M>& m, const T& val);
    template <class M>
    inline void ConjugateSelf(BaseMatrix_Band_Mutable<M>& m);
    template <class M, class F>
    inline void ApplyToAll(BaseMatrix_Band_Mutable<M>& m, const F& f);


    // A helper class for returning views without necessarily
    // making a new object.
    template <bool ref, class type, class view_type>
    struct MakeBandView_Helper;

    template <class type, class view_type>
    struct MakeBandView_Helper<true,type,view_type>
    {
        typedef type& ret_type;
        typedef const type& const_ret_type;
        static TMV_INLINE ret_type call(type& m) { return m; }
        static TMV_INLINE const_ret_type call(const type& m) { return m; }
    };

    template <class type, class view_type>
    struct MakeBandView_Helper<false,type,view_type>
    {
        typedef view_type ret_type;
        typedef view_type const_ret_type;
        static TMV_INLINE ret_type call(type& m) 
        {
            return view_type(
                m.ptr(),m.colsize(),m.rowsize(),m.nlo(),m.nhi(),
                m.stepi(),m.stepj()); 
        }
        static TMV_INLINE const_ret_type call(const type& m) 
        {
            return view_type(
                m.cptr(),m.colsize(),m.rowsize(),m.nlo(),m.nhi(),
                m.stepi(),m.stepj()); 
        }
    };

    template <class type, class view_type>
    struct MakeBandView
    {
        enum { ref = Traits2<type,view_type>::sametype };
        typedef MakeBandView_Helper<ref,type,view_type> helper;

        static TMV_INLINE typename helper::ret_type call(type& m)
        { return helper::call(m); }
        static TMV_INLINE typename helper::const_ret_type call(const type& m)
        { return helper::call(m); }
    };

 
    template <class M>
    class BaseMatrix_Band : 
        public BaseMatrix_Calc<M>
    {
    public:
        enum { _colsize = Traits<M>::_colsize };
        enum { _rowsize = Traits<M>::_rowsize };
        enum { _nlo = Traits<M>::_nlo };
        enum { _nhi = Traits<M>::_nhi };
        enum { _shape = Traits<M>::_shape };
        enum { _fort = Traits<M>::_fort };
        enum { _rowmajor = Traits<M>::_rowmajor }; 
        enum { _colmajor = Traits<M>::_colmajor }; 
        enum { _diagmajor = Traits<M>::_diagmajor }; 
        enum { _calc = Traits<M>::_calc };
        enum { _stepi = Traits<M>::_stepi };
        enum { _stepj = Traits<M>::_stepj };
        enum { _diagstep = Traits<M>::_diagstep };
        enum { _conj = Traits<M>::_conj };
        enum { _canlin = Traits<M>::_canlin };

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

        typedef typename Traits<M>::const_row_sub_type const_row_sub_type;
        typedef typename Traits<M>::const_col_sub_type const_col_sub_type;
        typedef typename Traits<M>::const_diag_type const_diag_type;
        typedef typename Traits<M>::const_diag_sub_type const_diag_sub_type;

        typedef typename Traits<M>::const_cmview_type const_cmview_type;
        typedef typename Traits<M>::const_rmview_type const_rmview_type;
        typedef typename Traits<M>::const_dmview_type const_dmview_type;

        typedef typename Traits<M>::const_submatrix_type const_submatrix_type;
        typedef typename Traits<M>::const_submatrix_step_type 
            const_submatrix_step_type;

        typedef typename Traits<M>::const_subvector_type const_subvector_type;

        typedef typename Traits<M>::const_subbandmatrix_type 
            const_subbandmatrix_type;
        typedef typename Traits<M>::const_subbandmatrix_step_type 
            const_subbandmatrix_step_type;

        typedef typename Traits<M>::const_colrange_type const_colrange_type;
        typedef typename Traits<M>::const_rowrange_type const_rowrange_type;
        typedef typename Traits<M>::const_diagrange_type const_diagrange_type;

        typedef typename Traits<M>::const_upperband_type const_upperband_type;
        typedef typename Traits<M>::const_lowerband_type const_lowerband_type;
        typedef typename Traits<M>::const_upperbandoff_type 
            const_upperbandoff_type;
        typedef typename Traits<M>::const_lowerbandoff_type 
            const_lowerbandoff_type;

        typedef typename Traits<M>::const_rowmajor_iterator
            const_rowmajor_iterator;
        typedef typename Traits<M>::const_colmajor_iterator
            const_colmajor_iterator;
        typedef typename Traits<M>::const_diagmajor_iterator
            const_diagmajor_iterator;


        //
        // Constructor
        //

    protected:
        TMV_INLINE BaseMatrix_Band() {}
        TMV_INLINE BaseMatrix_Band(const BaseMatrix_Band<M>&) {}
        TMV_INLINE ~BaseMatrix_Band() {}

    private:
        void operator=(const BaseMatrix_Band<M>&);
    public:


        //
        // Access 
        //

        // The get_ routines always use CStyle indexing.
        TMV_INLINE const_row_sub_type get_row(ptrdiff_t i, ptrdiff_t j1, ptrdiff_t j2) const
        {
            return const_row_sub_type(
                cptr()+i*stepi()+j1*stepj(),j2-j1,stepj()); 
        }

        TMV_INLINE const_col_sub_type get_col(ptrdiff_t j, ptrdiff_t i1, ptrdiff_t i2) const
        {
            return const_col_sub_type(
                cptr()+j*stepj()+i1*stepi(),i2-i1,stepi()); 
        }

        TMV_INLINE const_diag_sub_type get_diag(ptrdiff_t i) const
        {
            return const_diag_sub_type(
                cptr() + (i<0?(-i*stepi()):(i*stepj())),
                ( i<0 ? 
                  TMV_MIN(colsize()+i,rowsize()) :
                  TMV_MIN(colsize(),rowsize()-i)),
                diagstep());
        }

        TMV_INLINE const_diag_sub_type get_diag(ptrdiff_t i, ptrdiff_t j1, ptrdiff_t j2) const
        {
            return const_diag_sub_type(
                cptr() + (i<0?(-i*stepi()):(i*stepj())) + j1*diagstep(),
                j2-j1, diagstep());
        }


        // The regular versions respect the indexing style for i and j:
        TMV_INLINE_ND const_row_sub_type row(ptrdiff_t i, ptrdiff_t j1, ptrdiff_t j2) const
        {
            CheckRowIndex<_fort>(i,colsize());
            CheckColRange_Band<_fort>(i,j1,j2,rowsize(),nlo(),nhi());
            return get_row(i,j1,j2); 
        }

        TMV_INLINE_ND const_col_sub_type col(ptrdiff_t j, ptrdiff_t i1, ptrdiff_t i2) const
        {
            CheckColIndex<_fort>(j,rowsize());
            CheckRowRange_Band<_fort>(j,i1,i2,colsize(),nlo(),nhi());
            return get_col(j,i1,i2); 
        }

        // No need for a get_ routine for diag()
        TMV_INLINE const_diag_type diag() const
        {
            return const_diag_type(
                cptr(),TMV_MIN(colsize(),rowsize()),diagstep()); 
        }

        TMV_INLINE_ND const_diag_sub_type diag(ptrdiff_t i) const
        {
            CheckDiagIndex_Band<_fort>(i,colsize(),rowsize(),nlo(),nhi());
            return get_diag(i);
        }

        TMV_INLINE_ND const_diag_sub_type diag(ptrdiff_t i, ptrdiff_t j1, ptrdiff_t j2) const
        {
            CheckDiagIndex_Band<_fort>(i,j1,j2,colsize(),rowsize(),nlo(),nhi());
            return get_diag(i,j1,j2);
        }


        //
        // Functions
        //

        TMV_INLINE value_type sumElements() const
        { return tmv::DoSumElements(mat()); }

        TMV_INLINE float_type sumAbsElements() const
        { return tmv::DoSumAbsElements(mat()); }

        TMV_INLINE real_type sumAbs2Elements() const
        { return tmv::DoSumAbs2Elements(mat()); }

        TMV_INLINE float_type maxAbsElement() const
        { return tmv::DoMaxAbsElement(mat()); }

        TMV_INLINE real_type maxAbs2Element() const
        { return tmv::DoMaxAbs2Element(mat()); }

        TMV_INLINE real_type normSq() const
        { return tmv::DoNormSq(mat()); }

        TMV_INLINE float_type normSq(const float_type scale) const
        { return tmv::DoNormSq(mat(),scale); }

        TMV_INLINE float_type normF() const 
        { return tmv::DoNormF(mat()); }

        TMV_INLINE float_type norm() const
        { return normF(); }

        TMV_INLINE float_type norm1() const
        { return tmv::DoNorm1(mat()); }

        TMV_INLINE float_type normInf() const
        { return tmv::DoNormInf(mat()); }

        TMV_INLINE float_type norm2() const
        { return tmv::DoNorm2(mat()); }

        TMV_INLINE float_type condition() const
        { return tmv::DoCondition(mat()); }



        //
        // subMatrix, etc.
        //

        // These versions always uses CStyle
        TMV_INLINE const_submatrix_type cSubMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2) const
        {
            return const_submatrix_type(
                cptr()+i1*stepi()+j1*stepj(), i2-i1, j2-j1, stepi(), stepj()); 
        }

        TMV_INLINE const_submatrix_step_type cSubMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2, ptrdiff_t istep, ptrdiff_t jstep) const
        {
            return const_submatrix_step_type(
                cptr()+i1*stepi()+j1*stepj(), (i2-i1)/istep, (j2-j1)/jstep,
                istep*stepi(), jstep*stepj());
        }

        TMV_INLINE const_subvector_type cSubVector(
            ptrdiff_t i, ptrdiff_t j, ptrdiff_t istep, ptrdiff_t jstep, ptrdiff_t s) const
        {
            return const_subvector_type(
                cptr()+i*stepi()+j*stepj(), s, istep*stepi() + jstep*stepj());
        }

        TMV_INLINE const_subbandmatrix_type cSubBandMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2, ptrdiff_t newnlo, ptrdiff_t newnhi) const
        {
            return const_subbandmatrix_type(
                cptr()+i1*stepi()+j1*stepj(), i2-i1, j2-j1,
                newnlo, newnhi, stepi(), stepj()); 
        }

        const_subbandmatrix_type cSubBandMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2) const
        {
            const ptrdiff_t newnlo = TMV_MIN(nlo()+j1-i1,i2-i1-1);
            const ptrdiff_t newnhi = TMV_MIN(nhi()+i1-j1,j2-j1-1);
            CheckSubBandMatrix(i1,i2,j1,j2,newnlo,newnhi,nlo(),nhi());
            return cSubBandMatrix(i1,i2,j1,j2,newnlo,newnhi);
        }

        TMV_INLINE const_subbandmatrix_step_type cSubBandMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2, ptrdiff_t newnlo, ptrdiff_t newnhi,
            ptrdiff_t istep, ptrdiff_t jstep) const
        {
            return const_subbandmatrix_step_type(
                cptr()+i1*stepi()+j1*stepj(), (i2-i1)/istep, (j2-j1)/jstep,
                newnlo, newnhi, istep*stepi(), jstep*stepj());
        }

        const_colrange_type cColRange(ptrdiff_t j1, ptrdiff_t j2) const
        {
            const ptrdiff_t i1 = j1 > nhi() ? j1-nhi() : 0;
            const ptrdiff_t i2 = TMV_MIN(j2 + nlo(),colsize());
            const ptrdiff_t newnhi = j1 < nhi() ? TMV_MIN(nhi(),j2-1) - j1 : 0;
            const ptrdiff_t newnlo = i1==i2 ? 0 : TMV_MIN(nlo()+nhi()-newnhi,i2-i1-1);
            CheckSubBandMatrix(i1,i2,j1,j2,newnlo,newnhi,nlo(),nhi());
            return cSubBandMatrix(i1,i2,j1,j2,newnlo,newnhi);
        }

        const_rowrange_type cRowRange(ptrdiff_t i1, ptrdiff_t i2) const
        {
            const ptrdiff_t j1 = i1 > nlo() ? i1-nlo() : 0;
            const ptrdiff_t j2 = TMV_MIN(i2 + nhi(),rowsize());
            const ptrdiff_t newnlo = i1 < nlo() ? TMV_MIN(nlo(),i2-1) - i1 : 0;
            const ptrdiff_t newnhi = j1==j2 ? 0 : TMV_MIN(nlo()+nhi()-newnlo,j2-j1-1);
            CheckSubBandMatrix(i1,i2,j1,j2,newnlo,newnhi,nlo(),nhi());
            return cSubBandMatrix(i1,i2,j1,j2,newnlo,newnhi);
        }

        const_diagrange_type cDiagRange(ptrdiff_t k1, ptrdiff_t k2) const
        {
            const ptrdiff_t i1 = k2 <= 0 ? -k2+1 : 0;
            const ptrdiff_t i2 = TMV_MIN(rowsize()-k1,colsize());
            const ptrdiff_t j1 = k1 <= 0 ? 0 : k1;
            const ptrdiff_t j2 = TMV_MIN(rowsize(),colsize()+k2-1);
            const ptrdiff_t newnlo = k2 <= 0 ? k2-k1-1 : k1 < 0 ? -k1 : 0;
            const ptrdiff_t newnhi = k2 <= 0 ? 0 : k1 < 0 ? k2-1 : k2-k1-1;
            CheckSubBandMatrix(i1,i2,j1,j2,newnlo,newnhi,nlo(),nhi());
            return cSubBandMatrix(i1,i2,j1,j2,newnlo,newnhi);
        }



        // These check the indices according the the indexing style being
        // used, and then calls the above CStyle versions.
        TMV_INLINE_ND const_submatrix_type subMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2) const
        {
            CheckRowRange<_fort>(i1,i2,colsize());
            CheckColRange<_fort>(j1,j2,rowsize());
            CheckSubMatrix_Band(i1,i2,j1,j2,nlo(),nhi());
            return cSubMatrix(i1,i2,j1,j2);
        }

        TMV_INLINE_ND const_submatrix_step_type subMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2, ptrdiff_t istep, ptrdiff_t jstep) const
        {
            CheckRowRange<_fort>(i1,i2,colsize());
            CheckColRange<_fort>(j1,j2,rowsize());
            CheckSubMatrix_Band(i1,i2,j1,j2,istep,jstep,nlo(),nhi());
            return cSubMatrix(i1,i2,j1,j2,istep,jstep);
        }

        TMV_INLINE_ND const_subvector_type subVector(
            ptrdiff_t i, ptrdiff_t j, ptrdiff_t istep, ptrdiff_t jstep, ptrdiff_t s) const
        {
            CheckMatSubVector_Band<_fort>(
                i,j,istep,jstep,s,colsize(),rowsize(),nlo(),nhi());
            return cSubVector(i,j,istep,jstep,s); 
        }

        TMV_INLINE_ND const_subbandmatrix_type subBandMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2, ptrdiff_t newnlo, ptrdiff_t newnhi) const
        {
            CheckRowRange<_fort>(i1,i2,colsize());
            CheckColRange<_fort>(j1,j2,rowsize());
            CheckSubBandMatrix(i1,i2,j1,j2,newnlo,newnhi,nlo(),nhi());
            return cSubBandMatrix(i1,i2,j1,j2,newnlo,newnhi);
        }

        TMV_INLINE_ND const_subbandmatrix_type subBandMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2) const
        {
            CheckRowRange<_fort>(i1,i2,colsize());
            CheckColRange<_fort>(j1,j2,rowsize());
            return cSubBandMatrix(i1,i2,j1,j2);
        }

        TMV_INLINE_ND const_subbandmatrix_step_type subBandMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2, ptrdiff_t newnlo, ptrdiff_t newnhi,
            ptrdiff_t istep, ptrdiff_t jstep) const
        {
            CheckRowRange<_fort>(i1,i2,colsize());
            CheckColRange<_fort>(j1,j2,rowsize());
            CheckSubBandMatrix(
                i1,i2,j1,j2,newnlo,newnhi,istep,jstep,nlo(),nhi());
            return cSubBandMatrix(i1,i2,j1,j2,newnlo,newnhi,istep,jstep);
        }

        TMV_INLINE_ND const_colrange_type colRange(ptrdiff_t j1, ptrdiff_t j2) const
        {
            CheckColRange<_fort>(j1,j2,rowsize());
            return cColRange(j1,j2);
        }

        TMV_INLINE_ND const_rowrange_type rowRange(ptrdiff_t i1, ptrdiff_t i2) const
        {
            CheckRowRange<_fort>(i1,i2,colsize());
            return cRowRange(i1,i2);
        }

        TMV_INLINE_ND const_diagrange_type diagRange(ptrdiff_t k1, ptrdiff_t k2) const
        {
            CheckDiagRange_Band<_fort>(k1,k2,nlo(),nhi());
            return cDiagRange(k1,k2);
        }


        //
        // Views
        //

        TMV_INLINE TMV_MAYBE_CREF(type,const_view_type) view() const
        { return MakeBandView<type,const_view_type>::call(mat()); }

        TMV_INLINE TMV_MAYBE_CREF(type,const_cview_type) cView() const
        { return MakeBandView<type,const_cview_type>::call(mat()); }

        TMV_INLINE TMV_MAYBE_CREF(type,const_fview_type) fView() const
        { return MakeBandView<type,const_fview_type>::call(mat()); }

        TMV_INLINE TMV_MAYBE_CREF(type,const_xview_type) xView() const
        { return MakeBandView<type,const_xview_type>::call(mat()); }

        TMV_INLINE TMV_MAYBE_CREF(type,const_cmview_type) cmView() const
        { return MakeBandView<type,const_cmview_type>::call(mat()); }

        TMV_INLINE TMV_MAYBE_CREF(type,const_rmview_type) rmView() const
        { return MakeBandView<type,const_rmview_type>::call(mat()); }

        TMV_INLINE TMV_MAYBE_CREF(type,const_dmview_type) dmView() const
        { return MakeBandView<type,const_dmview_type>::call(mat()); }

        TMV_INLINE TMV_MAYBE_CREF(type,const_view_type) constView() const
        { return MakeBandView<type,const_view_type>::call(mat()); }

        TMV_INLINE const_transpose_type transpose() const
        {
            return const_transpose_type(
                cptr(),rowsize(),colsize(),nhi(),nlo(),stepj(),stepi()); 
        }

        TMV_INLINE TMV_MAYBE_CREF(type,const_conjugate_type) conjugate() const
        { return MakeBandView<type,const_conjugate_type>::call(mat()); }

        TMV_INLINE const_adjoint_type adjoint() const
        {
            return const_adjoint_type(
                cptr(),rowsize(),colsize(),nhi(),nlo(),stepj(),stepi()); 
        }

        TMV_INLINE const_upperband_type upperBand() const
        {
            return const_upperband_type(
                cptr(),TMV_MIN(colsize(),rowsize()),
                TMV_MIN(colsize()+nhi(),rowsize()),
                0,nhi(),stepi(),stepj());
        }

        TMV_INLINE const_lowerband_type lowerBand() const
        {
            return const_lowerband_type(
                cptr(),TMV_MIN(colsize(),rowsize()+nlo()),
                TMV_MIN(colsize(),rowsize()),
                nlo(),0,stepi(),stepj());
        }

        TMV_INLINE_ND const_upperbandoff_type upperBandOff() const
        {
            CheckUpperBandOff(nhi(),rowsize());
            return const_upperbandoff_type(
                cptr()+stepj(),TMV_MIN(colsize(),rowsize()-1),
                TMV_MIN(colsize()+nhi(),rowsize()-1),
                0,nhi()-1,stepi(),stepj());
        }

        TMV_INLINE_ND const_lowerbandoff_type lowerBandOff() const
        {
            CheckLowerBandOff(nlo(),colsize());
            return const_lowerbandoff_type(
                cptr()+stepi(),TMV_MIN(colsize()-1,rowsize()+nlo()),
                TMV_MIN(colsize()-1,rowsize()),
                nlo()-1,0,stepi(),stepj());
        }

        TMV_INLINE_ND const_upperbandoff_type upperBandOff(ptrdiff_t noff) const
        {
            CheckUpperBandOff(noff,nhi(),rowsize());
            return const_upperbandoff_type(
                cptr()+noff*stepj(),TMV_MIN(colsize(),rowsize()-noff),
                TMV_MIN(colsize()+nhi(),rowsize())-noff,
                0,nhi()-noff,stepi(),stepj());
        }

        TMV_INLINE_ND const_lowerbandoff_type lowerBandOff(ptrdiff_t noff) const
        {
            CheckLowerBandOff(noff,nlo(),colsize());
            return const_lowerbandoff_type(
                cptr()+noff*stepi(),TMV_MIN(colsize(),rowsize()+nlo())-noff,
                TMV_MIN(colsize()-noff,rowsize()),
                nlo()-noff,0,stepi(),stepj());
        }

        TMV_INLINE const_realpart_type realPart() const
        {
            const bool isreal = Traits<value_type>::isreal;
            return const_realpart_type(
                reinterpret_cast<const real_type*>(cptr()),
                colsize(), rowsize(), nlo(), nhi(),
                isreal ? stepi() : 2*stepi(), isreal ? stepj() : 2*stepj());
        }

        TMV_INLINE const_imagpart_type imagPart() const
        {
            const bool isreal = Traits<value_type>::isreal;
            TMVStaticAssert(Traits<value_type>::iscomplex);
            return const_imagpart_type(
                reinterpret_cast<const real_type*>(cptr())+1,
                colsize(), rowsize(), nlo(), nhi(),
                isreal ? stepi() : 2*stepi(), isreal ? stepj() : 2*stepj());
        }

        TMV_INLINE TMV_MAYBE_CREF(type,const_nonconj_type) nonConj() const
        { return MakeBandView<type,const_nonconj_type>::call(mat()); }

        TMV_INLINE nonconst_type nonConst() const
        {
            return nonconst_type(
                const_cast<value_type*>(cptr()),
                colsize(),rowsize(),nlo(),nhi(),stepi(),stepj());
        }


        //
        // Auxilliary routines
        //

        template <class M2>
        TMV_INLINE_ND void assignTo(BaseMatrix_Mutable<M2>& m2) const
        {
            TMVStaticAssert((ShapeTraits2<_shape,M2::_shape>::assignable));
            tmv::Copy(mat(),m2.mat()); 
        }

        TMV_INLINE const type& mat() const
        { return static_cast<const type&>(*this); }

        TMV_INLINE ptrdiff_t diagstep() const 
        { return _diagstep == Unknown ? stepi() + stepj() : _diagstep; }
        TMV_INLINE bool isconj() const { return _conj; }

        // Note that these last functions need to be defined in a more derived
        // class than this, or an infinite loop will result when compiling.
        // Also, cref from BaseMatrix.

        TMV_INLINE ptrdiff_t colsize() const { return mat().colsize(); }
        TMV_INLINE ptrdiff_t rowsize() const { return mat().rowsize(); }
        TMV_INLINE ptrdiff_t nlo() const { return mat().nlo(); }
        TMV_INLINE ptrdiff_t nhi() const { return mat().nhi(); }
        TMV_INLINE ptrdiff_t ls() const { return mat().ls(); }
        TMV_INLINE ptrdiff_t stepi() const { return mat().stepi(); }
        TMV_INLINE ptrdiff_t stepj() const { return mat().stepj(); }
        TMV_INLINE bool isrm() const { return mat().isrm(); }
        TMV_INLINE bool iscm() const { return mat().iscm(); }
        TMV_INLINE bool isdm() const { return mat().isdm(); }

        TMV_INLINE const value_type* cptr() const { return mat().cptr(); }
        TMV_INLINE const value_type* start_mem() const 
        { return mat().start_mem(); }

        ptrdiff_t nElements() const 
        { 
            const ptrdiff_t cs = colsize();
            const ptrdiff_t rs = rowsize();
            const ptrdiff_t lo = nlo();
            const ptrdiff_t hi = nhi();
            if (cs == 0 || rs == 0) return 0;
            else if (cs == rs) {
                return cs*(lo+hi+1) - (lo*(lo+1)/2) - (hi*(hi+1)/2);
            } else if (cs > rs) {
                // lox = number of subdiagonals that are clipped.
                ptrdiff_t lox = TMV_MAX(rs+lo-cs,ptrdiff_t(0));
                return rs*(lo+hi+1) - (lox*(lox+1)/2) - (hi*(hi+1)/2);
            } else {
                // hix = number of superdiagonals that are clipped.
                ptrdiff_t hix = TMV_MAX(cs+hi-rs,ptrdiff_t(0));
                return cs*(lo+hi+1) - (lo*(lo+1)/2) - (hix*(hix+1)/2);
            }
        }

        TMV_INLINE ptrdiff_t rowstart(ptrdiff_t i) const 
        { return TMV_MAX(ptrdiff_t(0),i-nlo()); }
        TMV_INLINE ptrdiff_t rowend(ptrdiff_t i) const 
        { return TMV_MIN(rowsize(),i+nhi()+1); }
        TMV_INLINE ptrdiff_t colstart(ptrdiff_t j) const 
        { return TMV_MAX(ptrdiff_t(0),j-nhi()); }
        TMV_INLINE ptrdiff_t colend(ptrdiff_t j) const 
        { return TMV_MIN(colsize(),j+nlo()+1); }

        TMV_INLINE const_rowmajor_iterator rowmajor_begin() const
        { return const_rowmajor_iterator(&mat(),0,0); }
        TMV_INLINE const_rowmajor_iterator rowmajor_end() const
        {
            return const_rowmajor_iterator(
                &mat(),TMV_MIN(colsize(),rowsize()+nlo()),
                rowstart(TMV_MIN(colsize(),rowsize()+nlo())));
        }
        TMV_INLINE const_colmajor_iterator colmajor_begin() const
        { return const_colmajor_iterator(&mat(),0,0); }
        TMV_INLINE const_colmajor_iterator colmajor_end() const
        {
            return const_colmajor_iterator(
                &mat(),colstart(TMV_MIN(rowsize(),colsize()+nhi())),
                TMV_MIN(rowsize(),colsize()+nhi()));
        }
        TMV_INLINE const_diagmajor_iterator diagmajor_begin() const
        { return const_diagmajor_iterator(&mat(),nlo(),0); }
        TMV_INLINE const_diagmajor_iterator diagmajor_end() const
        { return const_diagmajor_iterator(&mat(),0,nhi()+1); }

    }; // BaseMatrix_Band

    template <class M>
    class BaseMatrix_Band_Mutable : 
        public BaseMatrix_Band<M>,
        public BaseMatrix_Mutable<M>
    {
    public:
        enum { _colsize = Traits<M>::_colsize };
        enum { _rowsize = Traits<M>::_rowsize };
        enum { _nlo = Traits<M>::_nlo };
        enum { _nhi = Traits<M>::_nhi };
        enum { _shape = Traits<M>::_shape };
        enum { _fort = Traits<M>::_fort };
        enum { _calc = Traits<M>::_calc };
        enum { _rowmajor = Traits<M>::_rowmajor }; 
        enum { _colmajor = Traits<M>::_colmajor }; 
        enum { _diagmajor = Traits<M>::_diagmajor }; 
        enum { _stepi = Traits<M>::_stepi };
        enum { _stepj = Traits<M>::_stepj };
        enum { _diagstep = Traits<M>::_diagstep };
        enum { _conj = Traits<M>::_conj };
        enum { _canlin = Traits<M>::_canlin };

        typedef M type;
        typedef BaseMatrix_Band<M> base;
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

        typedef typename base::const_row_sub_type const_row_sub_type;
        typedef typename base::const_col_sub_type const_col_sub_type;
        typedef typename base::const_diag_type const_diag_type;
        typedef typename base::const_diag_sub_type const_diag_sub_type;
        typedef typename base::const_cmview_type const_cmview_type;
        typedef typename base::const_rmview_type const_rmview_type;
        typedef typename base::const_dmview_type const_dmview_type;
        typedef typename base::const_submatrix_type const_submatrix_type;
        typedef typename base::const_submatrix_step_type 
            const_submatrix_step_type;
        typedef typename base::const_subvector_type const_subvector_type;
        typedef typename base::const_subbandmatrix_type 
            const_subbandmatrix_type;
        typedef typename base::const_subbandmatrix_step_type 
            const_subbandmatrix_step_type;
        typedef typename base::const_colrange_type const_colrange_type;
        typedef typename base::const_rowrange_type const_rowrange_type;
        typedef typename base::const_diagrange_type const_diagrange_type;
        typedef typename base::const_upperband_type const_upperband_type;
        typedef typename base::const_lowerband_type const_lowerband_type;
        typedef typename base::const_upperbandoff_type const_upperbandoff_type;
        typedef typename base::const_lowerbandoff_type const_lowerbandoff_type;

        typedef typename base::const_rowmajor_iterator const_rowmajor_iterator;
        typedef typename base::const_colmajor_iterator const_colmajor_iterator;
        typedef typename base::const_diagmajor_iterator const_diagmajor_iterator;

        typedef typename Traits<M>::row_sub_type row_sub_type;
        typedef typename Traits<M>::col_sub_type col_sub_type;
        typedef typename Traits<M>::diag_type diag_type;
        typedef typename Traits<M>::diag_sub_type diag_sub_type;

        typedef typename Traits<M>::cmview_type cmview_type;
        typedef typename Traits<M>::rmview_type rmview_type;
        typedef typename Traits<M>::dmview_type dmview_type;

        typedef typename Traits<M>::submatrix_type submatrix_type;
        typedef typename Traits<M>::submatrix_step_type submatrix_step_type;

        typedef typename Traits<M>::subvector_type subvector_type;

        typedef typename Traits<M>::subbandmatrix_type subbandmatrix_type;
        typedef typename Traits<M>::subbandmatrix_step_type 
            subbandmatrix_step_type;

        typedef typename Traits<M>::colrange_type colrange_type;
        typedef typename Traits<M>::rowrange_type rowrange_type;
        typedef typename Traits<M>::diagrange_type diagrange_type;

        typedef typename Traits<M>::upperband_type upperband_type;
        typedef typename Traits<M>::lowerband_type lowerband_type;
        typedef typename Traits<M>::upperbandoff_type upperbandoff_type;
        typedef typename Traits<M>::lowerbandoff_type lowerbandoff_type;

        typedef typename Traits<M>::rowmajor_iterator rowmajor_iterator;
        typedef typename Traits<M>::colmajor_iterator colmajor_iterator;
        typedef typename Traits<M>::diagmajor_iterator diagmajor_iterator;


        //
        // Constructor
        //

    protected:
        TMV_INLINE BaseMatrix_Band_Mutable() {}
        TMV_INLINE BaseMatrix_Band_Mutable(const BaseMatrix_Band_Mutable<M>&) {}
        TMV_INLINE ~BaseMatrix_Band_Mutable() {}
    public:


        //
        // Access 
        //

        TMV_INLINE_ND reference operator()(ptrdiff_t i, ptrdiff_t j)
        {
            CheckRowIndex<_fort>(i,colsize());
            CheckColIndex<_fort>(j,rowsize());
            CheckInBand(i,j,nlo(),nhi());
            return ref(i,j);
        }

        // The get_ routines always use CStyle indexing.
        TMV_INLINE row_sub_type get_row(ptrdiff_t i, ptrdiff_t j1, ptrdiff_t j2) 
        { return row_sub_type(ptr()+i*stepi()+j1*stepj(),j2-j1,stepj()); }

        TMV_INLINE col_sub_type get_col(ptrdiff_t j, ptrdiff_t i1, ptrdiff_t i2) 
        { return col_sub_type(ptr()+j*stepj()+i1*stepi(),i2-i1,stepi()); }

        TMV_INLINE diag_sub_type get_diag(ptrdiff_t i) 
        {
            return diag_sub_type(
                ptr() + (i<0?(-i*stepi()):(i*stepj())),
                ( i<0 ? 
                  TMV_MIN(colsize()+i,rowsize()) :
                  TMV_MIN(colsize(),rowsize()-i) ),
                diagstep());
        }

        TMV_INLINE diag_sub_type get_diag(ptrdiff_t i, ptrdiff_t j1, ptrdiff_t j2) 
        {
            return diag_sub_type(
                ptr() + (i<0?(-i*stepi()):(i*stepj())) + j1*diagstep(),
                j2-j1, diagstep());
        }


        // The regular versions respect the indexing style for i and j:
        TMV_INLINE_ND row_sub_type row(ptrdiff_t i, ptrdiff_t j1, ptrdiff_t j2) 
        {
            CheckRowIndex<_fort>(i,colsize());
            CheckColRange_Band<_fort>(i,j1,j2,rowsize(),nlo(),nhi());
            return get_row(i,j1,j2); 
        }

        TMV_INLINE_ND col_sub_type col(ptrdiff_t j, ptrdiff_t i1, ptrdiff_t i2) 
        {
            CheckColIndex<_fort>(j,rowsize());
            CheckRowRange_Band<_fort>(j,i1,i2,colsize(),nlo(),nhi());
            return get_col(j,i1,i2); 
        }

        // No need for a get_ routine for diag()
        TMV_INLINE diag_type diag() 
        { return diag_type(ptr(),TMV_MIN(colsize(),rowsize()),diagstep()); }

        TMV_INLINE_ND diag_sub_type diag(ptrdiff_t i) 
        {
            CheckDiagIndex<_fort>(i,colsize(),rowsize());
            return get_diag(i);
        }

        TMV_INLINE_ND diag_sub_type diag(ptrdiff_t i, ptrdiff_t j1, ptrdiff_t j2) 
        {
            CheckDiagIndex<_fort>(i,j1,j2,colsize(),rowsize());
            return get_diag(i,j1,j2);
        }


        // We need to repeat the const versions so the non-const ones
        // don't clobber them.
        TMV_INLINE value_type operator()(ptrdiff_t i, ptrdiff_t j) const
        { return base::operator()(i,j); }

        TMV_INLINE const_diag_type diag() const
        { return base::diag(); }
        TMV_INLINE const_row_sub_type get_row(ptrdiff_t i, ptrdiff_t j1, ptrdiff_t j2) const
        { return base::get_row(i,j1,j2); }
        TMV_INLINE const_col_sub_type get_col(ptrdiff_t j, ptrdiff_t i1, ptrdiff_t i2) const
        { return base::get_col(j,i1,i2); }
        TMV_INLINE const_diag_sub_type get_diag(ptrdiff_t i) const
        { return base::get_diag(i); }
        TMV_INLINE const_diag_sub_type get_diag(ptrdiff_t i, ptrdiff_t j1, ptrdiff_t j2) const
        { return base::get_diag(i,j1,j2); }

        TMV_INLINE const_row_sub_type row(ptrdiff_t i, ptrdiff_t j1, ptrdiff_t j2) const
        { return base::row(i,j1,j2); }
        TMV_INLINE const_col_sub_type col(ptrdiff_t j, ptrdiff_t i1, ptrdiff_t i2) const
        { return base::col(j,i1,i2); }
        TMV_INLINE const_diag_sub_type diag(ptrdiff_t i) const
        { return base::diag(i); }
        TMV_INLINE const_diag_sub_type diag(ptrdiff_t i, ptrdiff_t j1, ptrdiff_t j2) const
        { return base::diag(i,j1,j2); }


        //
        // Op =
        //

        TMV_INLINE_ND type& operator=(const BaseMatrix_Band_Mutable<M>& m2) 
        {
            TMVAssert(colsize() == m2.colsize());
            TMVAssert(rowsize() == m2.rowsize());
            TMVAssert(nlo() >= m2.nlo());
            TMVAssert(nhi() >= m2.nhi());
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
            TMVStaticAssert((ShapeTraits2<M2::_shape,_shape>::assignable));
            m2.assignTo(mat());
            return mat(); 
        }

        template <class M2>
        type& operator=(const BaseMatrix_Tri<M2>& m2) 
        {
            TMVStaticAssert((Sizes<_colsize,M2::_colsize>::same));
            TMVStaticAssert((Sizes<_rowsize,M2::_rowsize>::same));
            TMVAssert(colsize() == m2.colsize());
            TMVAssert(rowsize() == m2.rowsize());
            TMVAssert((m2.isupper() ? nhi() : nlo()) == m2.size()-1);
            const int shape = M2::_upper ? UpperTri : LowerTri;
            const ptrdiff_t cs = _colsize;
            const ptrdiff_t rs = _rowsize;
            const ptrdiff_t si = _stepi;
            const ptrdiff_t sj = _stepj;
            const int A2 = 
                (_conj ? Conj : NonConj) |
                (_rowmajor ? RowMajor : _colmajor ? ColMajor : NonMajor);
            typename MViewHelper<value_type,shape,cs,rs,si,sj,A2>::type U(
                ptr(),rowsize(),stepi(),stepj(),NonUnitDiag);
            m2.assignTo(U);
            if (Maybe<M2::_upper>::select(nlo(),nhi()) > 0) 
                Maybe<!M2::_upper>::upperbandoff(mat()).setZero();
            return mat(); 
        }

        template <class M2>
        type& operator=(const BaseMatrix_Diag<M2>& m2) 
        {
            TMVStaticAssert((Sizes<_colsize,M2::_colsize>::same));
            TMVStaticAssert((Sizes<_rowsize,M2::_rowsize>::same));
            TMVAssert(colsize() == m2.colsize());
            TMVAssert(rowsize() == m2.rowsize());
            diag() = m2.diag();
            if (nhi() > 0) upperBandOff().setZero();
            if (nlo() > 0) lowerBandOff().setZero();
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

        // First the ones from BaseMatrix_Mutable:
        TMV_INLINE type& setZero()
        { tmv::SetZero(mat()); return mat(); }

        TMV_INLINE type& setAllTo(value_type val)
        { tmv::SetAllTo(mat(),val); return mat(); }

        TMV_INLINE type& addToAll(value_type val)
        { tmv::AddToAll(mat(),val); return mat(); }

        TMV_INLINE type& clip(float_type thresh)
        { tmv::Clip(mat(),thresh); return mat(); }

        template <class F>
        TMV_INLINE type& applyToAll(const F& f)
        { tmv::ApplyToAll(mat(),f); return mat(); }

        TMV_INLINE type& conjugateSelf()
        { tmv::ConjugateSelf(mat()); return mat(); }

        // Some more that are added for Band shape:
        TMV_INLINE type& transposeSelf() 
        { tmv::TransposeSelf(mat()); return mat(); }

        type& setToIdentity(const value_type x=value_type(1))
        {
            TMVAssert(colsize() == rowsize());
            this->setZero(); diag().setAllTo(x);
            return mat();
        }


        //
        // subMatrix, etc.
        //

        // These versions always uses CStyle
        TMV_INLINE submatrix_type cSubMatrix(ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2) 
        {
            return submatrix_type(
                ptr()+i1*stepi()+j1*stepj(), i2-i1, j2-j1, stepi(), stepj()); 
        }

        TMV_INLINE submatrix_step_type cSubMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2, ptrdiff_t istep, ptrdiff_t jstep) 
        {
            return submatrix_step_type(
                ptr()+i1*stepi()+j1*stepj(), (i2-i1)/istep, (j2-j1)/jstep,
                istep*stepi(), jstep*stepj());
        }

        TMV_INLINE subvector_type cSubVector(
            ptrdiff_t i, ptrdiff_t j, ptrdiff_t istep, ptrdiff_t jstep, ptrdiff_t s) 
        {
            return subvector_type(
                ptr()+i*stepi()+j*stepj(), s, istep*stepi() + jstep*stepj());
        }

        TMV_INLINE subbandmatrix_type cSubBandMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2, ptrdiff_t newnlo, ptrdiff_t newnhi)
        {
            return subbandmatrix_type(
                ptr()+i1*stepi()+j1*stepj(), i2-i1, j2-j1,
                newnlo, newnhi, stepi(), stepj()); 
        }

        subbandmatrix_type cSubBandMatrix(ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2) 
        {
            const ptrdiff_t newnlo = TMV_MIN(nlo()+j1-i1,i2-i1-1);
            const ptrdiff_t newnhi = TMV_MIN(nhi()+i1-j1,j2-j1-1);
            CheckSubBandMatrix(i1,i2,j1,j2,newnlo,newnhi,nlo(),nhi());
            return cSubBandMatrix(i1,i2,j1,j2,newnlo,newnhi);
        }

        TMV_INLINE subbandmatrix_step_type cSubBandMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2, ptrdiff_t newnlo, ptrdiff_t newnhi,
            ptrdiff_t istep, ptrdiff_t jstep) 
        {
            return subbandmatrix_step_type(
                ptr()+i1*stepi()+j1*stepj(), (i2-i1)/istep, (j2-j1)/jstep,
                newnlo, newnhi, istep*stepi(), jstep*stepj());
        }

        colrange_type cColRange(ptrdiff_t j1, ptrdiff_t j2) 
        {
            const ptrdiff_t i1 = j1 > nhi() ? j1-nhi() : 0;
            const ptrdiff_t i2 = TMV_MIN(j2 + nlo(),colsize());
            const ptrdiff_t newnhi = j1 < nhi() ? TMV_MIN(nhi(),j2-1) - j1 : 0;
            const ptrdiff_t newnlo = i1==i2 ? 0 : TMV_MIN(nlo()+nhi()-newnhi,i2-i1-1);
            CheckSubBandMatrix(i1,i2,j1,j2,newnlo,newnhi,nlo(),nhi());
            return cSubBandMatrix(i1,i2,j1,j2,newnlo,newnhi);
        }

        rowrange_type cRowRange(ptrdiff_t i1, ptrdiff_t i2) 
        {
            const ptrdiff_t j1 = i1 > nlo() ? i1-nlo() : 0;
            const ptrdiff_t j2 = TMV_MIN(i2 + nhi(),rowsize());
            const ptrdiff_t newnlo = i1 < nlo() ? TMV_MIN(nlo(),i2-1) - i1 : 0;
            const ptrdiff_t newnhi = j1==j2 ? 0 : TMV_MIN(nlo()+nhi()-newnlo,j2-j1-1);
            CheckSubBandMatrix(i1,i2,j1,j2,newnlo,newnhi,nlo(),nhi());
            return cSubBandMatrix(i1,i2,j1,j2,newnlo,newnhi);
        }

        diagrange_type cDiagRange(ptrdiff_t k1, ptrdiff_t k2) 
        {
            const ptrdiff_t i1 = k2 <= 0 ? -k2+1 : 0;
            const ptrdiff_t i2 = TMV_MIN(rowsize()-k1,colsize());
            const ptrdiff_t j1 = k1 <= 0 ? 0 : k1;
            const ptrdiff_t j2 = TMV_MIN(rowsize(),colsize()+k2-1);
            const ptrdiff_t newnlo = k2 <= 0 ? k2-k1-1 : k1 < 0 ? -k1 : 0;
            const ptrdiff_t newnhi = k2 <= 0 ? 0 : k1 < 0 ? k2-1 : k2-k1-1;
            CheckSubBandMatrix(i1,i2,j1,j2,newnlo,newnhi,nlo(),nhi());
            return cSubBandMatrix(i1,i2,j1,j2,newnlo,newnhi);
        }



        // These check the indices according the the indexing style being
        // used, and then calls the above CStyle versions.
        TMV_INLINE_ND submatrix_type subMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2) 
        {
            CheckRowRange<_fort>(i1,i2,colsize());
            CheckColRange<_fort>(j1,j2,rowsize());
            CheckSubMatrix_Band(i1,i2,j1,j2,nlo(),nhi());
            return cSubMatrix(i1,i2,j1,j2);
        }

        TMV_INLINE_ND submatrix_step_type subMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2, ptrdiff_t istep, ptrdiff_t jstep) 
        {
            CheckRowRange<_fort>(i1,i2,istep,colsize());
            CheckColRange<_fort>(j1,j2,jstep,rowsize());
            CheckSubMatrix_Band(i1,i2,j1,j2,istep,jstep,nlo(),nhi());
            return cSubMatrix(i1,i2,j1,j2,istep,jstep);
        }

        TMV_INLINE_ND subvector_type subVector(
            ptrdiff_t i, ptrdiff_t j, ptrdiff_t istep, ptrdiff_t jstep, ptrdiff_t s) 
        {
            CheckMatSubVector_Band<_fort>(
                i,j,istep,jstep,s,colsize(),rowsize());
            return cSubVector(i,j,istep,jstep,s); 
        }

        TMV_INLINE_ND subbandmatrix_type subBandMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2, ptrdiff_t newnlo, ptrdiff_t newnhi) 
        {
            CheckRowRange<_fort>(i1,i2,colsize());
            CheckColRange<_fort>(j1,j2,rowsize());
            CheckSubBandMatrix(i1,i2,j1,j2,newnlo,newnhi,nlo(),nhi());
            return cSubBandMatrix(i1,i2,j1,j2,newnlo,newnhi);
        }

        TMV_INLINE_ND subbandmatrix_type subBandMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2) 
        {
            CheckRowRange<_fort>(i1,i2,colsize());
            CheckColRange<_fort>(j1,j2,rowsize());
            return cSubBandMatrix(i1,i2,j1,j2);
        }

        TMV_INLINE_ND subbandmatrix_step_type subBandMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2, ptrdiff_t newnlo, ptrdiff_t newnhi,
            ptrdiff_t istep, ptrdiff_t jstep) 
        {
            CheckRowRange<_fort>(i1,i2,colsize());
            CheckColRange<_fort>(j1,j2,rowsize());
            CheckSubBandMatrix(
                i1,i2,j1,j2,newnlo,newnhi,istep,jstep,nlo(),nhi());
            return cSubBandMatrix(i1,i2,j1,j2,newnlo,newnhi,istep,jstep);
        }

        TMV_INLINE_ND colrange_type colRange(ptrdiff_t j1, ptrdiff_t j2) 
        {
            CheckColRange<_fort>(j1,j2,rowsize());
            return cColRange(j1,j2);
        }

        TMV_INLINE_ND rowrange_type rowRange(ptrdiff_t i1, ptrdiff_t i2) 
        {
            CheckRowRange<_fort>(i1,i2,colsize());
            return cRowRange(i1,i2);
        }

        TMV_INLINE_ND diagrange_type diagRange(ptrdiff_t k1, ptrdiff_t k2) 
        {
            CheckDiagRange_Band<_fort>(k1,k2,nlo(),nhi());
            return cDiagRange(k1,k2);
        }


        // Repeat the const versions:
        TMV_INLINE const_submatrix_type cSubMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2) const
        { return base::cSubMatrix(i1,i2,j1,j2); }
        TMV_INLINE const_submatrix_step_type cSubMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2, ptrdiff_t istep, ptrdiff_t jstep) const
        { return base::cSubMatrix(i1,i2,j1,j2,istep,jstep); }
        TMV_INLINE const_subvector_type cSubVector(
            ptrdiff_t i, ptrdiff_t j, ptrdiff_t istep, ptrdiff_t jstep, ptrdiff_t s) const
        { return base::cSubVector(i,j,istep,jstep,s); }
        TMV_INLINE const_subbandmatrix_type cSubBandMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2, ptrdiff_t newnlo, ptrdiff_t newnhi) const
        { return base::cSubBandMatrix(i1,i2,j1,j2,newnlo,newnhi); }
        TMV_INLINE const_subbandmatrix_step_type cSubBandMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2, ptrdiff_t newnlo, ptrdiff_t newnhi,
            ptrdiff_t istep, ptrdiff_t jstep) const
        { return base::cSubBandMatrix(i1,i2,j1,j2,newnlo,newnhi,istep,jstep); }

        TMV_INLINE const_submatrix_type subMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2) const
        { return base::subMatrix(i1,i2,j1,j2); }
        TMV_INLINE const_submatrix_step_type subMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2, ptrdiff_t istep, ptrdiff_t jstep) const
        { return base::subMatrix(i1,i2,j1,j2,istep,jstep); }
        TMV_INLINE const_subvector_type subVector(
            ptrdiff_t i, ptrdiff_t j, ptrdiff_t istep, ptrdiff_t jstep, ptrdiff_t s) const
        { return base::subVector(i,j,istep,jstep,s); }
        TMV_INLINE const_subbandmatrix_type subBandMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2, ptrdiff_t newnlo, ptrdiff_t newnhi) const
        { return base::subBandMatrix(i1,i2,j1,j2,newnlo,newnhi); }
        TMV_INLINE const_subbandmatrix_type subBandMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2) const
        { return base::subBandMatrix(i1,i2,j1,j2); }
        TMV_INLINE const_subbandmatrix_type subBandMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2, ptrdiff_t newnlo, ptrdiff_t newnhi,
            ptrdiff_t istep, ptrdiff_t jstep) const
        { return base::subBandMatrix(i1,i2,j1,j2,newnlo,newnhi,istep,jstep); }
        TMV_INLINE const_colrange_type colRange(ptrdiff_t j1, ptrdiff_t j2) const
        { return base::colRange(j1,j2); }
        TMV_INLINE const_rowrange_type rowRange(ptrdiff_t i1, ptrdiff_t i2) const
        { return base::rowRange(i1,i2); }
        TMV_INLINE const_diagrange_type diagRange(ptrdiff_t k1, ptrdiff_t k2) const
        { return base::diagRange(k1,k2); }


        //
        // Views
        //

        TMV_INLINE TMV_MAYBE_REF(type,view_type) view() 
        { return MakeBandView<type,view_type>::call(mat()); }

        TMV_INLINE TMV_MAYBE_REF(type,cview_type) cView() 
        { return MakeBandView<type,cview_type>::call(mat()); }

        TMV_INLINE TMV_MAYBE_REF(type,fview_type) fView() 
        { return MakeBandView<type,fview_type>::call(mat()); }

        TMV_INLINE TMV_MAYBE_REF(type,xview_type) xView() 
        { return MakeBandView<type,xview_type>::call(mat()); }

        TMV_INLINE_ND TMV_MAYBE_REF(type,cmview_type) cmView() 
        {
            TMVAssert(iscm() && "Called cmView on non-ColMajor matrix");
            return MakeBandView<type,cmview_type>::call(mat()); 
        }

        TMV_INLINE_ND TMV_MAYBE_REF(type,rmview_type) rmView() 
        {
            TMVAssert(isrm() && "Called rmView on non-RowMajor matrix");
            return MakeBandView<type,rmview_type>::call(mat()); 
        }

        TMV_INLINE_ND TMV_MAYBE_REF(type,dmview_type) dmView() 
        {
            TMVAssert(isdm() && "Called rmView on non-DiagMajor matrix");
            return MakeBandView<type,dmview_type>::call(mat()); 
        }

        TMV_INLINE transpose_type transpose() 
        {
            return transpose_type(
                ptr(),rowsize(),colsize(),nhi(),nlo(),stepj(),stepi()); 
        }

        TMV_INLINE TMV_MAYBE_REF(type,conjugate_type) conjugate() 
        { return MakeBandView<type,conjugate_type>::call(mat()); }

        TMV_INLINE adjoint_type adjoint() 
        {
            return adjoint_type(
                ptr(),rowsize(),colsize(),nhi(),nlo(),stepj(),stepi()); 
        }

        TMV_INLINE upperband_type upperBand() 
        {
            return upperband_type(
                ptr(),TMV_MIN(colsize(),rowsize()),
                TMV_MIN(colsize()+nhi(),rowsize()),
                0,nhi(),stepi(),stepj());
        }

        TMV_INLINE lowerband_type lowerBand() 
        {
            return lowerband_type(
                ptr(),TMV_MIN(colsize(),rowsize()+nlo()),
                TMV_MIN(colsize(),rowsize()),
                nlo(),0,stepi(),stepj());
        }

        TMV_INLINE_ND upperbandoff_type upperBandOff() 
        {
            CheckUpperBandOff(nhi(),rowsize());
            return upperbandoff_type(
                ptr()+stepj(),TMV_MIN(colsize(),rowsize()-1),
                TMV_MIN(colsize()+nhi(),rowsize()-1),
                0,nhi()-1,stepi(),stepj());
        }

        TMV_INLINE_ND lowerbandoff_type lowerBandOff() 
        {
            CheckLowerBandOff(nlo(),colsize());
            return lowerbandoff_type(
                ptr()+stepi(),TMV_MIN(colsize()-1,rowsize()+nlo()),
                TMV_MIN(colsize()-1,rowsize()),
                nlo()-1,0,stepi(),stepj());
        }

        TMV_INLINE_ND upperbandoff_type upperBandOff(ptrdiff_t noff) 
        {
            CheckUpperBandOff(noff,nhi(),rowsize());
            return upperbandoff_type(
                ptr()+noff*stepj(),TMV_MIN(colsize(),rowsize()-noff),
                TMV_MIN(colsize()+nhi(),rowsize())-noff,
                0,nhi()-noff,stepi(),stepj());
        }

        TMV_INLINE_ND lowerbandoff_type lowerBandOff(ptrdiff_t noff) 
        {
            CheckLowerBandOff(noff,nlo(),colsize());
            return lowerbandoff_type(
                ptr()+noff*stepi(),TMV_MIN(colsize(),rowsize()+nlo())-noff,
                TMV_MIN(colsize()-noff,rowsize()),
                nlo()-noff,0,stepi(),stepj());
        }

        TMV_INLINE realpart_type realPart() 
        {
            const bool isreal = Traits<value_type>::isreal;
            return realpart_type(
                reinterpret_cast<real_type*>(ptr()), 
                colsize(), rowsize(), nlo(), nhi(),
                isreal ? stepi() : 2*stepi(), isreal ? stepj() : 2*stepj());
        }

        TMV_INLINE imagpart_type imagPart() 
        {
            const bool isreal = Traits<value_type>::isreal;
            TMVStaticAssert(Traits<value_type>::iscomplex);
            return imagpart_type(
                reinterpret_cast<real_type*>(ptr())+1,
                colsize(), rowsize(), nlo(), nhi(),
                isreal ? stepi() : 2*stepi(), isreal ? stepj() : 2*stepj());
        }

        TMV_INLINE TMV_MAYBE_REF(type,nonconj_type) nonConj()
        { return MakeBandView<type,nonconj_type>::call(mat()); }

        TMV_INLINE TMV_MAYBE_REF(type,noalias_type) noAlias()
        { return MakeBandView<type,noalias_type>::call(mat()); }

        TMV_INLINE TMV_MAYBE_REF(type,alias_type) alias()
        { return MakeBandView<type,alias_type>::call(mat()); }


        // Repeat the const versions:
        TMV_INLINE TMV_MAYBE_CREF(type,const_view_type) view() const
        { return base::view(); }
        TMV_INLINE TMV_MAYBE_CREF(type,const_cview_type) cView() const
        { return base::cView(); }
        TMV_INLINE TMV_MAYBE_CREF(type,const_fview_type) fView() const
        { return base::fView(); }
        TMV_INLINE TMV_MAYBE_CREF(type,const_xview_type) xView() const
        { return base::xView(); }
        TMV_INLINE TMV_MAYBE_CREF(type,const_cmview_type) cmView() const
        { return base::cmView(); }
        TMV_INLINE TMV_MAYBE_CREF(type,const_rmview_type) rmView() const
        { return base::rmView(); }
        TMV_INLINE TMV_MAYBE_CREF(type,const_dmview_type) dmView() const
        { return base::dmView(); }
        TMV_INLINE const_transpose_type transpose() const
        { return base::transpose(); }
        TMV_INLINE TMV_MAYBE_CREF(type,const_conjugate_type) conjugate() const
        { return base::conjugate(); }
        TMV_INLINE const_adjoint_type adjoint() const
        { return base::adjoint(); }
        TMV_INLINE const_upperband_type upperBand() const
        { return base::upperBand(); }
        TMV_INLINE const_lowerband_type lowerBand() const
        { return base::lowerBand(); }
        TMV_INLINE const_upperbandoff_type upperBandOff() const
        { return base::upperBandOff(); }
        TMV_INLINE const_lowerbandoff_type lowerBandOff() const
        { return base::lowerBandOff(); }
        TMV_INLINE const_upperbandoff_type upperBandOff(ptrdiff_t noff) const
        { return base::upperBandOff(noff); }
        TMV_INLINE const_lowerbandoff_type lowerBandOff(ptrdiff_t noff) const
        { return base::lowerBandOff(noff); }
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

        TMV_INLINE bool isconj() const { return _conj; }
        TMV_INLINE ptrdiff_t diagstep() const 
        { return _diagstep == Unknown ? stepi() + stepj() : _diagstep; }

        // Note that these last functions need to be defined in a more derived
        // class than this, or an infinite loop will result when compiling.
        // Also, cref and cptr from above.

        TMV_INLINE ptrdiff_t colsize() const { return mat().colsize(); }
        TMV_INLINE ptrdiff_t rowsize() const { return mat().rowsize(); }
        TMV_INLINE ptrdiff_t nlo() const { return mat().nlo(); }
        TMV_INLINE ptrdiff_t nhi() const { return mat().nhi(); }
        TMV_INLINE ptrdiff_t ls() const { return mat().ls(); }
        TMV_INLINE ptrdiff_t stepi() const { return mat().stepi(); }
        TMV_INLINE ptrdiff_t stepj() const { return mat().stepj(); }
        TMV_INLINE bool isrm() const { return mat().isrm(); }
        TMV_INLINE bool iscm() const { return mat().iscm(); }
        TMV_INLINE bool isdm() const { return mat().isdm(); }

        TMV_INLINE value_type* ptr() { return mat().ptr(); }
        TMV_INLINE reference ref(ptrdiff_t i, ptrdiff_t j) { return mat().ref(i,j); }
        TMV_INLINE value_type* start_mem() { return mat().start_mem(); }

        TMV_INLINE rowmajor_iterator rowmajor_begin() 
        { return rowmajor_iterator(&mat(),0,0); }
        TMV_INLINE rowmajor_iterator rowmajor_end() 
        {
            return rowmajor_iterator(
                &mat(),TMV_MIN(colsize(),rowsize()+nlo()),
                this->rowstart(TMV_MIN(colsize(),rowsize()+nlo())));
        }
        TMV_INLINE colmajor_iterator colmajor_begin() 
        { return colmajor_iterator(&mat(),0,0); }
        TMV_INLINE colmajor_iterator colmajor_end() 
        {
            return colmajor_iterator(
                &mat(),this->colstart(TMV_MIN(rowsize(),colsize()+nhi())),
                TMV_MIN(rowsize(),colsize()+nhi()));
        }
        TMV_INLINE diagmajor_iterator diagmajor_begin() 
        { return diagmajor_iterator(&mat(),nlo(),0); }
        TMV_INLINE diagmajor_iterator diagmajor_end() 
        { return diagmajor_iterator(&mat(),0,nhi()+1); }

        TMV_INLINE const_rowmajor_iterator rowmajor_begin() const
        { return base::rowmajor_begin(); }
        TMV_INLINE const_rowmajor_iterator rowmajor_end() const
        { return base::rowmajor_end(); }
        TMV_INLINE const_colmajor_iterator colmajor_begin() const
        { return base::colmajor_begin(); }
        TMV_INLINE const_colmajor_iterator colmajor_end() const
        { return base::colmajor_end(); }
        TMV_INLINE const_diagmajor_iterator diagmajor_begin() const
        { return base::diagmajor_begin(); }
        TMV_INLINE const_diagmajor_iterator diagmajor_end() const
        { return base::diagmajor_end(); }

    }; // BaseMatrix_Band_Mutable


    //
    // setZero
    //

    template <int algo, class M>
    struct SetZeroB_Helper;

    // algo 1: Linearize to vector version
    template <class M>
    struct SetZeroB_Helper<1,M>
    {
        static inline void call(M& m) 
        { m.linearView().setZero(); } 
    };

    // algo 2: RowMajor
    template <class M>
    struct SetZeroB_Helper<2,M>
    {
        static void call(M& m1) 
        {
            const ptrdiff_t m = m1.colsize();
            const ptrdiff_t n = m1.rowsize();
            if (m > 0 && n > 0) {
                ptrdiff_t j1=0;
                ptrdiff_t j2=m1.nhi()+1;
                ptrdiff_t k=m1.nlo();
                for(ptrdiff_t i=0;i<m;++i) {
                    m1.get_row(i,j1,j2).setZero();
                    if (k>0) --k; else ++j1;
                    if (j2<n) ++j2;
                    else if (j1==n) break;
                }
            }
        }
    };

    // algo 3: ColMajor
    template <class M>
    struct SetZeroB_Helper<3,M>
    {
        static void call(M& m1) 
        {
            const ptrdiff_t m = m1.colsize();
            const ptrdiff_t n = m1.rowsize();
            if (m > 0 && n > 0) {
                ptrdiff_t i1=0;
                ptrdiff_t i2=m1.nlo()+1;
                ptrdiff_t k=m1.nhi();
                for(ptrdiff_t j=0;j<n;++j) {
                    m1.get_col(j,i1,i2).setZero();
                    if (k>0) --k; else ++i1;
                    if (i2<m) ++i2;
                    else if (i1==m) break;
                }
            }
        } 
    };

    // algo 4: DiagMajor
    template <class M>
    struct SetZeroB_Helper<4,M>
    {
        static void call(M& m1) 
        {
            for (ptrdiff_t k=-m1.nlo();k<=m1.nhi();++k) 
                m1.diag(k).setZero();
        } 
    };

    template <class M>
    TMV_INLINE void SetZero(BaseMatrix_Band_Mutable<M>& m)
    {
        const int algo = (
            M::_canlin ? 1 :
            M::_rowmajor ? 2 : 
            M::_colmajor ? 3 : 
            4 );
        SetZeroB_Helper<algo,M>::call(m.mat());
    }

    //
    // setAllTo
    //

    template <int algo, class M, class T>
    struct SetAllToB_Helper;

    // algo 1: Linearize to vector version
    template <class M, class T>
    struct SetAllToB_Helper<1,M,T> // algo 1, linearize
    {
        static inline void call(M& m, const T& val) 
        { m.linearView().setAllTo(val); } 
    };

    // algo 2: RowMajor
    template <class M, class T>
    struct SetAllToB_Helper<2,M,T>
    {
        static void call(M& m1, const T& val) 
        {
            const ptrdiff_t m = m1.colsize();
            const ptrdiff_t n = m1.rowsize();
            if (m > 0 && n > 0) {
                ptrdiff_t j1=0;
                ptrdiff_t j2=m1.nhi()+1;
                ptrdiff_t k=m1.nlo();
                for(ptrdiff_t i=0;i<m;++i) {
                    m1.get_row(i,j1,j2).setAllTo(val);
                    if (k>0) --k; else ++j1;
                    if (j2<n) ++j2;
                    else if (j1==n) break;
                }
            }
        }
    };

    // algo 3: ColMajor
    template <class M, class T>
    struct SetAllToB_Helper<3,M,T>
    {
        static void call(M& m1, const T& val) 
        {
            const ptrdiff_t m = m1.colsize();
            const ptrdiff_t n = m1.rowsize();
            if (m > 0 && n > 0) {
                ptrdiff_t i1=0;
                ptrdiff_t i2=m1.nlo()+1;
                ptrdiff_t k=m1.nhi();
                for(ptrdiff_t j=0;j<n;++j) {
                    m1.get_col(j,i1,i2).setAllTo(val);
                    if (k>0) --k; else ++i1;
                    if (i2<m) ++i2;
                    else if (i1==m) break;
                }
            }
        }
    };

    // algo 4: DiagMajor
    template <class M, class T>
    struct SetAllToB_Helper<4,M,T>
    {
        static void call(M& m1, const T& val) 
        {
            for (ptrdiff_t k=-m1.nlo();k<=m1.nhi();++k) 
                m1.diag(k).setAllTo(val);
        } 
    };

    template <class M, class T>
    TMV_INLINE void SetAllTo(BaseMatrix_Band_Mutable<M>& m, const T& val)
    {
        const int algo = (
            M::_canlin ? 1 :
            M::_rowmajor ? 2 : 
            M::_colmajor ? 3 : 
            4 );
        SetAllToB_Helper<algo,M,T>::call(m.mat(),val);
    }

    //
    // addToAll
    //

    template <int algo, class M, class T>
    struct AddToAllB_Helper;

    // algo 1: Linearize to vector version
    template <class M, class T>
    struct AddToAllB_Helper<1,M,T> // algo 1, linearize
    {
        static inline void call(M& m, const T& val) 
        { m.linearView().addToAll(val); } 
    };

    // algo 2: RowMajor
    template <class M, class T>
    struct AddToAllB_Helper<2,M,T>
    {
        static void call(M& m1, const T& val) 
        {
            const ptrdiff_t m = m1.colsize();
            const ptrdiff_t n = m1.rowsize();
            if (m > 0 && n > 0) {
                ptrdiff_t j1=0;
                ptrdiff_t j2=m1.nhi()+1;
                ptrdiff_t k=m1.nlo();
                for(ptrdiff_t i=0;i<m;++i) {
                    m1.get_row(i,j1,j2).addToAll(val);
                    if (k>0) --k; else ++j1;
                    if (j2<n) ++j2;
                    else if (j1==n) break;
                }
            }
        }
    };

    // algo 3: ColMajor
    template <class M, class T>
    struct AddToAllB_Helper<3,M,T>
    {
        static void call(M& m1, const T& val) 
        {
            const ptrdiff_t m = m1.colsize();
            const ptrdiff_t n = m1.rowsize();
            if (m > 0 && n > 0) {
                ptrdiff_t i1=0;
                ptrdiff_t i2=m1.nlo()+1;
                ptrdiff_t k=m1.nhi();
                for(ptrdiff_t j=0;j<n;++j) {
                    m1.get_col(j,i1,i2).addToAll(val);
                    if (k>0) --k; else ++i1;
                    if (i2<m) ++i2;
                    else if (i1==m) break;
                }
            }
        }
    };

    // algo 4: DiagMajor
    template <class M, class T>
    struct AddToAllB_Helper<4,M,T>
    {
        static void call(M& m1, const T& val) 
        {
            for (ptrdiff_t k=-m1.nlo();k<=m1.nhi();++k) 
                m1.diag(k).AddToAll(val);
        } 
    };

    template <class M, class T>
    TMV_INLINE void AddToAll(BaseMatrix_Band_Mutable<M>& m, const T& val)
    {
        const int algo = (
            M::_canlin ? 1 :
            M::_rowmajor ? 2 : 
            M::_colmajor ? 3 : 
            4 );
        AddToAllB_Helper<algo,M,T>::call(m.mat(),val);
    }

    //
    // Clip
    //

    template <int algo, class M, class RT>
    struct ClipB_Helper;

    // algo 1: Linearize to vector version
    template <class M, class RT>
    struct ClipB_Helper<1,M,RT> // algo 1, linearize
    {
        static inline void call(M& m, const RT& thresh) 
        { m.linearView().clip(thresh); } 
    };

    // algo 2: RowMajor
    template <class M, class RT>
    struct ClipB_Helper<2,M,RT>
    {
        static void call(M& m1, const RT& thresh) 
        {
            const ptrdiff_t m = m1.colsize();
            const ptrdiff_t n = m1.rowsize();
            if (m > 0 && n > 0) {
                ptrdiff_t j1=0;
                ptrdiff_t j2=m1.nhi()+1;
                ptrdiff_t k=m1.nlo();
                for(ptrdiff_t i=0;i<m;++i) {
                    m1.get_row(i,j1,j2).clip(thresh);
                    if (k>0) --k; else ++j1;
                    if (j2<n) ++j2;
                    else if (j1==n) break;
                }
            }
        }
    };

    // algo 3: ColMajor
    template <class M, class RT>
    struct ClipB_Helper<3,M,RT>
    {
        static void call(M& m1, const RT& thresh) 
        {
            const ptrdiff_t m = m1.colsize();
            const ptrdiff_t n = m1.rowsize();
            if (m > 0 && n > 0) {
                ptrdiff_t i1=0;
                ptrdiff_t i2=m1.nlo()+1;
                ptrdiff_t k=m1.nhi();
                for(ptrdiff_t j=0;j<n;++j) {
                    m1.get_col(j,i1,i2).clip(thresh);
                    if (k>0) --k; else ++i1;
                    if (i2<m) ++i2;
                    else if (i1==m) break;
                }
            }
        }
    };

    // algo 4: DiagMajor
    template <class M, class RT>
    struct ClipB_Helper<4,M,RT>
    {
        static void call(M& m1, const RT& thresh) 
        {
            for (ptrdiff_t k=-m1.nlo();k<=m1.nhi();++k) 
                m1.diag(k).clip(thresh);
        } 
    };

    template <class M, class RT>
    TMV_INLINE void Clip(BaseMatrix_Band_Mutable<M>& m, const RT& thresh)
    {
        const int algo = (
            M::_canlin ? 1 :
            M::_rowmajor ? 2 : 
            M::_colmajor ? 3 : 
            4 );
        ClipB_Helper<algo,M,RT>::call(m.mat(),thresh);
    }

    //
    // applyToAll
    //

    template <int algo, class M, class F>
    struct ApplyToAllB_Helper;

    // algo 1: Linearize to vector version
    template <class M, class F>
    struct ApplyToAllB_Helper<1,M,F> // algo 1, linearize
    {
        static inline void call(M& m, const F& f) 
        { m.linearView().applyToAll(f); } 
    };

    // algo 2: RowMajor
    template <class M, class F>
    struct ApplyToAllB_Helper<2,M,F>
    {
        static inline void call(M& m1, const F& f) 
        {
            const ptrdiff_t m = m1.colsize();
            const ptrdiff_t n = m1.rowsize();
            if (m > 0 && n > 0) {
                ptrdiff_t j1=0;
                ptrdiff_t j2=m1.nhi()+1;
                ptrdiff_t k=m1.nlo();
                for(ptrdiff_t i=0;i<m;++i) {
                    m1.get_row(i,j1,j2).applyToAll(f);
                    if (k>0) --k; else ++j1;
                    if (j2<n) ++j2;
                    else if (j1==n) break;
                }
            }
        }
    };

    // algo 3: ColMajor
    template <class M, class F>
    struct ApplyToAllB_Helper<3,M,F>
    {
        static inline void call(M& m1, const F& f) 
        {
            const ptrdiff_t m = m1.colsize();
            const ptrdiff_t n = m1.rowsize();
            if (m > 0 && n > 0) {
                ptrdiff_t i1=0;
                ptrdiff_t i2=m1.nlo()+1;
                ptrdiff_t k=m1.nhi();
                for(ptrdiff_t j=0;j<n;++j) {
                    m1.get_col(j,i1,i2).applyToAll(f);
                    if (k>0) --k; else ++i1;
                    if (i2<m) ++i2;
                    else if (i1==m) break;
                }
            }
        }
    };

    // algo 4: DiagMajor
    template <class M, class F>
    struct ApplyToAllB_Helper<4,M,F>
    {
        static void call(M& m1, const F& f) 
        {
            for (ptrdiff_t k=-m1.nlo();k<=m1.nhi();++k) 
                m1.diag(k).applyToAll(f);
        } 
    };

    template <class M, class F>
    TMV_INLINE void ApplyToAll(BaseMatrix_Band_Mutable<M>& m, const F& f)
    {
        const int algo = (
            M::_canlin ? 1 :
            M::_rowmajor ? 2 : 
            M::_colmajor ? 3 : 
            4 );
        ApplyToAllB_Helper<algo,M,F>::call(m.mat(),f);
    }

    //
    // ConjugateSelf
    //

    template <int algo, class M>
    struct ConjugateB_Helper;

    // algo 0: Not complex, nothing to do
    template <class M>
    struct ConjugateB_Helper<0,M>
    { static inline void call(M& ) {} };

    // algo 1: Linearize to vector version
    template <class M>
    struct ConjugateB_Helper<1,M>
    {
        static inline void call(M& m)
        { m.linearView().conjugateSelf(); }
    };

    // In TMV_ScaleB.h
    template <int ix, class T, class M>
    inline void Scale(
        const Scaling<ix,T>& x, BaseMatrix_Band_Mutable<M>& m);

    // algo 2: m.imagPart() *= -1
    template <class M>
    struct ConjugateB_Helper<2,M>
    {
        static inline void call(M& m)
        {
            typedef typename M::real_type RT;
            typedef typename M::imagpart_type Mi;
            const Scaling<-1,RT> mone;
            Mi mi = m.imagPart();
            Scale(mone,mi);
        }
    };

    template <class M>
    TMV_INLINE void ConjugateSelf(BaseMatrix_Band_Mutable<M>& m)
    {
        const bool isreal = Traits<typename M::value_type>::isreal;
        const int algo = isreal ? 0 : M::_canlin ? 1 : 2;
        ConjugateB_Helper<algo,M>::call(m.mat());
    }



    //
    // TMV_Text 
    //

    template <class M>
    inline std::string TMV_Text(const BaseMatrix_Band<M>& m)
    {
        std::ostringstream s;
        s << "BaseMatrix_Band< "<<TMV_Text(m.mat())<<" >";
        return s.str();
    }

    template <class M>
    inline std::string TMV_Text(const BaseMatrix_Band_Mutable<M>& m)
    {
        std::ostringstream s;
        s << "BaseMatrix_Band_Mutable< "<<TMV_Text(m.mat())<<" >";
        return s.str();
    }

} // namespace tmv

#endif
