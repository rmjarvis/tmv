

//---------------------------------------------------------------------------
//
// This file defines the TMV SmallBandMatrix and ThinBandMatrix classes.
//
// Constructors:
//   
//    SmallBandMatrix<T,M,N,LO,HI,A>()
//        Makes a banded matrix with column size = M and row size = N
//        with HI non-zero superdiagonals and LO non-zero subdiagonals
//        with _uninitialized_ values.  
//   
//    ThinBandMatrix<T,LO,HI,A>(size_t colsize, size_t rowsize)
//        Makes a banded matrix with column size = colsize and 
//        row size = rowsize, with HI non-zero superdiagonals and 
//        LO non-zero subdiagonals with _uninitialized_ values.
//   
//    SmallBandMatrix<T,M,N,LO,HI,A>(T x)
//    ThinBandMatrix<T,LO,HI,A>(size_t colsize, size_t rowsize, T x)
//        The same as above, but all values are initialized to x.
//    
//    SmallBandMatrix<T,M,N,LO,HI,A>(const BaseMatrix<M>& m)
//    ThinBandMatrix<T,LO,HI,A>(const BaseMatrix<M><& m)
//        Makes a banded matrix which copies the corresponding elements of m.
//    
// Special Constructors
//
//    ThinBandMatrix<T,0,1,A>(const Vector& v1, const Vector& v2)
//        Makes an upper bidiagonal matrix with the vectors
//        v1 on the main diagonal, and v2 on the superdiagonal
//
//    ThinBandMatrix<T,1,0,A>(const Vector& v1, const Vector& v2)
//        Makes an lower bidiagonal matrix with the vectors
//        v1 on the subdiagonal, and v2 on the main diagonal
//
//    ThinBandMatrix<T,1,1,A>(const Vector& v1, const Vector& v2,
//                            const Vector& v3)
//        Makes a tridiagonal matrix with the vectors
//        v1 on the subdiagonal, v2 on the main diagonal, and 
//        v3 on the superdiagonal
//


#ifndef TMV_SmallBandMatrix_H
#define TMV_SmallBandMatrix_H

#include "TMV_BaseMatrix_Band.h"
#include "TMV_BaseMatrix_Rec.h"
#include "TMV_BaseMatrix_Tri.h"
#include "TMV_Divider.h"
#include "TMV_Array.h"
#include "TMV_VIt.h"

namespace tmv {

    template <class T, int M, int N, int LO, int HI, int A0, int A1>
    struct Traits<SmallBandMatrix<T,M,N,LO,HI,A0,A1> >
    {
        typedef typename Traits<T>::real_type real_type;

        enum { A01 = A0 | A1 };
        enum { A = (A01 & ~NoDivider) | (
                ((Attrib<A01>::rowmajor || Attrib<A01>::diagmajor) ? 
                 0 : ColMajor) |
                ( Attrib<A01>::withdivider ? 0 : NoDivider ) ) };
        enum { okA = (
                !Attrib<A>::conj &&
                (Attrib<A>::rowmajor || Attrib<A>::colmajor || 
                 Attrib<A>::diagmajor) &&
                !(Attrib<A>::rowmajor && Attrib<A>::colmajor) &&
                !(Attrib<A>::rowmajor && Attrib<A>::diagmajor) &&
                !(Attrib<A>::colmajor && Attrib<A>::diagmajor) &&
                !Attrib<A>::nonunitdiag &&
                !Attrib<A>::unitdiag &&
                !Attrib<A>::zerodiag &&
                !Attrib<A>::packed &&
                !Attrib<A>::lower &&
                !Attrib<A>::upper &&
                (Attrib<A>::nodivider || Attrib<A>::withdivider) &&
                (Attrib<A>::nodivider != int(Attrib<A>::withdivider) ) &&
                ( !Attrib<A>::withdivider ||
                  ( Traits<real_type>::isinst && 
                    !Traits<real_type>::isinteger) ) )};
        enum { _attrib = A };

        typedef T value_type;
        typedef typename Traits<T>::complex_type complex_type;
        enum { isreal = Traits<T>::isreal };
        enum { iscomplex = Traits<T>::iscomplex };

        typedef SmallBandMatrix<T,M,N,LO,HI,A0,A1> type;
        typedef const type& calc_type;
        typedef const type& eval_type;
        typedef SmallBandMatrix<T,M,N,LO,HI,A01> copy_type;

        enum { _colsize = M };
        enum { _rowsize = N };
        enum { _nlo = LO };
        enum { _nhi = HI };
        enum { _shape = Band };
        enum { _fort = Attrib<A>::fort };
        enum { _calc = true };
        enum { _rowmajor = Attrib<A>::rowmajor };
        enum { _colmajor = Attrib<A>::colmajor };
        enum { _diagmajor = Attrib<A>::diagmajor };
        enum { _stepi = _colmajor ? 1 : _rowmajor ? LO+HI : N>=M ? 1-M : -N };
        enum { _stepj = _rowmajor ? 1 : _colmajor ? LO+HI : N>=M ? M : 1+N };
        enum { _diagstep = _diagmajor ? 1 : LO+HI+1 };
        enum { _conj = false };
        enum { _checkalias = Attrib<A>::checkalias };
        enum { _canlin = true };

        enum { Mm1 = IntTraits<M>::Sm1 };
        enum { Nm1 = IntTraits<N>::Sm1 };
        enum { LOpHI = IntTraits2<LO,HI>::sum };
        enum { Mx = IntTraits2<M,N+LO>::min };
        enum { Nx = IntTraits2<N,M+HI>::min };
        enum { _linsize = (
                (M == 0 || N == 0) ? 0 :
                (M == N) ? IntTraits2<M,IntTraits2<Mm1,LOpHI>::prod>::sum :
                _rowmajor ? IntTraits2<N,IntTraits2<Mm1,LOpHI>::prod>::sum :
                _colmajor ? IntTraits2<M,IntTraits2<Nm1,LOpHI>::prod>::sum :
                (M > N) ? IntTraits2<N,IntTraits<LOpHI>::Sp1>::prod :
                IntTraits2<N,IntTraits2<Mm1,LOpHI>::prod>::sum ) };

        enum { twoSi = isreal ? int(_stepi) : IntTraits<_stepi>::twoS };
        enum { twoSj = isreal ? int(_stepj) : IntTraits<_stepj>::twoS };
        enum { minMN = IntTraits2<M,N>::min };

        enum { _hasdivider = Attrib<A>::withdivider };
#if 0
        typedef QuotXB<1,real_type,type> inverse_type;
        typedef typename TypeSelect<_hasdivider ,
                const InstBandLUD<T>& , LUD<copy_type> >::type lud_type;
        typedef typename TypeSelect<_hasdivider ,
                const InstBandQRD<T>& , QRD<copy_type> >::type qrd_type;
        typedef typename TypeSelect<_hasdivider ,
                const InstBandSVD<T>& , SVD<copy_type> >::type svd_type;
#else
        typedef InvalidType inverse_type;
        typedef InvalidType lud_type;
        typedef InvalidType qrd_type;
        typedef InvalidType svd_type;
#endif

        enum { vecA = (
                (_fort ? FortranStyle : 0) |
                (_checkalias ? CheckAlias : 0) )};
        enum { colA = vecA | (_colmajor ? Unit : 0) };
        enum { rowA = vecA | (_rowmajor ? Unit : 0) };
        enum { diagA = vecA | (_diagmajor ? Unit : 0) };
        enum { ndA = (A & ~AllDivStatus) };
        enum { cstyleA = ndA & ~FortranStyle };
        enum { fstyleA = ndA | FortranStyle };
        enum { nmA = (ndA & ~AllStorageType) };
        enum { cmA = nmA | ColMajor };
        enum { rmA = nmA | RowMajor };
        enum { dmA = nmA | DiagMajor };
        enum { conjA = iscomplex ? (ndA ^ Conj) : int(ndA) };
        enum { trA = _colmajor ? rmA : _rowmajor ? cmA : int(ndA) };
        enum { adjA = iscomplex ? (trA ^ Conj) : int(trA) };
        enum { nonconjA = ndA };
        enum { twosA = isreal ? int(nonconjA) : (nonconjA & ~AllStorageType) };
        enum { An = _checkalias ? (ndA & ~CheckAlias) : (ndA | NoAlias) };
        enum { vecAn = _checkalias ? (vecA & ~CheckAlias) : (vecA | NoAlias) };
        enum { colAn = vecAn | (_colmajor ? Unit : 0) };
        enum { rowAn = vecAn | (_rowmajor ? Unit : 0) };
        enum { diagAn = vecAn | (_diagstep == 1 ? Unit : 0) };
        enum { nmAn = (An & ~AllStorageType) };

        enum { xx = TMV_UNKNOWN }; // For brevity.
        typedef ConstVectorView<T,colAn> const_col_sub_type;
        typedef ConstVectorView<T,rowAn> const_row_sub_type;
        typedef ConstSmallVectorView<T,minMN,_diagstep,diagA> const_diag_type;
        typedef ConstSmallVectorView<T,xx,_diagstep,diagA> const_diag_sub_type;

        typedef ConstMatrixView<T,An> const_submatrix_type;
        typedef ConstMatrixView<T,nmAn> const_submatrix_step_type;
        typedef ConstVectorView<T,vecAn> const_subvector_type;
        typedef ConstBandMatrixView<T,An> const_subbandmatrix_type;
        typedef ConstBandMatrixView<T,nmAn> const_subbandmatrix_step_type;
        typedef ConstBandMatrixView<T,An> const_colrange_type;
        typedef ConstBandMatrixView<T,An> const_rowrange_type;
        typedef ConstBandMatrixView<T,An> const_diagrange_type;

        typedef ConstSmallBandMatrixView<T,M,N,LO,HI,_stepi,_stepj,ndA> 
            const_view_type;
        typedef ConstSmallBandMatrixView<T,M,N,LO,HI,_stepi,_stepj,cstyleA> 
            const_cview_type;
        typedef ConstSmallBandMatrixView<T,M,N,LO,HI,_stepi,_stepj,fstyleA> 
            const_fview_type;
        typedef ConstBandMatrixView<T> const_xview_type;
        typedef ConstSmallBandMatrixView<T,M,N,LO,HI,1,_stepj,cmA> 
            const_cmview_type;
        typedef ConstSmallBandMatrixView<T,M,N,LO,HI,_stepi,1,rmA> 
            const_rmview_type;
        typedef ConstSmallBandMatrixView<T,M,N,LO,HI,_stepi,_stepj,dmA> 
            const_dmview_type;
        typedef ConstSmallBandMatrixView<T,M,N,LO,HI,_stepi,_stepj,conjA> 
            const_conjugate_type;
        typedef ConstSmallBandMatrixView<T,N,M,HI,LO,_stepj,_stepi,trA> 
            const_transpose_type;
        typedef ConstSmallBandMatrixView<T,N,M,HI,LO,_stepj,_stepi,adjA> 
            const_adjoint_type;

        typedef ConstSmallBandMatrixView<T,N,N,0,HI,_stepi,_stepj,ndA> 
            const_upperband_type;
        typedef ConstSmallBandMatrixView<T,M,M,LO,0,_stepi,_stepj,ndA> 
            const_lowerband_type;
        typedef ConstSmallBandMatrixView<T,xx,xx,0,xx,_stepi,_stepj,ndA> 
            const_upperbandoff_type;
        typedef ConstSmallBandMatrixView<T,xx,xx,xx,0,_stepi,_stepj,ndA> 
            const_lowerbandoff_type;

        typedef ConstSmallBandMatrixView<real_type,M,N,LO,HI,twoSi,twoSj,twosA> 
            const_realpart_type;
        typedef const_realpart_type const_imagpart_type;
        typedef ConstSmallVectorView<T,_linsize,1,(vecAn|Unit)> 
            const_linearview_type;
        typedef ConstSmallBandMatrixView<T,M,N,LO,HI,_stepi,_stepj,nonconjA> 
            const_nonconj_type;
        typedef SmallBandMatrixView<T,M,N,LO,HI,_stepi,_stepj,ndA> 
            nonconst_type;

        typedef T& reference;
        typedef VIt<T,1,false> linear_iterator;

        typedef VectorView<T,colAn> col_sub_type;
        typedef VectorView<T,rowAn> row_sub_type;
        typedef SmallVectorView<T,minMN,_diagstep,diagA> diag_type;
        typedef SmallVectorView<T,xx,_diagstep,diagA> diag_sub_type;

        typedef MatrixView<T,An> submatrix_type;
        typedef MatrixView<T,nmAn> submatrix_step_type;
        typedef VectorView<T,vecAn> subvector_type;
        typedef BandMatrixView<T,An> subbandmatrix_type;
        typedef BandMatrixView<T,nmAn> subbandmatrix_step_type;
        typedef BandMatrixView<T,An> colrange_type;
        typedef BandMatrixView<T,An> rowrange_type;
        typedef BandMatrixView<T,An> diagrange_type;

        typedef SmallBandMatrixView<T,M,N,LO,HI,_stepi,_stepj,ndA> view_type;
        typedef SmallBandMatrixView<T,M,N,LO,HI,_stepi,_stepj,cstyleA> 
            cview_type;
        typedef SmallBandMatrixView<T,M,N,LO,HI,_stepi,_stepj,fstyleA> 
            fview_type;
        typedef BandMatrixView<T> xview_type;
        typedef SmallBandMatrixView<T,M,N,LO,HI,1,_stepj,cmA> cmview_type;
        typedef SmallBandMatrixView<T,M,N,LO,HI,_stepi,1,rmA> rmview_type;
        typedef SmallBandMatrixView<T,M,N,LO,HI,_stepi,_stepj,dmA> dmview_type;
        typedef SmallBandMatrixView<T,M,N,LO,HI,_stepi,_stepj,conjA> 
            conjugate_type;
        typedef SmallBandMatrixView<T,N,M,HI,LO,_stepj,_stepi,trA> 
            transpose_type;
        typedef SmallBandMatrixView<T,N,M,HI,LO,_stepj,_stepi,adjA> 
            adjoint_type;

        typedef SmallBandMatrixView<T,N,N,0,HI,_stepi,_stepj,ndA> 
            upperband_type;
        typedef SmallBandMatrixView<T,M,M,LO,0,_stepi,_stepj,ndA> 
            lowerband_type;
        typedef SmallBandMatrixView<T,xx,xx,0,xx,_stepi,_stepj,ndA> 
            upperbandoff_type;
        typedef SmallBandMatrixView<T,xx,xx,xx,0,_stepi,_stepj,ndA> 
            lowerbandoff_type;

        typedef SmallBandMatrixView<real_type,M,N,LO,HI,twoSi,twoSj,twosA> 
            realpart_type;
        typedef realpart_type imagpart_type;
        typedef SmallVectorView<T,_linsize,1,(vecAn|Unit)> linearview_type;
        typedef SmallBandMatrixView<T,M,N,LO,HI,_stepi,_stepj,nonconjA> 
            nonconj_type;
    };

    template <class T, int M, int N, int LO, int HI, int A0, int A1>
    class SmallBandMatrix : 
        public BaseMatrix_Band_Mutable<SmallBandMatrix<T,M,N,LO,HI,A0,A1> >,
        public BandMatrixDivHelper<SmallBandMatrix<T,M,N,LO,HI,A0,A1> >
    {
    public:

        typedef SmallBandMatrix<T,M,N,LO,HI,A0,A1> type;
        typedef BaseMatrix_Band_Mutable<type> base_mut;

        typedef typename Traits<T>::real_type real_type;

        enum { _colsize = Traits<type>::_colsize };
        enum { _rowsize = Traits<type>::_rowsize };
        enum { _nlo = Traits<type>::_nlo };
        enum { _nhi = Traits<type>::_nhi };
        enum { _shape = Traits<type>::_shape };
        enum { _fort = Traits<type>::_fort };
        enum { _calc = Traits<type>::_calc };
        enum { _rowmajor = Traits<type>::_rowmajor };
        enum { _colmajor = Traits<type>::_colmajor };
        enum { _diagmajor = Traits<type>::_diagmajor };
        enum { _stepi = Traits<type>::_stepi };
        enum { _stepj = Traits<type>::_stepj };
        enum { _diagstep = Traits<type>::_diagstep };
        enum { _conj = Traits<type>::_conj };
        enum { _checkalias = Traits<type>::_checkalias };
        enum { _canlin = Traits<type>::_canlin };
        enum { _attrib = Traits<type>::_attrib };
        enum { _linsize = Traits<type>::_linsize };
        enum { _stor = Attrib<_attrib>::stor };

        typedef typename Traits<type>::linearview_type linearview_type;
        typedef typename Traits<type>::const_linearview_type 
            const_linearview_type;
        typedef typename Traits<type>::linear_iterator linear_iterator;

        //
        // Constructors
        //

        SmallBandMatrix() :
            itsm(itsm1.get() - _diagmajor ? LO*_stepi : 0)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(M>=0);
            TMVStaticAssert(N>=0);
            TMVStaticAssert(LO >= 0 && LO < M);
            TMVStaticAssert(HI >= 0 && HI < N);
#ifdef TMV_DEBUG
            this->setAllTo(Traits<T>::constr_value());
#endif
        }

        explicit SmallBandMatrix(T x) :
            itsm(itsm1.get() - _diagmajor ? LO*_stepi : 0)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(M>=0);
            TMVStaticAssert(N>=0);
            TMVStaticAssert(LO >= 0 && LO < M);
            TMVStaticAssert(HI >= 0 && HI < N);
            this->setAllTo(x);
        }

        SmallBandMatrix(const type& m2) :
            itsm(itsm1.get() - _diagmajor ? LO*_stepi : 0)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(M>=0);
            TMVStaticAssert(N>=0);
            TMVStaticAssert(LO >= 0 && LO < M);
            TMVStaticAssert(HI >= 0 && HI < N);
            m2.newAssignTo(*this);
        }

        template <class M2>
        SmallBandMatrix(const BaseMatrix<M2>& m2) :
            itsm(itsm1.get() - _diagmajor ? LO*_stepi : 0)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(M>=0);
            TMVStaticAssert(N>=0);
            TMVStaticAssert(LO >= 0 && LO < M);
            TMVStaticAssert(HI >= 0 && HI < N);
            TMVStaticAssert((ShapeTraits2<M2::_shape,_shape>::assignable));
            TMVAssert(m2.colsize() == M);
            TMVAssert(m2.rowsize() == N);
            m2.newAssignTo(*this);
        }

        template <class M2>
        SmallBandMatrix(const BaseMatrix_Band<M2>& m2) :
            itsm(itsm1.get() - _diagmajor ? LO*_stepi : 0)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(M>=0);
            TMVStaticAssert(N>=0);
            TMVStaticAssert(LO >= 0 && LO < M);
            TMVStaticAssert(HI >= 0 && HI < N);
            TMVAssert(LO <= m2.nlo());
            TMVAssert(HI <= m2.nhi());
            TMVAssert(m2.colsize() == M);
            TMVAssert(m2.rowsize() == N);
            m2.cDiagRange(-LO,HI+1).newAssignTo(*this);
        }

        template <class M2>
        SmallBandMatrix(const BaseMatrix_Rec<M2>& m2) :
            itsm(itsm1.get() - _diagmajor ? LO*_stepi : 0)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(M>=0);
            TMVStaticAssert(N>=0);
            TMVStaticAssert(LO >= 0 && LO < M);
            TMVStaticAssert(HI >= 0 && HI < N);
            TMVAssert(m2.colsize() == M);
            TMVAssert(m2.rowsize() == N);
            typedef typename M2::value_type T2;
            ConstBandMatrixView<T2> m2b(
                m2.cptr(),M,N,LO,HI,m2.stepi(),m2.stepj());
            m2b.newAssignTo(*this);
        }

        template <class M2>
        SmallBandMatrix(const BaseMatrix_Tri<M2>& m2) :
            itsm(itsm1.get() - _diagmajor ? LO*_stepi : 0)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(M>=0);
            TMVStaticAssert(N>=0);
            TMVStaticAssert(LO >= 0 && LO < M);
            TMVStaticAssert(HI >= 0 && HI < N);
            TMVStaticAssert((M2::_upper ? LO : HI) == 0);
            TMVStaticAssert(M == N);
            TMVAssert(m2.size() == M);
            typedef typename M2::value_type T2;
            ConstBandMatrixView<T2> m2b(
                m2.cptr(),M,N,LO,HI,m2.stepi(),m2.stepj());
            m2b.newAssignTo(*this);
        }

        template <class M2>
        SmallBandMatrix(const BaseMatrix_Diag<M2>& m2) :
            itsm(itsm1.get() - _diagmajor ? LO*_stepi : 0)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(M>=0);
            TMVStaticAssert(N>=0);
            TMVStaticAssert(LO == 0);
            TMVStaticAssert(HI == 0);
            TMVStaticAssert(M == N);
            TMVAssert(m2.size() == M);
            m2.diag().newAssignTo(this->diag());
        }

        ~SmallBandMatrix() 
        {
#ifdef TMV_DEBUG
            this->setAllTo(Traits<T>::destr_value());
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

        template <class M2>
        TMV_INLINE type& operator=(const BaseMatrix_Band<M2>& m2)
        { base_mut::operator=(m2); return *this; }

        template <class M2>
        TMV_INLINE type& operator=(const BaseMatrix_Tri<M2>& m2)
        { base_mut::operator=(m2); return *this; }

        template <class M2>
        TMV_INLINE type& operator=(const BaseMatrix_Diag<M2>& m2)
        { base_mut::operator=(m2); return *this; }

        TMV_INLINE type& operator=(T x)
        { base_mut::operator=(x); return *this; }


        //
        // LinearView is only allowed for an actual BandMatrix
        // Views cannot be linearized.
        //
        
        TMV_INLINE_ND linearview_type linearView() 
        {
#ifdef TMV_USE_VALGRIND
            AlignedArray<int> ok(ls());
            for (size_t i=0;i<ls();++i) ok[i] = 0;
            for (size_t i=0;i<colsize();++i) for (size_t j=0;j<rowsize();++j) {
                if (okij(i,j)) {
                    int k = i*stepi() + j*stepj();
                    ok[k] = true;
                }
            }
            for (size_t k=0;k<ls();++k) if (!ok[k]) {
                ptr()[k] = T(-777);
            }
#endif
            return linearview_type(start_mem());
        }

        TMV_INLINE_ND const_linearview_type linearView() const
        {
#ifdef TMV_USE_VALGRIND
            AlignedArray<int> ok(ls());
            for (size_t i=0;i<ls();++i) ok[i] = 0;
            for (size_t i=0;i<colsize();++i) for (size_t j=0;j<rowsize();++j) {
                if (okij(i,j)) {
                    int k = i*stepi() + j*stepj();
                    ok[k] = true;
                }
            }
            for (size_t k=0;k<ls();++k) if (!ok[k]) {
                ptr()[k] = T(-777);
            }
#endif
            return const_linearview_type(start_mem(),ls(),1);
        }

        //
        // Auxilliary Functions
        //

        TMV_INLINE const T* cptr() const { return itsm; }
        TMV_INLINE T* ptr() { return itsm; }
        TMV_INLINE const T* start_mem() const { return itsm1; }
        TMV_INLINE T* start_mem() { return itsm1; }

        T cref(int i, int j) const
        { return itsm[i*stepi() + j*stepj()]; }

        T& ref(int i, int j)
        { return itsm[i*stepi() + j*stepj()]; }

        TMV_INLINE size_t ls() const { return _linsize; }
        TMV_INLINE size_t colsize() const { return M; }
        TMV_INLINE size_t rowsize() const { return N; }
        TMV_INLINE int nlo() const { return LO; }
        TMV_INLINE int nhi() const { return HI; }
        TMV_INLINE int nElements() const { return _linsize; }
        TMV_INLINE int stepi() const { return _stepi; }
        TMV_INLINE int stepj() const { return _stepj; }
        TMV_INLINE int diagstep() const { return _diagstep; }
        TMV_INLINE bool isconj() const { return false; }
        TMV_INLINE bool isrm() const { return _rowmajor; }
        TMV_INLINE bool iscm() const { return _colmajor; }
        TMV_INLINE bool isdm() const { return _diagmajor; }

    private:

        StackArray<T,_linsize> itsm1;
        T* itsm;

    }; // SmallBandMatrix

    template <class T, int LO, int HI, int A0, int A1>
    struct Traits<ThinBandMatrix<T,LO,HI,A0,A1> >
    {
        typedef typename Traits<T>::real_type real_type;

        enum { A01 = A0 | A1 };
        enum { A = A01 | (
                ((Attrib<A01>::rowmajor || Attrib<A01>::diagmajor) ? 
                 0 : ColMajor) |
                ( !Traits<real_type>::isinst ? NoDivider :
                  Traits<real_type>::isinteger ? NoDivider :
                  Attrib<A01>::nodivider ? 0 : WithDivider ) )};
        enum { okA = (
                !Attrib<A>::conj &&
                (Attrib<A>::rowmajor || Attrib<A>::colmajor || 
                 Attrib<A>::diagmajor) &&
                !(Attrib<A>::rowmajor && Attrib<A>::colmajor) &&
                !(Attrib<A>::rowmajor && Attrib<A>::diagmajor) &&
                !(Attrib<A>::colmajor && Attrib<A>::diagmajor) &&
                !Attrib<A>::nonunitdiag &&
                !Attrib<A>::unitdiag &&
                !Attrib<A>::zerodiag &&
                !Attrib<A>::packed &&
                !Attrib<A>::lower &&
                !Attrib<A>::upper &&
                (Attrib<A>::nodivider || Attrib<A>::withdivider) &&
                (Attrib<A>::nodivider != int(Attrib<A>::withdivider) ) &&
                ( !Attrib<A>::withdivider ||
                  ( Traits<real_type>::isinst && 
                    !Traits<real_type>::isinteger) ) )};
        enum { _attrib = A };

        typedef T value_type;
        typedef typename Traits<T>::complex_type complex_type;
        enum { isreal = Traits<T>::isreal };
        enum { iscomplex = Traits<T>::iscomplex };

        typedef ThinBandMatrix<T,LO,HI,A0,A1> type;
        typedef const type& calc_type;
        typedef const type& eval_type;
        typedef ThinBandMatrix<T,LO,HI,A01> copy_type;

        enum { _colsize = TMV_UNKNOWN };
        enum { _rowsize = TMV_UNKNOWN };
        enum { _nlo = LO };
        enum { _nhi = HI };
        enum { _shape = Band };
        enum { _fort = Attrib<A>::fort };
        enum { _calc = true };
        enum { _rowmajor = Attrib<A>::rowmajor };
        enum { _colmajor = Attrib<A>::colmajor };
        enum { _diagmajor = Attrib<A>::diagmajor };
        enum { _stepi = _colmajor ? 1 : _rowmajor ? LO+HI : TMV_UNKNOWN };
        enum { _stepj = _rowmajor ? 1 : _colmajor ? LO+HI : TMV_UNKNOWN };
        enum { _diagstep = _diagmajor ? 1 : LO+HI+1 };
        enum { _conj = false };
        enum { _checkalias = !Attrib<A>::noalias };
        enum { _canlin = true };

        enum { twoSi = isreal ? int(_stepi) : IntTraits<_stepi>::twoS };
        enum { twoSj = isreal ? int(_stepj) : IntTraits<_stepj>::twoS };

        enum { _hasdivider = Attrib<A>::withdivider };

#if 0
        typedef QuotXB<1,real_type,type> inverse_type;
        typedef typename TypeSelect<_hasdivider ,
                const InstBandLUD<T>& , LUD<copy_type> >::type lud_type;
        typedef typename TypeSelect<_hasdivider ,
                const InstBandQRD<T>& , QRD<copy_type> >::type qrd_type;
        typedef typename TypeSelect<_hasdivider ,
                const InstBandSVD<T>& , SVD<copy_type> >::type svd_type;
#else
        typedef InvalidType inverse_type;
        typedef InvalidType lud_type;
        typedef InvalidType qrd_type;
        typedef InvalidType svd_type;
#endif

        enum { vecA = (
                (_fort ? FortranStyle : 0) |
                (_checkalias ? CheckAlias : 0) )};
        enum { colA = vecA | (_colmajor ? Unit : 0) };
        enum { rowA = vecA | (_rowmajor ? Unit : 0) };
        enum { diagA = vecA | (_diagmajor ? Unit : 0) };
        enum { ndA = (A & ~AllDivStatus) };
        enum { cstyleA = ndA & ~FortranStyle };
        enum { fstyleA = ndA | FortranStyle };
        enum { nmA = (ndA & ~AllStorageType) };
        enum { cmA = nmA | ColMajor };
        enum { rmA = nmA | RowMajor };
        enum { dmA = nmA | DiagMajor };
        enum { conjA = iscomplex ? (ndA ^ Conj) : int(ndA) };
        enum { trA = _colmajor ? rmA : _rowmajor ? cmA : int(ndA) };
        enum { adjA = iscomplex ? (trA ^ Conj) : int(trA) };
        enum { nonconjA = ndA };
        enum { twosA = isreal ? int(nonconjA) : (nonconjA & ~AllStorageType) };
        enum { An = _checkalias ? (ndA & ~CheckAlias) : (ndA | NoAlias) };
        enum { vecAn = _checkalias ? (vecA & ~CheckAlias) : (vecA | NoAlias) };
        enum { colAn = vecAn | (_colmajor ? Unit : 0) };
        enum { rowAn = vecAn | (_rowmajor ? Unit : 0) };
        enum { diagAn = vecAn | (_diagstep == 1 ? Unit : 0) };
        enum { nmAn = (An & ~AllStorageType) };

        enum { xx = TMV_UNKNOWN }; // For brevity.
        typedef ConstVectorView<T,colAn> const_col_sub_type;
        typedef ConstVectorView<T,rowAn> const_row_sub_type;
        typedef ConstSmallVectorView<T,xx,_diagstep,diagA> const_diag_type;
        typedef ConstSmallVectorView<T,xx,_diagstep,diagA> const_diag_sub_type;

        typedef ConstMatrixView<T,An> const_submatrix_type;
        typedef ConstMatrixView<T,nmAn> const_submatrix_step_type;
        typedef ConstVectorView<T,vecAn> const_subvector_type;
        typedef ConstBandMatrixView<T,An> const_subbandmatrix_type;
        typedef ConstBandMatrixView<T,nmAn> const_subbandmatrix_step_type;
        typedef ConstBandMatrixView<T,An> const_colrange_type;
        typedef ConstBandMatrixView<T,An> const_rowrange_type;
        typedef ConstBandMatrixView<T,An> const_diagrange_type;

        typedef ConstSmallBandMatrixView<T,xx,xx,LO,HI,_stepi,_stepj,ndA> 
            const_view_type;
        typedef ConstSmallBandMatrixView<T,xx,xx,LO,HI,_stepi,_stepj,cstyleA> 
            const_cview_type;
        typedef ConstSmallBandMatrixView<T,xx,xx,LO,HI,_stepi,_stepj,fstyleA> 
            const_fview_type;
        typedef ConstBandMatrixView<T> const_xview_type;
        typedef ConstSmallBandMatrixView<T,xx,xx,LO,HI,1,_stepj,cmA> 
            const_cmview_type;
        typedef ConstSmallBandMatrixView<T,xx,xx,LO,HI,_stepi,1,rmA> 
            const_rmview_type;
        typedef ConstSmallBandMatrixView<T,xx,xx,LO,HI,_stepi,_stepj,dmA> 
            const_dmview_type;
        typedef ConstSmallBandMatrixView<T,xx,xx,LO,HI,_stepi,_stepj,conjA> 
            const_conjugate_type;
        typedef ConstSmallBandMatrixView<T,xx,xx,HI,LO,_stepj,_stepi,trA> 
            const_transpose_type;
        typedef ConstSmallBandMatrixView<T,xx,xx,HI,LO,_stepj,_stepi,adjA> 
            const_adjoint_type;

        typedef ConstSmallBandMatrixView<T,xx,xx,0,HI,_stepi,_stepj,ndA> 
            const_upperband_type;
        typedef ConstSmallBandMatrixView<T,xx,xx,LO,0,_stepi,_stepj,ndA> 
            const_lowerband_type;
        typedef ConstSmallBandMatrixView<T,xx,xx,0,xx,_stepi,_stepj,ndA> 
            const_upperbandoff_type;
        typedef ConstSmallBandMatrixView<T,xx,xx,xx,0,_stepi,_stepj,ndA> 
            const_lowerbandoff_type;

        typedef ConstSmallBandMatrixView<real_type,xx,xx,LO,HI,twoSi,twoSj,twosA> 
            const_realpart_type;
        typedef const_realpart_type const_imagpart_type;
        typedef ConstVectorView<T,(vecAn|Unit)> const_linearview_type;
        typedef ConstSmallBandMatrixView<T,xx,xx,LO,HI,_stepi,_stepj,nonconjA> 
            const_nonconj_type;
        typedef SmallBandMatrixView<T,xx,xx,LO,HI,_stepi,_stepj,ndA> 
            nonconst_type;

        typedef T& reference;
        typedef VIt<T,1,false> linear_iterator;

        typedef VectorView<T,colA> col_sub_type;
        typedef VectorView<T,rowA> row_sub_type;
        typedef SmallVectorView<T,xx,_diagstep,diagA> diag_type;
        typedef SmallVectorView<T,xx,_diagstep,diagA> diag_sub_type;

        typedef MatrixView<T,ndA> submatrix_type;
        typedef MatrixView<T,nmA> submatrix_step_type;
        typedef VectorView<T,vecA> subvector_type;
        typedef BandMatrixView<T,ndA> subbandmatrix_type;
        typedef BandMatrixView<T,nmA> subbandmatrix_step_type;
        typedef BandMatrixView<T,ndA> colrange_type;
        typedef BandMatrixView<T,ndA> rowrange_type;
        typedef BandMatrixView<T,ndA> diagrange_type;

        typedef SmallBandMatrixView<T,xx,xx,LO,HI,_stepi,_stepj,ndA> 
            view_type;
        typedef SmallBandMatrixView<T,xx,xx,LO,HI,_stepi,_stepj,cstyleA> 
            cview_type;
        typedef SmallBandMatrixView<T,xx,xx,LO,HI,_stepi,_stepj,fstyleA> 
            fview_type;
        typedef BandMatrixView<T> xview_type;
        typedef SmallBandMatrixView<T,xx,xx,LO,HI,1,_stepj,cmA> 
            cmview_type;
        typedef SmallBandMatrixView<T,xx,xx,LO,HI,_stepi,1,rmA> 
            rmview_type;
        typedef SmallBandMatrixView<T,xx,xx,LO,HI,_stepi,_stepj,dmA> 
            dmview_type;
        typedef SmallBandMatrixView<T,xx,xx,LO,HI,_stepi,_stepj,conjA> 
            conjugate_type;
        typedef SmallBandMatrixView<T,xx,xx,HI,LO,_stepj,_stepi,trA> 
            transpose_type;
        typedef SmallBandMatrixView<T,xx,xx,HI,LO,_stepj,_stepi,adjA> 
            adjoint_type;

        typedef SmallBandMatrixView<T,xx,xx,0,HI,_stepi,_stepj,ndA> 
            upperband_type;
        typedef SmallBandMatrixView<T,xx,xx,LO,0,_stepi,_stepj,ndA> 
            lowerband_type;
        typedef SmallBandMatrixView<T,xx,xx,0,xx,_stepi,_stepj,ndA> 
            upperbandoff_type;
        typedef SmallBandMatrixView<T,xx,xx,xx,0,_stepi,_stepj,ndA> 
            lowerbandoff_type;

        typedef SmallBandMatrixView<real_type,xx,xx,LO,HI,twoSi,twoSj,twosA> 
            realpart_type;
        typedef realpart_type imagpart_type;
        typedef SmallBandMatrixView<T,xx,xx,LO,HI,_stepi,_stepj,nonconjA> 
            nonconj_type;
    };


    template <class T, int LO, int HI, int A0, int A1>
    class ThinBandMatrix : 
        public BaseMatrix_Band_Mutable<ThinBandMatrix<T,LO,HI,A0,A1> >,
        public BandMatrixDivHelper<ThinBandMatrix<T,LO,HI,A0,A1> >
    {
    public:

        typedef ThinBandMatrix<T,LO,HI<A0,A1> type;
        typedef BaseMatrix_Band_Mutable<type> base_mut;
        typedef BandMatrixDivHelper<type> divhelper;

        typedef typename Traits<T>::real_type real_type;

        enum { _colsize = Traits<type>::_colsize };
        enum { _rowsize = Traits<type>::_rowsize };
        enum { _nlo = Traits<type>::_nlo };
        enum { _nhi = Traits<type>::_nhi };
        enum { _shape = Traits<type>::_shape };
        enum { _fort = Traits<type>::_fort };
        enum { _calc = Traits<type>::_calc };
        enum { _rowmajor = Traits<type>::_rowmajor };
        enum { _colmajor = Traits<type>::_colmajor };
        enum { _diagmajor = Traits<type>::_diagmajor };
        enum { _stepi = Traits<type>::_stepi };
        enum { _stepj = Traits<type>::_stepj };
        enum { _diagstep = Traits<type>::_diagstep };
        enum { _conj = Traits<type>::_conj };
        enum { _checkalias = Traits<type>::_checkalias };
        enum { _canlin = Traits<type>::_canlin };
        enum { _attrib = Traits<type>::_attrib };
        enum { _stor = Attrib<_attrib>::stor };

        typedef typename Traits<type>::linearview_type linearview_type;
        typedef typename Traits<type>::const_linearview_type 
            const_linearview_type;
        typedef typename Traits<type>::linear_iterator linear_iterator;

        //
        // Constructors
        //

        ThinBandMatrix() :
            itscs(0), itsrs(0), linsize(0), itssi(0), itssj(0), 
            itsm1(0), itsm(0)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(LO >= 0);
            TMVStaticAssert(HI >= 0);
#ifdef TMV_DEBUG
            this->setAllTo(Traits<T>::constr_value());
#endif
        }

        ThinBandMatrix(size_t cs, size_t rs) :
            itscs(cs), itsrs(rs), 
            linsize(BandStorageLength(_stor,itscs,itsrs,LO,HI)),
            itssi(_rowmajor ? LO+HI : _colmajor ? 1 :
                  rs >= cs ? 1-int(cs) : -int(rs)),
            itssj(_rowmajor ? 1 : _colmajor ? LO+HI : -itssi+1),
            itsm1(linsize), itsm(itsm1.get() - _diagmajor ? LO*itssi : 0)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(LO >= 0);
            TMVStaticAssert(HI >= 0);
            TMVAssert(LO < cs);
            TMVAssert(HI < rs);
#ifdef TMV_DEBUG
            this->setAllTo(Traits<T>::constr_value());
#endif
        }

        ThinBandMatrix(size_t cs, size_t rs, T x) :
            itscs(cs), itsrs(rs), 
            linsize(BandStorageLength(_stor,itscs,itsrs,LO,HI)),
            itssi(_rowmajor ? LO+HI : _colmajor ? 1 :
                  rs >= cs ? 1-int(cs) : -int(rs)),
            itssj(_rowmajor ? 1 : _colmajor ? LO+HI : -itssi+1),
            itsm1(linsize), itsm(itsm1.get() - _diagmajor ? LO*itssi : 0)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(LO >= 0);
            TMVStaticAssert(HI >= 0);
            TMVAssert(LO < cs);
            TMVAssert(HI < rs);
            this->setAllTo(x);
        }

        ThinBandMatrix(const type& m2) :
            itscs(m2.itscs), itsrs(m2.itsrs), 
            linsize(BandStorageLength(_stor,itscs,itsrs,LO,HI)),
            itssi(_rowmajor ? LO+HI : _colmajor ? 1 :
                  itsrs >= itscs ? 1-int(itscs) : -int(itsrs)),
            itssj(_rowmajor ? 1 : _colmajor ? LO+HI : -itssi+1),
            itsm1(linsize), itsm(itsm1.get() - _diagmajor ? LO*itssi : 0)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(LO >= 0);
            TMVStaticAssert(HI >= 0);
            m2.newAssignTo(*this);
        }

        template <class M2>
        ThinBandMatrix(const BaseMatrix<M2>& m2) :
            itscs(m2.colsize()), itsrs(m2.rowsize()), 
            linsize(BandStorageLength(_stor,itscs,itsrs,LO,HI)),
            itssi(_rowmajor ? LO+HI : _colmajor ? 1 :
                  itsrs >= itscs ? 1-int(itscs) : -int(itsrs)),
            itssj(_rowmajor ? 1 : _colmajor ? LO+HI : -itssi+1),
            itsm1(linsize), itsm(itsm1.get() - _diagmajor ? LO*itssi : 0)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(LO >= 0);
            TMVStaticAssert(HI >= 0);
            TMVAssert(LO < m2.colsize());
            TMVAssert(HI < m2.rowsize());
            TMVStaticAssert((ShapeTraits2<M2::_shape,_shape>::assignable));
            m2.newAssignTo(*this);
        }

        template <class M2>
        ThinBandMatrix(const BaseMatrix_Band<M2>& m2) :
            itscs(m2.colsize()), itsrs(m2.rowsize()), 
            linsize(BandStorageLength(_stor,itscs,itsrs,LO,HI)),
            itssi(_rowmajor ? LO+HI : _colmajor ? 1 :
                  itsrs >= itscs ? 1-int(itscs) : -int(itsrs)),
            itssj(_rowmajor ? 1 : _colmajor ? LO+HI : -itssi+1),
            itsm1(linsize), itsm(itsm1.get() - _diagmajor ? LO*itssi : 0)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(LO >= 0);
            TMVStaticAssert(HI >= 0);
            TMVAssert(LO <= m2.nlo());
            TMVAssert(HI <= m2.nhi());
            m2.cDiagRange(-LO,HI+1).newAssignTo(*this);
        }

        template <class M2>
        ThinBandMatrix(const BaseMatrix_Rec<M2>& m2) :
            itscs(m2.colsize()), itsrs(m2.rowsize()), 
            linsize(BandStorageLength(_stor,itscs,itsrs,LO,HI)),
            itssi(_rowmajor ? LO+HI : _colmajor ? 1 :
                  itsrs >= itscs ? 1-int(itscs) : -int(itsrs)),
            itssj(_rowmajor ? 1 : _colmajor ? LO+HI : -itssi+1),
            itsm1(linsize), itsm(itsm1.get() - _diagmajor ? LO*itssi : 0)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(LO >= 0);
            TMVStaticAssert(HI >= 0);
            TMVAssert(LO < m2.colsize());
            TMVAssert(HI < m2.rowsize());
            typedef typename M2::value_type T2;
            ConstBandMatrixView<T2> m2b(
                m2.cptr(),m2.colsize(),m2.rowsize(),LO,HI,
                m2.stepi(),m2.stepj());
            m2b.newAssignTo(*this);
        }

        template <class M2>
        ThinBandMatrix(const BaseMatrix_Tri<M2>& m2) :
            itscs(m2.size()), itsrs(m2.size()), 
            linsize(BandStorageLength(_stor,itscs,itsrs,LO,HI)),
            itssi(_rowmajor ? LO+HI : _colmajor ? 1 : 1-int(itscs)),
            itssj(_rowmajor ? 1 : _colmajor ? LO+HI : int(itscs)),
            itsm1(linsize), itsm(itsm1.get() - _diagmajor ? LO*itssi : 0)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(LO >= 0);
            TMVStaticAssert(HI >= 0);
            TMVStaticAssert((M2::_upper ? LO : HI) == 0);
            TMVAssert((M2::_upper ? HI : LO) < m2.size());
            typedef typename M2::value_type T2;
            ConstBandMatrixView<T2> m2b(
                m2.cptr(),m2.size(),m2.size(),LO,HI,
                m2.stepi(),m2.stepj());
            m2b.newAssignTo(*this);
        }

        template <class M2>
        ThinBandMatrix(const BaseMatrix_Diag<M2>& m2) :
            itscs(m2.size()), itsrs(m2.size()), 
            linsize(BandStorageLength(_stor,itscs,itsrs,LO,HI)),
            itssi(_rowmajor ? 0 : _colmajor ? 1 : 1-int(itscs)),
            itssj(_rowmajor ? 1 : _colmajor ? 0 : int(itscs)),
            itsm1(linsize), itsm(itsm1.get())
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(LO == 0);
            TMVStaticAssert(HI == 0);
            m2.diag().newAssignTo(this->diag());
        }

        template <class V0, class V1>
        ThinBandMatrix(const BaseVector<V0>& v0, const BaseVector<V1>& v1) :
            itscs(LO==0?v0.size():v0.size()+1),
            itsrs(HI==0?v1.size():v1.size()+1),
            linsize(BandStorageLength(_stor,itscs,itsrs,LO,HI)),
            itssi(_rowmajor ? 0 : _colmajor ? 1 : 1-int(itscs)),
            itssj(_rowmajor ? 1 : _colmajor ? 0 : int(itscs)),
            itsm1(linsize), itsm(itsm1.get())
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert((LO == 0 && HI == 1) || (LO == 1 && HI == 0) );
            TMVAssert( (v0.size() == v1.size()) ||
                       (v0.size() == v1.size()-1 && LO==1 && HI==0) ||
                       (v0.size() == v1.size()+1 && LO==0 && HI==1) );
            v0.newAssignTo(this->diag(LO==1 ? -1 : 0));
            v1.newAssignTo(this->diag(LO==1 ? 0 : 1));
        }

        template <class V0, class V1, class V2>
        ThinBandMatrix(
            const BaseVector<V0>& v0, const BaseVector<V1>& v1,
            const BaseVector<V2>& v2) :
            itscs(v0.size()+1), itsrs(v2.size()+1),
            linsize(BandStorageLength(_stor,itscs,itsrs,LO,HI)),
            itssi(_rowmajor ? 0 : _colmajor ? 1 : 1-int(itscs)),
            itssj(_rowmajor ? 1 : _colmajor ? 0 : int(itscs)),
            itsm1(linsize), itsm(itsm1.get())
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(LO == 1 && HI == 1);
            TMVAssert( (v0.size() == v1.size()) || (v0.size() == v1.size()-1) );
            TMVAssert( (v2.size() == v1.size()) || (v2.size() == v1.size()-1) );
            v0.newAssignTo(this->diag(-1));
            v1.newAssignTo(this->diag());
            v2.newAssignTo(this->diag(1));
        }

        ~ThinBandMatrix() 
        {
#ifdef TMV_DEBUG
            this->setAllTo(Traits<T>::destr_value());
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

        template <class M2>
        TMV_INLINE type& operator=(const BaseMatrix_Band<M2>& m2)
        { base_mut::operator=(m2); return *this; }

        template <class M2>
        TMV_INLINE type& operator=(const BaseMatrix_Tri<M2>& m2)
        { base_mut::operator=(m2); return *this; }

        template <class M2>
        TMV_INLINE type& operator=(const BaseMatrix_Diag<M2>& m2)
        { base_mut::operator=(m2); return *this; }

        TMV_INLINE type& operator=(T x)
        { base_mut::operator=(x); return *this; }


        //
        // LinearView is only allowed for an actual BandMatrix
        // Views cannot be linearized.
        //
        
        TMV_INLINE_ND linearview_type linearView() 
        {
#ifdef TMV_USE_VALGRIND
            AlignedArray<int> ok(ls());
            for (size_t i=0;i<ls();++i) ok[i] = 0;
            for (size_t i=0;i<colsize();++i) for (size_t j=0;j<rowsize();++j) {
                if (okij(i,j)) {
                    int k = i*stepi() + j*stepj();
                    ok[k] = true;
                }
            }
            for (size_t k=0;k<ls();++k) if (!ok[k]) {
                ptr()[k] = T(-777);
            }
#endif
            return linearview_type(start_mem(),ls(),1);
        }

        TMV_INLINE_ND const_linearview_type linearView() const
        {
#ifdef TMV_USE_VALGRIND
            AlignedArray<int> ok(ls());
            for (size_t i=0;i<ls();++i) ok[i] = 0;
            for (size_t i=0;i<colsize();++i) for (size_t j=0;j<rowsize();++j) {
                if (okij(i,j)) {
                    int k = i*stepi() + j*stepj();
                    ok[k] = true;
                }
            }
            for (size_t k=0;k<ls();++k) if (!ok[k]) {
                ptr()[k] = T(-777);
            }
#endif
            return const_linearview_type(start_mem(),ls(),1);
        }

        //
        // Auxilliary Functions
        //

        TMV_INLINE const T* cptr() const { return itsm; }
        TMV_INLINE T* ptr() { return itsm; }
        TMV_INLINE const T* start_mem() const { return itsm1; }
        TMV_INLINE T* start_mem() { return itsm1; }

        T cref(int i, int j) const
        { return itsm[i*stepi() + j*stepj()]; }

        T& ref(int i, int j)
        { return itsm[i*stepi() + j*stepj()]; }

        void swapWith(type& m2)
        {
            TMVAssert(m2.colsize() == colsize());
            TMVAssert(m2.rowsize() == rowsize());
            TMVAssert(m2.nlo() == nlo());
            TMVAssert(m2.nhi() == nhi());
            if (itsm1.get() == m2.itsm1.get()) return;
            itsm1.swapWith(m2.itsm1);
            TMV_SWAP(itsm,m2.itsm);
        }

        void resize(size_t cs, size_t rs)
        {
#ifdef TMV_DEBUG
            this->setAllTo(Traits<T>::destr_value());
#endif
            itscs = cs;
            itsrs = rs;
            divhelper::resetDivType();
            linsize = BandStorageLength(_stor,itscs,itsrs,LO,HI);
            itsm1.resize(linsize);
            itsm = itsm1.get() - _diagmajor ? LO*itssi : 0;
            itssi = _rowmajor ? LO+HI : _colmajor ? 1 :
                rs >= cs ? 1-int(cs) : -int(rs);
            itssj = _rowmajor ? 1 : _colmajor ? LO+HI : -itssi+1;
#ifdef TMV_DEBUG
            this->setAllTo(Traits<T>::constr_value());
#endif
        }


        TMV_INLINE size_t ls() const { return linsize; }
        TMV_INLINE size_t colsize() const { return itscs; }
        TMV_INLINE size_t rowsize() const { return itsrs; }
        TMV_INLINE int nlo() const { return LO; }
        TMV_INLINE int nhi() const { return HI; }
        TMV_INLINE int nElements() const { return linsize; }
        TMV_INLINE int stepi() const { return itssi; }
        TMV_INLINE int stepj() const { return itssj; }
        TMV_INLINE int diagstep() const 
        { return _diagmajor ? 1 : stepi()+stepj(); }
        TMV_INLINE bool isconj() const { return false; }
        TMV_INLINE bool isrm() const { return _rowmajor; }
        TMV_INLINE bool iscm() const { return _colmajor; }
        TMV_INLINE bool isdm() const { return _diagmajor; }

    private:

        size_t itscs;
        size_t itsrs;
        size_t linsize;
        CheckedInt<_stepi> itssi;
        CheckedInt<_stepj> itssj;
        AlignedArray<T> itsm1;
        T* itsm;

    }; // ThinBandMatrix

    template <class T, int M, int N, int LO, int HI, int Si, int Sj, int A0>
    struct Traits<ConstSmallBandMatrixView<T,M,N,LO,HI,Si,Sj,A0> >
    {
        typedef typename Traits<T>::real_type real_type;

        enum { SipSj = IntTraits2<Si,Sj>::sum };
        enum { A = (A0 & ~NoDivider) | (
                ( Attrib<A0>::colmajor ? 0 :
                  Attrib<A0>::rowmajor ? 0 :
                  Attrib<A0>::diagmajor ? 0 :
                  Si == 1 ? ColMajor :
                  Sj == 1 ? RowMajor : 
                  SipSj == 1 ? DiagMajor : 0 ) |
                ( Attrib<A0>::withdivider ? 0 : NoDivider ) ) };
        enum { okA = (
                !Attrib<A>::nonunitdiag &&
                !Attrib<A>::unitdiag &&
                !Attrib<A>::zerodiag &&
                !Attrib<A>::packed &&
                !Attrib<A>::lower &&
                !Attrib<A>::upper &&
                (Attrib<A>::nodivider || Attrib<A>::withdivider) &&
                (Attrib<A>::nodivider != int(Attrib<A>::withdivider)) &&
                ( !Attrib<A>::withdivider ||
                  ( Traits<real_type>::isinst && 
                    !Traits<real_type>::isinteger ) ) &&
                ( Traits<T>::iscomplex || !Attrib<A>::conj ) &&
                !( Si!=TMV_UNKNOWN && Si!=1 && Attrib<A>::colmajor ) &&
                !( Sj!=TMV_UNKNOWN && Sj!=1 && Attrib<A>::rowmajor ) &&
                !( SipSj!=TMV_UNKNOWN && SipSj!=1 && Attrib<A>::diagmajor ) )};
        enum { _attrib = A };

        typedef T value_type;
        typedef typename Traits<T>::complex_type complex_type;
        enum { isreal = Traits<T>::isreal };
        enum { iscomplex = Traits<T>::iscomplex };

        typedef ConstSmallBandMatrixView<T,M,N,LO,HI,Si,Sj,A0> type;
        typedef const type& calc_type;
        typedef const type& eval_type;

        enum { _colsize = M };
        enum { _rowsize = N };
        enum { _nlo = LO };
        enum { _nhi = HI };
        enum { _shape = Band };
        enum { _fort = Attrib<A>::fort };
        enum { _calc = true };
        enum { _rowmajor = Attrib<A>::rowmajor };
        enum { _colmajor = Attrib<A>::colmajor };
        enum { _diagmajor = Attrib<A>::diagmajor };
        enum { _stepi = Si };
        enum { _stepj = Sj };
        enum { _diagstep = SipSj };
        enum { _conj = Attrib<A>::conj };
        enum { _checkalias = Attrib<A>::checkalias };
        enum { _canlin = false };
        enum { twoSi = isreal ? int(_stepi) : IntTraits<_stepi>::twoS };
        enum { twoSj = isreal ? int(_stepj) : IntTraits<_stepj>::twoS };
        enum { minMN = IntTraits2<M,N>::min };

        enum { bandknown = LO != TMV_UNKNOWN && HI != TMV_UNKNOWN };
        enum { allknown = M != TMV_UNKNOWN && N != TMV_UNKNOWN && bandknown };
        enum { copyA = (
                (_rowmajor ? RowMajor : _colmajor ? ColMajor : DiagMajor) |
                (_fort ? FortranStyle : CStyle) |
                NoDivider ) };
        typedef typename TypeSelect<allknown,
                SmallBandMatrix<T,M,N,LO,HI,copyA>,
                typename TypeSelect<bandknown,
                ThinBandMatrix<T,LO,HI,copyA>,
                BandMatrix<T,copyA> >::type>::type copy_type;

        enum { _hasdivider = Attrib<A>::withdivider };

#if 0
        typedef QuotXB<1,real_type,type> inverse_type;
        typedef typename TypeSelect<_hasdivider ,
                const InstBandLUD<T>& , BandLUD<copy_type> >::type lud_type;
        typedef typename TypeSelect<_hasdivider ,
                const InstBandQRD<T>& , BandQRD<copy_type> >::type qrd_type;
        typedef typename TypeSelect<_hasdivider ,
                const InstBandSVD<T>& , BandSVD<copy_type> >::type svd_type;
#else
        typedef InvalidType inverse_type;
        typedef InvalidType lud_type;
        typedef InvalidType qrd_type;
        typedef InvalidType svd_type;
#endif

        enum { vecA = (
                (_fort ? FortranStyle : 0) |
                (_conj ? Conj : 0 ) |
                (_checkalias ? CheckAlias : 0) )};
        enum { colA = vecA | (_colmajor ? Unit : 0) };
        enum { rowA = vecA | (_rowmajor ? Unit : 0) };
        enum { diagA = vecA | (_diagmajor ? Unit : 0) };
        enum { ndA = (A & ~AllDivStatus) };
        enum { cstyleA = ndA & ~FortranStyle };
        enum { fstyleA = ndA | FortranStyle };
        enum { nmA = (ndA & ~AllStorageType) };
        enum { cmA = nmA | ColMajor };
        enum { rmA = nmA | RowMajor };
        enum { dmA = nmA | DiagMajor };
        enum { conjA = iscomplex ? (ndA ^ Conj) : int(ndA) };
        enum { trA = _colmajor ? rmA : _rowmajor ? cmA : int(ndA) };
        enum { adjA = iscomplex ? (trA ^ Conj) : int(trA) };
        enum { nonconjA = ndA & ~Conj };
        enum { twosA = isreal ? int(ndA) : (ndA & ~Conj & ~AllStorageType) };
        enum { An = _checkalias ? (ndA & ~CheckAlias) : (ndA | NoAlias) };
        enum { vecAn = _checkalias ? (vecA & ~CheckAlias) : (vecA | NoAlias) };
        enum { colAn = vecAn | (_colmajor ? Unit : 0) };
        enum { rowAn = vecAn | (_rowmajor ? Unit : 0) };
        enum { diagAn = vecAn | (_diagstep == 1 ? Unit : 0) };
        enum { nmAn = (An & ~AllStorageType) };

        enum { xx = TMV_UNKNOWN }; // For brevity.
        typedef ConstVectorView<T,colAn> const_col_sub_type;
        typedef ConstVectorView<T,rowAn> const_row_sub_type;
        typedef ConstSmallVectorView<T,minMN,_diagstep,diagA> const_diag_type;
        typedef ConstSmallVectorView<T,xx,_diagstep,diagA> const_diag_sub_type;

        typedef ConstMatrixView<T,An> const_submatrix_type;
        typedef ConstMatrixView<T,nmAn> const_submatrix_step_type;
        typedef ConstVectorView<T,vecAn> const_subvector_type;
        typedef ConstBandMatrixView<T,An> const_subbandmatrix_type;
        typedef ConstBandMatrixView<T,nmAn> const_subbandmatrix_step_type;
        typedef ConstBandMatrixView<T,An> const_colrange_type;
        typedef ConstBandMatrixView<T,An> const_rowrange_type;
        typedef ConstBandMatrixView<T,An> const_diagrange_type;

        typedef ConstSmallBandMatrixView<T,M,N,LO,HI,_stepi,_stepj,ndA> 
            const_view_type;
        typedef ConstSmallBandMatrixView<T,M,N,LO,HI,_stepi,_stepj,cstyleA> 
            const_cview_type;
        typedef ConstSmallBandMatrixView<T,M,N,LO,HI,_stepi,_stepj,fstyleA> 
            const_fview_type;
        typedef ConstBandMatrixView<T,(_conj ? Conj : NonConj)> const_xview_type;
        typedef ConstSmallBandMatrixView<T,M,N,LO,HI,1,_stepj,cmA> 
            const_cmview_type;
        typedef ConstSmallBandMatrixView<T,M,N,LO,HI,_stepi,1,rmA> 
            const_rmview_type;
        typedef ConstSmallBandMatrixView<T,M,N,LO,HI,_stepi,_stepj,dmA> 
            const_dmview_type;
        typedef ConstSmallBandMatrixView<T,M,N,LO,HI,_stepi,_stepj,conjA> 
            const_conjugate_type;
        typedef ConstSmallBandMatrixView<T,N,M,HI,LO,_stepj,_stepi,trA> 
            const_transpose_type;
        typedef ConstSmallBandMatrixView<T,N,M,HI,LO,_stepj,_stepi,adjA> 
            const_adjoint_type;

        typedef ConstSmallBandMatrixView<T,N,N,0,HI,_stepi,_stepj,ndA> 
            const_upperband_type;
        typedef ConstSmallBandMatrixView<T,M,M,LO,0,_stepi,_stepj,ndA> 
            const_lowerband_type;
        typedef ConstSmallBandMatrixView<T,xx,xx,0,xx,_stepi,_stepj,ndA> 
            const_upperbandoff_type;
        typedef ConstSmallBandMatrixView<T,xx,xx,xx,0,_stepi,_stepj,ndA> 
            const_lowerbandoff_type;

        typedef ConstSmallBandMatrixView<real_type,M,N,LO,HI,twoSi,twoSj,twosA> 
            const_realpart_type;
        typedef const_realpart_type const_imagpart_type;
        typedef ConstVectorView<T,(vecAn|Unit)> const_linearview_type;
        typedef ConstSmallBandMatrixView<T,M,N,LO,HI,_stepi,_stepj,nonconjA> 
            const_nonconj_type;
        typedef SmallBandMatrixView<T,M,N,LO,HI,_stepi,_stepj,ndA> 
            nonconst_type;
    };

    template <class T, int M, int N, int LO, int HI, int Si, int Sj, int A>
    class ConstSmallBandMatrixView : 
        public BaseMatrix_Band<ConstSmallBandMatrixView<T,M,N,LO,HI,Si,Sj,A> >,
        public BandMatrixDivHelper<ConstSmallBandMatrixView<T,M,N,LO,HI,Si,Sj,A> >
    {
    public:
        typedef ConstSmallBandMatrixView<T,M,N,LO,HI,Si,Sj,A> type;

        enum { _colsize = Traits<type>::_colsize };
        enum { _rowsize = Traits<type>::_rowsize };
        enum { _nlo = Traits<type>::_nlo };
        enum { _nhi = Traits<type>::_nhi };
        enum { _shape = Traits<type>::_shape };
        enum { _fort = Traits<type>::_fort };
        enum { _calc = Traits<type>::_calc };
        enum { _rowmajor = Traits<type>::_rowmajor };
        enum { _colmajor = Traits<type>::_colmajor };
        enum { _diagmajor = Traits<type>::_diagmajor };
        enum { _stepi = Traits<type>::_stepi };
        enum { _stepj = Traits<type>::_stepj };
        enum { _diagstep = Traits<type>::_diagstep };
        enum { _conj = Traits<type>::_conj };
        enum { _checkalias = Traits<type>::_checkalias };
        enum { _canlin = Traits<type>::_canlin };
        enum { _attrib = Traits<type>::_attrib };

        //
        // Constructors
        //

        TMV_INLINE ConstSmallBandMatrixView(
            const T* m, size_t cs, size_t rs, int lo, int hi, int si, int sj) :
            itsm(m), itscs(cs), itsrs(rs), itsnlo(lo), itsnhi(hi),
            itssi(si), itssj(sj) 
        { 
            TMVStaticAssert(Traits<type>::okA);
        }

        TMV_INLINE ConstSmallBandMatrixView(
            const T* m, size_t cs, size_t rs, int lo, int hi, int si) :
            itsm(m), itscs(cs), itsrs(rs), itsnlo(lo), itsnhi(hi),
            itssi(si), itssj(Sj) 
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Sj != TMV_UNKNOWN); 
        }

        TMV_INLINE ConstSmallBandMatrixView(
            const T* m, size_t cs, size_t rs, int lo, int hi) :
            itsm(m), itscs(cs), itsrs(rs), itsnlo(lo), itsnhi(hi),
            itssi(Si), itssj(Sj)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Si != TMV_UNKNOWN);
            TMVStaticAssert(Sj != TMV_UNKNOWN); 
        }

        TMV_INLINE ConstSmallBandMatrixView(
            const T* m, size_t cs, size_t rs, int lo) :
            itsm(m), itscs(cs), itsrs(rs), itsnlo(lo), itsnhi(HI),
            itssi(Si), itssj(Sj)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(HI != TMV_UNKNOWN); 
            TMVStaticAssert(Si != TMV_UNKNOWN);
            TMVStaticAssert(Sj != TMV_UNKNOWN); 
        }

        TMV_INLINE ConstSmallBandMatrixView(const T* m, size_t cs, size_t rs) :
            itsm(m), itscs(cs), itsrs(rs), itsnlo(LO), itsnhi(HI),
            itssi(Si), itssj(Sj)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(LO != TMV_UNKNOWN);
            TMVStaticAssert(HI != TMV_UNKNOWN); 
            TMVStaticAssert(Si != TMV_UNKNOWN);
            TMVStaticAssert(Sj != TMV_UNKNOWN); 
        }

        TMV_INLINE ConstSmallBandMatrixView(const T* m, size_t cs) :
            itsm(m), itscs(cs), itsrs(N), itsnlo(LO), itsnhi(HI),
            itssi(Si), itssj(Sj)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(M != TMV_UNKNOWN);
            TMVStaticAssert(LO != TMV_UNKNOWN);
            TMVStaticAssert(HI != TMV_UNKNOWN); 
            TMVStaticAssert(Si != TMV_UNKNOWN);
            TMVStaticAssert(Sj != TMV_UNKNOWN); 
        }

        TMV_INLINE ConstSmallBandMatrixView(const T* m) :
            itsm(m), itscs(M), itsrs(N), itsnlo(LO), itsnhi(HI),
            itssi(Si), itssj(Sj)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(M != TMV_UNKNOWN);
            TMVStaticAssert(N != TMV_UNKNOWN);
            TMVStaticAssert(LO != TMV_UNKNOWN);
            TMVStaticAssert(HI != TMV_UNKNOWN); 
            TMVStaticAssert(Si != TMV_UNKNOWN);
            TMVStaticAssert(Sj != TMV_UNKNOWN); 
        }

        TMV_INLINE ConstSmallBandMatrixView(const type& m2) :
            itsm(m2.cptr()), itscs(m2.colsize()), itsrs(m2.rowsize()),
            itsnlo(m2.nlo()), itsnhi(m2.nhi()),
            itssi(m2.stepi()), itssj(m2.stepj()) 
        { 
            TMVStaticAssert(Traits<type>::okA);
        }

        template <int A2>
        TMV_INLINE ConstSmallBandMatrixView(const ConstBandMatrixView<T,A2>& m2) :
            itsm(m2.cptr()), itscs(m2.colsize()), itsrs(m2.rowsize()),
            itsnlo(m2.nlo()), itsnhi(m2.nhi()),
            itssi(m2.stepi()), itssj(m2.stepj()) 
        { 
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Attrib<A>::conj == int(Attrib<A2>::conj)); 
        }

        template <int A2>
        TMV_INLINE ConstSmallBandMatrixView(const BandMatrixView<T,A2>& m2) :
            itsm(m2.cptr()), itscs(m2.colsize()), itsrs(m2.rowsize()),
            itsnlo(m2.nlo()), itsnhi(m2.nhi()),
            itssi(m2.stepi()), itssj(m2.stepj()) 
        { 
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Attrib<A>::conj == int(Attrib<A2>::conj)); 
        }

        template <int M2, int N2, int LO2, int HI2, int Si2, int Sj2, int A2>
        TMV_INLINE ConstSmallBandMatrixView(
            const ConstSmallBandMatrixView<T,M2,N2,LO2,HI2,Si2,Sj2,A2>& m2) :
            itsm(m2.cptr()), itscs(m2.colsize()), itsrs(m2.rowsize()),
            itsnlo(m2.nlo()), itsnhi(m2.nhi()),
            itssi(m2.stepi()), itssj(m2.stepj()) 
        { 
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Attrib<A>::conj == int(Attrib<A2>::conj)); 
        }

        template <int M2, int N2, int LO2, int HI2, int Si2, int Sj2, int A2>
        TMV_INLINE ConstSmallBandMatrixView(
            const SmallBandMatrixView<T,M2,N2,LO2,HI2,Si2,Sj2,A2>& m2) :
            itsm(m2.cptr()), itscs(m2.colsize()), itsrs(m2.rowsize()),
            itsnlo(m2.nlo()), itsnhi(m2.nhi()),
            itssi(m2.stepi()), itssj(m2.stepj()) 
        { 
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Attrib<A>::conj == int(Attrib<A2>::conj)); 
        }

        TMV_INLINE ~ConstSmallBandMatrixView() {
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

        T cref(int i, int j) const
        { return DoConj<_conj>(itsm[i*stepi()+j*stepj()]); }

        TMV_INLINE size_t colsize() const { return itscs; }
        TMV_INLINE size_t rowsize() const { return itsrs; }
        TMV_INLINE int nlo() const { return itsnlo; }
        TMV_INLINE int nhi() const { return itsnhi; }
        int nElements() const 
        { return colsize()*nlo() + rowsize()*(nhi()+1); }
        TMV_INLINE int stepi() const { return itssi; }
        TMV_INLINE int stepj() const { return itssj; }
        TMV_INLINE int diagstep() const 
        { return _diagmajor ? 1 : stepi()+stepj(); }
        TMV_INLINE bool isconj() const { return _conj; }
        TMV_INLINE bool isrm() const 
        { return _rowmajor || (!_colmajor && !_diagmajor && stepj() == 1); }
        TMV_INLINE bool iscm() const 
        { return _colmajor || (!_rowmajor && !_diagmajor && stepi() == 1); }
        TMV_INLINE bool isdm() const 
        { return _diagmajor || (!_rowmajor && !_colmajor && diagstep() == 1); }

    private :

        const T* itsm;
        const CheckedInt<M> itscs;
        const CheckedInt<N> itsrs;
        const CheckedInt<LO> itsnlo;
        const CheckedInt<HI> itsnhi;
        const CheckedInt<Si> itssi;
        const CheckedInt<Sj> itssj;

    }; // ConstSmallBandMatrixView

    template <class T, int M, int N, int LO, int HI, int Si, int Sj, int A0>
    struct Traits<SmallBandMatrixView<T,M,N,LO,HI,Si,Sj,A0> >
    {
        typedef typename Traits<T>::real_type real_type;

        enum { SipSj = IntTraits2<Si,Sj>::sum };
        enum { A = (A0 & ~NoDivider) | (
                ( Attrib<A0>::colmajor ? 0 :
                  Attrib<A0>::rowmajor ? 0 :
                  Attrib<A0>::diagmajor ? 0 :
                  Si == 1 ? ColMajor :
                  Sj == 1 ? RowMajor : 
                  SipSj == 1 ? DiagMajor : 0 ) |
                ( Attrib<A0>::withdivider ? 0 : NoDivider ) ) };
        enum { okA = (
                !Attrib<A>::nonunitdiag &&
                !Attrib<A>::unitdiag &&
                !Attrib<A>::zerodiag &&
                !Attrib<A>::packed &&
                !Attrib<A>::lower &&
                !Attrib<A>::upper &&
                (Attrib<A>::nodivider || Attrib<A>::withdivider) &&
                (Attrib<A>::nodivider != int(Attrib<A>::withdivider)) &&
                ( !Attrib<A>::withdivider ||
                  ( Traits<real_type>::isinst && 
                    !Traits<real_type>::isinteger ) ) &&
                ( Traits<T>::iscomplex || !Attrib<A>::conj ) &&
                !( Si!=TMV_UNKNOWN && Si!=1 && Attrib<A>::colmajor ) &&
                !( Sj!=TMV_UNKNOWN && Sj!=1 && Attrib<A>::rowmajor ) &&
                !( SipSj!=TMV_UNKNOWN && SipSj!=1 && Attrib<A>::diagmajor ) )};
        enum { _attrib = A };

        typedef T value_type;
        typedef typename Traits<T>::complex_type complex_type;
        enum { isreal = Traits<T>::isreal };
        enum { iscomplex = Traits<T>::iscomplex };

        typedef SmallBandMatrixView<T,M,N,LO,HI,Si,Sj,A0> type;
        typedef const type& calc_type;
        typedef const type& eval_type;

        enum { _colsize = M };
        enum { _rowsize = N };
        enum { _nlo = LO };
        enum { _nhi = HI };
        enum { _shape = Band };
        enum { _fort = Attrib<A>::fort };
        enum { _calc = true };
        enum { _rowmajor = Attrib<A>::rowmajor };
        enum { _colmajor = Attrib<A>::colmajor };
        enum { _diagmajor = Attrib<A>::diagmajor };
        enum { _stepi = Si };
        enum { _stepj = Sj };
        enum { _diagstep = SipSj };
        enum { _conj = Attrib<A>::conj };
        enum { _checkalias = Attrib<A>::checkalias };
        enum { _canlin = false };
        enum { twoSi = isreal ? int(_stepi) : IntTraits<_stepi>::twoS };
        enum { twoSj = isreal ? int(_stepj) : IntTraits<_stepj>::twoS };
        enum { minMN = IntTraits2<M,N>::min };

        enum { bandknown = LO != TMV_UNKNOWN && HI != TMV_UNKNOWN };
        enum { allknown = M != TMV_UNKNOWN && N != TMV_UNKNOWN && bandknown };
        enum { copyA = (
                (_rowmajor ? RowMajor : _colmajor ? ColMajor : DiagMajor) |
                (_fort ? FortranStyle : CStyle) |
                NoDivider ) };
        typedef typename TypeSelect<allknown,
                SmallBandMatrix<T,M,N,LO,HI,copyA>,
                typename TypeSelect<bandknown,
                ThinBandMatrix<T,LO,HI,copyA>,
                BandMatrix<T,copyA> >::type>::type copy_type;

        enum { _hasdivider = Attrib<A>::withdivider };

#if 0
        typedef QuotXB<1,real_type,type> inverse_type;
        typedef typename TypeSelect<_hasdivider ,
                const InstBandLUD<T>& , BandLUD<copy_type> >::type lud_type;
        typedef typename TypeSelect<_hasdivider ,
                const InstBandQRD<T>& , BandQRD<copy_type> >::type qrd_type;
        typedef typename TypeSelect<_hasdivider ,
                const InstBandSVD<T>& , BandSVD<copy_type> >::type svd_type;
#else
        typedef InvalidType inverse_type;
        typedef InvalidType lud_type;
        typedef InvalidType qrd_type;
        typedef InvalidType svd_type;
#endif

        enum { vecA = (
                (_fort ? FortranStyle : 0) |
                (_conj ? Conj : 0 ) |
                (_checkalias ? CheckAlias : 0) )};
        enum { colA = vecA | (_colmajor ? Unit : 0) };
        enum { rowA = vecA | (_rowmajor ? Unit : 0) };
        enum { diagA = vecA | (_diagmajor ? Unit : 0) };
        enum { ndA = (A & ~AllDivStatus) };
        enum { cstyleA = ndA & ~FortranStyle };
        enum { fstyleA = ndA | FortranStyle };
        enum { nmA = (ndA & ~AllStorageType) };
        enum { cmA = nmA | ColMajor };
        enum { rmA = nmA | RowMajor };
        enum { dmA = nmA | DiagMajor };
        enum { conjA = iscomplex ? (ndA ^ Conj) : int(ndA) };
        enum { trA = _colmajor ? rmA : _rowmajor ? cmA : int(ndA) };
        enum { adjA = iscomplex ? (trA ^ Conj) : int(trA) };
        enum { nonconjA = ndA & ~Conj };
        enum { twosA = isreal ? int(ndA) : (ndA & ~Conj & ~AllStorageType) };
        enum { An = _checkalias ? (ndA & ~CheckAlias) : (ndA | NoAlias) };
        enum { vecAn = _checkalias ? (vecA & ~CheckAlias) : (vecA | NoAlias) };
        enum { colAn = vecAn | (_colmajor ? Unit : 0) };
        enum { rowAn = vecAn | (_rowmajor ? Unit : 0) };
        enum { diagAn = vecAn | (_diagstep == 1 ? Unit : 0) };
        enum { nmAn = (An & ~AllStorageType) };

        enum { xx = TMV_UNKNOWN }; // For brevity.
        typedef ConstVectorView<T,colAn> const_col_sub_type;
        typedef ConstVectorView<T,rowAn> const_row_sub_type;
        typedef ConstSmallVectorView<T,minMN,_diagstep,diagA> const_diag_type;
        typedef ConstSmallVectorView<T,xx,_diagstep,diagA> const_diag_sub_type;

        typedef ConstMatrixView<T,An> const_submatrix_type;
        typedef ConstMatrixView<T,nmAn> const_submatrix_step_type;
        typedef ConstVectorView<T,vecAn> const_subvector_type;
        typedef ConstBandMatrixView<T,An> const_subbandmatrix_type;
        typedef ConstBandMatrixView<T,nmAn> const_subbandmatrix_step_type;
        typedef ConstBandMatrixView<T,An> const_colrange_type;
        typedef ConstBandMatrixView<T,An> const_rowrange_type;
        typedef ConstBandMatrixView<T,An> const_diagrange_type;

        typedef ConstSmallBandMatrixView<T,M,N,LO,HI,_stepi,_stepj,ndA> 
            const_view_type;
        typedef ConstSmallBandMatrixView<T,M,N,LO,HI,_stepi,_stepj,cstyleA> 
            const_cview_type;
        typedef ConstSmallBandMatrixView<T,M,N,LO,HI,_stepi,_stepj,fstyleA> 
            const_fview_type;
        typedef ConstBandMatrixView<T,(_conj ? Conj : NonConj)> const_xview_type;
        typedef ConstSmallBandMatrixView<T,M,N,LO,HI,1,_stepj,cmA> 
            const_cmview_type;
        typedef ConstSmallBandMatrixView<T,M,N,LO,HI,_stepi,1,rmA> 
            const_rmview_type;
        typedef ConstSmallBandMatrixView<T,M,N,LO,HI,_stepi,_stepj,dmA> 
            const_dmview_type;
        typedef ConstSmallBandMatrixView<T,M,N,LO,HI,_stepi,_stepj,conjA> 
            const_conjugate_type;
        typedef ConstSmallBandMatrixView<T,N,M,HI,LO,_stepj,_stepi,trA> 
            const_transpose_type;
        typedef ConstSmallBandMatrixView<T,N,M,HI,LO,_stepj,_stepi,adjA> 
            const_adjoint_type;

        typedef ConstSmallBandMatrixView<T,N,N,0,HI,_stepi,_stepj,ndA> 
            const_upperband_type;
        typedef ConstSmallBandMatrixView<T,M,M,LO,0,_stepi,_stepj,ndA> 
            const_lowerband_type;
        typedef ConstSmallBandMatrixView<T,xx,xx,0,xx,_stepi,_stepj,ndA> 
            const_upperbandoff_type;
        typedef ConstSmallBandMatrixView<T,xx,xx,xx,0,_stepi,_stepj,ndA> 
            const_lowerbandoff_type;

        typedef ConstSmallBandMatrixView<real_type,M,N,LO,HI,twoSi,twoSj,twosA> 
            const_realpart_type;
        typedef const_realpart_type const_imagpart_type;
        typedef ConstVectorView<T,(vecAn|Unit)> const_linearview_type;
        typedef ConstSmallBandMatrixView<T,M,N,LO,HI,_stepi,_stepj,nonconjA> 
            const_nonconj_type;
        typedef SmallBandMatrixView<T,M,N,LO,HI,_stepi,_stepj,ndA> 
            nonconst_type;

        typedef VectorView<T,colAn> col_sub_type;
        typedef VectorView<T,rowAn> row_sub_type;
        typedef SmallVectorView<T,minMN,_diagstep,diagA> diag_type;
        typedef SmallVectorView<T,xx,_diagstep,diagA> diag_sub_type;

        typedef MatrixView<T,An> submatrix_type;
        typedef MatrixView<T,nmAn> submatrix_step_type;
        typedef VectorView<T,vecAn> subvector_type;
        typedef BandMatrixView<T,An> subbandmatrix_type;
        typedef BandMatrixView<T,nmAn> subbandmatrix_step_type;
        typedef BandMatrixView<T,An> colrange_type;
        typedef BandMatrixView<T,An> rowrange_type;
        typedef BandMatrixView<T,An> diagrange_type;

        typedef SmallBandMatrixView<T,M,N,LO,HI,_stepi,_stepj,ndA> 
            view_type;
        typedef SmallBandMatrixView<T,M,N,LO,HI,_stepi,_stepj,cstyleA> 
            cview_type;
        typedef SmallBandMatrixView<T,M,N,LO,HI,_stepi,_stepj,fstyleA> 
            fview_type;
        typedef BandMatrixView<T,(_conj ? Conj : NonConj)> xview_type;
        typedef SmallBandMatrixView<T,M,N,LO,HI,1,_stepj,cmA> 
            cmview_type;
        typedef SmallBandMatrixView<T,M,N,LO,HI,_stepi,1,rmA> 
            rmview_type;
        typedef SmallBandMatrixView<T,M,N,LO,HI,_stepi,_stepj,dmA> 
            dmview_type;
        typedef SmallBandMatrixView<T,M,N,LO,HI,_stepi,_stepj,conjA> 
            conjugate_type;
        typedef SmallBandMatrixView<T,N,M,HI,LO,_stepj,_stepi,trA> 
            transpose_type;
        typedef SmallBandMatrixView<T,N,M,HI,LO,_stepj,_stepi,adjA> 
            adjoint_type;

        typedef SmallBandMatrixView<T,N,N,0,HI,_stepi,_stepj,ndA> 
            upperband_type;
        typedef SmallBandMatrixView<T,M,M,LO,0,_stepi,_stepj,ndA> 
            lowerband_type;
        typedef SmallBandMatrixView<T,xx,xx,0,xx,_stepi,_stepj,ndA> 
            upperbandoff_type;
        typedef SmallBandMatrixView<T,xx,xx,xx,0,_stepi,_stepj,ndA> 
            lowerbandoff_type;

        typedef SmallBandMatrixView<real_type,M,N,LO,HI,twoSi,twoSj,twosA> 
            realpart_type;
        typedef realpart_type imagpart_type;
        typedef VectorView<T,(vecAn|Unit)> linearview_type;
        typedef SmallBandMatrixView<T,M,N,LO,HI,_stepi,_stepj,nonconjA> 
            nonconj_type;
    };

    template <class T, int M, int N, int LO, int HI, int Si, int Sj, int A>
    class SmallBandMatrixView : 
        public BaseMatrix_Band_Mutable<SmallBandMatrixView<T,M,N,LO,HI,Si,Sj,A> >,
        public BandMatrixDivHelper<SmallBandMatrixView<T,M,N,LO,HI,Si,Sj,A> >
    {
    public:

        typedef SmallBandMatrixView<T,M,N,LO,HI,Si,Sj,A> type;
        typedef BaseMatrix_Band_Mutable<type> base_mut;
        typedef typename base_mut::reference reference;

        enum { _colsize = Traits<type>::_colsize };
        enum { _rowsize = Traits<type>::_rowsize };
        enum { _nlo = Traits<type>::_nlo };
        enum { _nhi = Traits<type>::_nhi };
        enum { _shape = Traits<type>::_shape };
        enum { _fort = Traits<type>::_fort };
        enum { _calc = Traits<type>::_calc };
        enum { _rowmajor = Traits<type>::_rowmajor };
        enum { _colmajor = Traits<type>::_colmajor };
        enum { _diagmajor = Traits<type>::_diagmajor };
        enum { _stepi = Traits<type>::_stepi };
        enum { _stepj = Traits<type>::_stepj };
        enum { _diagstep = Traits<type>::_diagstep };
        enum { _conj = Traits<type>::_conj };
        enum { _checkalias = Traits<type>::_checkalias };
        enum { _canlin = Traits<type>::_canlin };
        enum { _attrib = Traits<type>::_attrib };

        //
        // Constructors
        //

        TMV_INLINE SmallBandMatrixView(
            T* m, size_t cs, size_t rs, int lo, int hi, int si, int sj) :
            itsm(m), itscs(cs), itsrs(rs), itsnlo(lo), itsnhi(hi),
            itssi(si), itssj(sj) 
        {
            TMVStaticAssert(Traits<type>::okA);
        }

        TMV_INLINE SmallBandMatrixView(
            T* m, size_t cs, size_t rs, int lo, int hi, int si) :
            itsm(m), itscs(cs), itsrs(rs), itsnlo(lo), itsnhi(hi),
            itssi(si), itssj(Sj) 
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Sj != TMV_UNKNOWN); 
        }

        TMV_INLINE SmallBandMatrixView(
            T* m, size_t cs, size_t rs, int lo, int hi) :
            itsm(m), itscs(cs), itsrs(rs), itsnlo(lo), itsnhi(hi),
            itssi(Si), itssj(Sj)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Si != TMV_UNKNOWN);
            TMVStaticAssert(Sj != TMV_UNKNOWN); 
        }

        TMV_INLINE SmallBandMatrixView(T* m, size_t cs, size_t rs, int lo) :
            itsm(m), itscs(cs), itsrs(rs), itsnlo(lo), itsnhi(HI),
            itssi(Si), itssj(Sj)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(HI != TMV_UNKNOWN); 
            TMVStaticAssert(Si != TMV_UNKNOWN);
            TMVStaticAssert(Sj != TMV_UNKNOWN); 
        }

        TMV_INLINE SmallBandMatrixView(T* m, size_t cs, size_t rs) :
            itsm(m), itscs(cs), itsrs(rs), itsnlo(LO), itsnhi(HI),
            itssi(Si), itssj(Sj)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(LO != TMV_UNKNOWN);
            TMVStaticAssert(HI != TMV_UNKNOWN); 
            TMVStaticAssert(Si != TMV_UNKNOWN);
            TMVStaticAssert(Sj != TMV_UNKNOWN); 
        }

        TMV_INLINE SmallBandMatrixView(T* m, size_t cs) :
            itsm(m), itscs(cs), itsrs(N), itsnlo(LO), itsnhi(HI),
            itssi(Si), itssj(Sj)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(M != TMV_UNKNOWN);
            TMVStaticAssert(LO != TMV_UNKNOWN);
            TMVStaticAssert(HI != TMV_UNKNOWN); 
            TMVStaticAssert(Si != TMV_UNKNOWN);
            TMVStaticAssert(Sj != TMV_UNKNOWN); 
        }

        TMV_INLINE SmallBandMatrixView(T* m) :
            itsm(m), itscs(M), itsrs(N), itsnlo(LO), itsnhi(HI),
            itssi(Si), itssj(Sj)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(M != TMV_UNKNOWN);
            TMVStaticAssert(N != TMV_UNKNOWN);
            TMVStaticAssert(LO != TMV_UNKNOWN);
            TMVStaticAssert(HI != TMV_UNKNOWN); 
            TMVStaticAssert(Si != TMV_UNKNOWN);
            TMVStaticAssert(Sj != TMV_UNKNOWN); 
        }


        TMV_INLINE SmallBandMatrixView(const type& m2) :
            itsm(m2.itsm), itscs(m2.colsize()), itsrs(m2.rowsize()),
            itsnlo(m2.nlo()), itsnhi(m2.nhi()),
            itssi(m2.stepi()), itssj(m2.stepj())
        {
            TMVStaticAssert(Traits<type>::okA);
        }

        template <int A2>
        TMV_INLINE SmallBandMatrixView(BandMatrixView<T,A2> m2) :
            itsm(m2.ptr()), itscs(m2.colsize()), itsrs(m2.rowsize()),
            itsnlo(m2.nlo()), itsnhi(m2.nhi()),
            itssi(m2.stepi()), itssj(m2.stepj())
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Attrib<A>::conj == int(Attrib<A2>::conj)); 
        }

        template <int M2, int N2, int LO2, int HI2, int Si2, int Sj2, int A2>
        TMV_INLINE SmallBandMatrixView(
            SmallBandMatrixView<T,M2,N2,LO2,HI2,Si2,Sj2,A2> m2) :
            itsm(m2.ptr()), itscs(m2.colsize()), itsrs(m2.rowsize()),
            itsnlo(m2.nlo()), itsnhi(m2.nhi()),
            itssi(m2.stepi()), itssj(m2.stepj())
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Attrib<A>::conj == int(Attrib<A2>::conj)); 
        }

        TMV_INLINE ~SmallBandMatrixView() {
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

        template <class M2>
        TMV_INLINE type& operator=(const BaseMatrix_Band<M2>& m2)
        { base_mut::operator=(m2); return *this; }

        template <class M2>
        TMV_INLINE type& operator=(const BaseMatrix_Tri<M2>& m2)
        { base_mut::operator=(m2); return *this; }

        template <class M2>
        TMV_INLINE type& operator=(const BaseMatrix_Diag<M2>& m2)
        { base_mut::operator=(m2); return *this; }

        TMV_INLINE type& operator=(T x)
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

        TMV_INLINE size_t colsize() const { return itscs; }
        TMV_INLINE size_t rowsize() const { return itsrs; }
        TMV_INLINE int nlo() const { return itsnlo; }
        TMV_INLINE int nhi() const { return itsnhi; }
        int nElements() const 
        { return colsize()*nlo() + rowsize()*(nhi()+1); }
        TMV_INLINE int stepi() const { return itssi; }
        TMV_INLINE int stepj() const { return itssj; }
        TMV_INLINE int diagstep() const
        { return _diagmajor ? 1 : stepi()+stepj(); }
        TMV_INLINE bool isconj() const { return _conj; }
        TMV_INLINE bool isrm() const 
        { return _rowmajor || (!_colmajor && !_diagmajor && stepj() == 1); }
        TMV_INLINE bool iscm() const 
        { return _colmajor || (!_rowmajor && !_diagmajor && stepi() == 1); }
        TMV_INLINE bool isdm() const 
        { return _diagmajor || (!_rowmajor && !_colmajor && diagstep() == 1); }

    private :

        T* itsm;
        const CheckedInt<M> itscs;
        const CheckedInt<N> itsrs;
        const CheckedInt<LO> itsnlo;
        const CheckedInt<HI> itsnhi;
        const CheckedInt<Si> itssi;
        const CheckedInt<Sj> itssj;

    }; // SmallBandMatrixView


    //-------------------------------------------------------------------------

    //
    // Special Creators:
    //
    // View parts of a Matrix, another BandMatrix, a DiagMatrix,
    // or a TriMatrix as a BandMatrix:
    //
    //   BandMatrixViewOf(Matrix m, nlo, nhi)
    //   BandMatrixViewOf(BandMatrix m, nlo, nhi)
    //   BandMatrixViewOf(DiagMatrix m)
    //   BandMatrixViewOf(UpperTriMatrix m)
    //   BandMatrixViewOf(UpperTriMatrix m, nhi)
    //   BandMatrixViewOf(LowerTriMatrix m)
    //   BandMatrixViewOf(LowerTriMatrix m, nlo)
    //

    template <class M1>
    struct BMVO
    {
        typedef typename M1::value_type T;
        enum { M = M1::_colsize };
        enum { N = M1::_rowoize };
        enum { Si = M1::_stepi };
        enum { Sj = M1::_stepj };
        enum { A = (
                ( M1::_conj ? Conj : NonConj ) |
                ( M1::_fort ? FortranStyle : CStyle ) |
                ( M1::_colmajor ? ColMajor : M1::_rowmajor ? RowMajor : 
                  Si == 1 ? ColMajor : Sj == 1 ? RowMajor : 0 ) |
                NoDivider )};
        enum { xx = TMV_UNKNOWN };
        enum { small = (
                M != TMV_UNKNOWN || 
                N != TMV_UNKNOWN ||
                (Si != TMV_UNKNOWN && Si != 1) ||
                (Sj != TMV_UNKNOWN && Sj != 1) )};

        typedef typename TypeSelect<small,
            ConstBandMatrixView<T,A> ,
            ConstSmallBandMatrixView<T,M,N,xx,xx,Si,Sj,A> >::type cb;
        typedef typename TypeSelect<small,
            BandMatrixView<T,A> ,
            SmallBandMatrixView<T,M,N,xx,xx,Si,Sj,A> >::type b;
        typedef ConstSmallBandMatrixView<T,M,N,0,xx,Si,Sj,A> cbu;
        typedef SmallBandMatrixView<T,M,N,0,xx,Si,Sj,A> bu;
        typedef ConstSmallBandMatrixView<T,M,N,xx,0,Si,Sj,A> cbl;
        typedef SmallBandMatrixView<T,M,N,xx,0,Si,Sj,A> bl;
        typedef ConstSmallBandMatrixView<T,M,N,0,0,Si,Sj,A> cbd;
        typedef SmallBandMatrixView<T,M,N,0,0,Si,Sj,A> bd;
    };

    template <class M>
    static TMV_INLINE typename BMVO<M>::cb BandMatrixViewOf(
        const BaseMatrix_Rec<M>& m, int lo, int hi)
    {
        TMVAssert(lo < m.colsize());
        TMVAssert(hi < m.rowsize());
        return typename BMVO<M>::cb(
            m.cptr(),m.colsize(),m.rowsize(),lo,hi,m.stepi(),m.stepj());
    }
    template <class M>
    static TMV_INLINE typename BMVO<M>::b BandMatrixViewOf(
        BaseMatrix_Rec_Mutable<M>& m, int lo, int hi)
    {
        TMVAssert(lo < m.colsize());
        TMVAssert(hi < m.rowsize());
        return typename BMVO<M>::b(
            m.ptr(),m.colsize(),m.rowsize(),lo,hi,m.stepi(),m.stepj());
    }
    template <class T, int A>
    static TMV_INLINE typename BMVO<MatrixView<T,A> >::b BandMatrixViewOf(
        MatrixView<T,A> m, int lo, int hi)
    {
        TMVAssert(lo < m.colsize());
        TMVAssert(hi < m.rowsize());
        return typename BMVO<MatrixView<T,A> >::b(
            m.ptr(),m.colsize(),m.rowsize(),lo,hi,m.stepi(),m.stepj());
    }
    template <class T, int M, int N, int Si, int Sj, int A>
    static TMV_INLINE typename BMVO<SmallMatrixView<T,M,N,Si,Sj,A> >::b BandMatrixViewOf(
        SmallMatrixView<T,M,N,Si,Sj,A> m, int lo, int hi)
    {
        TMVAssert(lo < m.colsize());
        TMVAssert(hi < m.rowsize());
        return typename BMVO<SmallMatrixView<T,M,N,Si,Sj,A> >::b(
            m.ptr(),m.colsize(),m.rowsize(),lo,hi,m.stepi(),m.stepj());
    }

    template <class M>
    static TMV_INLINE typename BMVO<M>::cb BandMatrixViewOf(
        const BaseMatrix_Band<M>& m, int lo, int hi)
    {
        TMVAssert(lo <= m.nlo());
        TMVAssert(hi <= m.nhi());
        return typename BMVO<M>::cb(
            m.cptr(),m.colsize(),m.rowsize(),lo,hi,m.stepi(),m.stepj());
    }
    template <class M>
    static TMV_INLINE typename BMVO<M>::b BandMatrixViewOf(
        BaseMatrix_Band_Mutable<M>& m, int lo, int hi)
    {
        TMVAssert(lo <= m.nlo());
        TMVAssert(hi <= m.nhi());
        return typename BMVO<M>::b(
            m.ptr(),m.colsize(),m.rowsize(),lo,hi,m.stepi(),m.stepj());
    }
    template <class T, int A>
    static TMV_INLINE typename BMVO<BandMatrixView<T,A> >::b BandMatrixViewOf(
        BandMatrixView<T,A> m, int lo, int hi)
    {
        TMVAssert(lo <= m.nlo());
        TMVAssert(hi <= m.nhi());
        return typename BMVO<BandMatrixView<T,A> >::b(
            m.ptr(),m.colsize(),m.rowsize(),lo,hi,m.stepi(),m.stepj());
    }
    template <class T, int M, int N, int LO, int HI, int Si, int Sj, int A>
    static TMV_INLINE typename BMVO<SmallBandMatrixView<T,M,N,LO,HI,Si,Sj,A> >::b BandMatrixViewOf(
        SmallBandMatrixView<T,M,N,LO,HI,Si,Sj,A> m, int lo, int hi)
    {
        TMVAssert(lo <= m.nlo());
        TMVAssert(hi <= m.nhi());
        return typename BMVO<SmallBandMatrixView<T,M,N,LO,HI,Si,Sj,A> >::b(
            m.ptr(),m.colsize(),m.rowsize(),lo,hi,m.stepi(),m.stepj());
    }

    // MJ
    template <class M>
    static TMV_INLINE typename BMVOTri<M>::cb BandMatrixViewOf(
        const BaseMatrix_Tri<M>& m, int lohi)
    {
        TMVAssert(lohi < m.size());
        return typename BMVOTri<M>::cb(
            m.cptr(),m.colsize(),m.rowsize(),
            M::_upper ? 0 : lohi, M::_upper ? lohi : 0,
            m.stepi(),m.stepj());
    }
    template <class M>
    static TMV_INLINE typename BMVOTri<M>::b BandMatrixViewOf(
        BaseMatrix_Tri_Mutable<M>& m, int lohi)
    {
        TMVAssert(lohi < m.size());
        return typename BMVOTri<M>::b(
            m.ptr(),m.colsize(),m.rowsize(),
            M::_upper ? 0 : lohi, M::_upper ? lohi : 0,
            m.stepi(),m.stepj());
    }
    template <class T, int A>
    static TMV_INLINE typename BMVO<MatrixView<T,A> >::bu BandMatrixViewOf(
        UpperTriMatrixView<T,A> m, int hi)
    {
        return typename BMVO<MatrixView<T,A> >::b(
            m.ptr(),m.colsize(),m.rowsize(),lo,hi,m.stepi(),m.stepj());
    }
    template <class T, int N, int Si, int Sj, int A>
    static TMV_INLINE typename BMVO<SmallUpperTriMatrixView<T,N,Si,Sj,A> >::bu BandMatrixViewOf(
        SmallUpperTriMatrixView<T,N,Si,Sj,A> m, int lo, int hi)
    {
        return typename BMVO<SmallMatrixView<T,M,N,A> >::b(
            m.ptr(),m.colsize(),m.rowsize(),lo,hi,m.stepi(),m.stepj());
    }

    template <class M>
    static TMV_INLINE typename BMVO<M>::cb BandMatrixViewOf(
        const BaseMatrix_Rec<M>& m, int lo, int hi)
    {
        return typename BMVO<M>::cb(
            m.cptr(),m.colsize(),m.rowsize(),lo,hi,m.stepi(),m.stepj());
    }
    template <class M>
    static TMV_INLINE typename BMVO<M>::b BandMatrixViewOf(
        BaseMatrix_Rec_Mutable<M>& m, int lo, int hi)
    {
        return typename BMVO<M>::b(
            m.ptr(),m.colsize(),m.rowsize(),lo,hi,m.stepi(),m.stepj());
    }
    template <class T, int A>
    static TMV_INLINE typename BMVO<MatrixView<T,A> >::b BandMatrixViewOf(
        MatrixView<T,A> m, int lo, int hi)
    {
        return typename BMVO<MatrixView<T,A> >::b(
            m.ptr(),m.colsize(),m.rowsize(),lo,hi,m.stepi(),m.stepj());
    }
    template <class T, int M, int N, int A>
    static TMV_INLINE typename BMVO<SmallMatrixView<T,M,N,A> >::b BandMatrixViewOf(
        SmallMatrixView<T,M,N,A> m, int lo, int hi)
    {
        return typename BMVO<SmallMatrixView<T,M,N,A> >::b(
            m.ptr(),m.colsize(),m.rowsize(),lo,hi,m.stepi(),m.stepj());
    }

    template <class M>
    static TMV_INLINE typename BMVO<M>::cb BandMatrixViewOf(
        const BaseMatrix_Rec<M>& m, int lo, int hi)
    {
        return typename BMVO<M>::cb(
            m.cptr(),m.colsize(),m.rowsize(),lo,hi,m.stepi(),m.stepj());
    }
    template <class M>
    static TMV_INLINE typename BMVO<M>::b BandMatrixViewOf(
        BaseMatrix_Rec_Mutable<M>& m, int lo, int hi)
    {
        return typename BMVO<M>::b(
            m.ptr(),m.colsize(),m.rowsize(),lo,hi,m.stepi(),m.stepj());
    }
    template <class T, int A>
    static TMV_INLINE typename BMVO<MatrixView<T,A> >::b BandMatrixViewOf(
        MatrixView<T,A> m, int lo, int hi)
    {
        return typename BMVO<MatrixView<T,A> >::b(
            m.ptr(),m.colsize(),m.rowsize(),lo,hi,m.stepi(),m.stepj());
    }
    template <class T, int M, int N, int A>
    static TMV_INLINE typename BMVO<SmallMatrixView<T,M,N,A> >::b BandMatrixViewOf(
        SmallMatrixView<T,M,N,A> m, int lo, int hi)
    {
        return typename BMVO<SmallMatrixView<T,M,N,A> >::b(
            m.ptr(),m.colsize(),m.rowsize(),lo,hi,m.stepi(),m.stepj());
    }

    template <class M>
    static TMV_INLINE typename BMVO<M>::cb BandMatrixViewOf(
        const BaseMatrix_Rec<M>& m, int lo, int hi)
    {
        return typename BMVO<M>::cb(
            m.cptr(),m.colsize(),m.rowsize(),lo,hi,m.stepi(),m.stepj());
    }
    template <class M>
    static TMV_INLINE typename BMVO<M>::b BandMatrixViewOf(
        BaseMatrix_Rec_Mutable<M>& m, int lo, int hi)
    {
        return typename BMVO<M>::b(
            m.ptr(),m.colsize(),m.rowsize(),lo,hi,m.stepi(),m.stepj());
    }
    template <class T, int A>
    static TMV_INLINE typename BMVO<MatrixView<T,A> >::b BandMatrixViewOf(
        MatrixView<T,A> m, int lo, int hi)
    {
        return typename BMVO<MatrixView<T,A> >::b(
            m.ptr(),m.colsize(),m.rowsize(),lo,hi,m.stepi(),m.stepj());
    }
    template <class T, int M, int N, int A>
    static TMV_INLINE typename BMVO<SmallMatrixView<T,M,N,A> >::b BandMatrixViewOf(
        SmallMatrixView<T,M,N,A> m, int lo, int hi)
    {
        return typename BMVO<SmallMatrixView<T,M,N,A> >::b(
            m.ptr(),m.colsize(),m.rowsize(),lo,hi,m.stepi(),m.stepj());
    }

    template <class M>
    static TMV_INLINE typename BMVO<M>::cb BandMatrixViewOf(
        const BaseMatrix_Rec<M>& m, int lo, int hi)
    {
        return typename BMVO<M>::cb(
            m.cptr(),m.colsize(),m.rowsize(),lo,hi,m.stepi(),m.stepj());
    }
    template <class M>
    static TMV_INLINE typename BMVO<M>::b BandMatrixViewOf(
        BaseMatrix_Rec_Mutable<M>& m, int lo, int hi)
    {
        return typename BMVO<M>::b(
            m.ptr(),m.colsize(),m.rowsize(),lo,hi,m.stepi(),m.stepj());
    }
    template <class T, int A>
    static TMV_INLINE typename BMVO<MatrixView<T,A> >::b BandMatrixViewOf(
        MatrixView<T,A> m, int lo, int hi)
    {
        return typename BMVO<MatrixView<T,A> >::b(
            m.ptr(),m.colsize(),m.rowsize(),lo,hi,m.stepi(),m.stepj());
    }
    template <class T, int M, int N, int A>
    static TMV_INLINE typename BMVO<SmallMatrixView<T,M,N,A> >::b BandMatrixViewOf(
        SmallMatrixView<T,M,N,A> m, int lo, int hi)
    {
        return typename BMVO<SmallMatrixView<T,M,N,A> >::b(
            m.ptr(),m.colsize(),m.rowsize(),lo,hi,m.stepi(),m.stepj());
    }

    template <class M>
    static TMV_INLINE typename BMVO<M>::cb BandMatrixViewOf(
        const BaseMatrix_Rec<M>& m, int lo, int hi)
    {
        return typename BMVO<M>::cb(
            m.cptr(),m.colsize(),m.rowsize(),lo,hi,m.stepi(),m.stepj());
    }
    template <class M>
    static TMV_INLINE typename BMVO<M>::b BandMatrixViewOf(
        BaseMatrix_Rec_Mutable<M>& m, int lo, int hi)
    {
        return typename BMVO<M>::b(
            m.ptr(),m.colsize(),m.rowsize(),lo,hi,m.stepi(),m.stepj());
    }
    template <class T, int A>
    static TMV_INLINE typename BMVO<MatrixView<T,A> >::b BandMatrixViewOf(
        MatrixView<T,A> m, int lo, int hi)
    {
        return typename BMVO<MatrixView<T,A> >::b(
            m.ptr(),m.colsize(),m.rowsize(),lo,hi,m.stepi(),m.stepj());
    }
    template <class T, int M, int N, int A>
    static TMV_INLINE typename BMVO<SmallMatrixView<T,M,N,A> >::b BandMatrixViewOf(
        SmallMatrixView<T,M,N,A> m, int lo, int hi)
    {
        return typename BMVO<SmallMatrixView<T,M,N,A> >::b(
            m.ptr(),m.colsize(),m.rowsize(),lo,hi,m.stepi(),m.stepj());
    }

    template <class M>
    static TMV_INLINE typename BMVO<M>::cb BandMatrixViewOf(
        const BaseMatrix_Rec<M>& m, int lo, int hi)
    {
        return typename BMVO<M>::cb(
            m.cptr(),m.colsize(),m.rowsize(),lo,hi,m.stepi(),m.stepj());
    }
    template <class M>
    static TMV_INLINE typename BMVO<M>::b BandMatrixViewOf(
        BaseMatrix_Rec_Mutable<M>& m, int lo, int hi)
    {
        return typename BMVO<M>::b(
            m.ptr(),m.colsize(),m.rowsize(),lo,hi,m.stepi(),m.stepj());
    }
    template <class T, int A>
    static TMV_INLINE typename BMVO<MatrixView<T,A> >::b BandMatrixViewOf(
        MatrixView<T,A> m, int lo, int hi)
    {
        return typename BMVO<MatrixView<T,A> >::b(
            m.ptr(),m.colsize(),m.rowsize(),lo,hi,m.stepi(),m.stepj());
    }
    template <class T, int M, int N, int A>
    static TMV_INLINE typename BMVO<SmallMatrixView<T,M,N,A> >::b BandMatrixViewOf(
        SmallMatrixView<T,M,N,A> m, int lo, int hi)
    {
        return typename BMVO<SmallMatrixView<T,M,N,A> >::b(
            m.ptr(),m.colsize(),m.rowsize(),lo,hi,m.stepi(),m.stepj());
    }



    //
    // Swap
    //

    template <class T, int LO, int HI, int A0, int A1>
    static TMV_INLINE void Swap(
        ThinBandMatrix<T,LO,HI,A0,A1>& m1, ThinBandMatrix<T,LO,HI,A0,A1>& m2)
    { m1.swapWith(m2); }
    template <class T, int M, int N, int LO, int HI, int Si, int Sj, int A, class MM>
    static TMV_INLINE void Swap(
        BaseMatrix_Band_Mutable<MM>& m1,
        SmallBandMatrixView<T,M,N,LO,HI,Si,Sj,A> m2)
    { DoSwap(m1,m2); }
    template <class T, int M, int N, int LO, int HI, int Si, int Sj, int A, class MM>
    static TMV_INLINE void Swap(
        SmallBandMatrixView<T,M,N,LO,HI,Si,Sj,A> m1,
        BaseMatrix_Band_Mutable<MM>& m2)
    { DoSwap(m1,m2); }
    template <class T, int M, int N, int LO, int HI, int Si1, int Sj1, int A1, int Si2, int Sj2, int A2>
    static TMV_INLINE void Swap(
        SmallBandMatrixView<T,M,N,LO,HI,Si1,Sj1,A1> m1,
        SmallBandMatrixView<T,M,N,LO,HI,Si2,Sj2,A2> m2)
    { DoSwap(m1,m2); }
    template <class T, int M, int N, int LO, int HI, int Si1, int Sj1, int A1, int A2>
    static TMV_INLINE void Swap(
        SmallBandMatrixView<T,M,N,LO,HI,Si1,Sj1,A1> m1,
        BandMatrixView<T,A2> m2)
    { DoSwap(m1,m2); }
    template <class T, int M, int N, int LO, int HI, int A1, int Si2, int Sj2, int A2>
    static TMV_INLINE void Swap(
        BandMatrixView<T,A1> m1,
        SmallBandMatrixView<T,M,N,LO,HI,Si2,Sj2,A2> m2)
    { DoSwap(m1,m2); }


    //
    // Conjugate, Transpose, Adjoint
    //

    template <class T, int M, int N, int LO, int HI, int A0, int A1>
    static TMV_INLINE typename BandMatrix<T,M,N,LO,HI,A0,A1>::conjugate_type Conjugate(
        BandMatrix<T,M,N,LO,HI,A0,A1>& m)
    { return m.conjugate(); }
    template <class T, int M, int N, int LO, int HI, int Si, int Sj, int A>
    static TMV_INLINE typename BandMatrixView<T,M,N,LO,HI,Si,Sj,A>::conjugate_type Conjugate(
        BandMatrixView<T,M,N,LO,HI,Si,Sj,A> m)
    { return m.conjugate(); }

    template <class T, int M, int N, int LO, int HI, int A0, int A1>
    static TMV_INLINE typename BandMatrix<T,M,N,LO,HI,A0,A1>::transpose_type Transpose(
        BandMatrix<T,M,N,LO,HI,A0,A1>& m)
    { return m.transpose(); }
    template <class T, int M, int N, int LO, int HI, int Si, int Sj, int A>
    static TMV_INLINE typename BandMatrixView<T,M,N,LO,HI,Si,Sj,A>::transpose_type Transpose(
        BandMatrixView<T,M,N,LO,HI,Si,Sj,A> m)
    { return m.transpose(); }

    template <class T, int M, int N, int LO, int HI, int A0, int A1>
    static TMV_INLINE typename BandMatrix<T,M,N,LO,HI,A0,A1>::adjoint_type Adjoint(
        BandMatrix<T,M,N,LO,HI,A0,A1>& m)
    { return m.adjoint(); }
    template <class T, int M, int N, int LO, int HI, int Si, int Sj, int A>
    static TMV_INLINE typename BandMatrixView<T,M,N,LO,HI,Si,Sj,A>::adjoint_type Adjoint(
        BandMatrixView<T,M,N,LO,HI,Si,Sj,A> m)
    { return m.adjoint(); }


    //
    // TMV_Text 
    //

#ifdef TMV_TEXT
    template <class T, int M, int N, int LO, int HI, int A0, int A1>
    static inline std::string TMV_Text(
        const SmallBandMatrix<T,M,N,LO,HI,A0,A1>& m)
    {
        const int A = A0 | A1;
        std::ostringstream s;
        s << "SmallBandMatrix<"<<TMV_Text(T());
        s << ','<<M<<','<<N<<','<<LO<<','<<HI;
        s << ","<<Attrib<A>::text()<<">";
        s << "("<<m.colsize()<<","<<m.rowsize()<<",";
        s << m.nlo()<<","<<m.nhi()<<",";
        s << m.stepi()<<","<<m.stepj()<<")";
        return s.str();
    }

    template <class T, int LO, int HI, int A0, int A1>
    static inline std::string TMV_Text(
        const ThinBandMatrix<T,LO,HI,A0,A1>& m)
    {
        const int A = A0 | A1;
        std::ostringstream s;
        s << "ThinBandMatrix<"<<TMV_Text(T());
        s << ','<<LO<<','<<HI;
        s << ","<<Attrib<A>::text()<<">";
        s << "("<<m.colsize()<<","<<m.rowsize()<<",";
        s << m.nlo()<<","<<m.nhi()<<",";
        s << m.stepi()<<","<<m.stepj()<<")";
        return s.str();
    }

    template <class T, int M, int N, int LO, int HI, int Si, int Sj, int A>
    static inline std::string TMV_Text(
        const ConstBandMatrixView<T,M,N,LO,HI,Si,Sj,A>& m)
    {
        std::ostringstream s;
        s << "ConstSmallBandMatrixView<"<<TMV_Text(T());
        s << ','<<IntTraits<M>::text()<<','<<IntTraits<N>::text();
        s << ','<<IntTraits<LO>::text()<<','<<IntTraits<HI>::text();
        s << ','<<IntTraits<Si>::text()<<','<<IntTraits<Sj>::text();
        s << ","<<Attrib<A>::text()<<">";
        s << "("<<m.colsize()<<","<<m.rowsize()<<",";
        s << m.nlo()<<","<<m.nhi()<<",";
        s << m.stepi()<<","<<m.stepj()<<")";
        return s.str();
    }

    template <class T, int M, int N, int LO, int HI, int Si, int Sj, int A>
    static inline std::string TMV_Text(
        const SmallBandMatrixView<T,M,N,LO,HI,Si,Sj,A>& m)
    {
        std::ostringstream s;
        s << "SmallBandMatrixView<"<<TMV_Text(T());
        s << ','<<IntTraits<M>::text()<<','<<IntTraits<N>::text();
        s << ','<<IntTraits<LO>::text()<<','<<IntTraits<HI>::text();
        s << ','<<IntTraits<Si>::text()<<','<<IntTraits<Sj>::text();
        s << ","<<Attrib<A>::text()<<">";
        s << "("<<m.colsize()<<","<<m.rowsize()<<",";
        s << m.nlo()<<","<<m.nhi()<<",";
        s << m.stepi()<<","<<m.stepj()<<")";
        return s.str();
    }
#endif

} // namespace tmv

#endif
