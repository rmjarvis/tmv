

#ifndef TMV_SVDecompose_Bidiag_H
#define TMV_SVDecompose_Bidiag_H

#include "TMV_SVDecompose.h"
#include "TMV_Householder.h"

#ifdef PRINTALGO_SVD
#include <iostream>
#include "TMV_MatrixIO.h"
#include "TMV_VectorIO.h"
#define TMV_USE_WRITER
#endif

#ifdef XDEBUG_SVD
#include "TMV_UnpackQ.h"
#include "TMV_Norm.h"
#endif

#define dbgcout SafeWriter()
#include "TMV_SafeWriter.h"

// BLOCKSIZE is the block size to use in algo 21, etc.
#define TMV_BIDIAG_BLOCKSIZE 16

namespace tmv {

    // Defined in TMV_SVDecompose_Bidiag.cpp
    template <class T, class RT>
    void InstBidiagonalize(
        MatrixView<T> A, VectorView<RT> Ubeta, VectorView<RT> Vbeta,
        VectorView<T> D, VectorView<T> E);

    // Decompose A into U B V
    // The Bidiagonal Matrix B is stored as two vectors: D, E
    // D is the diagonal, E is the super-diagonal
    // A along with Ubeta and Vbeta hold the U and V matrices.
    template <int algo, int cs, int rs, class M1, class V1, class V2, class V3, class V4>
    struct Bidiagonalize_Helper;

    // algo 0: Trivial, nothing to do (M == 0, or N == 0)
    template <int cs, int rs, class M1, class V1, class V2, class V3, class V4>
    struct Bidiagonalize_Helper<0,cs,rs,M1,V1,V2,V3,V4>
    { static TMV_INLINE void call(M1& , V1& , V2& , V3&, V4& ) {} };

    // algo 11: Non-blocked algorithm
    // We use Householder reflections to reduce A to the bidiagonal form:
    template <int cs, int rs, class M1, class V1, class V2, class V3, class V4>
    struct Bidiagonalize_Helper<11,cs,rs,M1,V1,V2,V3,V4>
    {
        static void call(M1& A, V1& Ubeta, V2& Vbeta, V3& D, V4& E)
        {
            typedef typename M1::value_type T;

            const int M = cs==TMV_UNKNOWN ? int(A.colsize()) : cs;
            const int N = rs==TMV_UNKNOWN ? int(A.rowsize()) : rs;
            if (N == 0) return;
#ifdef PRINTALGO_SVD
            std::cout<<"Bidiagonalize algo 11: M,N,cs,rs = "<<
                M<<','<<N<<','<<cs<<','<<rs<<std::endl;
#endif
            typedef typename V1::iterator IT1;
            typedef typename V2::iterator IT2;
            typedef typename V3::iterator IT3;
            typedef typename V4::iterator IT4;
            typedef typename M1::col_sub_type M1c;
            typedef typename M1::row_sub_type M1r;
            typedef typename M1::submatrix_type M1s;
            typedef typename M1s::transpose_type M1st;
            typedef typename VCopyHelper<T,cs>::type Vx;
            typedef typename Vx::subvector_type Vxs;
            Vx tempBase = VectorSizer<T>(M);

            IT1 Ubj = Ubeta.begin();
            IT2 Vbj = Vbeta.begin();
            for(int j=0;j<N-1;++j,++Ubj,++Vbj) {
                M1c u1 = A.col(j,j+1,M);
                M1r A10 = A.row(j,j+1,N);
                M1s A1x = A.subMatrix(j+1,M,j+1,N);
                M1r u2 = A.row(j,j+2,N);
                M1c A20 = A.col(j+1,j+1,M);
                M1st A2x = A.subMatrix(j+1,M,j+2,N).transpose();

                HouseholderReflect(A.ref(j,j),u1,*Ubj);
                Vxs temp1 = tempBase.subVector(0,N-j-1);
                HouseholderMultEq(u1,*Ubj,A10,A1x,temp1);

                HouseholderReflect(A.ref(j,j+1),u2,*Vbj);
                Vxs temp2 = tempBase.subVector(0,M-j-1);
                HouseholderMultEq(u2,*Vbj,A20,A2x,temp2);
            }
            M1c u1 = A.col(N-1,N,M);
            HouseholderReflect(A.ref(N-1,N-1),u1,*Ubj);

            // The bidiagonal of A is the bidiagonal we want, so copy it to D,E
            D = A.diag();
            E = A.diag(1);
        }
    };

    // algo 21: Blocked algorithm
    template <int cs, int rs, class M1, class V1, class V2, class V3, class V4>
    struct Bidiagonalize_Helper<21,cs,rs,M1,V1,V2,V3,V4>
    {
        static void call(M1& A, V1& Ubeta, V2& Vbeta, V3& D, V4& E)
        {
            typedef typename M1::value_type T;
            typedef typename M1::real_type RT;

            const int M = cs==TMV_UNKNOWN ? int(A.colsize()) : cs;
            const int N = rs==TMV_UNKNOWN ? int(A.rowsize()) : rs;
#ifdef PRINTALGO_SVD
            std::cout<<"Bidiagonalize algo 21: M,N,cs,rs = "<<
                M<<','<<N<<','<<cs<<','<<rs<<std::endl;
#endif
            // Normally we keep the Z matrix for block Householder matrices 
            // where the block Householder is I - YZYt 
            // (and Z is upper triangular).
            //
            // However, since the bidiagonalizing process proceeds along both 
            // the rows and columns, we have a problem.  If we just kept the Z
            // matrix for each set, then we would not be able to update the 
            // resulting submatrix at the end of the block.
            // m' = (I-YZYt) m (I-XtWX)
            //
            // For the left multiply, the m needs to have full height for the 
            // Yt m product.  Likewise, for the right multiply, it needs full 
            // width.
            // Since we update the first K rows and columns with the Y matrix, 
            // this doesn't work.  So instead of keeping Z,W we are forced to 
            // use a bit more temporary storage and store the products 
            // ZYtm and mXtW.
            //
            // Furthermore, the m in these products is maintained such that the
            // it already has the appropriate multiplies from the other side.
            // Then, when we are done with the block, the update becomes just:
            //
            // m' = m' - Y (ZYtm) - (mXtW) X
            //

            typedef typename V1::iterator IT1;
            typedef typename V2::iterator IT2;
            typedef typename M1::col_sub_type M1c;
            typedef typename M1c::conjugate_type M1cc;
            typedef typename M1::row_sub_type M1r;
            typedef typename M1r::conjugate_type M1rc;
            typedef typename M1::submatrix_type M1s;
            typedef typename M1s::transpose_type M1st;
            typedef typename VCopyHelper<T,cs>::type Vx;
            typedef typename Vx::subvector_type Vxs;
            Vx tempBase = VectorSizer<T>(M);

            const int NB = TMV_BIDIAG_BLOCKSIZE;
            typedef typename MCopyHelper<T,Rec,NB,rs,true>::type M3;
            typedef typename MCopyHelper<T,Rec,cs,NB>::type M4;

            typedef typename M3::row_sub_type M3r;
            typedef typename M3::submatrix_type M3s;
            typedef typename M4::col_sub_type M4c;
            typedef typename M4::submatrix_type M4s;

            IT1 Ubj = Ubeta.begin();
            IT2 Vbj = Vbeta.begin();
            int j1 = 0;
            RT bu, bv;

            M3 ZYtm = MatrixSizer<T>(NB,N);
            M4 mXtW = MatrixSizer<T>(M,NB);

            for(int j2=j1+NB;j2<N-1;j1=j2,j2+=NB) {
                for(int j=j1,jj=0;j<j2;++j,++jj,++Ubj,++Vbj) {
                    // jj = j-j1

                    // Update current column:
                    // A(j:M,j) -= Y(j:M,0:j) ZYtm(0:j,j) 
                    //            + mXtW(j:M,0:j) X(0:j,j)

                    M1c u = A.col(j,j,M);
                    M1c u1 = A.col(j,j+1,M);
                    M1s Y0 = A.subMatrix(j,M,j1,j);
                    M3s ZYtm0 = ZYtm.subMatrix(0,jj,j,N);
                    M4s mXtW0 = mXtW.subMatrix(j,M,0,jj);
                    M1s X0 = A.subMatrix(j1,j,j,N);
                    if (jj > 0) {
                        u -= Y0 * ZYtm0.col(0);
                        u -= mXtW0 * X0.col(0);
                    }

                    // Do the Householder reflection for U
                    // Copy the reflection into D(j), and set the top of 
                    // the Householder vector to be explicitly 1.
                    // (It makes life easier if it's actually 1 rather 
                    // than dealing with it implicitly.)
                    HouseholderReflect(u.ref(0),u1,bu);
                    *Ubj = bu;
                    D.ref(j) = u.cref(0);
                    u.ref(0) = T(1);

                    // Update ZYtm:
                    //
                    // ZYtm(j,j+1:N) = 
                    //     Z(j,0:j+1) Yt(0:j+1,0:M) m'(0:M,j+1:N)
                    // m' is the matrix after taking into account the 
                    // Householder multiplies that we have already done:
                    //
                    // m' = m' - Y ZYtm - mXtW X
                    //
                    // The new ZYtm(j,j+1:N) = (Z Yt m')(j,j+1:N)
                    // = Z(j,0:j+1) Yt(0:j+1,0:M) m'(0:M,j+1:N)
                    //
                    // Z is Upper Triangular, so Z(j,0:j+1) = 0
                    // Also Z(j,j) = bu, so:
                    // ZYtm(j,j+1:N) = bu Yt(j,0:M) m'(0:M,j+1:N)
                    //
                    // Y is Lower Unit Trapezoidal, so Yt(j,0:j) = 0
                    // The rest, Yt(j,j:M) is just ut, so:
                    // ZYtm(j,j+1:N) = bu ut m'(j:M,j+1:N)
                    //
                    // Finally, expand out the m':
                    // m'(j:M,j+1:N) = A(j:M,j+1:N) 
                    //                 - Y(j:M,0:j) ZYtm(0:j,j+1:N)
                    //                 - mXtW(j:M,0:j) X(0:j,j+1:N)
                    //
                    M3r ZYtmj = ZYtm.row(jj,j+1,N);
                    M3r temp = ZYtm.row(jj,j1,j);
                    M1cc ut = u.conjugate();
                    M3s ZYtm0a = ZYtm0.colRange(1,ZYtm0.rowsize());
                    M1s X0a = X0.colRange(1,X0.rowsize());
                    ZYtmj = ut * A.subMatrix(j,M,j+1,N);
                    if (jj > 0) {
                        ZYtmj -= (temp = ut*Y0) * ZYtm0a;
                        ZYtmj -= (temp = ut*mXtW0) * X0a;
                    }
                    ZYtmj *= bu;

                    // Update the current row:
                    // A(j,j+1:N) -= Y(j,0:j+1) ZYtm(0:j+1,j+1:N) 
                    //              + mXtW(j,0:j) X(0:j:j+1,N)
                    //
                    M1s Y1 = A.subMatrix(j,M,j1,j+1);
                    M3s ZYtm1 = ZYtm.subMatrix(0,jj+1,j+1,N);
                    M1r v = A.row(j,j+1,N);
                    M1r v1 = A.row(j,j+2,N);
                    v -= Y1.row(0) * ZYtm1;
                    if (jj > 0) {
                        v -= mXtW0.row(0) * X0a;
                    }

                    // Do the Householder reflection for V
                    //
                    HouseholderReflect(v.ref(0),v1,bv);
                    *Vbj = bv;
                    E.ref(j) = v.cref(0);
                    v.ref(0) = T(1);

                    // Update mXtW:
                    //
                    // mXtW(j+1:M,j) = 
                    //       m'(j+1:M,0:N) Xt(0:N,0:j+1) W(0:j+1,j)
                    // = bv m'(j+1:M,j+1:N) vt
                    //
                    // And m' is:
                    //
                    // m'(j+1:M,j+1:N) = A(j+1:M,j+1:N) 
                    //                   - Y(j+1:M,0:j+1) ZYtm(0:j+1,j+1:N) 
                    //                   - mXtW(j+1:M,0:j) X(0:j,j+1:N)
                    //
                    M4c mXtWj = mXtW.col(jj,j+1,M);
                    M4c temp1 = mXtW.col(jj,j1,j+1);
                    M4c temp2 = mXtW.col(jj,j1,j);
                    M1rc vt = v.conjugate();
                    M1s Y1a = Y1.rowRange(1,Y1.colsize());
                    M4s mXtW0a = mXtW0.rowRange(1,mXtW0.colsize());
                    mXtWj = A.subMatrix(j+1,M,j+1,N)*vt;
                    mXtWj -= Y1a * (temp1 = ZYtm1*vt);
                    mXtWj -= mXtW0a * (temp2 = X0a*vt);
                    mXtWj *= bv;
                }

                // Update the rest of the matrix:
                A.subMatrix(j2,M,j2,N) -= A.subMatrix(j2,M,j1,j2) *
                    ZYtm.colRange(j2,N);
                A.subMatrix(j2,M,j2,N) -= mXtW.rowRange(j2,M) *
                    A.subMatrix(j1,j2,j2,N);
            }

            // Do the ones that don't fit into blocks.
            for(int j=j1;j<N-1;++j,++Ubj,++Vbj) {
                M1c u1 = A.col(j,j+1,M);
                M1r A10 = A.row(j,j+1,N);
                M1s A1x = A.subMatrix(j+1,M,j+1,N);
                M1r u2 = A.row(j,j+2,N);
                M1c A20 = A.col(j+1,j+1,M);
                M1st A2x = A.subMatrix(j+1,M,j+2,N).transpose();

                HouseholderReflect(A.ref(j,j),u1,*Ubj);
                Vxs temp1 = tempBase.subVector(0,N-j-1);
                HouseholderMultEq(u1,*Ubj,A10,A1x,temp1);

                HouseholderReflect(A.ref(j,j+1),u2,*Vbj);
                Vxs temp2 = tempBase.subVector(0,M-j-1);
                HouseholderMultEq(u2,*Vbj,A20,A2x,temp2);
            }

            // Do the last U Householder vector:
            M1c u1 = A.col(N-1,N,M);
            HouseholderReflect(A.ref(N-1,N-1),u1,*Ubj);

            D.subVector(j1,N) = A.diag(0,j1,N);
            E.subVector(j1,N-1) = A.diag(1,j1,N-1);
        }
    };

    // algo 31: Decide which algorithm to use based on runtime size
    template <int cs, int rs, class M1, class V1, class V2, class V3, class V4>
    struct Bidiagonalize_Helper<31,cs,rs,M1,V1,V2,V3,V4>
    {
        static void call(M1& A, V1& Ubeta, V2& Vbeta, V3& D, V4& E)
        {
            typedef typename M1::value_type T;

            const int M = cs==TMV_UNKNOWN ? int(A.colsize()) : cs;
            const int N = rs==TMV_UNKNOWN ? int(A.rowsize()) : rs;
#ifdef PRINTALGO_SVD
            std::cout<<"Bidiagonalize algo 31: M,N,cs,rs = "<<
                M<<','<<N<<','<<cs<<','<<rs<<std::endl;
#endif
            const int l2cache = TMV_L2_CACHE*1024/sizeof(T);
            if (M*N <= l2cache)
                Bidiagonalize_Helper<11,cs,rs,M1,V1,V2,V3,V4>::call(
                    A,Ubeta,Vbeta,D,E);
            else
                Bidiagonalize_Helper<21,cs,rs,M1,V1,V2,V3,V4>::call(
                    A,Ubeta,Vbeta,D,E);
        }
    };

    // algo 90: call InstSV_DecomposeFromBidiagonal
    template <int cs, int rs, class M1, class V1, class V2, class V3, class V4>
    struct Bidiagonalize_Helper<90,cs,rs,M1,V1,V2,V3,V4>
    {
        static TMV_INLINE void call(M1& A, V1& Ubeta, V2& Vbeta, V3& D, V4& E)
        { 
            InstBidiagonalize(
                A.xView(),Ubeta.xView(),Vbeta.xView(),D.xView(),E.xView()); 
        }
    };

    // algo -4: No copies or branches
    template <int cs, int rs, class M1, class V1, class V2, class V3, class V4>
    struct Bidiagonalize_Helper<-4,cs,rs,M1,V1,V2,V3,V4>
    {
        static TMV_INLINE void call(M1& A, V1& Ubeta, V2& Vbeta, V3& D, V4& E)
        {
#if 0
            const int algo = 21;
#else
            typedef typename M1::value_type T;
            const int csrs = IntTraits2<cs,rs>::prod;
            const int rsrs = IntTraits2<rs,rs>::prod;
            const int l2cache = TMV_L2_CACHE*1024/sizeof(T);
            const int algo = 
                cs == 0 || rs == 0 ? 0 :
                TMV_OPT == 0 ? 11 :
                rs == TMV_UNKNOWN ? 31 :
                rsrs > l2cache ? 21 :
                cs == TMV_UNKNOWN ? 31 : 
                csrs <= l2cache ? 11 :
                21;
#endif
#ifdef PRINTALGO_SVD
            std::cout<<"Inline Bidiagonalize: \n";
            std::cout<<"A = "<<TMV_Text(A)<<std::endl;
            std::cout<<"Ubeta = "<<TMV_Text(Ubeta)<<std::endl;
            std::cout<<"Vbeta = "<<TMV_Text(Vbeta)<<std::endl;
            std::cout<<"D = "<<TMV_Text(D)<<std::endl;
            std::cout<<"E = "<<TMV_Text(E)<<std::endl;
            std::cout<<"cs = "<<cs<<"  rs = "<<rs<<std::endl;
            std::cout<<"sizes = "<<A.colsize()<<"  "<<A.rowsize()<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
            //std::cout<<"m = "<<U<<std::endl;
#endif
#ifdef XDEBUG_SVD
            const int N = rs==TMV_UNKNOWN ? int(D.size()) : rs;
            dbgcout<<"Start Bidiagonalize:\n";
            dbgcout<<"A = "<<TMV_Text(A)<<"   "<<A<<std::endl;
            typedef typename M1::value_type T;
            Matrix<T> A0(A);
#endif
            Bidiagonalize_Helper<algo,cs,rs,M1,V1,V2,V3,V4>::call(
                A,Ubeta,Vbeta,D,E);
            //std::cout<<"U = "<<U<<std::endl;
            //std::cout<<"S = "<<S<<std::endl;
            //std::cout<<"V = "<<V<<std::endl;
#ifdef XDEBUG_SVD
            Matrix<T> U(A);
            UnpackQ(U,Ubeta);
            Matrix<T> V(N,N);
            V.setToIdentity();
            V.subMatrix(1,N,1,N) = A.subMatrix(0,N-1,1,N);
            MatrixView<T,RowMajor> V1t = V.subMatrix(1,N,1,N).transpose();
            UnpackQ(V1t,Vbeta);
            Matrix<T> B(N,N,T(0));
            B.diag() = D;
            B.diag(1) = E;
            Matrix<T> AA = U*B*V;
            if (!(Norm(A0-AA) < 0.001*Norm(A0))) {
                std::cerr<<"Bidiagonalize: A = "<<TMV_Text(A)<<"  "<<A0<<std::endl;
                std::cerr<<"A = "<<A<<std::endl;
                std::cerr<<"Ubeta = "<<Ubeta<<std::endl;
                std::cerr<<"Vbeta = "<<Vbeta<<std::endl;
                std::cerr<<"U = "<<U<<std::endl;
                std::cerr<<"B = "<<B<<std::endl;
                std::cerr<<"V = "<<V<<std::endl;
                std::cerr<<"UBV = "<<AA<<std::endl;
                std::cerr<<"Norm(diff) = "<<Norm(A0-AA)<<std::endl;
                std::cerr<<"cf. "<<0.001*Norm(A0)<<std::endl;
                abort();
            }
#endif
        }
    };

    // algo -3: Determine which algorithm to use
    template <int cs, int rs, class M1, class V1, class V2, class V3, class V4>
    struct Bidiagonalize_Helper<-3,cs,rs,M1,V1,V2,V3,V4>
    {
        static TMV_INLINE void call(M1& A, V1& Ubeta, V2& Vbeta, V3& D, V4& E)
        {
            const int algo = -4;
            Bidiagonalize_Helper<algo,cs,rs,M1,V1,V2,V3,V4>::call(
                A,Ubeta,Vbeta,D,E);
        }
    };

    // algo -2: Check for inst
    template <int cs, int rs, class M1, class V1, class V2, class V3, class V4>
    struct Bidiagonalize_Helper<-2,cs,rs,M1,V1,V2,V3,V4>
    {
        static TMV_INLINE void call(M1& A, V1& Ubeta, V2& Vbeta, V3& D, V4& E)
        {
            typedef typename M1::value_type T;
            const bool inst = 
                (cs == TMV_UNKNOWN || cs > 16) &&
                (rs == TMV_UNKNOWN || rs > 16) &&
                Traits<T>::isinst;
            const int algo = 
                inst ? 90 :
                -3;
            Bidiagonalize_Helper<algo,cs,rs,M1,V1,V2,V3,V4>::call(
                A,Ubeta,Vbeta,D,E);
        }
    };

    template <int cs, int rs, class M1, class V1, class V2, class V3, class V4>
    struct Bidiagonalize_Helper<-1,cs,rs,M1,V1,V2,V3,V4>
    {
        static TMV_INLINE void call(M1& A, V1& Ubeta, V2& Vbeta, V3& D, V4& E)
        {
            Bidiagonalize_Helper<-2,cs,rs,M1,V1,V2,V3,V4>::call(
                A,Ubeta,Vbeta,D,E);
        }
    };

    template <class M1, class V1, class V2, class V3, class V4>
    static inline void InlineBidiagonalize(
        BaseMatrix_Rec_Mutable<M1>& A,
        BaseVector_Mutable<V1>& Ubeta, BaseVector_Mutable<V2>& Vbeta, 
        BaseVector_Mutable<V3>& D, BaseVector_Mutable<V4>& E)
    {
        TMVStaticAssert((Traits2<
                         typename M1::value_type,
                         typename V1::value_type>::samebase));
        TMVStaticAssert(Traits<typename V1::value_type>::isreal);
        TMVStaticAssert((Traits2<
                         typename V1::value_type,
                         typename V2::value_type>::sametype));
        TMVStaticAssert((Traits2<
                         typename M1::value_type,
                         typename V3::value_type>::sametype));
        TMVStaticAssert((Traits2<
                         typename M1::value_type,
                         typename V4::value_type>::sametype));
        TMVStaticAssert((Sizes<M1::_rowsize,V1::_size>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,V3::_size>::same));
        TMVStaticAssert((Sizes<V3::_size,IntTraits<V4::_size>::Sp1>::same));
        TMVStaticAssert((Sizes<V2::_size,V4::_size>::same));
        TMVAssert(A.colsize() >= A.rowsize());
        TMVAssert(A.rowsize() == D.size());
        TMVAssert(A.rowsize() == Ubeta.size());
        TMVAssert(A.rowsize() == E.size()+1);
        TMVAssert(A.rowsize() == Vbeta.size()+1);
        const int cs = M1::_colsize;
        const int rs1 = Sizes<M1::_rowsize,V1::_size>::size;
        const int rs = Sizes<rs1,V3::_size>::size;
        typedef typename M1::cview_type M1v;
        typedef typename V1::cview_type V1v;
        typedef typename V2::cview_type V2v;
        typedef typename V3::cview_type V3v;
        typedef typename V4::cview_type V4v;
        TMV_MAYBE_REF(M1,M1v) Av = A.cView();
        TMV_MAYBE_REF(V1,V1v) Ubv = Ubeta.cView();
        TMV_MAYBE_REF(V2,V2v) Vbv = Vbeta.cView();
        TMV_MAYBE_REF(V3,V3v) Dv = D.cView();
        TMV_MAYBE_REF(V4,V4v) Ev = E.cView();
        Bidiagonalize_Helper<-3,cs,rs,M1v,V1v,V2v,V3v,V4v>::call(
            Av,Ubv,Vbv,Dv,Ev);
    }

    template <class M1, class V1, class V2, class V3, class V4>
    static inline void Bidiagonalize(
        BaseMatrix_Rec_Mutable<M1>& A,
        BaseVector_Mutable<V1>& Ubeta, BaseVector_Mutable<V2>& Vbeta, 
        BaseVector_Mutable<V3>& D, BaseVector_Mutable<V4>& E)
    {
        TMVStaticAssert((Traits2<
                         typename M1::value_type,
                         typename V1::value_type>::samebase));
        TMVStaticAssert(Traits<typename V1::value_type>::isreal);
        TMVStaticAssert((Traits2<
                         typename V1::value_type,
                         typename V2::value_type>::sametype));
        TMVStaticAssert((Traits2<
                         typename M1::value_type,
                         typename V3::value_type>::sametype));
        TMVStaticAssert((Traits2<
                         typename M1::value_type,
                         typename V4::value_type>::sametype));
        TMVStaticAssert((Sizes<M1::_rowsize,V1::_size>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,V3::_size>::same));
        TMVStaticAssert((Sizes<V3::_size,IntTraits<V4::_size>::Sp1>::same));
        TMVStaticAssert((Sizes<V2::_size,V4::_size>::same));
        TMVAssert(A.colsize() >= A.rowsize());
        TMVAssert(A.rowsize() == D.size());
        TMVAssert(A.rowsize() == Ubeta.size());
        TMVAssert(A.rowsize() == E.size()+1);
        TMVAssert(A.rowsize() == Vbeta.size()+1);
        const int cs = M1::_colsize;
        const int rs1 = Sizes<M1::_rowsize,V1::_size>::size;
        const int rs = Sizes<rs1,V3::_size>::size;
        typedef typename M1::cview_type M1v;
        typedef typename V1::cview_type V1v;
        typedef typename V2::cview_type V2v;
        typedef typename V3::cview_type V3v;
        typedef typename V4::cview_type V4v;
        TMV_MAYBE_REF(M1,M1v) Av = A.cView();
        TMV_MAYBE_REF(V1,V1v) Ubv = Ubeta.cView();
        TMV_MAYBE_REF(V2,V2v) Vbv = Vbeta.cView();
        TMV_MAYBE_REF(V3,V3v) Dv = D.cView();
        TMV_MAYBE_REF(V4,V4v) Ev = E.cView();
        Bidiagonalize_Helper<-2,cs,rs,M1v,V1v,V2v,V3v,V4v>::call(
            Av,Ubv,Vbv,Dv,Ev);
    }


} // namespace tmv

#endif

