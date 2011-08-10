

#ifndef TMV_SVDecompose_H
#define TMV_SVDecompose_H

#include "TMV_BaseMatrix_Rec.h"
#include "TMV_BaseMatrix_Diag.h"
#include "TMV_Det.h"
#include "TMV_Givens.h"
#include "TMV_Permutation.h"

#ifdef XDEBUG_SVD
#define THRESH 1.e-4
#include "TMV_SmallDiagMatrix.h"
#include "TMV_ProdMM.h"
#include "TMV_MultMM.h"
#include "TMV_MultMD.h"
#include "TMV_AddMM.h"
#include "TMV_MultXM.h"
#endif

#ifdef PRINTALGO_SVD
#include <iostream>
#include "TMV_MatrixIO.h"
#include "TMV_VectorIO.h"
#define TMV_USE_WRITER
#endif

#define dbgcout SafeWriter()
#include "TMV_SafeWriter.h"

namespace tmv {

    //
    // First some helper functions that we will need.
    // The first few are substantial enough that they are 
    // in other files.
    //
    
    // Defined in TMV_SVDecompose_Bidiag.h
    template <class M, class V1, class V2, class V3, class V4>
    inline void Bidiagonalize(
        BaseMatrix_Rec_Mutable<M>& A,
        BaseVector_Mutable<V1>& Ubeta, BaseVector_Mutable<V2>& Vbeta, 
        BaseVector_Mutable<V3>& D, BaseVector_Mutable<V4>& E);

    // Defined in TMV_SVDecompose_QR.h
    template <class M1, class V1, class V2, class M2>
    inline void SV_DecomposeFromBidiagonal_QR(
        BaseMatrix_Rec_Mutable<M1>& U,
        BaseVector_Mutable<V1>& D, BaseVector_Mutable<V2>& E,
        BaseMatrix_Rec_Mutable<M2>& V, bool UisI, bool VisI);

    // Defined in TMV_SVDecompose_DC.h
    template <class M1, class V1, class V2, class M2>
    inline void SV_DecomposeFromBidiagonal_DC(
        BaseMatrix_Rec_Mutable<M1>& U,
        BaseVector_Mutable<V1>& D, BaseVector_Mutable<V2>& E,
        BaseMatrix_Rec_Mutable<M2>& V, bool UisI, bool VisI);

    
    // This routines sets to 0 any elements in D,E which
    // are essentially 0, given the machine precision.
    // if |D(i)|^2*Epsilon == 0, the D(i) = 0
    // if |E(i)| < Epsilon * (|D(i)| + |D(i+1)|), then E(i) = 0
    template <class V1, class V2>
    static inline void BidiagonalChopSmallElements(
        BaseVector_Mutable<V1>& D, BaseVector_Mutable<V2>& E, bool& zd)
    {
        TMVAssert(E.size() == D.size()-1);
        TMVAssert(D.step() == 1);
        TMVAssert(E.step() == 1);
        typedef typename V1::real_type RT;
        const RT eps = TMV_Epsilon<RT>();

        typename V1::iterator Di = D.begin();
        typename V1::iterator Ei = E.begin();

        if (TMV_Underflow(TMV_NORM(*Di))) {
            *Di = RT(0);
            zd = true;
        }
        ++Di;
        for(int k=E.size();k>0;--k,++Di,++Ei) {

            // if |D(i)|^2 underflows, the set D(i) = 0
            if (TMV_Underflow(TMV_NORM(*Di))) {
                *Di = RT(0);
                zd = true;
            }

            // if |E(i)| < Epsilon * (|D(i)| + |D(i+1)|), then E(i) = 0
            // Do it as !(|e| > eps (|d1| + |d2|)) so nan's will get set to 
            // zero too and not iterate forever.
            if ( !(TMV_ABS2(*Ei) > eps*(TMV_ABS2(*Di)+TMV_ABS2(*(Di-1))))
                 || TMV_Underflow(*Ei) ) {
                *Ei = RT(0);
            }

            // This last one checks that when the product of a D*E underflows
            // that at least one of the two values have been set to zero.
            // Otherwise, set the smaller one to 0.
            // e.g. if underflow happens at 1.e-310, and you have
            // Di = 1.e-154, Ei = 1.e-157, then none of the above tests
            // will trigger a zero.  But the product will underflow, leading
            // to problems with Reduce not reducing.
            if (TMV_Underflow(*Di * *Ei) && *Di!=RT(0) && *Ei!=RT(0)) {
                if (TMV_ABS2(*Ei) <= TMV_ABS2(*Di) ) *Ei = RT(0);
                else *Di = RT(0);
            }
            if (TMV_Underflow(*(Di-1) * *Ei) && *(Di-1)!=RT(0) && *Ei!=RT(0)) {
                if (TMV_ABS2(*Ei) <= TMV_ABS2(*(Di-1)) ) *Ei = RT(0);
                else *(Di-1) = RT(0);
            }
        }
    }

    // 
    // Now some auxilliary functions for zeroing out a column or row
    // in particular special cases.
    //

    template <class M1, class V1, class V2>
    static void BidiagonalZeroFirstRow(
        BaseMatrix_Rec_Mutable<M1>& U,
        BaseVector_Mutable<V1>& D, BaseVector_Mutable<V2>& E)
    {
        // Input D,E form an N+1 x N bidiagonal matrix
        // (eg. for N = 4)
        //     [ x 0 0 0 ]
        //     [ x x 0 0 ]
        // B = [ 0 x x 0 ]
        //     [ 0 0 x x ]
        //     [ 0 0 0 x ]
        // Zero out the first row maintaining the constancy of U B
        // using Givens transformations.
        TMVAssert(E.size() == D.size());
        if (U.ptr()) TMVAssert(U.rowsize() == D.size()+1);

        typedef typename V1::real_type RT;
        const int N = D.size();
        typedef typename V1::iterator Dit;
        typedef typename V2::iterator Eit;
        Dit Di = D.begin();
        Eit Ei = E.begin();

        RT x = *Ei;
        if (x != RT(0)) {
            *Ei = RT(0);
            ++Ei;
            // Loop Invariant: x = B(0,i)
            for(int i=0; i<N; ++i,++Di,++Ei) {
                Givens<RT> G = GivensRotate(*Di,x);
                // Make new B = G B
                if (i<N) {
                    G.mult(*Ei,x);
                }
                // Make new U = U Gt
                // Ut = G Ut
                if (U.ptr()) {
                    typename M1::colpair_type::adjoint_type Ucpt =
                        U.colPair(i+1,0).adjoint();
                    G.mult(Ucpt);
                }
            }
        }
    }

    template <class V1, class V2, class M2>
    static void BidiagonalZeroLastCol(
        BaseVector_Mutable<V1>& D, BaseVector_Mutable<V2>& E,
        BaseMatrix_Rec_Mutable<M2>& V)
    {
        // Input D,E form an N x N+1 bidiagonal matrix
        // (eg. for N = 4)
        //     [ x x 0 0 0 ]
        // B = [ 0 x x 0 0 ]
        //     [ 0 0 x x 0 ]
        //     [ 0 0 0 x x ]
        // Zero out the last col maintaining the constancy of B V
        // using Givens transformations.
        TMVAssert(E.size() == D.size());
        if (V.ptr()) TMVAssert(V.colsize() == D.size()+1);

        typedef typename V1::real_type RT;
        const int N = D.size();
        typedef typename V1::reverse_iterator Dit;
        typedef typename V2::reverse_iterator Eit;
        Dit Di = D.rbegin();
        Eit Ei = E.rbegin();

        RT x = *Ei;
        if (x != RT(0)) {
            *Ei = RT(0);
            // Loop Invariant: x = B(i,N-1)
            for(int i=N-1; i>=0; --i,++Di) {
                Givens<RT> G = GivensRotate(*Di,x);
                // Make new B = B GT
                if (i>0) {
                    G.mult(*(++Ei),x);
                }
                // Make new V = G* V 
                if (V.ptr()) {
                    typename M2::rowpair_type Vrp = V.rowPair(i,N);
                    G.conjMult(Vrp);
                }
            }
        }
    }

    //
    // Next the routine for solving the SVD after it has been
    // bidiagonalized.
    //
 
    // Defined in TMV_SVDecompose.cpp
    template <class Tu, class RT, class Tv>
    void InstSV_DecomposeFromBidiagonal(
        MatrixView<Tu> U, VectorView<RT> D, VectorView<RT> E,
        MatrixView<Tv> V, bool setUV);

    template <int algo, int cs, int rs, class Mu, class Vd, class Ve, class Mv>
    struct SVDecomposeFromBidiagonal_Helper;

    // algo 0: Trivial, nothing to do (M == 0, or N == 0)
    template <int cs, int rs, class Mu, class Vd, class Ve, class Mv>
    struct SVDecomposeFromBidiagonal_Helper<0,cs,rs,Mu,Vd,Ve,Mv>
    { static TMV_INLINE void call(Mu& , Vd& , Ve& , Mv& , bool ) {} };

    // algo 1: D,E are complex.  Make them real first.
    template <int cs, int rs, class Mu, class Vd, class Ve, class Mv>
    struct SVDecomposeFromBidiagonal_Helper<1,cs,rs,Mu,Vd,Ve,Mv>
    {
        static void call(Mu& U, Vd& D, Ve& E, Mv& V, bool setUV)
        {
            typedef typename Vd::value_type T;
            typedef typename Vd::real_type RT;
            TMVStaticAssert(Traits<T>::iscomplex);

            const int N = rs==TMV_UNKNOWN ? int(D.size()) : rs;
            if (N == 0) return;
#ifdef PRINTALGO_SVD
            const int M = cs==TMV_UNKNOWN ? int(U.colsize()) : cs;
            std::cout<<"SVDecomposeFromBidiagonal algo 1: M,N,cs,rs = "<<
                M<<','<<N<<','<<cs<<','<<rs<<std::endl;
#endif
            const int rsm1 = IntTraits<rs>::Sm1;
            typedef typename VCopyHelper<RT,rs>::type rVd;
            typedef typename VCopyHelper<RT,rsm1>::type rVe;
            rVd rD = VectorSizer<T>(N);
            rVe rE = VectorSizer<T>(N-1);

            if (setUV) {
                TMVAssert(U.ptr() && V.ptr());
                U.setToIdentity();
                V.setToIdentity();
            }

            for(int j=0;j<N-1;++j) {
                RT absDj = TMV_ABS(D.cref(j));
                T signDj = TMV_SIGN(D.cref(j),absDj);
                rD.ref(j) = absDj;
                E.ref(j) *= TMV_CONJ(signDj);
                if (U.ptr()) {
                    if (setUV) U.ref(j,j) = signDj;
                    else U.col(j) *= signDj;
                }

                RT absEj = TMV_ABS(E.cref(j));
                T signEj = TMV_SIGN(E.cref(j),absEj);
                rE.ref(j) = absEj;
                D.ref(j+1) *= TMV_CONJ(signEj);
                if (V.ptr()) {
                    if (setUV) V.ref(j+1,j+1) = signEj;
                    else V.row(j+1) *= signEj;
                }
            }
            RT absDj = TMV_ABS(D.cref(N-1));
            T signDj = TMV_SIGN(D.cref(N-1),absDj);
            rD.ref(N-1) = absDj;
            if (U.ptr()) {
                if (setUV) U.ref(N-1,N-1) = signDj;
                else U.col(N-1) *= signDj;
            }

            SVDecomposeFromBidiagonal_Helper<-2,cs,rs,Mu,rVd,rVe,Mv>::call(
                U,rD,rE,V,false);
            E.setZero();
            D = rD;
        }
    };

    // algo 11: Normal algorithm, D,E are real
    template <int cs, int rs, class Mu, class Vd, class Ve, class Mv>
    struct SVDecomposeFromBidiagonal_Helper<11,cs,rs,Mu,Vd,Ve,Mv>
    {
        static void call(Mu& U, Vd& D, Ve& E, Mv& V, bool setUV)
        {
            typedef typename Vd::value_type RT;
            TMVStaticAssert(Traits<RT>::isreal);

            const int N = rs==TMV_UNKNOWN ? int(D.size()) : rs;
            if (N == 0) return;
#ifdef PRINTALGO_SVD
            const int M = cs==TMV_UNKNOWN ? int(U.colsize()) : cs;
            std::cout<<"SVDecomposeFromBidiagonal algo 11: M,N,cs,rs = "<<
                M<<','<<N<<','<<cs<<','<<rs<<std::endl;
#endif
            if (setUV) {
                TMVAssert(U.ptr() && V.ptr());
                U.setToIdentity();
                V.setToIdentity();
            }

            // Before running the normal algorithms, rescale D,E by the maximum
            // value to help avoid overflow and underflow.
            RT scale = TMV_MAX(D.maxAbsElement(),E.maxAbsElement());
            dbgcout<<"scale = "<<scale<<std::endl;
            dbgcout<<"1/scale = "<<RT(1)/scale<<std::endl;
            if (TMV_Underflow(scale)) {
                // Hopeless case.  Just zero out D,E and call it done.
                D.setZero();
                E.setZero();
                return;
            }
            D /= scale;
            E /= scale;
            dbgcout<<"After scale: \n";
            dbgcout<<"D = "<<D<<std::endl;
            dbgcout<<"E = "<<E<<std::endl;

            dbgcout<<"Before Reduction: D = "<<D<<std::endl;
            dbgcout<<"U = "<<U<<std::endl;
            SVDecomposeFromBidiagonal_Helper<21,cs,rs,Mu,Vd,Ve,Mv>::call(
                U,D,E,V,setUV,setUV);
            dbgcout<<"After Reduction: D = "<<D<<std::endl;
            dbgcout<<"U = "<<U<<std::endl;

            // Make all of the singular values positive
            typename Vd::iterator Di = D.begin();
            for(int i=0;i<N;++i,++Di) if (*Di < 0) {
                *Di = -(*Di);
                if (V.ptr()) V.row(i) = -V.row(i);
            }
            dbgcout<<"After make positive: \n";
            dbgcout<<"D = "<<D<<std::endl;

            // Now A = U * S * V
            // Sort output singular values 
            Permutation P(N);
            D.sort(P,Descend);
            if (U.ptr()) U = U * P.transpose();
            if (V.ptr()) V = P * V;

            // Undo the scaling
            D *= scale;
            dbgcout<<"After undo scale: \n";
            dbgcout<<"D = "<<D<<std::endl;

        }
    };

    // algo 21: After basic prep work, this routine looks for sub-problems
    // that have no zeros in D or E. 
    template <int cs, int rs, class Mu, class Vd, class Ve, class Mv>
    struct SVDecomposeFromBidiagonal_Helper<21,cs,rs,Mu,Vd,Ve,Mv>
    {
        static void call(Mu& U, Vd& D, Ve& E, Mv& V, bool UisI, bool VisI)
        {
            typedef typename Vd::value_type RT;
            TMVStaticAssert(Traits<RT>::isreal);

            const int N = rs==TMV_UNKNOWN ? int(D.size()) : rs;
            if (N == 0) return;
#ifdef PRINTALGO_SVD
            const int M = cs==TMV_UNKNOWN ? int(U.colsize()) : cs;
            std::cout<<"SVDecomposeFromBidiagonal algo 21: M,N,cs,rs = "<<
                M<<','<<N<<','<<cs<<','<<rs<<std::endl;
#endif
            // First chop any small elements in D,E
            bool zd;
            BidiagonalChopSmallElements(D,E,zd);
            dbgcout<<"After Chop: D = "<<D<<std::endl;
            dbgcout<<"After Chop: E = "<<E<<std::endl;

            const int xx = TMV_UNKNOWN;
            typedef typename Mu::colrange_type Mus;
            typedef typename Vd::subvector_type Vds;
            typedef typename Ve::subvector_type Ves;
            typedef typename Mv::rowrange_type Mvs;

            // Find sub-problems to solve:
            for(int q = N-1; q>0; ) {
                dbgcout<<"Looking for sub-problem:\n";
                dbgcout<<"q = "<<q<<std::endl;
                if (E.cref(q-1) == RT(0)) --q;
                else if (D.cref(q) == RT(0)) {
                    dbgcout<<"D(q) == 0, so do ZeroLastCol\n";
                    // We have the end looking like:
                    //   ? ?
                    //     ? x
                    //       0
                    // So we need to find a p where all E(i) with p<=i<q are 
                    // non-zero.
                    int p = q-1;
                    while (p>0 && !(E.cref(p-1) == RT(0))) --p;
                    // Now Zero out the last column:
                    Vds D1 = D.subVector(p,q);
                    Ves E1 = E.subVector(p,q);
                    Mvs V1 = V.ptr() ? V.rowRange(p,q+1) : V.rowRange(0,0);
                    BidiagonalZeroLastCol(D1,E1,V1);
                    VisI = false;
                    --q;
                } else {
                    // Find first p before q with either E(p) = 0 or D(p) = 0
                    int p=q-1;
                    while (p>0 && 
                           !(E.cref(p-1)==RT(0)) && 
                           !(D.cref(p)==RT(0))) --p;
                    dbgcout<<"p = "<<p<<std::endl;
                    if (D.cref(p) == RT(0)) {
                        dbgcout<<"D(p) == 0, so do ZeroFirstRow\n";
                        // We have a block looking like:
                        //   0 x
                        //     x x 
                        //       x x
                        //         x
                        Mus U1 = U.ptr() ? U.colRange(p,q+1) : U.colRange(0,0);
                        Vds D1 = D.subVector(p+1,q+1);
                        Ves E1 = E.subVector(p,q);
                        BidiagonalZeroFirstRow(U1,D1,E1);
                        UisI = false;
                        ++p;
                    }
                    if (q > p) {
                        dbgcout<<"No zeros in D,E:\n";
                        dbgcout<<"D = "<<D.subVector(p,q+1)<<std::endl;
                        dbgcout<<"E = "<<E.subVector(p,q)<<std::endl;
                        Mus U1 = U.ptr() ? U.colRange(p,q+1) : U.colRange(0,0);
                        Vds D1 = D.subVector(p,q+1);
                        Ves E1 = E.subVector(p,q);
                        Mvs V1 = V.ptr() ? V.rowRange(p,q+1) : V.rowRange(0,0);

                        SVDecomposeFromBidiagonal_Helper<
                            31,xx,xx,Mus,Vds,Ves,Mvs>::call(
                                U1,D1,E1,V1,
                                U.ptr() && UisI && p==0 && q+1==N,
                                V.ptr() && VisI && p==0 && q+1==N);
                    }
                    q = p;
                }
            }
        }
    };

    // algo 31: Do a sub-problem with no zeros in D or E.
    template <int cs, int rs, class Mu, class Vd, class Ve, class Mv>
    struct SVDecomposeFromBidiagonal_Helper<31,cs,rs,Mu,Vd,Ve,Mv>
    {
        static TMV_INLINE void call(
            Mu& U, Vd& D, Ve& E, Mv& V, bool UisI, bool VisI)
        {
#ifdef PRINTALGO_SVD
            const int M = cs==TMV_UNKNOWN ? int(U.colsize()) : cs;
            const int N = rs==TMV_UNKNOWN ? int(D.size()) : rs;
            std::cout<<"SVDecomposeFromBidiagonal algo 31: M,N,cs,rs = "<<
                M<<','<<N<<','<<cs<<','<<rs<<std::endl;
#endif
#if 1
            SV_DecomposeFromBidiagonal_DC(U,D,E,V,UisI,VisI);
#else
            SV_DecomposeFromBidiagonal_QR(U,D,E,V,UisI,VisI);
#endif
        }
    };

    // algo 90: call InstSV_DecomposeFromBidiagonal
    template <int cs, int rs, class Mu, class Vd, class Ve, class Mv>
    struct SVDecomposeFromBidiagonal_Helper<90,cs,rs,Mu,Vd,Ve,Mv>
    {
        static TMV_INLINE void call(Mu& U, Vd& D, Ve& E, Mv& V, bool setUV)
        { 
            InstSV_DecomposeFromBidiagonal(
                U.xView(),D.xView(),E.xView(),V.xView(),setUV); 
        }
    };

    // algo -3: Determine which algorithm to use
    template <int cs, int rs, class Mu, class Vd, class Ve, class Mv>
    struct SVDecomposeFromBidiagonal_Helper<-3,cs,rs,Mu,Vd,Ve,Mv>
    {
        static TMV_INLINE void call(Mu& U, Vd& D, Ve& E, Mv& V, bool setUV)
        {
            const int algo = 11;
#ifdef PRINTALGO_SVD
            std::cout<<"Inline SVDecomposeFromBidiagonal: \n";
            std::cout<<"U = "<<TMV_Text(U)<<std::endl;
            std::cout<<"D = "<<TMV_Text(D)<<std::endl;
            std::cout<<"E = "<<TMV_Text(E)<<std::endl;
            std::cout<<"V = "<<TMV_Text(V)<<std::endl;
            std::cout<<"cs = "<<cs<<"  rs = "<<rs<<std::endl;
            std::cout<<"sizes = "<<U.colsize()<<"  "<<U.rowsize()<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
#endif
#ifdef XDEBUG_SVD
            const int N = rs==TMV_UNKNOWN ? int(D.size()) : rs;
            dbgcout<<"Start Decompose from Bidiag:\n";
            if (U.ptr()) dbgcout<<"U = "<<TMV_Text(U)<<std::endl;
            if (V.ptr()) dbgcout<<"V = "<<TMV_Text(V)<<std::endl;
            dbgcout<<"D = "<<TMV_Text(D)<<"  step "<<D.step()<<"  "<<D<<std::endl;
            dbgcout<<"E = "<<TMV_Text(E)<<"  step "<<E.step()<<"  "<<E<<std::endl;
            //if (U.ptr()) dbgcout<<"U = "<<U<<std::endl;
            //if (V.ptr()) dbgcout<<"V = "<<V<<std::endl;

            dbgcout<<"setUV = "<<setUV<<std::endl;
            typedef typename Traits2<typename Mu::value_type, typename Mv::value_type>::type T;
            typedef typename Mu::real_type RT;
            Matrix<RT> B(N,N,RT(0));
            B.diag() = D;
            B.diag(1) = E;
            const int M1 = U.ptr() && V.ptr() ? int(U.colsize()) : N;
            Matrix<T> A0(M1,N);
            if (U.ptr() && V.ptr() && !setUV) A0 = U * B * V;
            else A0 = B;
            //dbgcout<<"A0 = "<<A0<<std::endl;
#endif
            SVDecomposeFromBidiagonal_Helper<algo,cs,rs,Mu,Vd,Ve,Mv>::call( 
                U,D,E,V,setUV);
            //std::cout<<"U = "<<U<<std::endl;
            //std::cout<<"S = "<<S<<std::endl;
            //std::cout<<"V = "<<V<<std::endl;
#ifdef XDEBUG_SVD
            if (U.ptr() && V.ptr()) {
                Matrix<T> AA = U * DiagMatrixViewOf(D) * V;
                if (!(Norm(A0-AA) < THRESH*Norm(A0))) {
                    std::cerr<<"SV_DecomposeFromBidiagonal: \n";
                    std::cerr<<"input B = "<<B<<std::endl;
                    std::cerr<<"U => "<<U<<std::endl;
                    std::cerr<<"S => "<<D<<std::endl;
                    std::cerr<<"V => "<<V<<std::endl;
                    std::cerr<<"UBV = "<<A0<<std::endl;
                    std::cerr<<"USV = "<<AA<<std::endl;
                    std::cerr<<"diff = ";
                    (A0-AA).write(std::cerr,(A0-AA).maxAbsElement()*1.e-3);
                    std::cerr<<std::endl;
                    std::cerr<<"Norm(diff) = "<<Norm(A0-AA)<<std::endl;
                    std::cerr<<"THRESH*Norm(A0) = "<<THRESH<<"*"<<Norm(A0)<<
                        " = "<<THRESH*Norm(A0)<<std::endl;
                    abort();
                }
            }
#endif
        }
    };

    // algo -2: Check for inst
    template <int cs, int rs, class Mu, class Vd, class Ve, class Mv>
    struct SVDecomposeFromBidiagonal_Helper<-2,cs,rs,Mu,Vd,Ve,Mv>
    {
        static TMV_INLINE void call(Mu& U, Vd& D, Ve& E, Mv& V, bool setUV)
        {
            typedef typename Mu::value_type T;
            const bool inst = 
                (cs == TMV_UNKNOWN || cs > 16) &&
                (rs == TMV_UNKNOWN || rs > 16) &&
                Traits<T>::isinst;
            const int algo = 
                Traits<typename Vd::value_type>::iscomplex ? 1 :
                inst ? 90 :
                -3;
            SVDecomposeFromBidiagonal_Helper<algo,cs,rs,Mu,Vd,Ve,Mv>::call(
                U,D,E,V,setUV);
        }
    };

    template <int cs, int rs, class Mu, class Vd, class Ve, class Mv>
    struct SVDecomposeFromBidiagonal_Helper<-1,cs,rs,Mu,Vd,Ve,Mv>
    {
        static TMV_INLINE void call(Mu& U, Vd& D, Ve& E, Mv& V, bool setUV)
        {
            SVDecomposeFromBidiagonal_Helper<-2,cs,rs,Mu,Vd,Ve,Mv>::call(
                U,D,E,V,setUV); 
        }
    };

    template <class Mu, class Vd, class Ve, class Mv>
    static inline void InlineSV_DecomposeFromBidiagonal(
        BaseMatrix_Rec_Mutable<Mu>& U,
        BaseVector_Mutable<Vd>& D, BaseVector_Mutable<Ve>& E, 
        BaseMatrix_Rec_Mutable<Mv>& V, bool setUV)
    {
        TMVStaticAssert((Traits2<
                         typename Mu::value_type,
                         typename Vd::value_type>::samebase));
        TMVStaticAssert((Traits2<
                         typename Vd::value_type,
                         typename Ve::value_type>::sametype));
        TMVStaticAssert((Traits2<
                         typename Mu::value_type,
                         typename Mv::value_type>::samebase));
        TMVStaticAssert((Sizes<Mu::_rowsize,Vd::_size>::same));
        TMVStaticAssert((Sizes<Vd::_size,IntTraits<Ve::_size>::Sp1>::same));
        TMVStaticAssert((Sizes<Mu::_rowsize,Mv::_colsize>::same));
        TMVStaticAssert((Sizes<Mu::_rowsize,Mv::_rowsize>::same));
        TMVAssert(D.size() == E.size()+1);
        if (U.ptr()) {
            TMVAssert(U.colsize() >= U.rowsize());
            TMVAssert(U.rowsize() == D.size());
        }
        if (V.ptr()) {
            TMVAssert(V.colsize() == D.size());
            TMVAssert(V.rowsize() == D.size());
        }
        const int cs = Mu::_colsize;
        const int rs1 = Sizes<Mu::_rowsize,Vd::_size>::size;
        const int rs2 = Sizes<Mv::_rowsize,Mv::_colsize>::size;
        const int rs = Sizes<rs1,rs2>::size;
        typedef typename Mu::cview_type Muv;
        typedef typename Vd::cview_type Vdv;
        typedef typename Ve::cview_type Vev;
        typedef typename Mv::cview_type Mvv;
        TMV_MAYBE_REF(Mu,Muv) Uv = U.cView();
        TMV_MAYBE_REF(Vd,Vdv) Dv = D.cView();
        TMV_MAYBE_REF(Ve,Vev) Ev = E.cView();
        TMV_MAYBE_REF(Mv,Mvv) Vv = V.cView();
        SVDecomposeFromBidiagonal_Helper<-3,cs,rs,Muv,Vdv,Vev,Mvv>::call(
            Uv,Dv,Ev,Vv,setUV);
    }

    template <class Mu, class Vd, class Ve, class Mv>
    static inline void SV_DecomposeFromBidiagonal(
        BaseMatrix_Rec_Mutable<Mu>& U,
        BaseVector_Mutable<Vd>& D, BaseVector_Mutable<Ve>& E, 
        BaseMatrix_Rec_Mutable<Mv>& V, bool setUV)
    {
        TMVStaticAssert((Traits2<
                         typename Mu::value_type,
                         typename Vd::value_type>::samebase));
        TMVStaticAssert((Traits2<
                         typename Vd::value_type,
                         typename Ve::value_type>::sametype));
        TMVStaticAssert((Traits2<
                         typename Mu::value_type,
                         typename Mv::value_type>::samebase));
        TMVStaticAssert((Sizes<Mu::_rowsize,Vd::_size>::same));
        TMVStaticAssert((Sizes<Vd::_size,IntTraits<Ve::_size>::Sp1>::same));
        TMVStaticAssert((Sizes<Mu::_rowsize,Mv::_colsize>::same));
        TMVStaticAssert((Sizes<Mu::_rowsize,Mv::_rowsize>::same));
        TMVAssert(D.size() == E.size()+1);
        if (U.ptr()) {
            TMVAssert(U.colsize() >= U.rowsize());
            TMVAssert(U.rowsize() == D.size());
        }
        if (V.ptr()) {
            TMVAssert(V.colsize() == D.size());
            TMVAssert(V.rowsize() == D.size());
        }
        const int cs = Mu::_colsize;
        const int rs1 = Sizes<Mu::_rowsize,Vd::_size>::size;
        const int rs2 = Sizes<Mv::_rowsize,Mv::_colsize>::size;
        const int rs = Sizes<rs1,rs2>::size;
        typedef typename Mu::cview_type Muv;
        typedef typename Vd::cview_type Vdv;
        typedef typename Ve::cview_type Vev;
        typedef typename Mv::cview_type Mvv;
        TMV_MAYBE_REF(Mu,Muv) Uv = U.cView();
        TMV_MAYBE_REF(Vd,Vdv) Dv = D.cView();
        TMV_MAYBE_REF(Ve,Vev) Ev = E.cView();
        TMV_MAYBE_REF(Mv,Mvv) Vv = V.cView();
        SVDecomposeFromBidiagonal_Helper<-2,cs,rs,Muv,Vdv,Vev,Mvv>::call(
            Uv,Dv,Ev,Vv,setUV);
    }


    //
    // And finally the main driver routine.
    // Decompose A (input as U) into U S V
    // where S is a diagonal real matrix, and U,V are unitary matrices.
    // A,U are M x N (M >= N)
    // S,V are N x N
    // The determinant of UV is returned in signdet, logdet.
    //
    
    // Defined in TMV_SVDecompose.cpp
    template <class T, class RT>
    void InstSV_Decompose(
        MatrixView<T> U, DiagMatrixView<RT> S,
        MatrixView<T> V, T& signuv, RT& logdet, bool StoreU);

    template <int algo, int cs, int rs, class Mu, class Ms, class Mv>
    struct SVDecompose_Helper;

    // algo 0: Trivial, nothing to do (M == 0, or N == 0)
    template <int cs, int rs, class Mu, class Ms, class Mv>
    struct SVDecompose_Helper<0,cs,rs,Mu,Ms,Mv>
    {
        typedef typename Mu::float_type FT;
        typedef typename Mu::zfloat_type ZT;
        static TMV_INLINE void call(Mu& , Ms& , Mv& , ZT& , FT& , bool ) {} 
    };

    // algo 11: Normal algorithm
    template <int cs, int rs, class Mu, class Ms, class Mv>
    struct SVDecompose_Helper<11,cs,rs,Mu,Ms,Mv>
    {
        template <bool isreal, int dummy>
        struct ImagIsZeroHelper;

        template <int dummy> 
        struct ImagIsZeroHelper<true,dummy>
        {
            template <class M>
            static inline bool call(M& ) { return true; }
        };
        template <int dummy> 
        struct ImagIsZeroHelper<false,dummy>
        {
            template <class M>
            static inline bool call(M& m)
            { return Norm1(m.imagPart()) == typename M::real_type(0); }
        };


        typedef typename Mu::float_type FT;
        typedef typename Mu::zfloat_type ZT;
        static void call(
            Mu& U, Ms& S, Mv& V, ZT& signdet, FT& logdet, bool StoreU)
        {
            const int M = cs==TMV_UNKNOWN ? int(U.colsize()) : cs;
            const int N = rs==TMV_UNKNOWN ? int(U.rowsize()) : rs;
            if (N == 0) return;
#ifdef PRINTALGO_SVD
            std::cout<<"SVDecompose algo 11: M,N,cs,rs = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<std::endl;
#endif
            typedef typename Mu::real_type RT;
            typedef typename Mu::value_type T;

            if (M <= 5*N/3) {
                // First we reduce A to bidiagonal form: A = U * B * V
                // using a series of Householder transformations.
                // The diagonal of the Bidiagonal Matrix B is stored in D.
                // The superdiagonal is stored in E.
                const int rsm1 = IntTraits<rs>::Sm1;
                typename VCopyHelper<T,rs>::type D = VectorSizer<T>(N);
                typename VCopyHelper<T,rsm1>::type E = VectorSizer<T>(N-1);
                typename VCopyHelper<RT,rs>::type Ubeta = VectorSizer<T>(N);
                typename VCopyHelper<RT,rsm1>::type Vbeta = VectorSizer<T>(N-1);
                dbgcout<<"Before Bidiagonalize:\n";
                dbgcout<<"U.maxAbs = "<<U.maxAbsElement()<<std::endl;
                Bidiagonalize(U,Ubeta,Vbeta,D,E);
                dbgcout<<"After Bidiagonalize:\n";
                dbgcout<<"U.maxAbs = "<<U.maxAbsElement()<<std::endl;
                dbgcout<<"D.maxAbs = "<<D.maxAbsElement()<<std::endl;
                dbgcout<<"E.maxAbs = "<<E.maxAbsElement()<<std::endl;
                dbgcout<<"D = "<<D<<std::endl;
                dbgcout<<"E = "<<E<<std::endl;

                // The determinant of B is just the product of the 
                // diagonal elements.
                // Calculate it now so we don't have to keep track of
                // it during the reduction stage.
                if (signdet != ZT(0)) {
                    logdet += D.logProdElements(&signdet);
                    signdet *= CalculateDetQ(Ubeta) * CalculateDetQ(Vbeta);
                }

                // Now UV stores Householder vectors for U in lower 
                // diagonal columns (HLi) and Householder vectors for V in 
                // upper diagonal rows (HRi).
                // The Householder matrices for U are actually the adjoints 
                // of the matrices that bidiagonalize A, and for V are the 
                // transposes:
                // U = HLn-1t ... HL1t HL0t A HR0T HR1T ... HRn-2T
                // Using the fact that H Ht = I, we get A = U B V with:
                // U = HL0 ... HLn-1 
                if (V.ptr()) {
                    V.row(0).makeBasis(0);
                    V.rowRange(1,N) = U.rowRange(0,N-1);
                    V.col(0,1,N).setZero();
                    typename Mv::submatrix_type::transpose_type V1t =
                        V.subMatrix(1,N,1,N).transpose();
                    UnpackQ(V1t,Vbeta);
                    dbgcout<<"V => "<<U<<std::endl;
                    dbgcout<<"Norm(VtV-1) = "<<
                        Norm(V.adjoint()*V-T(1))<<std::endl;
                    dbgcout<<"Norm(VVt-1) = "<<
                        Norm(V*V.adjoint()-T(1))<<std::endl;
                }
                if (StoreU) {
                    UnpackQ(U,Ubeta);
                    dbgcout<<"U => "<<U<<std::endl;
                    dbgcout<<"Norm(UtU-1) = "<<
                        Norm(U.adjoint()*U-T(1))<<std::endl;
                }

                if (StoreU) {
                    SV_DecomposeFromBidiagonal(U,D,E,V,false);
                    dbgcout<<"After DecomposeFromBidiag: Norm(UtU-1) = "<<
                        Norm(U.adjoint()*U-T(1))<<std::endl;
                    dbgcout<<"U = "<<U<<std::endl;
                    dbgcout<<"D = "<<D<<std::endl;
                    dbgcout<<"E = "<<E<<std::endl;
                    dbgcout<<"V = "<<V<<std::endl;
                } else {
                    MatrixView<T> U0(0,0,0,1,1);
                    SV_DecomposeFromBidiagonal(U0,D,E,V,false);
                }
                
                // Now S is the real part of D:
                S.diag() = D.realPart();
                dbgcout<<"S = "<<S.diag()<<std::endl;
                TMVAssert((ImagIsZeroHelper<Traits<T>::isreal,1>::call(D)));
            } else {
                // If M is much larger than N (technically M > 5/3 N),
                // then it is quicker to start by doing a QR decomposition 
                // and then do SVD on the square R matrix.  
                // Thus, the final U of the SVD is Q (from the QR decomp)
                // times U from R's SVD.
                if (StoreU) {
                    dbgcout<<"Large M, StoreU\n";
                    typename VCopyHelper<RT,rs>::type Qbeta(N);
                    QR_Decompose(U,Qbeta);
                    dbgcout<<"After QR_Decompose: U = "<<U<<std::endl;
                    typename MCopyHelper<T,Rec,rs,rs>::type R = U.upperTri();
                    dbgcout<<"R = "<<R<<std::endl;
                    UnpackQ(U,Qbeta);
                    dbgcout<<"U = "<<U<<std::endl;
                    SV_Decompose(R,S,V,signdet,logdet,StoreU);
                    dbgcout<<"U1 = "<<R<<std::endl;
                    dbgcout<<"S = "<<S<<std::endl;
                    dbgcout<<"V = "<<V<<std::endl;
                    signdet *= CalculateDetQ(Qbeta);
                    // Now R is a Unitary Matrix U'.  Need to multiply U' by Q
                    U = U*R;
                    dbgcout<<"U*U1 = "<<U<<std::endl;
                } else {
                    typename VCopyHelper<RT,rs>::type Qbeta(N);
                    QR_Decompose(U,Qbeta);
                    typename Mu::rowrange_type R = U.rowRange(0,N);
                    if (N > 1) R.lowerTri().offDiag().setZero();
                    SV_Decompose(R,S,V,signdet,logdet,StoreU);
                    signdet *= CalculateDetQ(Qbeta);
                }
            }
        }
    };

    // algo 81: Copy U to colmajor
    template <int cs, int rs, class Mu, class Ms, class Mv>
    struct SVDecompose_Helper<81,cs,rs,Mu,Ms,Mv>
    {
        typedef typename Mu::float_type FT;
        typedef typename Mu::zfloat_type ZT;
        static inline void call(
            Mu& U, Ms& S, Mv& V, ZT& signdet, FT& logdet, bool StoreU)
        {
#ifdef PRINTALGO_SVD
            std::cout<<"SVDecompose algo 81: cs,rs = "<<cs<<','<<rs<<std::endl;
#endif
            typedef typename Mu::value_type T;
            typedef typename MCopyHelper<T,Rec,cs,rs>::type Mucm;
            Mucm Ucm = U;
            SVDecompose_Helper<-2,cs,rs,Mucm,Ms,Mv>::call(
                Ucm,S,V,signdet,logdet,StoreU);
            if (StoreU) NoAliasCopy(Ucm,U);
        }
    };

    // algo 90: call InstSV_Decompose
    template <int cs, int rs, class Mu, class Ms, class Mv>
    struct SVDecompose_Helper<90,cs,rs,Mu,Ms,Mv>
    {
        typedef typename Mu::float_type FT;
        typedef typename Mu::zfloat_type ZT;
        static TMV_INLINE void call(
            Mu& U, Ms& S, Mv& V, ZT& signdet, FT& logdet, bool StoreU)
        {
            InstSV_Decompose(
                U.xView(),S.xView(),V.xView(),signdet,logdet,StoreU); 
        }
    };

    // algo 95: Conjugate V
    template <int cs, int rs, class Mu, class Ms, class Mv>
    struct SVDecompose_Helper<95,cs,rs,Mu,Ms,Mv>
    {
        typedef typename Mu::float_type FT;
        typedef typename Mu::zfloat_type ZT;
        static TMV_INLINE void call(
            Mu& U, Ms& S, Mv& V, ZT& signdet, FT& logdet, bool StoreU)
        {
            typedef typename Mv::conjugate_type Mvc;
            Mvc Vc = V.conjugate();
            SVDecompose_Helper<-2,cs,rs,Mu,Ms,Mvc>::call(
                U,S,Vc,signdet,logdet,StoreU);
            V.conjugateSelf();
        }
    };

    // algo 97: Conjugate
    template <int cs, int rs, class Mu, class Ms, class Mv>
    struct SVDecompose_Helper<97,cs,rs,Mu,Ms,Mv>
    {
        typedef typename Mu::float_type FT;
        typedef typename Mu::zfloat_type ZT;
        static TMV_INLINE void call(
            Mu& U, Ms& S, Mv& V, ZT& signdet, FT& logdet, bool StoreU)
        {
            typedef typename Mu::conjugate_type Muc;
            typedef typename Mv::conjugate_type Mvc;
            Muc Uc = U.conjugate();
            Mvc Vc = V.conjugate();
            SVDecompose_Helper<-2,cs,rs,Muc,Ms,Mvc>::call(
                Uc,S,Vc,signdet,logdet,StoreU);
        }
    };

    // algo -3: Determine which algorithm to use
    template <int cs, int rs, class Mu, class Ms, class Mv>
    struct SVDecompose_Helper<-3,cs,rs,Mu,Ms,Mv>
    {
        typedef typename Mu::float_type FT;
        typedef typename Mu::zfloat_type ZT;
        static TMV_INLINE void call(
            Mu& U, Ms& S, Mv& V, ZT& signdet, FT& logdet, bool StoreU)
        {
            const int algo = (
                ( cs != TMV_UNKNOWN && rs != TMV_UNKNOWN &&
                  cs <= 16 && rs <= 16 ) ? 11 :
                ( TMV_OPT >= 2 && !Mu::_colmajor ) ? 81 :
                11 );
#ifdef PRINTALGO_SVD
            std::cout<<"Inline SVDecompose: \n";
            std::cout<<"U = "<<TMV_Text(U)<<std::endl;
            std::cout<<"S = "<<TMV_Text(S)<<std::endl;
            std::cout<<"V = "<<TMV_Text(V)<<std::endl;
            std::cout<<"cs = "<<cs<<"  rs = "<<rs<<std::endl;
            std::cout<<"sizes = "<<U.colsize()<<"  "<<U.rowsize()<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
#endif
#ifdef XDEBUG_SVD
            typedef typename Mu::value_type T;
            Matrix<T> A0(U);
            dbgcout<<"SVDecompose:\n";
            //dbgcout<<"A0 = "<<A0<<std::endl;
            dbgcout<<"StoreU = "<<StoreU<<std::endl;
#endif
            SVDecompose_Helper<algo,cs,rs,Mu,Ms,Mv>::call(
                U,S,V,signdet,logdet,StoreU);
#ifdef XDEBUG_SVD
            dbgcout<<"S = "<<S.diag()<<std::endl;
            if (StoreU && V.ptr() && S.size()>0) {
                Matrix<T> A2 = U * S * V;
                dbgcout<<"SVDecompose: Norm(A0-A2) = "<<Norm(A0-A2)<<std::endl;
                dbgcout<<"cf "<<THRESH*Norm(U)*Norm(S)*Norm(V)<<std::endl;
                if (!(Norm(A0-A2) < THRESH * Norm(U) * Norm(S) * Norm(V))) {
                    std::cerr<<"SV_Decompose:\n";
                    std::cerr<<"A = "<<A0<<std::endl;
                    std::cerr<<"U = "<<U<<std::endl;
                    std::cerr<<"S = "<<S.diag()<<std::endl;
                    std::cerr<<"V = "<<V<<std::endl;
                    std::cerr<<"USV = "<<A2<<std::endl;
                    abort();
                }
            }
#endif
#ifdef PRINTALGO_SVD
            std::cout<<"Done SVDecompose\n";
#endif
        }
    };

    // algo -2: Check for inst
    template <int cs, int rs, class Mu, class Ms, class Mv>
    struct SVDecompose_Helper<-2,cs,rs,Mu,Ms,Mv>
    {
        typedef typename Mu::float_type FT;
        typedef typename Mu::zfloat_type ZT;
        static TMV_INLINE void call(
            Mu& U, Ms& S, Mv& V, ZT& signdet, FT& logdet, bool StoreU)
        {
            typedef typename Mu::value_type T;
            const bool inst = 
                (cs == TMV_UNKNOWN || cs > 16) &&
                (rs == TMV_UNKNOWN || rs > 16) &&
                Traits<T>::isinst;
            const int algo = 
                cs == 0 || rs == 0 ? 0 :
                Mu::_conj ? 97 :
                inst ? (Mv::_conj ? 95 : 90) :
                -3;
            SVDecompose_Helper<algo,cs,rs,Mu,Ms,Mv>::call(
                U,S,V,signdet,logdet,StoreU);
        }
    };

    template <int cs, int rs, class Mu, class Ms, class Mv>
    struct SVDecompose_Helper<-1,cs,rs,Mu,Ms,Mv>
    {
        typedef typename Mu::float_type FT;
        typedef typename Mu::zfloat_type ZT;
        static TMV_INLINE void call(
            Mu& U, Ms& S, Mv& V, ZT& signdet, FT& logdet, bool StoreU)
        {
            SVDecompose_Helper<-2,cs,rs,Mu,Ms,Mv>::call(
                U,S,V,signdet,logdet,StoreU); 
        }
    };

    template <class Mu, class Ms, class Mv>
    static inline void InlineSV_Decompose(
        BaseMatrix_Rec_Mutable<Mu>& U, BaseMatrix_Diag_Mutable<Ms>& S,
        BaseMatrix_Rec_Mutable<Mv>& V, 
        typename Mu::zfloat_type& signdet, typename Mu::float_type& logdet,
        bool StoreU)
    {
        TMVStaticAssert(Ms::isreal);
        TMVStaticAssert((Traits2<
                         typename Mu::value_type,
                         typename Ms::value_type>::samebase));
        TMVStaticAssert((Traits2<
                         typename Mu::value_type,
                         typename Mv::value_type>::sametype));
        TMVStaticAssert((Sizes<Mu::_rowsize,Ms::_size>::same));
        TMVStaticAssert((Sizes<Mu::_rowsize,Mv::_colsize>::same));
        TMVStaticAssert((Sizes<Mu::_rowsize,Mv::_rowsize>::same));
        TMVAssert(U.colsize() >= U.rowsize());
        TMVAssert(S.size() == U.rowsize());
        if (V.ptr()) {
            TMVAssert(V.colsize() == U.rowsize());
            TMVAssert(V.rowsize() == U.rowsize());
        }
        const int cs = Mu::_colsize;
        const int rs1 = Sizes<Mu::_rowsize,Ms::_size>::size;
        const int rs2 = Sizes<Mv::_rowsize,Mv::_colsize>::size;
        const int rs = Sizes<rs1,rs2>::size;
        typedef typename Mu::cview_type Muv;
        typedef typename Ms::cview_type Msv;
        typedef typename Mv::cview_type Mvv;
        TMV_MAYBE_REF(Mu,Muv) Uv = U.cView();
        TMV_MAYBE_REF(Ms,Msv) Sv = S.cView();
        TMV_MAYBE_REF(Mv,Mvv) Vv = V.cView();
        SVDecompose_Helper<-3,cs,rs,Muv,Msv,Mvv>::call(
            Uv,Sv,Vv,signdet,logdet,StoreU);
    }

    // This is the basic functionality
    template <class Mu, class Ms, class Mv>
    static inline void SV_Decompose(
        BaseMatrix_Rec_Mutable<Mu>& U, BaseMatrix_Diag_Mutable<Ms>& S,
        BaseMatrix_Rec_Mutable<Mv>& V,
        typename Mu::zfloat_type& signdet, typename Mu::float_type& logdet,
        bool StoreU)
    {
        TMVStaticAssert(Ms::isreal);
        TMVStaticAssert((Traits2<
                         typename Mu::value_type,
                         typename Ms::value_type>::samebase));
        TMVStaticAssert((Traits2<
                         typename Mu::value_type,
                         typename Mv::value_type>::sametype));
        TMVStaticAssert((Sizes<Mu::_rowsize,Ms::_size>::same));
        TMVStaticAssert((Sizes<Mu::_rowsize,Mv::_colsize>::same));
        TMVStaticAssert((Sizes<Mu::_rowsize,Mv::_rowsize>::same));
        TMVAssert(U.colsize() >= U.rowsize());
        TMVAssert(S.size() == U.rowsize());
        if (V.ptr()) {
            TMVAssert(V.colsize() == U.rowsize());
            TMVAssert(V.rowsize() == U.rowsize());
        }
        const int cs = Mu::_colsize;
        const int rs1 = Sizes<Mu::_rowsize,Ms::_size>::size;
        const int rs2 = Sizes<Mv::_rowsize,Mv::_colsize>::size;
        const int rs = Sizes<rs1,rs2>::size;
        typedef typename Mu::cview_type Muv;
        typedef typename Ms::cview_type Msv;
        typedef typename Mv::cview_type Mvv;
        TMV_MAYBE_REF(Mu,Muv) Uv = U.cView();
        TMV_MAYBE_REF(Ms,Msv) Sv = S.cView();
        TMV_MAYBE_REF(Mv,Mvv) Vv = V.cView();
        SVDecompose_Helper<-2,cs,rs,Muv,Msv,Mvv>::call(
            Uv,Sv,Vv,signdet,logdet,StoreU);
    }

    // The rest of these below are basically convenience functions
    // to allow the user to provide fewer or different arguments.
    template <class Mu, class Ms, class Mv>
    static inline void SV_Decompose(
        BaseMatrix_Rec_Mutable<Mu>& U, BaseMatrix_Diag_Mutable<Ms>& S,
        BaseMatrix_Rec_Mutable<Mv>& V, bool StoreU=true)
    {
        typename Mu::zfloat_type signdet(0);
        typename Mu::float_type logdet;
        SV_Decompose(U,S,V,signdet,logdet,StoreU);
    }

    template <class Mu, class Ms>
    static inline void SV_Decompose(
        BaseMatrix_Rec_Mutable<Mu>& U, BaseMatrix_Diag_Mutable<Ms>& S,
        bool StoreU)
    {
        MatrixView<typename Mu::value_type> V(0,0,0,1,1);
        typename Mu::zfloat_type signdet(0);
        typename Mu::float_type logdet;
        SV_Decompose(U,S,V,signdet,logdet,StoreU);
    }

    // Allow views as an argument by value (for convenience)
    template <class T, int Au, int As, int Av>
    static inline void SV_Decompose(
        MatrixView<T,Au> U,
        DiagMatrixView<typename Traits<T>::real_type,As> S,
        MatrixView<T,Av> V, bool StoreU=true)
    {
        typedef typename Traits<T>::real_type RT;
        typedef MatrixView<T,Au> Mu;
        typedef DiagMatrixView<RT,As> Ms;
        typedef MatrixView<T,Av> Mv;
        SV_Decompose(
            static_cast<BaseMatrix_Rec_Mutable<Mu>&>(U),
            static_cast<BaseMatrix_Diag_Mutable<Ms>&>(S),
            static_cast<BaseMatrix_Rec_Mutable<Mv>&>(V), StoreU);
    }

    template <class T, int M, int N, int Siu, int Sju, int Au, int Ss, int As, int Siv, int Sjv, int Av>
    static inline void SV_Decompose(
        SmallMatrixView<T,M,N,Siu,Sju,Au> U,
        SmallDiagMatrixView<typename Traits<T>::real_type,N,Ss,As> S,
        SmallMatrixView<T,N,N,Siv,Sjv,Av> V, bool StoreU=true)
    {
        typedef typename Traits<T>::real_type RT;
        typedef SmallMatrixView<T,M,N,Siu,Sju,Au> Mu;
        typedef SmallDiagMatrixView<RT,N,Ss,As> Ms;
        typedef SmallMatrixView<T,N,N,Siv,Sjv,Av> Mv;
        SV_Decompose(
            static_cast<BaseMatrix_Rec_Mutable<Mu>&>(U),
            static_cast<BaseMatrix_Diag_Mutable<Ms>&>(S),
            static_cast<BaseMatrix_Rec_Mutable<Mv>&>(V), StoreU);
    }

    template <class T, int N, int Au, int As, int Siv, int Sjv, int Av>
    static inline void SV_Decompose(
        MatrixView<T,Au> U,
        DiagMatrixView<typename Traits<T>::real_type,As> S,
        SmallMatrixView<T,N,N,Siv,Sjv,Av> V, bool StoreU=true)
    {
        typedef typename Traits<T>::real_type RT;
        typedef MatrixView<T,Au> Mu;
        typedef DiagMatrixView<RT,As> Ms;
        typedef SmallMatrixView<T,N,N,Siv,Sjv,Av> Mv;
        SV_Decompose(
            static_cast<BaseMatrix_Rec_Mutable<Mu>&>(U),
            static_cast<BaseMatrix_Diag_Mutable<Ms>&>(S),
            static_cast<BaseMatrix_Rec_Mutable<Mv>&>(V), StoreU);
    }

    template <class T, int N, int Au, int Ss, int As, int Av>
    static inline void SV_Decompose(
        MatrixView<T,Au> U,
        SmallDiagMatrixView<typename Traits<T>::real_type,N,Ss,As> S,
        MatrixView<T,Av> V, bool StoreU=true)
    {
        typedef typename Traits<T>::real_type RT;
        typedef MatrixView<T,Au> Mu;
        typedef SmallDiagMatrixView<RT,N,Ss,As> Ms;
        typedef MatrixView<T,Av> Mv;
        SV_Decompose(
            static_cast<BaseMatrix_Rec_Mutable<Mu>&>(U),
            static_cast<BaseMatrix_Diag_Mutable<Ms>&>(S),
            static_cast<BaseMatrix_Rec_Mutable<Mv>&>(V), StoreU);
    }

    template <class T, int N, int Au, int Ss, int As, int Siv, int Sjv, int Av>
    static inline void SV_Decompose(
        MatrixView<T,Au> U,
        SmallDiagMatrixView<typename Traits<T>::real_type,N,Ss,As> S,
        SmallMatrixView<T,N,N,Siv,Sjv,Av> V, bool StoreU=true)
    {
        typedef typename Traits<T>::real_type RT;
        typedef MatrixView<T,Au> Mu;
        typedef SmallDiagMatrixView<RT,N,Ss,As> Ms;
        typedef SmallMatrixView<T,N,N,Siv,Sjv,Av> Mv;
        SV_Decompose(
            static_cast<BaseMatrix_Rec_Mutable<Mu>&>(U),
            static_cast<BaseMatrix_Diag_Mutable<Ms>&>(S),
            static_cast<BaseMatrix_Rec_Mutable<Mv>&>(V), StoreU);
    }

    template <class T, int M, int N, int Siu, int Sju, int Au, int As, int Av>
    static inline void SV_Decompose(
        SmallMatrixView<T,M,N,Siu,Sju,Au> U,
        DiagMatrixView<typename Traits<T>::real_type,As> S,
        MatrixView<T,Av> V, bool StoreU=true)
    {
        typedef typename Traits<T>::real_type RT;
        typedef SmallMatrixView<T,M,N,Siu,Sju,Au> Mu;
        typedef DiagMatrixView<RT,As> Ms;
        typedef MatrixView<T,Av> Mv;
        SV_Decompose(
            static_cast<BaseMatrix_Rec_Mutable<Mu>&>(U),
            static_cast<BaseMatrix_Diag_Mutable<Ms>&>(S),
            static_cast<BaseMatrix_Rec_Mutable<Mv>&>(V), StoreU);
    }

    template <class T, int M, int N, int Siu, int Sju, int Au, int As, int Siv, int Sjv, int Av>
    static inline void SV_Decompose(
        SmallMatrixView<T,M,N,Siu,Sju,Au> U,
        DiagMatrixView<typename Traits<T>::real_type,As> S,
        SmallMatrixView<T,N,N,Siv,Sjv,Av> V, bool StoreU=true)
    {
        typedef typename Traits<T>::real_type RT;
        typedef SmallMatrixView<T,M,N,Siu,Sju,Au> Mu;
        typedef DiagMatrixView<RT,As> Ms;
        typedef SmallMatrixView<T,N,N,Siv,Sjv,Av> Mv;
        SV_Decompose(
            static_cast<BaseMatrix_Rec_Mutable<Mu>&>(U),
            static_cast<BaseMatrix_Diag_Mutable<Ms>&>(S),
            static_cast<BaseMatrix_Rec_Mutable<Mv>&>(V), StoreU);
    }

    template <class T, int M, int N, int Siu, int Sju, int Au, int Ss, int As, int Av>
    static inline void SV_Decompose(
        SmallMatrixView<T,M,N,Siu,Sju,Au> U,
        SmallDiagMatrixView<typename Traits<T>::real_type,N,Ss,As> S,
        MatrixView<T,Av> V, bool StoreU=true)
    {
        typedef typename Traits<T>::real_type RT;
        typedef SmallMatrixView<T,M,N,Siu,Sju,Au> Mu;
        typedef SmallDiagMatrixView<RT,N,Ss,As> Ms;
        typedef MatrixView<T,Av> Mv;
        SV_Decompose(
            static_cast<BaseMatrix_Rec_Mutable<Mu>&>(U),
            static_cast<BaseMatrix_Diag_Mutable<Ms>&>(S),
            static_cast<BaseMatrix_Rec_Mutable<Mv>&>(V), StoreU);
    }

    template <class T, int Au, int As>
    static inline void SV_Decompose(
        MatrixView<T,Au> U,
        DiagMatrixView<typename Traits<T>::real_type,As> S, bool StoreU=true)
    {
        typedef typename Traits<T>::real_type RT;
        typedef MatrixView<T,Au> Mu;
        typedef DiagMatrixView<RT,As> Ms;
        SV_Decompose(
            static_cast<BaseMatrix_Rec_Mutable<Mu>&>(U),
            static_cast<BaseMatrix_Diag_Mutable<Ms>&>(S), StoreU);
    }

    template <class T, int N, int Au, int Ss, int As>
    static inline void SV_Decompose(
        MatrixView<T,Au> U,
        SmallDiagMatrixView<typename Traits<T>::real_type,N,Ss,As> S,
        bool StoreU=true)
    {
        typedef typename Traits<T>::real_type RT;
        typedef MatrixView<T,Au> Mu;
        typedef SmallDiagMatrixView<RT,N,Ss,As> Ms;
        SV_Decompose(
            static_cast<BaseMatrix_Rec_Mutable<Mu>&>(U),
            static_cast<BaseMatrix_Diag_Mutable<Ms>&>(S), StoreU);
    }

    template <class T, int M, int N, int Siu, int Sju, int Au, int As>
    static inline void SV_Decompose(
        SmallMatrixView<T,M,N,Siu,Sju,Au> U,
        DiagMatrixView<typename Traits<T>::real_type,As> S, bool StoreU=true)
    {
        typedef typename Traits<T>::real_type RT;
        typedef SmallMatrixView<T,M,N,Siu,Sju,Au> Mu;
        typedef DiagMatrixView<RT,As> Ms;
        SV_Decompose(
            static_cast<BaseMatrix_Rec_Mutable<Mu>&>(U),
            static_cast<BaseMatrix_Diag_Mutable<Ms>&>(S), StoreU);
    }

    template <class T, int M, int N, int Siu, int Sju, int Au, int Ss, int As>
    static inline void SV_Decompose(
        SmallMatrixView<T,M,N,Siu,Sju,Au> U,
        SmallDiagMatrixView<typename Traits<T>::real_type,N,Ss,As> S,
        bool StoreU=true)
    {
        typedef typename Traits<T>::real_type RT;
        typedef SmallMatrixView<T,M,N,Siu,Sju,Au> Mu;
        typedef SmallDiagMatrixView<RT,N,Ss,As> Ms;
        SV_Decompose(
            static_cast<BaseMatrix_Rec_Mutable<Mu>&>(U),
            static_cast<BaseMatrix_Diag_Mutable<Ms>&>(S), StoreU);
    }

    // Don't forget the ones that mix *MatrixView and BaseMatrix_*_Mutable:
    // U is _Mutable
    template <class Mu, class T, int As, int Av>
    static inline void SV_Decompose(
        BaseMatrix_Rec_Mutable<Mu>& U,
        DiagMatrixView<typename Traits<T>::real_type,As> S,
        MatrixView<T,Av> V, bool StoreU=true)
    {
        typedef typename Traits<T>::real_type RT;
        typedef DiagMatrixView<RT,As> Ms;
        typedef MatrixView<T,Av> Mv;
        SV_Decompose(
            U, static_cast<BaseMatrix_Diag_Mutable<Ms>&>(S),
            static_cast<BaseMatrix_Rec_Mutable<Mv>&>(V), StoreU);
    }

    template <class Mu, class T, int N, int As, int Siv, int Sjv, int Av>
    static inline void SV_Decompose(
        BaseMatrix_Rec_Mutable<Mu>& U,
        DiagMatrixView<typename Traits<T>::real_type,As> S,
        SmallMatrixView<T,N,N,Siv,Sjv,Av> V, bool StoreU=true)
    {
        typedef typename Traits<T>::real_type RT;
        typedef DiagMatrixView<RT,As> Ms;
        typedef SmallMatrixView<T,N,N,Siv,Sjv,Av> Mv;
        SV_Decompose(
            U, static_cast<BaseMatrix_Diag_Mutable<Ms>&>(S),
            static_cast<BaseMatrix_Rec_Mutable<Mv>&>(V), StoreU);
    }

    template <class Mu, class T, int N, int Ss, int As, int Av>
    static inline void SV_Decompose(
        BaseMatrix_Rec_Mutable<Mu>& U,
        SmallDiagMatrixView<typename Traits<T>::real_type,N,Ss,As> S,
        MatrixView<T,Av> V, bool StoreU=true)
    {
        typedef typename Traits<T>::real_type RT;
        typedef SmallDiagMatrixView<RT,N,Ss,As> Ms;
        typedef MatrixView<T,Av> Mv;
        SV_Decompose(
            U, static_cast<BaseMatrix_Diag_Mutable<Ms>&>(S),
            static_cast<BaseMatrix_Rec_Mutable<Mv>&>(V), StoreU);
    }

    template <class Mu, class T, int N, int Ss, int As, int Siv, int Sjv, int Av>
    static inline void SV_Decompose(
        BaseMatrix_Rec_Mutable<Mu>& U,
        SmallDiagMatrixView<typename Traits<T>::real_type,N,Ss,As> S,
        SmallMatrixView<T,N,N,Siv,Sjv,Av> V, bool StoreU=true)
    {
        typedef typename Traits<T>::real_type RT;
        typedef SmallDiagMatrixView<RT,N,Ss,As> Ms;
        typedef SmallMatrixView<T,N,N,Siv,Sjv,Av> Mv;
        SV_Decompose(
            U, static_cast<BaseMatrix_Diag_Mutable<Ms>&>(S),
            static_cast<BaseMatrix_Rec_Mutable<Mv>&>(V), StoreU);
    }

    template <class Mu, class T, int As>
    static inline void SV_Decompose(
        BaseMatrix_Rec_Mutable<Mu>& U,
        DiagMatrixView<typename Traits<T>::real_type,As> S, bool StoreU=true)
    {
        typedef typename Traits<T>::real_type RT;
        typedef DiagMatrixView<RT,As> Ms;
        SV_Decompose(U, static_cast<BaseMatrix_Diag_Mutable<Ms>&>(S), StoreU);
    }

    template <class Mu, class T, int N, int Ss, int As>
    static inline void SV_Decompose(
        BaseMatrix_Rec_Mutable<Mu>& U,
        SmallDiagMatrixView<typename Traits<T>::real_type,N,Ss,As> S,
        bool StoreU=true)
    {
        typedef typename Traits<T>::real_type RT;
        typedef SmallDiagMatrixView<RT,N,Ss,As> Ms;
        SV_Decompose(U, static_cast<BaseMatrix_Diag_Mutable<Ms>&>(S), StoreU);
    }

    // S is _Mutable
    template <class T, int Au, class Ms, int Av>
    static inline void SV_Decompose(
        MatrixView<T,Au> U, BaseMatrix_Diag_Mutable<Ms>& S,
        MatrixView<T,Av> V, bool StoreU=true)
    {
        typedef MatrixView<T,Au> Mu;
        typedef MatrixView<T,Av> Mv;
        SV_Decompose(
            static_cast<BaseMatrix_Rec_Mutable<Mu>&>(U), S,
            static_cast<BaseMatrix_Rec_Mutable<Mv>&>(V), StoreU);
    }

    template <class T, int N, int Au, class Ms, int Siv, int Sjv, int Av>
    static inline void SV_Decompose(
        MatrixView<T,Au> U, BaseMatrix_Diag_Mutable<Ms>& S,
        SmallMatrixView<T,N,N,Siv,Sjv,Av> V, bool StoreU=true)
    {
        typedef MatrixView<T,Au> Mu;
        typedef SmallMatrixView<T,N,N,Siv,Sjv,Av> Mv;
        SV_Decompose(
            static_cast<BaseMatrix_Rec_Mutable<Mu>&>(U), S,
            static_cast<BaseMatrix_Rec_Mutable<Mv>&>(V), StoreU);
    }

    template <class T, int M, int N, int Siu, int Sju, int Au, class Ms, int Av>
    static inline void SV_Decompose(
        SmallMatrixView<T,M,N,Siu,Sju,Au> U,
        BaseMatrix_Diag_Mutable<Ms>& S,
        MatrixView<T,Av> V, bool StoreU=true)
    {
        typedef SmallMatrixView<T,M,N,Siu,Sju,Au> Mu;
        typedef MatrixView<T,Av> Mv;
        SV_Decompose(
            static_cast<BaseMatrix_Rec_Mutable<Mu>&>(U), S,
            static_cast<BaseMatrix_Rec_Mutable<Mv>&>(V), StoreU);
    }

    template <class T, int M, int N, int Siu, int Sju, int Au, class Ms, int Siv, int Sjv, int Av>
    static inline void SV_Decompose(
        SmallMatrixView<T,M,N,Siu,Sju,Au> U,
        BaseMatrix_Diag_Mutable<Ms>& S,
        SmallMatrixView<T,N,N,Siv,Sjv,Av> V, bool StoreU=true)
    {
        typedef SmallMatrixView<T,M,N,Siu,Sju,Au> Mu;
        typedef SmallMatrixView<T,N,N,Siv,Sjv,Av> Mv;
        SV_Decompose(
            static_cast<BaseMatrix_Rec_Mutable<Mu>&>(U), S,
            static_cast<BaseMatrix_Rec_Mutable<Mv>&>(V), StoreU);
    }

    template <class T, int Au, class Ms>
    static inline void SV_Decompose(
        MatrixView<T,Au> U, BaseMatrix_Diag_Mutable<Ms>& S, bool StoreU=true)
    {
        typedef MatrixView<T,Au> Mu;
        SV_Decompose(static_cast<BaseMatrix_Rec_Mutable<Mu>&>(U), S, StoreU);
    }

    template <class T, int M, int N, int Siu, int Sju, int Au, class Ms>
    static inline void SV_Decompose(
        SmallMatrixView<T,M,N,Siu,Sju,Au> U,
        BaseMatrix_Diag_Mutable<Ms>& S, bool StoreU=true)
    {
        typedef SmallMatrixView<T,M,N,Siu,Sju,Au> Mu;
        SV_Decompose(static_cast<BaseMatrix_Rec_Mutable<Mu>&>(U), S, StoreU);
    }

    // V is _Mutable
    template <class T, int Au, int As, class Mv>
    static inline void SV_Decompose(
        MatrixView<T,Au> U,
        DiagMatrixView<typename Traits<T>::real_type,As> S,
        BaseMatrix_Rec_Mutable<Mv> V, bool StoreU=true)
    {
        typedef typename Traits<T>::real_type RT;
        typedef MatrixView<T,Au> Mu;
        typedef DiagMatrixView<RT,As> Ms;
        SV_Decompose(
            static_cast<BaseMatrix_Rec_Mutable<Mu>&>(U),
            static_cast<BaseMatrix_Diag_Mutable<Ms>&>(S), V, StoreU);
    }

    template <class T, int N, int Au, int Ss, int As, class Mv>
    static inline void SV_Decompose(
        MatrixView<T,Au> U,
        SmallDiagMatrixView<typename Traits<T>::real_type,N,Ss,As> S,
        BaseMatrix_Rec_Mutable<Mv> V, bool StoreU=true)
    {
        typedef typename Traits<T>::real_type RT;
        typedef MatrixView<T,Au> Mu;
        typedef SmallDiagMatrixView<RT,N,Ss,As> Ms;
        SV_Decompose(
            static_cast<BaseMatrix_Rec_Mutable<Mu>&>(U),
            static_cast<BaseMatrix_Diag_Mutable<Ms>&>(S), V, StoreU);
    }

    template <class T, int M, int N, int Siu, int Sju, int Au, int As, class Mv>
    static inline void SV_Decompose(
        SmallMatrixView<T,M,N,Siu,Sju,Au> U,
        DiagMatrixView<typename Traits<T>::real_type,As> S,
        BaseMatrix_Rec_Mutable<Mv> V, bool StoreU=true)
    {
        typedef typename Traits<T>::real_type RT;
        typedef SmallMatrixView<T,M,N,Siu,Sju,Au> Mu;
        typedef DiagMatrixView<RT,As> Ms;
        SV_Decompose(
            static_cast<BaseMatrix_Rec_Mutable<Mu>&>(U),
            static_cast<BaseMatrix_Diag_Mutable<Ms>&>(S), V, StoreU);
    }

    template <class T, int M, int N, int Siu, int Sju, int Au, int Ss, int As, class Mv>
    static inline void SV_Decompose(
        SmallMatrixView<T,M,N,Siu,Sju,Au> U,
        SmallDiagMatrixView<typename Traits<T>::real_type,N,Ss,As> S,
        BaseMatrix_Rec_Mutable<Mv> V, bool StoreU=true)
    {
        typedef typename Traits<T>::real_type RT;
        typedef SmallMatrixView<T,M,N,Siu,Sju,Au> Mu;
        typedef SmallDiagMatrixView<RT,N,Ss,As> Ms;
        SV_Decompose(
            static_cast<BaseMatrix_Rec_Mutable<Mu>&>(U),
            static_cast<BaseMatrix_Diag_Mutable<Ms>&>(S), V, StoreU);
    }

    // U,S are _Mutable
    template <class Mu, class Ms, class T, int Av>
    static inline void SV_Decompose(
        BaseMatrix_Rec_Mutable<Mu>& U,
        BaseMatrix_Diag_Mutable<Ms>& S,
        MatrixView<T,Av> V, bool StoreU=true)
    {
        typedef MatrixView<T,Av> Mv;
        SV_Decompose(
            U, S, static_cast<BaseMatrix_Rec_Mutable<Mv>&>(V), StoreU);
    }

    template <class Mu, class Ms, class T, int N, int Siv, int Sjv, int Av>
    static inline void SV_Decompose(
        BaseMatrix_Rec_Mutable<Mu>& U,
        BaseMatrix_Diag_Mutable<Ms>& S,
        SmallMatrixView<T,N,N,Siv,Sjv,Av> V, bool StoreU=true)
    {
        typedef SmallMatrixView<T,N,N,Siv,Sjv,Av> Mv;
        SV_Decompose(
            U, S, static_cast<BaseMatrix_Rec_Mutable<Mv>&>(V), StoreU);
    }

    // S,V are _Mutable
    template <class T, int Au, class Ms, class Mv>
    static inline void SV_Decompose(
        MatrixView<T,Au> U, BaseMatrix_Diag_Mutable<Ms>& S,
        BaseMatrix_Rec_Mutable<Mv>& V, bool StoreU=true)
    {
        typedef MatrixView<T,Au> Mu;
        SV_Decompose(
            static_cast<BaseMatrix_Rec_Mutable<Mu>&>(U), S, V, StoreU);
    }

    template <class T, int M, int N, int Siu, int Sju, int Au, class Ms, class Mv>
    static inline void SV_Decompose(
        SmallMatrixView<T,M,N,Siu,Sju,Au> U, BaseMatrix_Diag_Mutable<Ms>& S,
        BaseMatrix_Rec_Mutable<Mv>& V, bool StoreU=true)
    {
        typedef SmallMatrixView<T,M,N,Siu,Sju,Au> Mu;
        SV_Decompose(
            static_cast<BaseMatrix_Rec_Mutable<Mu>&>(U), S, V, StoreU);
    }

    // U,V are _Mutable
    template <class Mu, class T, int As, class Mv>
    static inline void SV_Decompose(
        BaseMatrix_Rec_Mutable<Mu>& U,
        DiagMatrixView<typename Traits<T>::real_type,As> S,
        BaseMatrix_Rec_Mutable<Mv>& V, bool StoreU=true)
    {
        typedef typename Traits<T>::real_type RT;
        typedef DiagMatrixView<RT,As> Ms;
        SV_Decompose(
            U, static_cast<BaseMatrix_Diag_Mutable<Ms>&>(S), V, StoreU);
    }

    template <class Mu, class T, int N, int Ss, int As, class Mv>
    static inline void SV_Decompose(
        BaseMatrix_Rec_Mutable<Mu>& U,
        SmallDiagMatrixView<typename Traits<T>::real_type,N,Ss,As> S,
        BaseMatrix_Rec_Mutable<Mv> V, bool StoreU=true)
    {
        typedef typename Traits<T>::real_type RT;
        typedef SmallDiagMatrixView<RT,N,Ss,As> Ms;
        SV_Decompose(
            U, static_cast<BaseMatrix_Diag_Mutable<Ms>&>(S), V, StoreU);
    }

} // namespace tmv

#endif

