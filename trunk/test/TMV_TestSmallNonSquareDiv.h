
#include "TMV_Test.h"
#include "TMV_Test_3.h"
#include "TMV.h"
#define NORDIVEQ
#define NOLDIVEQ
#include "TMV_TestMatrixDivArith.h"

template <class T> void TestSmallNonSquareDiv_a();
template <class T> void TestSmallNonSquareDiv_b();
template <class T> void TestSmallNonSquareDiv_c();
template <class T> void TestSmallNonSquareDiv_d();
template <class T> void TestSmallNonSquareDiv_e();
template <class T> void TestSmallNonSquareDiv_f();
template <class T> void TestSmallNonSquareDiv_g();

template <class T, tmv::StorageType stor, int M, int N> 
static void TestSmallNonSquareDiv_Basic(std::string label)
{
    typedef typename tmv::Traits<T>::float_type FT;

    tmv::SmallMatrix<T,M,N,stor> m;
    
    for(int i=0;i<M;++i) for(int j=0;j<N;++j) m(i,j) = T(2+4*i-5*j);
    m.diag() *= T(11);
    m /= T(7);
    if (M > 1) m(1,0) = T(-2);
    if (M > 2) m(2,0) = T(7);
    if (M > 3) m(3,0) = T(-10);
    if (M > 2 && N > 2) m(2,2) = T(30);

    tmv::SmallVector<T,N> x = 2.5*m.row(0);
    x(0) = T(2);
    if (N  > 1) x(1) = T(-10);
    if (N  > 2) x(2) = T(5);
    if (N  > 3) x(3) = T(-5);

    tmv::SmallMatrix<T,N,M> minv = m.inverse();
    T eps = EPS * Norm(m) * Norm(minv);
    tmv::SmallVector<T,M> b = m * x;
    tmv::SmallVector<T,N> x2 = b/m;
    if (showacc) {
        std::cout<<"m = "<<m<<std::endl;
        std::cout<<"minv = "<<minv<<std::endl;
        std::cout<<"eps = "<<eps<<std::endl;
        std::cout<<"b = m*x = "<<b<<std::endl;
        std::cout<<"x2 = b/m = "<<x2<<std::endl;
        std::cout<<"Norm(x-x2) = "<<Norm(x-x2)<<std::endl;
        std::cout<<"eps*Norm(x) = "<<eps*Norm(x)<<std::endl;
    }
    Assert(Norm(x2-x) < eps*Norm(x),label+" NonSquare exact b/m");

    tmv::SmallVector<T,M> b2 = x%m;
    x2 = b2*m;
    Assert(Norm(x2-x) < eps*Norm(x),label+" NonSquare x%m");

    tmv::SmallMatrix<T,N,N> id = minv*m;
    T kappa = Norm(m) * Norm(minv);
    eps *= kappa;
    tmv::SmallMatrix<T,M,M> nonid = m*minv;
    if (showacc) {
        std::cout<<"minv = "<<minv<<std::endl;
        std::cout<<"minv*m = "<<id<<std::endl;
        std::cout<<"m*minv = "<<nonid<<std::endl;
        std::cout<<"Norm(id-I) = "<<Norm(id-T(1))<<std::endl;
        std::cout<<"(m*minv)T = "<<nonid.transpose()<<std::endl;
        std::cout<<"(m*minv) - (m*minv)T = "<<
            (nonid-nonid.transpose())<<std::endl;
    }
    Assert(Norm(id-T(1)) < eps,label+" NonSquare Inverse");
    Assert(Norm(nonid-nonid.transpose()) < eps,label+" NonSquare Pseudo-Inverse");

    tmv::SmallMatrix<T,N,N> mata;
    m.makeInverseATA(mata);
    tmv::SmallMatrix<T,N,N> mtm = m.transpose()*m;
    tmv::SmallMatrix<T,N,N> mtminv = mtm.inverse();
    if (showacc) {
        std::cout<<"mata = "<<mata<<std::endl;
        std::cout<<"mtminv = "<<mtminv<<std::endl;
        std::cout<<"Norm(diff) = "<<Norm(mata-mtminv)<<std::endl;
        std::cout<<"  cf. "<<eps<<" * "<<Norm(mata)<<" = "<<eps*Norm(mata)<<std::endl;
        std::cout<<"kappa = "<<kappa<<std::endl;
    }
    Assert(Norm(mata-mtminv) < eps*Norm(mata),label+" NonSquare InverseATA");

    tmv::SmallMatrix<std::complex<T>,M,N,stor> c = m * std::complex<T>(1,2);
    c.diag() *= std::complex<T>(6,-9);
    if (M > 2 && N > 3) c(2,3) += std::complex<T>(2,3);
    if (M > 1) c(1,0) *= std::complex<T>(0,2);
    if (N > 1) c.col(1) *= std::complex<T>(-1,3);
    if (M > 3) 
        c.row(3) += tmv::SmallVector<std::complex<T>,N>(std::complex<T>(1,9));

    tmv::SmallVector<std::complex<T>,N> y = std::complex<T>(2.5,1.5)*c.row(0);
    y(0) = std::complex<T>(2,9);
    if (N > 1) y(1) = std::complex<T>(-10,4);
    if (N > 2) y(2) = std::complex<T>(5,-1);
    if (N > 3) y(3) = std::complex<T>(-5,-2);

    T ceps = EPS * Norm(c) * Norm(c.inverse());
    tmv::SmallVector<std::complex<T>,M> e = c * y;
    tmv::SmallVector<std::complex<T>,N> y2 = e/c;

    Assert(Norm(y2-y) < ceps*Norm(y),label+" NonSquare exact e/c");

    tmv::SmallVector<std::complex<T>,M> e2 = y%c;
    y2 = e2*c;
    Assert(Norm(y2-y) < ceps*Norm(y),label+" NonSquare e%c");

    tmv::SmallMatrix<std::complex<T>,N,M> cinv = c.inverse();
    T ckappa = Norm(c) * Norm(cinv);
    ceps *= ckappa;
    tmv::SmallMatrix<std::complex<T>,N,N> cid = cinv*c;
    Assert(Norm(cid-T(1)) < ceps,label+" NonSquare CInverse");
    tmv::SmallMatrix<std::complex<T>,M,M > cnonid = c*cinv;
    Assert(Norm(cnonid-cnonid.adjoint()) < ceps,label+" NonSquare CPseudo-Inverse");

    tmv::SmallMatrix<std::complex<T>,N,N> cata;
    c.makeInverseATA(cata);
    tmv::SmallMatrix<std::complex<T>,N,N> ctc = c.adjoint()*c;
    tmv::SmallMatrix<std::complex<T>,N,N> ctcinv = ctc.inverse();
    Assert(Norm(cata-ctcinv) < ceps*Norm(cata),label+" NonSquare CInverseATA");

    // Test short matrix (M < N)
    tmv::SmallMatrix<T,N,M,stor> ms = m.transpose();

    b = x * ms;
    x2 = b%ms;
    Assert(Norm(x2-x) < eps*Norm(x),label+" NonSquare exact b%ms");

    b2 = x/ms;
    x2 = ms*b2;
    Assert(Norm(x2-x) < eps*Norm(x),label+" NonSquare x/ms");
}


template <class T, tmv::StorageType stor, int M, int N> 
static void TestSmallNonSquareDiv_Arith(std::string label)
{
    tmv::SmallMatrix<T,M,N,stor> m;

    for(int i=0;i<M;++i) for(int j=0;j<N;++j) m(i,j) = T(2+4*i-5*j);
    m.diag() *= T(11);
    m /= T(7);
    if (M > 1) m(1,0) = -2;
    if (M > 2) m(2,0) = 7;
    if (M > 3) m(3,0) = -10;
    if (M > 2 && N > 2) m(2,2) = 30;

    tmv::SmallMatrix<std::complex<T>,M,N,stor> c = m;

    tmv::SmallMatrix<T,M,N,stor> a1 = m;
    tmv::SmallMatrix<std::complex<T>,M,N,stor> c1 = a1 * std::complex<T>(1,2);
    c1.diag().addToAll(std::complex<T>(3,1));

    tmv::SmallVector<T,M> b = m.col(0);
    tmv::SmallVector<std::complex<T>,M> e = c.col(0);
    tmv::SmallVector<T,N> x = m.row(0);
    tmv::SmallVector<std::complex<T>,N> y = c.row(0);

    tmv::SmallMatrix<T,N,M,stor> a2 = m.transpose();
    tmv::SmallMatrix<std::complex<T>,N,M,stor> c2 = a2 * std::complex<T>(-3,4);
    c2.diag().addToAll(std::complex<T>(-5,8));
    c2.row(0).addToAll(std::complex<T>(-2,-11));

    tmv::SmallMatrix<T,M,M,stor> a3 = a1*a1.transpose();
    tmv::SmallMatrix<std::complex<T>,M,M,tmv::ColMajor> c3 = c1*a1.transpose();

    tmv::SmallMatrix<T,N,N,stor> a4 = a1.transpose()*a1;
    tmv::SmallMatrix<std::complex<T>,N,N,tmv::ColMajor> c4 = a1.transpose()*c1;

    const int Q = 3*(M+N);

    tmv::SmallMatrix<T,N,Q,stor> a5;
    for(int i=0;i<N;++i) for(int j=0;j<Q;++j) a5(i,j) = T(1-3*i+2*j);
    a5.subMatrix(0,N,2,2+M) += a1.transpose();
    tmv::SmallMatrix<std::complex<T>,N,Q,stor> c5 = a5*std::complex<T>(1,2);
    c5.subMatrix(0,N,2,2+M) += c1.transpose();

    tmv::SmallMatrix<T,Q,N,stor> a6 = a5.transpose();
    a6.subMatrix(1,1+M,0,N) += a1;
    tmv::SmallMatrix<std::complex<T>,Q,N,stor> c6 = c5.adjoint();
    c6.subMatrix(1,1+M,0,N) += c1;

    tmv::SmallMatrix<T,M,Q,stor> a7;
    for(int i=0;i<M;++i) for(int j=0;j<Q;++j) a7(i,j) = T(1-3*i+2*j);
    a7.subMatrix(0,M,2,2+N) += a1;
    tmv::SmallMatrix<std::complex<T>,M,Q,stor> c7 = a7*std::complex<T>(1,2);
    c7.subMatrix(0,M,2,2+N) += c1;

    tmv::SmallMatrix<T,Q,M,stor> a8 = a7.transpose();
    a8.subMatrix(1,1+N,0,M) += a1.transpose();
    tmv::SmallMatrix<std::complex<T>,Q,M,stor> c8 = c7.adjoint();
    c8.subMatrix(1,1+N,0,M) += c1.transpose();

    TestMatrixDivArith3a<T>(tmv::QR,a1,c1,label+" NonSquare"); 
    TestMatrixDivArith3a<T>(tmv::QR,a2,c2,label+" NonSquare"); 
    TestMatrixDivArith3d<T>(tmv::QR,a1,b,x,c1,e,y,label+" V/NonSquare"); 
    TestMatrixDivArith3d<T>(tmv::QR,a2,x,b,c2,y,e,label+" V/NonSquare"); 
    TestMatrixDivArith3e<T>(tmv::QR,a1,x,b,c1,y,e,label+" V%NonSquare"); 
    TestMatrixDivArith3e<T>(tmv::QR,a2,b,x,c2,e,y,label+" V%NonSquare"); 
    TestMatrixDivArith3b<T>(tmv::QR,a1,a3,a2,c1,c3,c2,label+" Square/NonSquare"); 
    TestMatrixDivArith3b<T>(tmv::QR,a2,a4,a1,c2,c4,c1,label+" Square/NonSquare"); 
    TestMatrixDivArith3c<T>(tmv::QR,a1,a4,a2,c1,c4,c2,label+" Square%NonSquare"); 
    TestMatrixDivArith3c<T>(tmv::QR,a2,a3,a1,c2,c3,c1,label+" Square%NonSquare"); 
    TestMatrixDivArith3b<T>(tmv::QR,a1,a7,a5,c1,c7,c5,label+" NonSquare/NonSquare");
    TestMatrixDivArith3b<T>(tmv::QR,a2,a5,a7,c2,c5,c7,label+" NonSquare/NonSquare");
    TestMatrixDivArith3c<T>(tmv::QR,a1,a6,a8,c1,c6,c8,label+" NonSquare%NonSquare");
    TestMatrixDivArith3c<T>(tmv::QR,a2,a8,a6,c2,c8,c6,label+" NonSquare%NonSquare");
}

template <class T, tmv::StorageType stor, int M, int N> 
static void TestSmallNonSquareDiv(std::string label)
{
    TestSmallNonSquareDiv_Basic<T,stor,M,N>(label);
    TestSmallNonSquareDiv_Arith<T,stor,M,N>(label);
}
