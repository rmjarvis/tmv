
#include "TMV_Test.h"
#include "TMV_Test_3.h"
#include "TMV.h"
#include "TMV_TestMatrixDivArith.h"

template <class T> void TestSmallSquareDiv_a();
template <class T> void TestSmallSquareDiv_b();
template <class T> void TestSmallSquareDiv_c();
template <class T> void TestSmallSquareDiv_a();
template <class T> void TestSmallSquareDiv_d();
template <class T> void TestSmallSquareDiv_e();
template <class T> void TestSmallSquareDiv_f();
template <class T> void TestSmallSquareDiv_g();

template <class T, tmv::StorageType stor, int N> 
static void TestSmallSquareDiv_Basic(std::string label)
{
    typedef typename tmv::Traits<T>::float_type FT;

    tmv::SmallMatrix<T,N,N,stor> m;

    for(int i=0;i<N;++i) for(int j=0;j<N;++j) m(i,j) = T(2+4*i-5*j);
    m.diag() *= T(11);
    m /= T(7);
    if (N > 1) m(1,0) = -2;
    if (N > 2) m(2,0) = 7;
    if (N > 3) m(3,0) = -10;
    if (N > 2) m(2,2) = 30;

    tmv::SmallVector<T,N> b = 2.5*m.row(0);
    b(0) = 2;
    if (N > 1) b(1) = -10;
    if (N > 2) b(2) = 5;
    if (N > 3) b(3) = -5;

    tmv::SmallMatrix<T,N,N> minv = m.inverse();
    T kappa = Norm(m) * Norm(minv);
    FT eps = EPS * kappa;

    tmv::SmallMatrix<T,N,N> id = m*minv;
    if (showacc) {
        std::cout<<"minv = "<<minv<<std::endl;
        std::cout<<"m*minv = "<<id<<std::endl;
        std::cout<<"minv*m = "<<minv*m<<std::endl;
        std::cout<<"Norm(id-I) = "<<Norm(id-T(1))<<std::endl;
        std::cout<<"eps = "<<eps<<std::endl;
    }
    Assert(Norm(id-T(1)) < eps,label+" Square Inverse");

    tmv::SmallVector<T,N> x = b/m;
    tmv::SmallVector<T,N> b2 = m*x;
    if (showacc) {
        std::cout<<"b = "<<b<<std::endl;
        std::cout<<"x = b/m = "<<x<<std::endl;
        std::cout<<"b2 = "<<b2<<std::endl;
        std::cout<<"Norm(b-b2) = "<<Norm(b-b2)<<std::endl;
        std::cout<<"eps*Norm(b) = "<<eps*Norm(b)<<std::endl;
    }
    Assert(Norm(b2-b) < eps*Norm(b),label+" Square b/m");

    x = b%m;
    b2 = x*m;
    if (showacc) {
        std::cout<<"b = "<<b<<std::endl;
        std::cout<<"x = b%m = "<<x<<std::endl;
        std::cout<<"b2 = "<<b2<<std::endl;
        std::cout<<"Norm(b-b2) = "<<Norm(b-b2)<<std::endl;
        std::cout<<"eps*Norm(b) = "<<eps*Norm(b)<<std::endl;
    }
    Assert(Norm(b2-b) < eps*Norm(b),label+" Square b%m");

    tmv::SmallMatrix<T,N,N> mtm = m.adjoint() * m;
    tmv::SmallMatrix<T,N,N> mata;
    m.makeInverseATA(mata);
    if (showacc) {
        std::cout<<"mtm = "<<mtm<<std::endl;
        std::cout<<"mtm.inv = "<<mtm.inverse()<<std::endl;
        std::cout<<"m.invata = "<<mata<<std::endl;
        std::cout<<"minv*minvt = "<<(minv*minv.adjoint())<<std::endl;
        std::cout<<"Norm(diff) = "<<Norm(mata-mtm.inverse())<<std::endl;
        std::cout<<"c.f. eps*Norm(mata) = "<<eps*kappa<<" * "<<Norm(mata)<<" = "<<eps*Norm(mata)<<std::endl;
    }
    Assert(Norm(mata-mtm.inverse()) < eps*kappa*Norm(mata),label+" Square inverseATA");

    T mdet = Det(tmv::Matrix<T>(m));
    if (T(1)/tmv::TMV_ABS(mdet) > T(0)) {
        // If det is inf, then skip these tests
        // (Happens for the larger sized matrices we test on.)
        if (showacc) {
            std::cout<<"Det(m) = "<<Det(m)<<std::endl;
            std::cout<<"mdet = "<<mdet<<std::endl;
            std::cout<<"abs(det-mdet) = "<<std::abs(Det(m)-mdet);
            std::cout<<"  EPS*abs(mdet) = "<<eps*std::abs(mdet)<<std::endl;
            std::cout<<"abs(logdet-log(mdet)) = "<<
                std::abs(m.logDet()-std::log(std::abs(mdet)))<<std::endl;
        }
        Assert(std::abs(Det(m)-mdet) < eps*std::abs(mdet),label+" Square Det");
        T sdet;
        Assert(std::abs(m.logDet(&sdet)-std::log(std::abs(mdet))) < eps,
               "Square LogDet");
        Assert(std::abs(sdet-mdet/std::abs(mdet)) < eps,label+" Square LogDet - sign");
    }

    tmv::SmallMatrix<std::complex<T>,N,N,stor> c = m*std::complex<T>(2,3);
    c.diag() *= std::complex<T>(6,-9);
    if (N > 3) c(2,3) += std::complex<T>(2,3);
    if (N > 1) c(1,0) *= std::complex<T>(0,2);
    if (N > 1) c.col(1) *= std::complex<T>(-1,3);
    if (N > 3) c.row(3) += 
        tmv::SmallVector<std::complex<T>,N>(std::complex<T>(1,9));

    tmv::SmallMatrix<std::complex<T>,N,N> cinv = c.inverse();
    FT ckappa = Norm(c)*Norm(cinv);
    FT ceps = EPS * ckappa;

    std::complex<T> cdet = Det(tmv::Matrix<std::complex<T> >(c));
    if (T(1)/tmv::TMV_ABS(cdet) > T(0)) {
        // If det is inf or (nan,nan) then skip these tests
        if (showacc) {
            std::cout<<"cdet = "<<cdet<<std::endl;
            std::cout<<"Det(c) = "<<Det(c)<<std::endl;
            std::cout<<"abs(det-cdet) = "<<std::abs(Det(c)-cdet);
            std::cout<<"  EPS*abs(cdet) = "<<ceps*std::abs(cdet)<<std::endl;
            std::cout<<"abs(logdet-log(cdet)) = "<<
                std::abs(c.logDet()-std::log(cdet))<<std::endl;
        }
        Assert(std::abs(Det(c)-cdet) < ceps*std::abs(cdet),label+" Square CDet");
        std::complex<T> csdet;
        Assert(std::abs(c.logDet(&csdet)-std::log(std::abs(cdet))) < eps,
               "Square CLogDet");
        Assert(std::abs(csdet-cdet/std::abs(cdet)) < eps,label+" Square CLogDet - sign");
    }

    tmv::SmallMatrix<std::complex<T>,N,N> cid = c*cinv;
    Assert(Norm(cid-T(1)) < ceps,label+" Square CInverse");

    tmv::SmallMatrix<std::complex<T>,N,N> ctc = c.adjoint() * c;
    tmv::SmallMatrix<std::complex<T>,N,N> cata;
    tmv::SmallMatrix<std::complex<T>,N,N> ctcinv = ctc.inverse();
    c.makeInverseATA(cata);
    Assert(Norm(cata-ctcinv) < ceps*ckappa*Norm(cata),label+" Square CInverseATA");

    tmv::SmallVector<std::complex<T>,N> e = std::complex<T>(2.5,1.5)*c.row(0);
    e = b*std::complex<T>(1,2);
    if (N > 1) e(1) += std::complex<T>(-1,5);
    if (N > 2) e(2) -= std::complex<T>(-1,5);

    // test real / complex
    tmv::SmallVector<std::complex<T>,N> y = b/c;
    tmv::SmallVector<std::complex<T>,N> b3 = c*y;
    Assert(Norm(b3-b) < ceps*Norm(b),label+" Square b/c");
    y = b%c;
    b3 = y*c;
    Assert(Norm(b3-b) < ceps*Norm(b),label+" Square b%c");

    // test complex / real
    y = e/m;
    b3 = m*y;
    Assert(Norm(b3-e) < eps*Norm(e),label+" Square e/m");
    y = e%m;
    b3 = y*m;
    Assert(Norm(b3-e) < eps*Norm(e),label+" Square e%m");

    // test complex / complex
    y = e/c;
    b3 = c*y;
    Assert(Norm(b3-e) < ceps*Norm(e),label+" Square e/c");
    y = e%c;
    b3 = y*c;
    Assert(Norm(b3-e) < ceps*Norm(e),label+" Square e%c");
}

template <class T, tmv::StorageType stor, int N> 
static void TestSmallSquareDiv_Arith(std::string label)
{
    tmv::SmallMatrix<T,N,N,stor> m;

    for(int i=0;i<N;++i) for(int j=0;j<N;++j) m(i,j) = T(2+4*i-5*j);
    m.diag() *= T(11);
    m /= T(7);
    if (N > 1) m(1,0) = -2;
    if (N > 2) m(2,0) = 7;
    if (N > 3) m(3,0) = -10;
    if (N > 2) m(2,2) = 30;

    tmv::SmallMatrix<std::complex<T>,N,N,stor> c = m;

    tmv::SmallMatrix<T,N,N,stor> a1 = m;
    tmv::SmallMatrix<std::complex<T>,N,N,stor> c1 = a1 * std::complex<T>(1,2);
    c1.diag().addToAll(std::complex<T>(3,1));

    tmv::SmallVector<T,N> b = m.row(0);
    tmv::SmallVector<std::complex<T>,N> e = c.row(0);
    tmv::SmallVector<T,N> x = m.col(0);
    tmv::SmallVector<std::complex<T>,N> y = c.col(0);

    tmv::SmallMatrix<T,N,N,tmv::ColMajor> a2a = m.transpose();
    tmv::SmallMatrix<std::complex<T>,N,N,tmv::ColMajor> c2a = a2a * std::complex<T>(-3,4);
    c2a.diag().addToAll(std::complex<T>(-5,8));
    c2a.row(0).addToAll(std::complex<T>(-2,-11));
    tmv::SmallMatrix<T,N,N,tmv::RowMajor> a2b = a2a;
    tmv::SmallMatrix<std::complex<T>,N,N,tmv::RowMajor> c2b = c2a;

    tmv::SmallMatrix<T,N,N,stor> a3 = a1;
    tmv::SmallMatrix<std::complex<T>,N,N,stor> c3 = c1;

    tmv::SmallMatrix<T,N,3*N,stor> a4;
    for(int i=0;i<N;++i) for(int j=0;j<3*N;++j) a4(i,j) = T(1-3*i+2*j);
    a4.subMatrix(0,N,2,2+N) += a1;
    tmv::SmallMatrix<std::complex<T>,N,3*N,stor> c4 = a4*std::complex<T>(1,2);
    c4.subMatrix(0,N,2,2+N) += c1;

    tmv::SmallMatrix<T,3*N,N,stor> a5 = a4.transpose();
    a5.subMatrix(1,1+N,0,N) -= a1;
    tmv::SmallMatrix<std::complex<T>,3*N,N,stor> c5 = c4.adjoint();
    c5.subMatrix(1,1+N,0,N) -= c1;

    tmv::SmallMatrix<T,N,3*N,stor> a4b = a4;
    tmv::SmallMatrix<std::complex<T>,N,3*N,stor> c4b = c4;

    tmv::SmallMatrix<T,3*N,N,stor> a5b = a5;
    tmv::SmallMatrix<std::complex<T>,3*N,N,stor> c5b = c5;

    TestMatrixDivArith3a<T>(tmv::LU,a1,c1,label+" Square"); 
    TestMatrixDivArith3d<T>(tmv::LU,a1,b,x,c1,e,y,label+" V/Square"); 
    TestMatrixDivArith3e<T>(tmv::LU,a1,b,x,c1,e,y,label+" V/Square"); 
    TestMatrixDivArith3b<T>(tmv::LU,a1,a2a,a3,c1,c2a,c3,label+" Square/Square"); 
    TestMatrixDivArith3b<T>(tmv::LU,a1,a2b,a3,c1,c2b,c3,label+" Square/Square"); 
    TestMatrixDivArith3c<T>(tmv::LU,a1,a2a,a3,c1,c2a,c3,label+" Square/Square"); 
    TestMatrixDivArith3c<T>(tmv::LU,a1,a2b,a3,c1,c2b,c3,label+" Square/Square"); 
    TestMatrixDivArith3b<T>(tmv::LU,a1,a4,a4b,c1,c4,c4b,label+" NonSquare/Square");
    TestMatrixDivArith3c<T>(tmv::LU,a1,a5,a5b,c1,c5,c5b,label+" NonSquare/Square");
}

template <class T, tmv::StorageType stor, int N> 
static void TestSmallSquareDiv(std::string label)
{
    TestSmallSquareDiv_Basic<T,stor,N>(label);
    TestSmallSquareDiv_Arith<T,stor,N>(label);
}
