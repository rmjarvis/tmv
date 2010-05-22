

#include "TMV_Test.h"
#include "TMV_Test3.h"
#include "TMV.h"
#include "TMV_TestMatrixDivArith.h"

template <class T, tmv::StorageType stor, int N> 
static void TestSmallSquareDiv_Basic()
{
    tmv::SmallMatrix<T,N,N,stor> m;

    for(int i=0;i<N;++i) for(int j=0;j<N;++j) m(i,j) = T(2+4*i-5*j);
    m.diag() *= T(11);
    m /= T(7);
    if (N > 1) m(1,0) = -2;
    if (N > 2) m(2,0) = 7;
    if (N > 3) m(3,0) = -10;
    if (N > 2) m(2,2) = 30;

    tmv::SmallVector<T,N> b(T(7));
    b(0) = 2;
    if (N > 1) b(1) = -10;
    if (N > 2) b(2) = 5;
    if (N > 3) b(3) = -5;

    if (showstartdone) {
        std::cout<<"Start TestSmallSquareDiv\n";
        std::cout<<"m = "<<TMV_Text(m)<<" "<<m<<std::endl;
    }

    tmv::SmallMatrix<T,N,N> minv = m.inverse();
    T kappa = Norm(m) * Norm(minv);
    T eps = EPS * kappa;

    tmv::SmallVector<T,N> x = b/m;
    tmv::SmallVector<T,N> b2 = m*x;
    if (showacc) {
        std::cout<<"b = "<<b<<std::endl;
        std::cout<<"x = b/m = "<<x<<std::endl;
        std::cout<<"b2 = "<<b2<<std::endl;
        std::cout<<"Norm(b-b2) = "<<Norm(b-b2)<<std::endl;
        std::cout<<"eps*Norm(b) = "<<eps*Norm(b)<<std::endl;
    }
    Assert(Norm(b2-b) < eps*Norm(b),"Square b/m");

    x = b%m;
    b2 = x*m;
    if (showacc) {
        std::cout<<"b = "<<b<<std::endl;
        std::cout<<"x = b%m = "<<x<<std::endl;
        std::cout<<"b2 = "<<b2<<std::endl;
        std::cout<<"Norm(b-b2) = "<<Norm(b-b2)<<std::endl;
        std::cout<<"eps*Norm(b) = "<<eps*Norm(b)<<std::endl;
    }
    Assert(Norm(b2-b) < eps*Norm(b),"Square b%m");

    tmv::SmallMatrix<T,N,N> id = m*minv;
    if (showacc) {
        std::cout<<"minv = "<<minv<<std::endl;
        std::cout<<"m*minv = "<<id<<std::endl;
        std::cout<<"minv*m = "<<minv*m<<std::endl;
        std::cout<<"Norm(id-I) = "<<Norm(id-T(1))<<std::endl;
        std::cout<<"eps = "<<eps<<std::endl;
    }
    Assert(Norm(id-T(1)) < eps,"Square Inverse");

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
    Assert(Norm(mata-mtm.inverse()) < eps*kappa*Norm(mata),"Square inverseATA");

    T mdet = Det(tmv::Matrix<T>(m));
    if (T(1)/mdet != T(0)) {
        // If det is inf, then skip these tests
        // (Happens for the larger sized matrices we test on.
        if (showacc) {
            std::cout<<"Det(m) = "<<Det(m)<<std::endl;
            std::cout<<"mdet = "<<mdet<<std::endl;
            std::cout<<"abs(det-mdet) = "<<std::abs(Det(m)-mdet);
            std::cout<<"  EPS*abs(mdet) = "<<eps*std::abs(mdet)<<std::endl;
            std::cout<<"abs(logdet-log(mdet)) = "<<
                std::abs(m.logDet()-std::log(std::abs(mdet)))<<std::endl;
        }
        Assert(std::abs(Det(m)-mdet) < eps*std::abs(mdet),"Square Det");
        T sdet;
        Assert(std::abs(m.logDet(&sdet)-std::log(std::abs(mdet))) < eps,
               "Square LogDet");
        Assert(std::abs(sdet-mdet/std::abs(mdet)) < eps,"Square LogDet - sign");
    }

    tmv::SmallMatrix<std::complex<T>,N,N,stor> c = m*std::complex<T>(2,3);
    c.diag() *= std::complex<T>(6,-9);
    if (N > 3) c(2,3) += std::complex<T>(2,3);
    if (N > 1) c(1,0) *= std::complex<T>(0,2);
    if (N > 1) c.col(1) *= std::complex<T>(-1,3);
    if (N > 3) c.row(3) += 
        tmv::SmallVector<std::complex<T>,N>(std::complex<T>(1,9));

    tmv::SmallMatrix<std::complex<T>,N,N> cinv = c.inverse();
    T ckappa = Norm(c)*Norm(cinv);
    T ceps = EPS * ckappa;

    std::complex<T> cdet = Det(tmv::Matrix<std::complex<T> >(c));
    if (T(1)/cdet != T(0)) {
        // If det is inf, then skip these tests
        if (showacc) {
            std::cout<<"cdet = "<<cdet<<std::endl;
            std::cout<<"Det(c) = "<<Det(c)<<std::endl;
            std::cout<<"abs(det-cdet) = "<<std::abs(Det(c)-cdet);
            std::cout<<"  EPS*abs(cdet) = "<<ceps*std::abs(cdet)<<std::endl;
            std::cout<<"abs(logdet-log(cdet)) = "<<
                std::abs(c.logDet()-std::log(cdet))<<std::endl;
        }
        Assert(std::abs(Det(c)-cdet) < ceps*std::abs(cdet),"Square CDet");
        std::complex<T> csdet;
        Assert(std::abs(c.logDet(&csdet)-std::log(std::abs(cdet))) < eps,
               "Square CLogDet");
        Assert(std::abs(csdet-cdet/std::abs(cdet)) < eps,"Square CLogDet - sign");
    }

    tmv::SmallMatrix<std::complex<T>,N,N> cid = c*cinv;
    Assert(Norm(cid-T(1)) < ceps,"Square CInverse");

    tmv::SmallMatrix<std::complex<T>,N,N> ctc = c.adjoint() * c;
    tmv::SmallMatrix<std::complex<T>,N,N> cata;
    tmv::SmallMatrix<std::complex<T>,N,N> ctcinv = ctc.inverse();
    c.makeInverseATA(cata);
    Assert(Norm(cata-ctcinv) < ceps*ckappa*Norm(cata),"Square CInverseATA");

    tmv::SmallVector<std::complex<T>,N> e;
    e = b*std::complex<T>(1,2);
    if (N > 1) e(1) += std::complex<T>(-1,5);
    if (N > 2) e(2) -= std::complex<T>(-1,5);

    // test real / complex
    tmv::SmallVector<std::complex<T>,N> y = b/c;
    tmv::SmallVector<std::complex<T>,N> b3 = c*y;
    Assert(Norm(b3-b) < ceps*Norm(b),"Square b/c");
    y = b%c;
    b3 = y*c;
    Assert(Norm(b3-b) < ceps*Norm(b),"Square b%c");

    // test complex / real
    y = e/m;
    b3 = m*y;
    Assert(Norm(b3-e) < eps*Norm(e),"Square e/m");
    y = e%m;
    b3 = y*m;
    Assert(Norm(b3-e) < eps*Norm(e),"Square e%m");

    // test complex / complex
    y = e/c;
    b3 = c*y;
    Assert(Norm(b3-e) < ceps*Norm(e),"Square e/c");
    y = e%c;
    b3 = y*c;
    Assert(Norm(b3-e) < ceps*Norm(e),"Square e%c");
}

template <class T, tmv::StorageType stor, int N> 
static void TestSmallSquareDiv_Arith()
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
    tmv::SmallVector<T,N> x;
    tmv::SmallVector<std::complex<T>,N> y;

    tmv::SmallMatrix<T,N,N,tmv::ColMajor> a2a = m.transpose();
    tmv::SmallMatrix<std::complex<T>,N,N,tmv::ColMajor> c2a = a2a * std::complex<T>(-3,4);
    c2a.diag().addToAll(std::complex<T>(-5,8));
    c2a.row(0).addToAll(std::complex<T>(-2,-11));
    tmv::SmallMatrix<T,N,N,tmv::RowMajor> a2b = a2a;
    tmv::SmallMatrix<std::complex<T>,N,N,tmv::RowMajor> c2b = c2a;

    tmv::SmallMatrix<T,N,N,stor> a0;
    tmv::SmallMatrix<std::complex<T>,N,N,stor> c0;
    tmv::SmallMatrix<T,N,N,stor> a3;
    tmv::SmallMatrix<std::complex<T>,N,N,stor> c3;

    tmv::SmallMatrix<T,N,3*N,stor> a4;
    for(int i=0;i<N;++i) for(int j=0;j<3*N;++j) a4(i,j) = T(1-3*i+2*j);
    a4.subMatrix(0,N,2,2+N) += a1;
    tmv::SmallMatrix<std::complex<T>,N,3*N,stor> c4 = a4*std::complex<T>(1,2);
    c4.subMatrix(0,N,2,2+N) += c1;

    tmv::SmallMatrix<T,3*N,N,stor> a5 = a4.transpose();
    a5.subMatrix(1,1+N,0,N) -= a1;
    tmv::SmallMatrix<std::complex<T>,3*N,N,stor> c5 = c4.adjoint();
    c5.subMatrix(1,1+N,0,N) -= c1;

    tmv::SmallMatrix<T,N,3*N,stor> a4x;
    tmv::SmallMatrix<T,N,3*N,stor> a4b;
    tmv::SmallMatrix<std::complex<T>,N,3*N,stor> c4x;
    tmv::SmallMatrix<std::complex<T>,N,3*N,stor> c4b;

    tmv::SmallMatrix<T,3*N,N,stor> a5x;
    tmv::SmallMatrix<T,3*N,N,stor> a5b;
    tmv::SmallMatrix<std::complex<T>,3*N,N,stor> c5x;
    tmv::SmallMatrix<std::complex<T>,3*N,N,stor> c5b;

    TestMatrixDivArith3a<T>(tmv::LU,a1,c1,"Square"); 
    TestMatrixDivArith3d<T>(tmv::LU,a1,b,x,c1,e,y,"V/Square"); 
    TestMatrixDivArith3e<T>(tmv::LU,a1,b,x,c1,e,y,"V/Square"); 
    TestMatrixDivArith3b<T>(tmv::LU,a0,c0,a1,a2a,a3,c1,c2a,c3,"Square/Square"); 
    TestMatrixDivArith3b<T>(tmv::LU,a0,c0,a1,a2b,a3,c1,c2b,c3,"Square/Square"); 
    TestMatrixDivArith3c<T>(tmv::LU,a0,c0,a1,a2a,a3,c1,c2a,c3,"Square/Square"); 
    TestMatrixDivArith3c<T>(tmv::LU,a0,c0,a1,a2b,a3,c1,c2b,c3,"Square/Square"); 
    TestMatrixDivArith3b<T>(tmv::LU,a4x,c4x,a1,a4,a4b,c1,c4,c4b,"NonSquare/Square");
    TestMatrixDivArith3c<T>(tmv::LU,a5x,c5x,a1,a5,a5b,c1,c5,c5b,"NonSquare/Square");

#if XTEST & 32
    tmv::SmallMatrix<T,N,N,stor,tmv::FortranStyle> a1f = a1;
    tmv::SmallMatrix<std::complex<T>,N,N,stor,tmv::FortranStyle> c1f = c1;

    tmv::SmallVector<T,N,tmv::FortranStyle> bf = b;
    tmv::SmallVector<std::complex<T>,N,tmv::FortranStyle> ef = e;
    tmv::SmallVector<T,N,tmv::FortranStyle> xf = x;
    tmv::SmallVector<std::complex<T>,N,tmv::FortranStyle> yf = y;

    tmv::SmallMatrix<T,N,N,tmv::ColMajor,tmv::FortranStyle> a2fa = a2a;
    tmv::SmallMatrix<std::complex<T>,N,N,tmv::ColMajor,tmv::FortranStyle> c2fa = c2a;
    tmv::SmallMatrix<T,N,N,tmv::RowMajor,tmv::FortranStyle> a2fb = a2b;
    tmv::SmallMatrix<std::complex<T>,N,N,tmv::RowMajor,tmv::FortranStyle> c2fb = c2b;
    tmv::SmallMatrix<T,N,N,stor,tmv::FortranStyle> a3f = a3;
    tmv::SmallMatrix<std::complex<T>,N,N,stor,tmv::FortranStyle> c3f = c3;

    tmv::SmallMatrix<T,N,3*N,stor,tmv::FortranStyle> a4f = a4;
    tmv::SmallMatrix<std::complex<T>,N,3*N,stor,tmv::FortranStyle> c4f = c4;
    tmv::SmallMatrix<T,N,3*N,stor,tmv::FortranStyle> a4fb = a4;
    tmv::SmallMatrix<std::complex<T>,N,3*N,stor,tmv::FortranStyle> c4fb = c4;

    tmv::SmallMatrix<T,3*N,N,stor,tmv::FortranStyle> a5f = a5;
    tmv::SmallMatrix<std::complex<T>,3*N,N,stor,tmv::FortranStyle> c5f = c5;
    tmv::SmallMatrix<T,3*N,N,stor,tmv::FortranStyle> a5fb = a5;
    tmv::SmallMatrix<std::complex<T>,3*N,N,stor,tmv::FortranStyle> c5fb = c5;

    TestMatrixDivArith3a<T>(tmv::LU,a1f,c1f,"Square"); 

    TestMatrixDivArith3d<T>(tmv::LU,a1f,b,x,c1f,e,y,"V/Square"); 
    TestMatrixDivArith3d<T>(tmv::LU,a1f,bf,x,c1f,ef,y,"V/Square"); 
    TestMatrixDivArith3d<T>(tmv::LU,a1f,bf,xf,c1f,ef,yf,"V/Square"); 

    TestMatrixDivArith3e<T>(tmv::LU,a1f,b,x,c1f,e,y,"V/Square"); 
    TestMatrixDivArith3e<T>(tmv::LU,a1f,bf,x,c1f,ef,y,"V/Square"); 
    TestMatrixDivArith3e<T>(tmv::LU,a1f,bf,xf,c1f,ef,yf,"V/Square"); 

    TestMatrixDivArith3b<T>(tmv::LU,a0,c0,a1f,a2a,a3,c1f,c2a,c3,"Square/Square"); 
    TestMatrixDivArith3b<T>(tmv::LU,a0,c0,a1f,a2fa,a3,c1f,c2fa,c3,"Square/Square"); 
    TestMatrixDivArith3b<T>(tmv::LU,a0,c0,a1f,a2fa,a3f,c1f,c2fa,c3f,"Square/Square"); 
    TestMatrixDivArith3b<T>(tmv::LU,a0,c0,a1f,a2b,a3,c1f,c2b,c3,"Square/Square"); 
    TestMatrixDivArith3b<T>(tmv::LU,a0,c0,a1f,a2fb,a3,c1f,c2fb,c3,"Square/Square"); 
    TestMatrixDivArith3b<T>(tmv::LU,a0,c0,a1f,a2fb,a3f,c1f,c2fb,c3f,"Square/Square"); 

    TestMatrixDivArith3c<T>(tmv::LU,a0,c0,a1f,a2a,a3,c1f,c2a,c3,"Square/Square"); 
    TestMatrixDivArith3c<T>(tmv::LU,a0,c0,a1f,a2fa,a3,c1f,c2fa,c3,"Square/Square"); 
    TestMatrixDivArith3c<T>(tmv::LU,a0,c0,a1f,a2fa,a3f,c1f,c2fa,c3f,"Square/Square"); 
    TestMatrixDivArith3c<T>(tmv::LU,a0,c0,a1f,a2b,a3,c1f,c2b,c3,"Square/Square"); 
    TestMatrixDivArith3c<T>(tmv::LU,a0,c0,a1f,a2fb,a3,c1f,c2fb,c3,"Square/Square"); 
    TestMatrixDivArith3c<T>(tmv::LU,a0,c0,a1f,a2fb,a3f,c1f,c2fb,c3f,"Square/Square"); 

    TestMatrixDivArith3b<T>(tmv::LU,a4x,c4x,a1f,a4,a4b,c1f,c4,c4b,"NonSquare/Square");
    TestMatrixDivArith3b<T>(tmv::LU,a4x,c4x,a1f,a4f,a4b,c1f,c4f,c4b,"NonSquare/Square");
    TestMatrixDivArith3b<T>(tmv::LU,a4x,c4x,a1f,a4f,a4fb,c1f,c4f,c4fb,"NonSquare/Square");

    TestMatrixDivArith3c<T>(tmv::LU,a5x,c5x,a1f,a5,a5b,c1f,c5,c5b,"NonSquare/Square");
    TestMatrixDivArith3c<T>(tmv::LU,a5x,c5x,a1f,a5f,a5b,c1f,c5f,c5b,"NonSquare/Square");
    TestMatrixDivArith3c<T>(tmv::LU,a5x,c5x,a1f,a5f,a5fb,c1f,c5f,c5fb,"NonSquare/Square");
#endif
}

template <class T, tmv::StorageType stor, int N> 
static void TestSmallSquareDiv()
{
    TestSmallSquareDiv_Basic<T,stor,N>();
    TestSmallSquareDiv_Arith<T,stor,N>();
}
