
#include "TMV.h"
#include "TMV_Small.h"
#include "TMV_Test.h"
#include "TMV_Test_3.h"

template <class T, tmv::StorageType stor, int N> 
static void TestSmallSquareDiv()
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

    tmv::SmallVector<T,N> vtemp;
    tmv::SmallMatrix<T,N,N> minv = m.inverse();
    T eps = EPS * Norm(m) * Norm(minv);

    tmv::SmallVector<T,N> x = b/m;
    tmv::SmallVector<T,N> b2 = m*x;
    if (showacc) {
        std::cout<<"eps = "<<eps<<std::endl;
        std::cout<<"b = "<<b<<std::endl;
        std::cout<<"x = b/m = "<<x<<std::endl;
        std::cout<<"b2 = "<<b2<<std::endl;
        std::cout<<"Norm(b-b2) = "<<Norm(b-b2)<<std::endl;
        std::cout<<"eps*Norm(b) = "<<eps*Norm(b)<<std::endl;
    }
    Assert(Norm(vtemp=b2-b) < eps*Norm(b),"Square b/m");

    x = b%m;
    b2 = x*m;
    if (showacc) {
        std::cout<<"b = "<<b<<std::endl;
        std::cout<<"x = b%m = "<<x<<std::endl;
        std::cout<<"b2 = "<<b2<<std::endl;
        std::cout<<"Norm(b-b2) = "<<Norm(vtemp=b-b2)<<std::endl;
        std::cout<<"eps*Norm(b) = "<<eps*Norm(b)<<std::endl;
    }
    Assert(Norm(vtemp=b2-b) < eps*Norm(b),"Square b%m");

    tmv::SmallMatrix<T,N,N> id = m*minv;
    tmv::SmallMatrix<T,N,N> mtemp;
    if (showacc) {
        std::cout<<"minv = "<<minv<<std::endl;
        std::cout<<"m*minv = "<<id<<std::endl;
        std::cout<<"minv*m = "<<(mtemp=minv*m)<<std::endl;
        std::cout<<"Norm(id-I) = "<<Norm(mtemp=id-T(1))<<std::endl;
        std::cout<<"eps = "<<eps<<std::endl;
    }
    Assert(Norm(mtemp=id-T(1)) < eps,"Square Inverse");

    tmv::SmallMatrix<T,N,N> mtm = m.adjoint() * m;
    tmv::SmallMatrix<T,N,N> mata;
    tmv::SmallMatrix<T,N,N> mtminv = mtm.inverse();
    m.makeInverseATA(mata);
    if (showacc) {
        std::cout<<"mtm = "<<mtm<<std::endl;
        std::cout<<"mtm.inv = "<<mtminv<<std::endl;
        std::cout<<"m.invata = "<<mata<<std::endl;
        std::cout<<"minv*minvt = "<<(mtemp=minv*minv.adjoint())<<std::endl;
        std::cout<<"Norm(diff) = "<<Norm(mtemp=mata-mtminv)<<std::endl;
    }
    Assert(Norm(mtemp=mata-mtminv) < eps*Norm(mata),"Square InverseATA");

    T mdet = Det(tmv::Matrix<T>(m));
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

    tmv::SmallMatrix<std::complex<T>,N,N,stor> c = m*std::complex<T>(2,3);
    c.diag() *= std::complex<T>(6,-9);
    if (N > 3) c(2,3) += std::complex<T>(2,3);
    if (N > 1) c(1,0) *= std::complex<T>(0,2);
    if (N > 1) c.col(1) *= std::complex<T>(-1,3);
    if (N > 3) c.row(3) += 
        tmv::SmallVector<std::complex<T>,N>(std::complex<T>(1,9));

    tmv::SmallMatrix<std::complex<T>,N,N> cinv = c.inverse();
    T ceps = EPS * Norm(c) * Norm(cinv);

    std::complex<T> cdet = Det(tmv::Matrix<std::complex<T> >(c));
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

    tmv::SmallMatrix<std::complex<T>,N,N> cid = c*cinv;
    tmv::SmallMatrix<std::complex<T>,N,N> ctemp;
    Assert(Norm(ctemp=cid-T(1)) < ceps,"Square CInverse");

    tmv::SmallMatrix<std::complex<T>,N,N> ctc = c.adjoint() * c;
    tmv::SmallMatrix<std::complex<T>,N,N> cata;
    tmv::SmallMatrix<std::complex<T>,N,N> ctcinv = ctc.inverse();
    c.makeInverseATA(cata);
    Assert(Norm(ctemp=cata-ctcinv) < ceps*Norm(cata),"Square CInverseATA");

    tmv::SmallVector<std::complex<T>,N> e;
    e = b*std::complex<T>(1,2);
    if (N > 1) e(1) += std::complex<T>(-1,5);
    if (N > 2) e(2) -= std::complex<T>(-1,5);

    // test real / complex
    tmv::SmallVector<std::complex<T>,N> y = b/c;
    tmv::SmallVector<std::complex<T>,N> b3 = c*y;
    tmv::SmallVector<std::complex<T>,N> cvtemp;
    Assert(Norm(cvtemp=b3-b) < ceps*Norm(b),"Square b/c");
    y = b%c;
    b3 = y*c;
    Assert(Norm(cvtemp=b3-b) < ceps*Norm(b),"Square b%c");

    // test complex / real
    y = e/m;
    b3 = m*y;
    Assert(Norm(cvtemp=b3-e) < eps*Norm(e),"Square e/m");
    y = e%m;
    b3 = y*m;
    Assert(Norm(cvtemp=b3-e) < eps*Norm(e),"Square e%m");

    // test complex / complex
    y = e/c;
    b3 = c*y;
    Assert(Norm(cvtemp=b3-e) < ceps*Norm(e),"Square e/c");
    y = e%c;
    b3 = y*c;
    Assert(Norm(cvtemp=b3-e) < ceps*Norm(e),"Square e%c");
}

template <class T, tmv::StorageType stor, int N> 
static void TestSmallNonSquareDiv()
{
    tmv::SmallMatrix<T,6,N,stor> m;
    for(int i=0;i<6;++i) for(int j=0;j<N;++j) m(i,j) = T(2+4*i-5*j);
    m.diag() *= T(11);
    m /= T(7);
    if (N > 1) m(1,0) = T(-2);
    if (N > 2) m(2,0) = T(7);
    if (N > 3) m(3,0) = T(-10);
    if (N > 2) m(2,2) = T(30);

    tmv::SmallVector<T,N> x(T(7));
    x(0) = T(2);
    if (N  > 1) x(1) = T(-10);
    if (N  > 2) x(2) = T(5);
    if (N  > 3) x(3) = T(-5);

    if (showstartdone) {
        std::cout<<"Start TestSmallNonSquareDiv\n";
        std::cout<<"m = "<<TMV_Text(m)<<" "<<m<<std::endl;
    }

    tmv::SmallMatrix<T,N,6> mn6temp;

    T eps = EPS * Norm(m) * Norm(mn6temp=m.inverse());
    tmv::SmallVector<T,6> b = m * x;
    tmv::SmallVector<T,N> x2 = b/m;
    tmv::SmallVector<T,N> vntemp;
    if (showacc) {
        std::cout<<"x = "<<x<<std::endl;
        std::cout<<"b = m*x "<<b<<std::endl;
        std::cout<<"x2 = b/m "<<x2<<std::endl;
        std::cout<<"Norm(diff) = "<<Norm(vntemp=x2-x)<<std::endl;
        std::cout<<"eps*Norm(x) = "<<eps<<" * "<<Norm(x)<<std::endl;
    }
    Assert(Norm(vntemp=x2-x) < eps*Norm(x),"NonSquare exact b/m");

    tmv::SmallVector<T,6> b2 = x%m;
    x2 = b2*m;
    Assert(Norm(vntemp=x2-x) < eps*Norm(x),"NonSquare x%m");

    tmv::SmallMatrix<T,N,6> minv = m.inverse();
    tmv::SmallMatrix<T,N,N> id = minv*m;
    tmv::SmallMatrix<T,6,6> nonid = m*minv;
    tmv::SmallMatrix<T,N,N> mnntemp;
    if (showacc) {
        std::cout<<"minv = "<<minv<<std::endl;
        std::cout<<"minv*m = "<<id<<std::endl;
        std::cout<<"m*minv = "<<nonid<<std::endl;
        std::cout<<"Norm(id-I) = "<<Norm(mnntemp=id-T(1))<<std::endl;
        std::cout<<"(m*minv)T = "<<nonid.transpose()<<std::endl;
        std::cout<<"(m*minv) - (m*minv)T = "<<
            (nonid-nonid.transpose())<<std::endl;
    }
    Assert(Norm(mnntemp=id-T(1)) < eps,"NonSquare Inverse");
    Assert(Norm(nonid-nonid.transpose()) < eps,"NonSquare Pseudo-Inverse");

    tmv::SmallMatrix<T,N,N> mata;
    m.makeInverseATA(mata);
    tmv::SmallMatrix<T,N,N> mtm = m.transpose()*m;
    tmv::SmallMatrix<T,N,N> mtminv = mtm.inverse();
    Assert(Norm(mnntemp=mata-mtminv) < eps*Norm(mata),"NonSquare InverseATA");

    tmv::SmallMatrix<std::complex<T>,6,N,stor> c = m * std::complex<T>(1,2);
    c.diag() *= std::complex<T>(6,-9);
    if (N > 3) c(2,3) += std::complex<T>(2,3);
    if (N > 1) c(1,0) *= std::complex<T>(0,2);
    if (N > 1) c.col(1) *= std::complex<T>(-1,3);
    c.row(3) += tmv::SmallVector<std::complex<T>,N>(std::complex<T>(1,9));

    tmv::SmallVector<std::complex<T>,N> y(std::complex<T>(7,3));
    y(0) = std::complex<T>(2,9);
    if (N > 1) y(1) = std::complex<T>(-10,4);
    if (N > 2) y(2) = std::complex<T>(5,-1);
    if (N > 3) y(3) = std::complex<T>(-5,-2);

    tmv::SmallMatrix<std::complex<T>,N,6> cn6temp;
    T ceps = EPS * Norm(c) * Norm(cn6temp=c.inverse());
    tmv::SmallVector<std::complex<T>,6> e = c * y;
    tmv::SmallVector<std::complex<T>,N> y2 = e/c;
    tmv::SmallVector<std::complex<T>,N> cvntemp;

    Assert(Norm(cvntemp=y2-y) < ceps*Norm(y),"NonSquare exact e/c");

    tmv::SmallVector<std::complex<T>,6> e2 = y%c;
    y2 = e2*c;
    Assert(Norm(cvntemp=y2-y) < ceps*Norm(y),"NonSquare e%c");

    tmv::SmallMatrix<std::complex<T>,N,6> cinv = c.inverse();
    tmv::SmallMatrix<std::complex<T>,N,N> cid = cinv*c;
    tmv::SmallMatrix<std::complex<T>,N,N> cnntemp;
    Assert(Norm(cnntemp=cid-T(1)) < ceps,"NonSquare CInverse");
    tmv::SmallMatrix<std::complex<T>,6,6 > cnonid = c*cinv;
    Assert(Norm(cnonid-cnonid.adjoint()) < ceps,"NonSquare CPseudo-Inverse");

    tmv::SmallMatrix<std::complex<T>,N,N> cata;
    c.makeInverseATA(cata);
    tmv::SmallMatrix<std::complex<T>,N,N> ctc = c.adjoint()*c;
    tmv::SmallMatrix<std::complex<T>,N,N> ctcinv = ctc.inverse();
    Assert(Norm(cnntemp=cata-ctcinv) < ceps*Norm(cata),
           "NonSquare CInverseATA");

    // Test short matrix (M < N)
    tmv::SmallMatrix<T,N,6,stor> ms = m.transpose();

    b = x * ms;
    x2 = b%ms;
    Assert(Norm(vntemp=x2-x) < eps*Norm(x),"NonSquare exact b%ms");

    b2 = x/ms;
    x2 = ms*b2;
    Assert(Norm(vntemp=x2-x) < eps*Norm(x),"NonSquare x/ms");
}

template <class T> 
void TestAllSmallMatrixDiv()
{
    TestSmallSquareDiv<T,tmv::ColMajor,2>();
    TestSmallSquareDiv<T,tmv::ColMajor,5>();
    TestSmallNonSquareDiv<T,tmv::ColMajor,2>();
    TestSmallNonSquareDiv<T,tmv::ColMajor,5>();

#if (XTEST & 2)
    TestSmallSquareDiv<T,tmv::ColMajor,1>();
    TestSmallSquareDiv<T,tmv::ColMajor,3>();
    TestSmallSquareDiv<T,tmv::ColMajor,4>();
    TestSmallSquareDiv<T,tmv::RowMajor,1>();
    TestSmallSquareDiv<T,tmv::RowMajor,2>();
    TestSmallSquareDiv<T,tmv::RowMajor,3>();
    TestSmallSquareDiv<T,tmv::RowMajor,4>();
    TestSmallSquareDiv<T,tmv::RowMajor,5>();

    TestSmallNonSquareDiv<T,tmv::ColMajor,1>();
    TestSmallNonSquareDiv<T,tmv::ColMajor,3>();
    TestSmallNonSquareDiv<T,tmv::ColMajor,4>();
    TestSmallNonSquareDiv<T,tmv::RowMajor,1>();
    TestSmallNonSquareDiv<T,tmv::RowMajor,2>();
    TestSmallNonSquareDiv<T,tmv::RowMajor,3>();
    TestSmallNonSquareDiv<T,tmv::RowMajor,4>();
    TestSmallNonSquareDiv<T,tmv::RowMajor,5>();
#endif

    std::cout<<"SmallMatrix<"<<tmv::TMV_Text(T())<<"> Division ";
    std::cout<<" passed all basic tests\n";
}

#ifdef TEST_DOUBLE
template void TestAllSmallMatrixDiv<double>();
#endif
#ifdef TEST_FLOAT
template void TestAllSmallMatrixDiv<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestAllSmallMatrixDiv<long double>();
#endif
