
#include "TMV_Test.h"
#include "TMV_Test_3.h"

#if 0
template <class T, tmv::StorageType stor, int N> 
static void TestSmallNonSquareDiv()
{
    typedef typename tmv::Traits<T>::float_type FT;
    tmv::SmallMatrix<T,6,N,stor> m;
    for(int i=0;i<6;++i) for(int j=0;j<N;++j) m(i,j) = T(2+4*i-5*j);
    m.diag() *= T(11);
    m /= T(7);
    if (N > 1) m(1,0) = T(-2);
    if (N > 2) m(2,0) = T(7);
    if (N > 3) m(3,0) = T(-10);
    if (N > 2) m(2,2) = T(30);

    tmv::SmallVector<T,N> x;
    x(0) = T(2);
    if (N  > 1) x(1) = T(-10);
    if (N  > 2) x(2) = T(5);
    if (N  > 3) x(3) = T(-5);

    if (showstartdone) {
        std::cout<<"Start TestSmallSquareDiv\n";
        std::cout<<"m = "<<TMV_Text(m)<<" "<<m<<std::endl;
    }

    tmv::SmallMatrix<T,N,6> mn6temp;

    T eps = EPS * Norm(m) * Norm(mn6temp=m.inverse());
    tmv::SmallVector<T,6> b = m * x;
    tmv::SmallVector<T,N> x2 = b/m;
    tmv::SmallVector<T,N> vntemp;
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

    tmv::SmallVector<std::complex<T>,N> y;
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
#endif

template <class T> 
void TestAllSmallMatrixDiv()
{
    TestSmallSquareDiv_a<T>();
    TestSmallSquareDiv_b<T>();
    TestSmallSquareDiv_c<T>();
    TestSmallSquareDiv_d<T>();
    TestSmallSquareDiv_e<T>();
    TestSmallSquareDiv_f<T>();
    TestSmallSquareDiv_g<T>();
    TestSmallSquareDiv_h<T>();
    TestSmallSquareDiv_i<T>();
    TestSmallSquareDiv_j<T>();

    std::cout<<"Square SmallMatrix<"<<tmv::TMV_Text(T())<<"> Division ";
    std::cout<<"passed all tests\n";

#if 0
    TestSmallNonSquareDiv<T,tmv::ColMajor,1>();
    TestSmallNonSquareDiv<T,tmv::ColMajor,3>();
    TestSmallNonSquareDiv<T,tmv::ColMajor,4>();
    TestSmallNonSquareDiv<T,tmv::RowMajor,1>();
    TestSmallNonSquareDiv<T,tmv::RowMajor,2>();
    TestSmallNonSquareDiv<T,tmv::RowMajor,3>();
    TestSmallNonSquareDiv<T,tmv::RowMajor,4>();
    TestSmallNonSquareDiv<T,tmv::RowMajor,5>();
    std::cout<<"NonSquare SmallMatrix<"<<tmv::TMV_Text(T())<<"> Division ";
    std::cout<<"passed all tests\n";
#endif

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
