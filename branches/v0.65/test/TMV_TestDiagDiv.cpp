
#include "TMV_Test.h"
#include "TMV_Test_1.h"
#include "TMV.h"
#include <fstream>
#include <cstdio>

template <class T> 
void TestDiagDiv()
{
    tmv::DiagMatrix<T> m(10);

    for(int i=0;i<10;++i) m(i) = T(2+4*i);

    tmv::DiagMatrix<T> minv = m;
    minv.invertSelf();

    for(int i=0;i<10;++i) {
        Assert(std::abs(m(i)*minv(i) - T(1)) < EPS,"Diag invertSelf");
    }

    tmv::Vector<T> b(10,4.);
    b(0) = 2;
    b(3) = -10;
    b(5) = 5;
    b(9) = -5;

    T eps = EPS * Norm(m) * Norm(minv);

    tmv::Vector<T> x = b/m;
    tmv::Vector<T> b2 = m*x;
    if (showacc) {
        std::cout<<"b = "<<b<<std::endl;
        std::cout<<"x = b/m = "<<x<<std::endl;
        std::cout<<"b2 = "<<b2<<std::endl;
        std::cout<<"Norm(b-b2) = "<<Norm(b-b2)<<std::endl;
    }
    Assert(Norm(b2-b) < eps*Norm(b),"Diag b/m");

    x = b%m;
    b2 = x*m;
    if (showacc) {
        std::cout<<"b = "<<b<<std::endl;
        std::cout<<"x = b%m = "<<x<<std::endl;
        std::cout<<"b2 = "<<b2<<std::endl;
        std::cout<<"Norm(b-b2) = "<<Norm(b-b2)<<std::endl;
    }
    Assert(Norm(b2-b) < eps*Norm(b),"Diag b%m");

    tmv::DiagMatrix<T> minv2 = m.inverse();
    tmv::DiagMatrix<T> id = m*minv2;
    if (showacc) {
        std::cout<<"minv = "<<minv<<std::endl;
        std::cout<<"minv2 = "<<minv2<<std::endl;
        std::cout<<"m*minv2 = "<<id<<std::endl;
        std::cout<<"minv2*m = "<<minv2*m<<std::endl;
        std::cout<<"Norm(id-I) = "<<Norm(id-T(1))<<std::endl;
    }
    Assert(Norm(id-T(1)) < eps,"Diag inverse (id)");
    Assert(Norm(minv-minv2) < eps*Norm(minv),"Diag inverse");

    tmv::Matrix<T> minv3 = m.inverse();
    Assert(Norm(minv-minv3) < eps*Norm(minv),"Matrix = Diag inverse");

    tmv::DiagMatrix<T> mtm = m.adjoint() * m;
    tmv::DiagMatrix<T> mata(10);
    m.makeInverseATA(mata);
    if (showacc) {
        std::cout<<"mtm = "<<mtm<<std::endl;
        std::cout<<"mtm.inv = "<<mtm.inverse()<<std::endl;
        std::cout<<"m.invata = "<<mata<<std::endl;
        std::cout<<"minv*minvt = "<<minv*minv.adjoint()<<std::endl;
        std::cout<<"Norm(diff) = "<<Norm(mata-mtm.inverse())<<std::endl;
    }
    Assert(Norm(mata-mtm.inverse()) < eps*Norm(mata),"Diag inverseATA");

    tmv::Matrix<T> mata2(10,10);
    m.makeInverseATA(mata2);
    Assert(Norm(mata2-mata) < eps*Norm(mata),"Matrix = Diag inverseATA");

    T det1(1);
    for(int i=0;i<10;++i) det1 *= m(i);
    T det = m.det();
    if (showacc) {
        std::cout<<"det1 = "<<det1<<std::endl;
        std::cout<<"m.det() = "<<det<<std::endl;
        std::cout<<"abs(det-det1) = "<<std::abs(det-det1);
        std::cout<<"  EPS*abs(det1) = "<<eps*std::abs(det1)<<std::endl;
    }
    Assert(std::abs(det-det1) < eps*std::abs(det1),"Diag Det");
    T signdet;
    Assert(std::abs(m.logDet(&signdet)-std::log(det1)) < eps,"Diag logDet");
    Assert(std::abs(signdet-1.) < eps,"Diag logDet - sign");

    tmv::DiagMatrix<std::complex<T> > c = m * std::complex<T>(1,2);;
    c(3) += std::complex<T>(2,3);
    c.subDiagMatrix(5,9) *= std::complex<T>(0,2);

    T ceps = EPS * Norm(c) * Norm(c.inverse());

    std::complex<T> cdet1(1);
    for(int i=0;i<10;++i) cdet1 *= c(i);
    std::complex<T> cdet = c.det();
    if (showacc) {
        std::cout<<"cdet1 = "<<cdet1<<std::endl;
        std::cout<<"c.det = "<<cdet<<std::endl;
        std::cout<<"abs(cdet-cdet1) = "<<std::abs(cdet-cdet1);
        std::cout<<"  EPS*abs(cdet) = "<<ceps*std::abs(cdet1)<<std::endl;
    }
    Assert(std::abs(cdet-cdet1) < ceps*std::abs(cdet1),"Diag Cdet");
    std::complex<T> csigndet;
    Assert(std::abs(c.logDet(&csigndet)-std::log(std::abs(cdet1))) < eps,
           "Diag ClogDet");
    Assert(std::abs(csigndet-cdet1/std::abs(cdet1)) < eps,"Diag ClogDet - sign");

    tmv::DiagMatrix<std::complex<T> > cinv = c.inverse();
    tmv::DiagMatrix<std::complex<T> > cid = c*cinv;
    Assert(Norm(cid-T(1)) < ceps,"Diag Cinverse");

    tmv::DiagMatrix<std::complex<T> > cinv2 = c;
    cinv2.invertSelf();
    Assert(Norm(cinv-cinv2) < eps*Norm(cinv),"Matrix = Diag CinverseSelf");

    tmv::Matrix<std::complex<T> > cinv3 = c.inverse();
    Assert(Norm(cinv-cinv3) < eps*Norm(cinv),"Matrix = Diag Cinverse");

    tmv::DiagMatrix<std::complex<T> > ctc = c.adjoint() * c;
    tmv::DiagMatrix<std::complex<T> > cata(10);
    c.makeInverseATA(cata);
    Assert(Norm(cata-ctc.inverse()) < ceps*Norm(cata),"Diag CinverseATA");

    tmv::Matrix<std::complex<T> > cata2(10,10);
    c.makeInverseATA(cata2);
    Assert(Norm(cata2-cata) < eps*Norm(cata),"Matrix = Diag CinverseATA");

    tmv::Vector<std::complex<T> > e(10);
    e = b*std::complex<T>(1,2);
    e(1) += std::complex<T>(-1,5);
    e(2) -= std::complex<T>(-1,5);

    // test real / complex
    tmv::Vector<std::complex<T> > y = b/c;
    tmv::Vector<std::complex<T> > b3 = c*y;
    Assert(Norm(b3-b) < ceps*Norm(b),"Diag b/c");
    y = b%c;
    b3 = y*c;
    Assert(Norm(b3-b) < ceps*Norm(b),"Diag b%c");

    // test complex / real
    y = e/m;
    b3 = m*y;
    Assert(Norm(b3-e) < eps*Norm(e),"Diag e/m");
    y = e%m;
    b3 = y*m;
    Assert(Norm(b3-e) < eps*Norm(e),"Diag e%m");

    // test complex / complex
    y = e/c;
    b3 = c*y;
    Assert(Norm(b3-e) < ceps*Norm(e),"Diag e/c");
    y = e%c;
    b3 = y*c;
    Assert(Norm(b3-e) < ceps*Norm(e),"Diag e%c");

    TestDiagDiv_A<T>();
    TestDiagDiv_B1<T>();
    TestDiagDiv_B2<T>();
    std::cout<<"DiagMatrix<"<<tmv::TMV_Text(T())<<
        "> Division passed all tests\n";
}


#ifdef TEST_DOUBLE
template void TestDiagDiv<double>();
#endif
#ifdef TEST_FLOAT
template void TestDiagDiv<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestDiagDiv<long double>();
#endif
