
#include "TMV_Test.h"
#include "TMV_Test_1.h"
#include "TMV.h"
#include <fstream>
#include <cstdio>

template <class T, class M, class CM> 
static void TestBasicTriDiv()
{
    tmv::Matrix<T,tmv::ColMajor> a(10,10);

    for(int i=0;i<10;++i) for(int j=0;j<10;++j) 
        a(i,j) = T(2+4*i-3*j);
    a.subMatrix(3,4,8,9) *= T(3);
    a /= T(10);
    a.diag().addToAll(10);

    M m(a);

    M minv = m;
    minv.invertSelf();

    T eps = EPS * Norm(m) * Norm(minv);

    M id1 = m*minv;
    M id2 = minv*m;
    if (showacc) {
        std::cout<<"m = "<<m<<std::endl;
        std::cout<<"minv = "<<minv<<std::endl;
        std::cout<<"m*minv = "<<id1<<std::endl;
        std::cout<<"minv*m = "<<id2<<std::endl;
        std::cout<<"Norm(id1-1) = "<<Norm(id1-T(1))<<std::endl;
        std::cout<<"Norm(id2-1) = "<<Norm(id2-T(1))<<std::endl;
    }
    Assert(Norm(id1-T(1)) < eps,"Tri inverse (id1)");
    Assert(Norm(id2-T(1)) < eps,"Tri inverse (id2)");

    tmv::Vector<T> b(10,4.);
    b(0) = 2;
    b(3) = -10;
    b(5) = 5;
    b(9) = -5;

    tmv::Vector<T> x = b/m;
    tmv::Vector<T> b2 = m*x;
    if (showacc) {
        std::cout<<"b = "<<b<<std::endl;
        std::cout<<"x = b/m = "<<x<<std::endl;
        std::cout<<"minv * b = "<<minv * b<<std::endl;
        std::cout<<"b2 = "<<b2<<std::endl;
        std::cout<<"Norm(b-b2) = "<<Norm(b-b2)<<std::endl;
    }
    Assert(Norm(b2-b) < eps*Norm(b),"Tri b/m");

    x = b%m;
    b2 = x*m;
    if (showacc) {
        std::cout<<"b = "<<b<<std::endl;
        std::cout<<"x = b%m = "<<x<<std::endl;
        std::cout<<"b * minv = "<<b * minv<<std::endl;
        std::cout<<"b2 = "<<b2<<std::endl;
        std::cout<<"Norm(b-b2) = "<<Norm(b-b2)<<std::endl;
    }
    Assert(Norm(b2-b) < eps*Norm(b),"Tri b%m");

    tmv::Matrix<T,tmv::ColMajor> P(10,10,4);
    P(0,3) = 2;
    P(3,1) = -10;
    P(5,9) = 5;
    P(9,4) = -5;
    P.row(4).addToAll(2);
    P.col(8).addToAll(-3);
    P.diag(4).setAllTo(9);
    P.subMatrix(2,5,3,8) *= T(3);

    tmv::Matrix<T,tmv::ColMajor> Q = P/m;
    tmv::Matrix<T,tmv::ColMajor> P2 = m*Q;
    if (showacc) {
        std::cout<<"P = "<<P<<std::endl;
        std::cout<<"Q = P/m = "<<Q<<std::endl;
        std::cout<<"minv * P = "<<minv * P<<std::endl;
        std::cout<<"P2 = "<<P2<<std::endl;
        std::cout<<"Norm(P-P2) = "<<Norm(P-P2)<<std::endl;
    }
    Assert(Norm(P2-P) < eps*Norm(P),"Tri P/m");

    Q = P%m;
    P2 = Q*m;
    if (showacc) {
        std::cout<<"P = "<<P<<std::endl;
        std::cout<<"Q = P%m = "<<Q<<std::endl;
        std::cout<<"P * minv = "<<P * minv<<std::endl;
        std::cout<<"P2 = "<<P2<<std::endl;
        std::cout<<"Norm(P-P2) = "<<Norm(P-P2)<<std::endl;
    }
    Assert(Norm(P2-P) < eps*Norm(P),"Tri P%m");

    tmv::Matrix<T,tmv::RowMajor> Pr = P;
    tmv::Matrix<T,tmv::RowMajor> Qr = Pr/m;
    tmv::Matrix<T,tmv::RowMajor> P2r = m*Qr;
    if (showacc) {
        std::cout<<"Pr = "<<Pr<<std::endl;
        std::cout<<"Qr = Pr/m = "<<Qr<<std::endl;
        std::cout<<"minv * Pr = "<<minv * Pr<<std::endl;
        std::cout<<"P2r = "<<P2r<<std::endl;
        std::cout<<"Norm(Pr-P2r) = "<<Norm(Pr-P2r)<<std::endl;
    }
    Assert(Norm(P2r-Pr) < eps*Norm(Pr),"Tri P/m rowmajor");

    Qr = Pr%m;
    P2r = Qr*m;
    if (showacc) {
        std::cout<<"Pr = "<<Pr<<std::endl;
        std::cout<<"Qr = Pr%m = "<<Qr<<std::endl;
        std::cout<<"Pr * minv = "<<Pr * minv<<std::endl;
        std::cout<<"P2r = "<<P2r<<std::endl;
        std::cout<<"Norm(Pr-P2r) = "<<Norm(Pr-P2r)<<std::endl;
    }
    Assert(Norm(P2r-Pr) < eps*Norm(Pr),"Tri P%m rowmajor");

    M U(P);
    M V = U/m;
    M U2 = m*V;
    if (showacc) {
        std::cout<<"U = "<<U<<std::endl;
        std::cout<<"V = U/m = "<<V<<std::endl;
        std::cout<<"minv * U = "<<minv * U<<std::endl;
        std::cout<<"U2 = "<<U2<<std::endl;
        std::cout<<"Norm(U-U2) = "<<Norm(U-U2)<<std::endl;
    }
    Assert(Norm(U2-U) < eps*Norm(U),"Tri U/m");

    V = U%m;
    U2 = V*m;
    if (showacc) {
        std::cout<<"U = "<<U<<std::endl;
        std::cout<<"V = U%m = "<<V<<std::endl;
        std::cout<<"U * minv = "<<U * minv<<std::endl;
        std::cout<<"U2 = "<<U2<<std::endl;
        std::cout<<"Norm(U-U2) = "<<Norm(U-U2)<<std::endl;
    }
    Assert(Norm(U2-U) < eps*Norm(U),"Tri U%m");

    M minv2 = m.inverse();
    Assert(Norm(minv-minv2) < eps*Norm(minv),"Tri inverse");

    tmv::Matrix<T> minv3 = m.inverse();
    Assert(Norm(minv-minv3) < eps*Norm(minv),"Matrix = Tri inverse");

    tmv::Matrix<T> mata1 = minv * minv.adjoint();
    tmv::Matrix<T> mata(10,10);
    m.makeInverseATA(mata);
    if (showacc) {
        std::cout<<"ata = "<<m.adjoint()*m<<std::endl;
        std::cout<<"m.invata = "<<mata<<std::endl;
        std::cout<<"minv*minvt = "<<mata1<<std::endl;
        std::cout<<"Norm(diff) = "<<Norm(mata-mata1)<<std::endl;
        std::cout<<"ata*mata = "<<(a.adjoint()*a) * mata<<std::endl;
    }
    Assert(Norm(mata-mata1) < eps*Norm(mata1),"Tri inverseATA");

    T det1(1);
    for(int i=0;i<10;++i) det1 *= m(i,i);
    T det = m.det();
    if (showacc) {
        std::cout<<"det1 = "<<det1<<std::endl;
        std::cout<<"m.det() = "<<det<<std::endl;
        std::cout<<"abs(det-det1) = "<<std::abs(det-det1);
        std::cout<<"  EPS*abs(det1) = "<<eps*std::abs(det1)<<std::endl;
    }
    Assert(std::abs(det-det1) < eps*std::abs(det1),"Tri Det");
    T signdet;
    Assert(std::abs(m.logDet(&signdet)-std::log(det1)) < eps,"Tri logDet");
    Assert(std::abs(signdet-1.) < eps,"Tri logDet - sign");

    tmv::Matrix<std::complex<T> > ca = a * std::complex<T>(1,2);
    CM c(ca);

    T ceps = EPS * Norm(c) * Norm(c.inverse());

    std::complex<T> cdet1(1);
    for(int i=0;i<10;++i) cdet1 *= c(i,i);
    std::complex<T> cdet = c.det();
    if (showacc) {
        std::cout<<"cdet1 = "<<cdet1<<std::endl;
        std::cout<<"c.det = "<<cdet<<std::endl;
        std::cout<<"abs(cdet-cdet1) = "<<std::abs(cdet-cdet1);
        std::cout<<"  EPS*abs(cdet) = "<<ceps*std::abs(cdet1)<<std::endl;
    }
    Assert(std::abs(cdet-cdet1) < ceps*std::abs(cdet1),"Tri Cdet");
    std::complex<T> csigndet;
    Assert(std::abs(c.logDet(&csigndet)-std::log(std::abs(cdet1))) < eps,
           "Tri ClogDet");
    Assert(std::abs(csigndet-cdet1/std::abs(cdet1)) < eps,"Tri ClogDet - sign");

    CM cinv = c.inverse();
    CM cid = c*cinv;
    Assert(Norm(cid-T(1)) < ceps,"Tri Cinverse");

    CM cinv2 = c;
    cinv2.invertSelf();
    Assert(Norm(cinv-cinv2) < eps*Norm(cinv),"Matrix = Tri CinverseSelf");

    tmv::Matrix<std::complex<T> > cinv3 = c.inverse();
    Assert(Norm(cinv-cinv3) < eps*Norm(cinv),"Matrix = Tri Cinverse");

    tmv::Matrix<std::complex<T> > cata1 = cinv * cinv.adjoint();
    tmv::Matrix<std::complex<T> > cata(10,10);
    c.makeInverseATA(cata);
    if (showacc) {
        std::cout<<"ctc = "<<c.adjoint()*c<<std::endl;
        std::cout<<"c.invata = "<<cata<<std::endl;
        std::cout<<"cinv*cinvt = "<<cata1<<std::endl;
        std::cout<<"Norm(diff) = "<<Norm(cata-cata1)<<std::endl;
        std::cout<<"ctc*cata = "<<(c.adjoint()*c) * cata<<std::endl;
    }
    Assert(Norm(cata-cata1) < ceps*Norm(cata1),"Tri CinverseATA");

    tmv::Vector<std::complex<T> > e(10);
    e = b*std::complex<T>(1,2);
    e(1) += std::complex<T>(-1,5);
    e(2) -= std::complex<T>(-1,5);

    // test complex / real
    tmv::Vector<std::complex<T> > y = e/m;
    tmv::Vector<std::complex<T> > e2 = m*y;
    Assert(Norm(e2-e) < eps*Norm(e),"Tri e/m");
    y = e%m;
    e2 = y*m;
    Assert(Norm(e2-e) < eps*Norm(e),"Tri e%m");

    // test real / complex
    y = b/c;
    tmv::Vector<std::complex<T> > b3 = c*y;
    Assert(Norm(b3-b) < ceps*Norm(b),"Tri b/c");
    y = b%c;
    b3 = y*c;
    Assert(Norm(b3-b) < ceps*Norm(b),"Tri b%c");

    // test complex / complex
    y = e/c;
    e2 = c*y;
    Assert(Norm(e2-e) < ceps*Norm(e),"Tri e/c");
    y = e%c;
    e2 = y*c;
    Assert(Norm(e2-e) < ceps*Norm(e),"Tri e%c");

    tmv::Matrix<std::complex<T>,tmv::ColMajor> R = std::complex<T>(1,2) * P;

    // matrix complex / real
    tmv::Matrix<std::complex<T>,tmv::ColMajor> S = R/m;
    tmv::Matrix<std::complex<T>,tmv::ColMajor> R2 = m*S;
    Assert(Norm(R2-R) < eps*Norm(R),"Tri R/m");
    S = R%m;
    R2 = S*m;
    Assert(Norm(R2-R) < eps*Norm(R),"Tri R%m");

    tmv::Matrix<std::complex<T>,tmv::RowMajor> Rr = R;
    tmv::Matrix<std::complex<T>,tmv::RowMajor> Sr = Rr/m;
    tmv::Matrix<std::complex<T>,tmv::RowMajor> R2r = m*Sr;
    Assert(Norm(R2r-Rr) < eps*Norm(Rr),"Tri R/m rowmajor");
    Sr = Rr%m;
    R2r = Sr*m;
    Assert(Norm(R2r-Rr) < eps*Norm(Rr),"Tri R%m rowmajor");

    // matrix real/complex
    S = P/c;
    tmv::Matrix<std::complex<T>,tmv::ColMajor> P3 = c*S;
    Assert(Norm(P3-P) < eps*Norm(P),"Tri P/c");
    S = P%c;
    P3 = S*c;
    Assert(Norm(P3-P) < eps*Norm(P),"Tri P%c");

    Sr = Pr/c;
    tmv::Matrix<std::complex<T>,tmv::RowMajor> P3r = c*Sr;
    Assert(Norm(P3r-Pr) < eps*Norm(Pr),"Tri P/c rowmajor");
    Sr = Pr%c;
    P3r = Sr*c;
    Assert(Norm(P3r-Pr) < eps*Norm(Pr),"Tri P%c rowmajor");
 
    // matrix complex / complex
    S = R/c;
    R2 = c*S;
    Assert(Norm(R2-R) < eps*Norm(R),"Tri R/c");
    S = R%c;
    R2 = S*c;
    Assert(Norm(R2-R) < eps*Norm(R),"Tri R%c");

    Rr = R;
    Sr = Rr/c;
    R2r = c*Sr;
    Assert(Norm(R2r-Rr) < eps*Norm(Rr),"Tri R/c rowmajor");
    Sr = Rr%c;
    R2r = Sr*c;
    Assert(Norm(R2r-Rr) < eps*Norm(Rr),"Tri R%c rowmajor");

    // trimatrix complex / real
    CM W(R);
    CM X = W/m;
    CM W2 = m*X;
    Assert(Norm(W2-W) < eps*Norm(W),"Tri W/m");
    X = W%m;
    W2 = X*m;
    Assert(Norm(W2-W) < eps*Norm(W),"Tri W%m");

    // trimatrix real/complex
    X = U/c;
    CM U3 = c*X;
    Assert(Norm(U3-U) < eps*Norm(U),"Tri U/c");
    X = U%c;
    U3 = X*c;
    Assert(Norm(U3-U) < eps*Norm(P),"Tri U%c");

    // trimatrix complex / complex
    X = W/c;
    W2 = c*X;
    Assert(Norm(W2-W) < eps*Norm(W),"Tri W/c");
    X = W%c;
    W2 = X*c;
    Assert(Norm(W2-W) < eps*Norm(W),"Tri W%c");
}

template <class T>
void TestTriDiv()
{
    typedef std::complex<T> CT;
    TestBasicTriDiv<T,
        tmv::UpperTriMatrix<T>,
        tmv::UpperTriMatrix<CT> >();
    TestBasicTriDiv<T,
        tmv::UpperTriMatrix<T,tmv::UnitDiag>,
        tmv::UpperTriMatrix<CT,tmv::UnitDiag> >();
    TestBasicTriDiv<T,
        tmv::LowerTriMatrix<T>,
        tmv::LowerTriMatrix<CT> >();
    TestBasicTriDiv<T,
        tmv::LowerTriMatrix<T,tmv::UnitDiag>,
        tmv::LowerTriMatrix<CT,tmv::UnitDiag> >();
#if (XTEST & 2)
    TestBasicTriDiv<T,
        tmv::UpperTriMatrix<T,tmv::NonUnitDiag|tmv::RowMajor>,
        tmv::UpperTriMatrix<CT,tmv::NonUnitDiag|tmv::RowMajor> >();
    TestBasicTriDiv<T,
        tmv::UpperTriMatrix<T,tmv::UnitDiag|tmv::RowMajor>,
        tmv::UpperTriMatrix<CT,tmv::UnitDiag|tmv::RowMajor> >();
    TestBasicTriDiv<T,
        tmv::LowerTriMatrix<T,tmv::NonUnitDiag|tmv::RowMajor>,
        tmv::LowerTriMatrix<CT,tmv::NonUnitDiag|tmv::RowMajor> >();
    TestBasicTriDiv<T,
        tmv::LowerTriMatrix<T,tmv::UnitDiag|tmv::RowMajor>,
        tmv::LowerTriMatrix<CT,tmv::UnitDiag|tmv::RowMajor> >();
#endif

    TestTriDiv_A1<T>();
    TestTriDiv_A2<T>();
    TestTriDiv_B1<T>();
    TestTriDiv_B2<T>();
    TestTriDiv_C1<T>();
    TestTriDiv_C2<T>();

    std::cout<<"TriMatrix<"<<tmv::TMV_Text(T())<<
        "> Division passed all tests\n";
}

#ifdef TEST_DOUBLE
template void TestTriDiv<double>();
#endif
#ifdef TEST_FLOAT
template void TestTriDiv<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestTriDiv<long double>();
#endif
