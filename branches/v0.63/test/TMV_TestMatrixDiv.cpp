
#include "TMV_Test.h"
#include "TMV_Test1.h"
#include "TMV.h"

#include "TMV_TestMatrixDivArith.h"

template <class T, tmv::StorageType stor> 
static void TestSquareDiv(tmv::DivType dt)
{
    tmv::Matrix<T,stor> m(4,4);

    for(int i=0;i<4;++i) for(int j=0;j<4;++j) m(i,j) = T(2+4*i-5*j);
    m(0,0) = 14;
    m(1,0) = -2;
    m(2,0) = 7;
    m(3,0) = -10;
    m(2,2) = 30;

    tmv::Vector<T> b(4);
    b(0) = 2;
    b(1) = -10;
    b(2) = 5;
    b(3) = -5;

    m.DivideUsing(dt);
    m.SaveDiv();
    m.SetDiv();
    std::ostream* dbgout = showdiv ? &std::cout : 0;
    Assert(m.CheckDecomp(dbgout),"CheckDecomp");

    T eps = EPS * Norm(m) * Norm(m.Inverse());

    tmv::Vector<T> x = b/m;
    tmv::Vector<T> b2 = m*x;
    if (showacc) {
        std::cout<<"b = "<<b<<std::endl;
        std::cout<<"x = b/m = "<<x<<std::endl;
        std::cout<<"b2 = "<<b2<<std::endl;
        std::cout<<"Norm(b-b2) = "<<Norm(b-b2)<<std::endl;
    }
    Assert(Norm(b2-b) < eps*Norm(b),"Square b/m");

    x = b%m;
    b2 = x*m;
    if (showacc) {
        std::cout<<"b = "<<b<<std::endl;
        std::cout<<"x = b%m = "<<x<<std::endl;
        std::cout<<"b2 = "<<b2<<std::endl;
        std::cout<<"Norm(b-b2) = "<<Norm(b-b2)<<std::endl;
    }
    Assert(Norm(b2-b) < eps*Norm(b),"Square b%m");

    tmv::Matrix<T> minv = m.Inverse();
    tmv::Matrix<T> id = m*minv;
    if (showacc) {
        std::cout<<"minv = "<<minv<<std::endl;
        std::cout<<"m*minv = "<<id<<std::endl;
        std::cout<<"minv*m = "<<minv*m<<std::endl;
        std::cout<<"Norm(id-I) = "<<Norm(id-T(1))<<std::endl;
    }
    Assert(Norm(id-T(1)) < eps,"Square Inverse");

    tmv::Matrix<T> mtm = m.Adjoint() * m;
    tmv::Matrix<T> mata(4,4);
    m.InverseATA(mata);
    if (showacc) {
        std::cout<<"mtm = "<<mtm<<std::endl;
        std::cout<<"mtm.inv = "<<mtm.Inverse()<<std::endl;
        std::cout<<"m.invata = "<<mata<<std::endl;
        std::cout<<"minv*minvt = "<<minv*minv.Adjoint()<<std::endl;
        std::cout<<"Norm(diff) = "<<Norm(mata-mtm.Inverse())<<std::endl;
    }
    Assert(Norm(mata-mtm.Inverse()) < eps*Norm(mata),"Square InverseATA");

    T mdet = 28800;
    if (showacc) {
        std::cout<<"abs(det-mdet) = "<<std::abs(m.Det()-mdet);
        std::cout<<"  EPS*abs(mdet) = "<<eps*std::abs(mdet)<<std::endl;
    }
    Assert(std::abs(m.Det()-mdet) < eps*std::abs(mdet),"Square Det");
    T sdet;
    Assert(std::abs(m.LogDet(&sdet)-std::log(mdet)) < eps,"Square LogDet");
    Assert(std::abs(sdet-1.) < eps,"Square LogDet - sign");

    tmv::Matrix<std::complex<T>,stor> c(4,4);
    c = m;
    c(2,3) += std::complex<T>(2,3);
    c(1,0) *= std::complex<T>(0,2);
    c.col(1) *= std::complex<T>(-1,3);
    c.row(3) += tmv::Vector<std::complex<T> >(4,std::complex<T>(1,9));

    c.DivideUsing(dt);
    c.SaveDiv();
    c.SetDiv();
    Assert(c.CheckDecomp(dbgout),"CheckDecomp");

    T ceps = EPS * Norm(m) * Norm(m.Inverse());

    std::complex<T> cdet(-103604,101272);
    if (showacc) {
        std::cout<<"cdet = "<<cdet<<std::endl;
        std::cout<<"C.Det = "<<c.Det()<<std::endl;
        std::cout<<"abs(det-cdet) = "<<std::abs(c.Det()-cdet);
        std::cout<<"  EPS*abs(cdet) = "<<ceps*std::abs(cdet)<<std::endl;
    }
    Assert(std::abs(c.Det()-cdet) < ceps*std::abs(cdet),"Square CDet");
    std::complex<T> csdet;
    Assert(std::abs(c.LogDet(&csdet)-std::log(std::abs(cdet))) < eps,
           "Square CLogDet");
    Assert(std::abs(csdet-cdet/std::abs(cdet)) < eps,"Square CLogDet - sign");

    tmv::Matrix<std::complex<T> > cinv = c.Inverse();
    tmv::Matrix<std::complex<T> > cid = c*cinv;
    Assert(Norm(cid-T(1)) < ceps,"Square CInverse");

    tmv::Matrix<std::complex<T> > ctc = c.Adjoint() * c;
    tmv::Matrix<std::complex<T> > cata(4,4);
    c.InverseATA(cata);
    Assert(Norm(cata-ctc.Inverse()) < ceps*Norm(cata),"Square CInverseATA");

    tmv::Vector<std::complex<T> > e(4);
    e = b*std::complex<T>(1,2);
    e(1) += std::complex<T>(-1,5);
    e(2) -= std::complex<T>(-1,5);

    // test real / complex
    tmv::Vector<std::complex<T> > y = b/c;
    tmv::Vector<std::complex<T> > b3 = c*y;
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

    // test really big one
    const size_t BIGN=100;
    T p[BIGN] = { 
        3,6,1,6,8,3,3,3,34,25, 
        76,4,67,52,3,2,1,2,4,6, 57,4,24,7,2,1,33,4,23,9, 3,6,1,6,8,3,3,3,34,5, 
        76,4,67,52,3,2,1,2,4,6, 57,4,24,7,2,1,33,4,23,9, 3,6,1,6,8,3,3,3,34,5, 
        76,4,67,52,3,2,1,2,4,6, 57,4,24,7,2,1,33,4,23,9, 3,6,1,6,8,3,3,3,34,5 
    };
    T q[BIGN] = { 
        12,3,5,34,52,4,4,23,42,68, 
        71,5,4,5,5,5,45,3,2,5, 6,32,2,53,6,5,2,43,1,1, 12,3,5,34,5,4,4,23,42,8,
        71,5,4,5,5,5,45,3,2,5, 6,32,2,53,6,5,2,43,1,1, 12,3,5,34,5,4,4,23,42,8,
        71,5,4,5,5,5,45,3,2,5, 6,32,2,53,6,5,2,43,1,1, 12,3,5,34,5,4,4,23,42,8
    };
    T r[BIGN] = {
        42,3,51,2,7,42,4,3,42,14, 
        4,24,2,8,4,24,6,1,6,7, 2,42,1,2,35,7,3,5,32,13, 42,3,1,2,7,42,4,3,42,1,
        4,24,2,8,4,24,6,1,6,7, 2,42,1,2,35,7,3,5,32,13, 42,3,1,2,7,42,4,3,42,1,
        4,24,2,8,4,24,6,1,6,7, 2,42,1,2,35,7,3,5,32,13, 42,3,1,2,7,42,4,3,42,1
    };
    tmv::Vector<T> P(BIGN,p);
    tmv::Vector<T> Q(BIGN,q);
    tmv::Vector<T> R(BIGN,r);
    tmv::Matrix<T,stor> M = P ^ Q;
    tmv::Matrix<std::complex<T>,stor> CM = P ^ (std::complex<T>(-4,10)*Q);
    M.diag().AddToAll(T(215));
    CM.diag().AddToAll(std::complex<T>(103,-53));
    M.row(size_t(floor(0.23*BIGN))) *= T(12);
    M(size_t(floor(0.12*BIGN)),1) -= T(142);
    CM(size_t(floor(0.06*BIGN)),2) += std::complex<T>(23,89);
    CM.col(size_t(floor(0.15*BIGN))) *= std::complex<T>(61,12);
    M.SubVector(size_t(floor(0.65*BIGN)),size_t(floor(0.05*BIGN)),1,3,
                size_t(floor(0.29*BIGN))) *= T(2);
    M.SubVector(size_t(floor(0.98*BIGN)),size_t(floor(0.12*BIGN)),-1,2,
                size_t(floor(0.18*BIGN))).AddToAll(T(197));
    CM.SubVector(size_t(floor(0.53*BIGN)),0,1,3,size_t(floor(0.31*BIGN))) *= 
        std::complex<T>(2,-1);
    CM.SubVector(size_t(floor(0.88*BIGN)),size_t(floor(0.18*BIGN)),-1,2,
                 size_t(floor(0.23*BIGN))).AddToAll(std::complex<T>(197,174));
    M.DivideUsing(dt);
    M.SaveDiv();
    CM.DivideUsing(dt);
    CM.SaveDiv();

    eps = EPS * Norm(M) * Norm(M.Inverse());
    ceps = EPS * Norm(CM) * Norm(CM.Inverse());

    tmv::Vector<T> S = R/M;
    tmv::Vector<T> R2 = M*S;
    Assert(M.CheckDecomp(),"CheckDecomp");
    //Assert(M.CheckDecomp(dbgout),"CheckDecomp");

    if (showacc) {
        std::cout<<"R/M Norm(R2-R) = "<<Norm(R2-R)<<std::endl;
        std::cout<<"EPS*Norm(R) = "<<eps*Norm(R)<<std::endl;
    }
    Assert(Norm(R2-R) < eps*Norm(R),"Square R/M");
    S = R%M;
    R2 = S*M;
    if (showacc) {
        std::cout<<"R%M Norm(R2-R) = "<<Norm(R2-R)<<std::endl;
        std::cout<<"EPS*Norm(R) = "<<eps*Norm(R)<<std::endl;
    }
    Assert(Norm(R2-R) < eps*Norm(R),"Square R%M");
    tmv::Vector<std::complex<T> > CS = R/CM;
    tmv::Vector<std::complex<T> > CR2 = CM*CS;
    Assert(CM.CheckDecomp(),"CheckDecomp");
    //Assert(CM.CheckDecomp(dbgout),"CheckDecomp");
    if (showacc) {
        std::cout<<"R/CM Norm(CR2-R) = "<<Norm(CR2-R)<<std::endl;
        std::cout<<"EPS*Norm(R) = "<<ceps*Norm(R)<<std::endl;
    }
    Assert(Norm(CR2-R) < ceps*Norm(R),"Square R/CM");
    CS = R%CM;
    CR2 = CS*CM;
    if (showacc) {
        std::cout<<"R%CM Norm(CR2-R) = "<<Norm(CR2-R)<<std::endl;
        std::cout<<"EPS*Norm(R) = "<<ceps*Norm(R)<<std::endl;
    }
    Assert(Norm(CR2-R) < ceps*Norm(R),"Square R%CM");
    tmv::Vector<std::complex<T> > CR = std::complex<T>(3,-4)*R;
    CS = CR/CM;
    CR2 = CM*CS;
    if (showacc) {
        std::cout<<"CR/CM Norm(CR2-CR) = "<<Norm(CR2-CR)<<std::endl;
        std::cout<<"EPS*Norm(CR) = "<<ceps*Norm(CR)<<std::endl;
    }
    Assert(Norm(CR2-CR) < ceps*Norm(CR),"Square CR/CM");
    CS = CR%CM;
    CR2 = CS*CM;
    if (showacc) {
        std::cout<<"CR%CM Norm(CR2-CR) = "<<Norm(CR2-CR)<<std::endl;
        std::cout<<"EPS*Norm(CR) = "<<ceps*Norm(CR)<<std::endl;
    }
    Assert(Norm(CR2-CR) < ceps*Norm(CR),"Square CR%CM");
    CS = CR/M;
    CR2 = M*CS;
    if (showacc) {
        std::cout<<"CR/M Norm(CR2-CR) = "<<Norm(CR2-CR)<<std::endl;
        std::cout<<"EPS*Norm(CR) = "<<eps*Norm(CR)<<std::endl;
    }
    Assert(Norm(CR2-CR) < eps*Norm(CR),"Square CR/M");
    CS = CR%M;
    CR2 = CS*M;
    if (showacc) {
        std::cout<<"CR%M Norm(CR2-CR) = "<<Norm(CR2-CR)<<std::endl;
        std::cout<<"EPS*Norm(CR) = "<<eps*Norm(CR)<<std::endl;
    }
    Assert(Norm(CR2-CR) < eps*Norm(CR),"Square CR%M");

    tmv::Matrix<T,stor> a1 = m;
    tmv::Matrix<T,stor> a2 = m.Transpose();
    a2.row(1) *= T(3);
    a2.col(2) -= tmv::Vector<T>(4,4.);

    tmv::Matrix<std::complex<T>,stor> c1 = a1 * std::complex<T>(1,2);
    tmv::Matrix<std::complex<T>,stor> c2 = a2 * std::complex<T>(-3,4);
    c1.diag().AddToAll(std::complex<T>(3,1));
    c2.diag().AddToAll(std::complex<T>(-5,8));
    c1.row(3).AddToAll(std::complex<T>(1,-6));
    c2.row(0).AddToAll(std::complex<T>(-2,-11));

    tmv::Matrix<T> a1x = a1;
    tmv::Matrix<std::complex<T> > c1x = c1;
    TestMatrixDivArith2<T>(
        dt,a1x,c1x,a1.View(),a2.View(),c1.View(),c2.View(), "Square"); 
#ifdef XTEST
    tmv::Matrix<T,stor,tmv::FortranStyle> a1f = a1;
    tmv::Matrix<T,stor,tmv::FortranStyle> a2f = a2;
    tmv::Matrix<std::complex<T>,stor,tmv::FortranStyle> c1f = c1;
    tmv::Matrix<std::complex<T>,stor,tmv::FortranStyle> c2f = c2;
    TestMatrixDivArith1<T>(
        dt,a1x,c1x,a1f.View(),a2.View(),c1f.View(),c2.View(), "Square"); 
    TestMatrixDivArith1<T>(
        dt,a1x,c1x,a1.View(),a2f.View(),c1.View(),c2f.View(), "Square"); 
    TestMatrixDivArith1<T>(
        dt,a1x,c1x,a1f.View(),a2f.View(),c1f.View(),c2f.View(), "Square"); 
#endif

    tmv::Matrix<T,stor> a3(7,4);
    for(int i=0;i<7;++i) for(int j=0;j<4;++j) a3(i,j) = T(1-3*i+2*j);
    tmv::Matrix<T,stor> a4 = a3.Transpose();
    a3.SubMatrix(2,6,0,4) += a1;
    a4.SubMatrix(0,4,1,5) -= a2;

    tmv::Matrix<std::complex<T>,stor> c3 = a3*std::complex<T>(1,2);
    tmv::Matrix<std::complex<T>,stor> c4 = c3.Adjoint();
    c3.SubMatrix(2,6,0,4) += c1;
    c4.SubMatrix(0,4,1,5) -= c2;
    c3.col(1) *= std::complex<T>(2,1);
    c3.row(2).AddToAll(std::complex<T>(-7,2));
    c4.col(3) *= std::complex<T>(-1,3);
    c4.row(0).AddToAll(std::complex<T>(1,9));

    tmv::Matrix<T,stor> a3x = a3;
    tmv::Matrix<T,stor> a4x = a4;
    tmv::Matrix<std::complex<T> > c3x = c3;
    tmv::Matrix<std::complex<T> > c4x = c4;
    TestMatrixDivArith1<T>(
        dt,a3x,c3x,a1.View(),a3.View(),c1.View(),c3.View(), "Square/NonSquare");
    TestMatrixDivArith1<T>(
        dt,a4x,c4x,a1.View(),a4.View(),c1.View(),c4.View(), "Square/NonSquare");

#ifdef XTEST
    tmv::Matrix<T,stor> a5(4,0);
    tmv::Matrix<T,stor> a6(0,4);
    tmv::Matrix<std::complex<T>,stor> c5 = a5;
    tmv::Matrix<std::complex<T>,stor> c6 = a6;

    tmv::Matrix<T,stor> a5x = a5;
    tmv::Matrix<T,stor> a6x = a6;
    tmv::Matrix<std::complex<T> > c5x = c5;
    tmv::Matrix<std::complex<T> > c6x = c6;
    TestMatrixDivArith1<T>(
        dt,a5,c5,a1.View(),a5.View(),c1.View(),c5.View(), "Square/Degenerate");
    TestMatrixDivArith1<T>(
        dt,a6,c6,a1.View(),a6.View(),c1.View(),c6.View(), "Square/Degenerate");
#endif

    if (stor == tmv::ColMajor) {
#ifdef XTEST
        TestSquareDiv<T,tmv::RowMajor>(dt);
    } else {
#endif
        std::cout<<"Square Matrix<"<<tmv::TMV_Text(T())<<"> Division using ";
        std::cout<<tmv::TMV_Text(dt)<<" passed all tests\n";
    }
}

template <class T, tmv::StorageType stor> 
static void TestNonSquareDiv(tmv::DivType dt)
{
    tmv::Matrix<T,stor> m(6,4);
    for(int i=0;i<6;++i) for(int j=0;j<4;++j) m(i,j) = T(2+4*i-5*j);
    m(0,0) = 14;
    m(1,0) = -2;
    m(2,0) = 7;
    m(3,0) = -10;
    m(2,2) = 30;

    tmv::Vector<T> x(4);
    x(0) = 2;
    x(1) = -10;
    x(2) = 5;
    x(3) = -5;

    m.DivideUsing(dt);
    m.SaveDiv();
    m.SetDiv();

    std::ostream* dbgout = showdiv ? &std::cout : 0;
    Assert(m.CheckDecomp(dbgout),"CheckDecomp");

    T eps = EPS * Norm(m) * Norm(m.Inverse());
    tmv::Vector<T> b = m * x;
    tmv::Vector<T> x2 = b/m;
    Assert(Norm(x2-x) < eps*Norm(x),"NonSquare exact b/m");

    tmv::Vector<T> b2 = x%m;
    x2 = b2*m;
    Assert(Norm(x2-x) < eps*Norm(x),"NonSquare x%m");

    b(0) += 100;
    x = b/m;
    b2 = m*x;
    T refnorm = Norm(b2-b);
    T dxval = T(sqrt(eps))*Norm(x);
    tmv::Vector<T> dx = dxval*tmv::BasisVector<T>(4,0);
    b2 = m*(x+dx);
    Assert(Norm(b2-b) >= refnorm,"NonSquare Least Squares b/m (1)");
    b2 = m*(x-dx);
    Assert(Norm(b2-b) >= refnorm,"NonSquare Least Squares b/m (2)");
    dx = dxval*tmv::BasisVector<T>(4,1);
    b2 = m*(x+dx);
    Assert(Norm(b2-b) >= refnorm,"NonSquare Least Squares b/m (3)");
    b2 = m*(x-dx);
    Assert(Norm(b2-b) >= refnorm,"NonSquare Least Squares b/m (4)");
    dx = dxval*tmv::BasisVector<T>(4,2);
    b2 = m*(x+dx);
    Assert(Norm(b2-b) >= refnorm,"NonSquare Least Squares b/m (5)");
    b2 = m*(x-dx);
    Assert(Norm(b2-b) >= refnorm,"NonSquare Least Squares b/m (6)");
    dx = dxval*tmv::BasisVector<T>(4,3);
    b2 = m*(x+dx);
    Assert(Norm(b2-b) >= refnorm,"NonSquare Least Squares b/m (7)");
    b2 = m*(x-dx);
    Assert(Norm(b2-b) >= refnorm,"NonSquare Least Squares b/m (8)");

    tmv::Matrix<T> minv = m.Inverse();
    tmv::Matrix<T> id = minv*m;
    tmv::Matrix<T> nonid = m*minv;
    if (showacc) {
        std::cout<<"minv = "<<minv<<std::endl;
        std::cout<<"minv*m = "<<id<<std::endl;
        std::cout<<"m*minv = "<<nonid<<std::endl;
        std::cout<<"Norm(id-I) = "<<Norm(id-T(1))<<std::endl;
    }
    Assert(Norm(id-T(1)) < eps,"NonSquare Inverse");
    Assert(Norm(nonid-nonid.Transpose()) < eps,"NonSquare Pseudo-Inverse");

    tmv::Matrix<T> mata(4,4);
    m.InverseATA(mata);
    tmv::Matrix<T> mtm = m.Transpose()*m;
    Assert(Norm(mata-mtm.Inverse()) < eps*Norm(mata),"NonSquare InverseATA");

    tmv::Matrix<std::complex<T>,stor> c(6,4);
    c = m;
    c(2,3) += std::complex<T>(2,3);
    c(1,0) *= std::complex<T>(0,2);
    c.col(1) *= std::complex<T>(-1,3);
    c.row(3) += tmv::Vector<std::complex<T> >(4,std::complex<T>(1,9));

    tmv::Vector<std::complex<T> > y(4);
    y(0) = std::complex<T>(2,9);
    y(1) = std::complex<T>(-10,4);
    y(2) = std::complex<T>(5,-1);
    y(3) = std::complex<T>(-5,-2);

    c.DivideUsing(dt);
    c.SaveDiv();
    c.SetDiv();

    T ceps = EPS * Norm(c) * Norm(c.Inverse());
    tmv::Vector<std::complex<T> > e = c * y;
    tmv::Vector<std::complex<T> > y2 = e/c;

    Assert(c.CheckDecomp(dbgout),"CheckDecomp");

    Assert(Norm(y2-y) < ceps*Norm(y),"NonSquare exact e/c");

    tmv::Vector<std::complex<T> > e2 = y%c;
    y2 = e2*c;
    Assert(Norm(y2-y) < ceps*Norm(y),"NonSquare e%c");

    tmv::Matrix<std::complex<T> > cinv = c.Inverse();
    tmv::Matrix<std::complex<T> > cid = cinv*c;
    Assert(Norm(cid-T(1)) < ceps,"NonSquare CInverse");
    tmv::Matrix<std::complex<T> > cnonid = c*cinv;
    Assert(Norm(cnonid-cnonid.Adjoint()) < ceps,"NonSquare CPseudo-Inverse");

    tmv::Matrix<std::complex<T> > cata(4,4);
    c.InverseATA(cata);
    tmv::Matrix<std::complex<T> > ctc = c.Adjoint()*c;
    Assert(Norm(cata-ctc.Inverse()) < ceps*Norm(cata),"NonSquare CInverseATA");

    // Test short matrix (M < N)
    tmv::Matrix<T,stor> ms = m.Transpose();
    ms.DivideUsing(dt);
    ms.SaveDiv();

    eps = EPS * Norm(ms) * Norm(ms.Inverse());
    b = x * ms;
    x2 = b%ms;
    Assert(ms.CheckDecomp(dbgout),"CheckDecomp");
    Assert(Norm(x2-x) < eps*Norm(x),"NonSquare exact b%ms");

    b2 = x/ms;
    x2 = ms*b2;
    Assert(Norm(x2-x) < eps*Norm(x),"NonSquare x/ms");

    // Test really long matrix
    tmv::Matrix<std::complex<T>,stor> a(30,10);
    a.DivideUsing(dt);
    a.SaveDiv();
    for(int i=0;i<30;++i) for(int j=0;j<10;++j) a(i,j) = T(7-13*i+11*j);
    a.SubMatrix(0,10,0,10) += std::complex<T>(30,20);
    a.SubMatrix(10,20,0,10) -= std::complex<T>(50,-123);
    a.SubMatrix(20,30,0,10) += std::complex<T>(10,-75);
    a.SubMatrix(1,10,1,10) += std::complex<T>(99,100);
    a.SubMatrix(2,10,2,10) -= std::complex<T>(51,37);

    tmv::Vector<std::complex<T> > s(10);
    for(int i=0;i<10;++i) s(i) = T(i+2);

    eps = EPS*Norm(a)*Norm(a.Inverse());

    tmv::Vector<std::complex<T> > t = a * s;
    tmv::Vector<std::complex<T> > s2 = t/a;
    Assert(Norm(s2-s) < eps*Norm(s),"NonSquare t/a");

    tmv::Vector<std::complex<T> > t2 = s%a;
    s2 = t2*a;
    Assert(Norm(s2-s) < eps*Norm(s),"NonSquare t%a");

    // Test QR Update/Downdate:

    tmv::Matrix<std::complex<T>,stor> q30 = a;
    tmv::UpperTriMatrix<std::complex<T>,tmv::NonUnitDiag,stor> r30(10,10);
    QR_Decompose(q30.View(),r30.View());
    Assert(Norm(q30*r30-a) < eps*Norm(a),"QR_Decompose");
    Assert(Norm(r30.Adjoint()*r30-a.Adjoint()*a) < eps*r30.NormSq(),
           "QR_Decompose (RtR)");

    tmv::Matrix<std::complex<T>,stor> q10 = a.Rows(0,10);
    tmv::UpperTriMatrix<std::complex<T>,tmv::NonUnitDiag,stor> r10(10,10);
    QR_Decompose(q10.View(),r10.View());
    tmv::UpperTriMatrix<std::complex<T>,tmv::NonUnitDiag,stor> r = r10;
    tmv::Matrix<std::complex<T>,stor> a1030 = a.Rows(10,30);
    QR_Update(r.View(),a1030.View());
    Assert(Norm(r.Adjoint()*r-r30.Adjoint()*r30) < eps*r.NormSq(),
           "QR_Update");
    r = r30;

    a1030 = a.Rows(10,30);
    QR_Downdate(r.View(),a1030.View());
    Assert(Norm(r.Adjoint()*r-r10.Adjoint()*r10) < eps*r.NormSq(),
           "QR_Downdate");
    r = r10;

    tmv::Matrix<std::complex<T>,stor> a1020 = a.Rows(10,20);
    tmv::Matrix<std::complex<T>,stor> a2030 = a.Rows(20,30);
    QR_Update(r.View(),a1020.View());
    QR_Update(r.View(),a2030.View());
    Assert(Norm(r.Adjoint()*r-r30.Adjoint()*r30) < eps*r.NormSq(),
           "QR_Update (double)");
    r = r30;

    a1020 = a.Rows(10,20);
    a2030 = a.Rows(20,30);
    QR_Downdate(r.View(),a1020.View());
    QR_Downdate(r.View(),a2030.View());
    Assert(Norm(r.Adjoint()*r-r10.Adjoint()*r10) < eps*r.NormSq(),
           "QR_Downdate (double)");
    r = r10;

    tmv::Matrix<std::complex<T>,stor> q29 = a.Rows(0,29);
    tmv::UpperTriMatrix<std::complex<T>,tmv::NonUnitDiag,stor> r29(10,10);
    QR_Decompose(q29.View(),r29.View());
    r = r30;
    tmv::Vector<std::complex<T> > a29 = a.row(29);
    QR_Downdate(r.View(),a29.View());
    Assert(Norm(r.Adjoint()*r-r29.Adjoint()*r29) < eps*r.NormSq(),
           "QR_Downdate (single row)");
    r = r29;

    a29 = a.row(29);
    QR_Update(r.View(),a29.View());
    Assert(Norm(r.Adjoint()*r-r30.Adjoint()*r30) < eps*r.NormSq(),
           "QR_Downdate (single row)");

    // Test with some identical eigenvalues.
    // First make an arbitrary unitary matrix:
    tmv::Matrix<std::complex<T>,stor> q = a;
    q.DivideUsing(dt);
    q.SaveDiv();
    QR_Decompose(q.View(),r.View());
    r.Zero();
    r(0,0) = 1;
    r(1,1) = 5;
    r(2,2) = 1;
    r(3,3) = 1;
    r(4,4) = -3;
    r(5,5) = 7;
    r(6,6) = 1;
    r(7,7) = 1;
    r(8,8) = -2;
    r(9,9) = 1;
    q = q*r;

    eps = EPS * Norm(q) * Norm(q.Inverse());

    t = q * s;
    s2 = t/q;
    Assert(q.CheckDecomp(),"CheckDecomp");
    //Assert(q.CheckDecomp(dbgout),"CheckDecomp");
    Assert(Norm(s2-s) < eps*Norm(s),"NonSquare t/q");

    t2 = s%q;
    s2 = t2*q;
    Assert(Norm(s2-s) < eps*Norm(s),"NonSquare t%q");

    tmv::Matrix<T,stor> a1 = m;
    tmv::Matrix<T,stor> a2 = m.Transpose() * m;
    tmv::Matrix<T,stor> a3 = m * m.Transpose();
    a2.row(1) *= T(3);
    a2.col(2).AddToAll(-4);
    a3.row(5) *= T(7);
    a3.col(3).AddToAll(7);
    tmv::Matrix<std::complex<T>,stor> c1 = a1 * std::complex<T>(1,2);
    tmv::Matrix<std::complex<T>,stor> c2 = a2 * std::complex<T>(-3,4);
    tmv::Matrix<std::complex<T>,stor> c3 = a3 * std::complex<T>(-4,8);

    tmv::Matrix<T> a2x = a2;
    tmv::Matrix<T> a3x = a3;
    tmv::Matrix<std::complex<T> > c2x = c2;
    tmv::Matrix<std::complex<T> > c3x = c3;
    TestMatrixDivArith2<T>(dt,a2x,c2,a1.View(),a2.View(),c1.View(),c2.View(),
                           "NonSquare/Square"); 
    TestMatrixDivArith2<T>(dt,a3x,c3x,a1.View(),a3.View(),c1.View(),c3.View(),
                           "NonSquare/Square"); 

    tmv::Matrix<T,stor> a4(7,4);
    for(int i=0;i<7;++i) for(int j=0;j<4;++j) a4(i,j) = T(1-3*i+2*j);
    tmv::Matrix<T,stor> a5 = a4.Transpose();
    a4.SubMatrix(0,6,0,4) += a1;
    a5.SubMatrix(0,4,1,5) -= a2;
    tmv::Matrix<std::complex<T>,stor> c4 = a4*std::complex<T>(1,2);
    tmv::Matrix<std::complex<T>,stor> c5 = c4.Adjoint();
    c4.SubMatrix(0,6,0,4) += c1;
    c5.SubMatrix(0,4,1,5) -= c2;
    c4.col(1) *= std::complex<T>(2,1);
    c4.row(2).AddToAll(std::complex<T>(-7,2));
    c5.col(3) *= std::complex<T>(-1,3);
    c5.row(0).AddToAll(std::complex<T>(1,9));

    tmv::Matrix<T,stor> a6(9,6);
    for(int i=0;i<9;++i) for(int j=0;j<6;++j) a6(i,j) = T(5+2*i-2*j);
    tmv::Matrix<T,stor> a7 = a6.Transpose();
    a6.SubMatrix(2,8,1,5) += a1;
    a7.SubMatrix(0,6,4,8) -= T(2)*a1;
    tmv::Matrix<std::complex<T>,stor> c6 = a6*std::complex<T>(1,2);
    tmv::Matrix<std::complex<T>,stor> c7 = c6.Adjoint();
    c6.SubMatrix(2,8,1,5) += c1;
    c7.SubMatrix(0,6,4,8) -= T(2)*c1;
    c6.col(1) *= std::complex<T>(2,1);
    c6.row(5).AddToAll(std::complex<T>(-7,2));
    c7.col(7) *= std::complex<T>(-1,3);
    c7.row(4).AddToAll(std::complex<T>(1,9));

    tmv::Matrix<T> a4x = a4;
    tmv::Matrix<T> a5x = a5;
    tmv::Matrix<T> a6x = a6;
    tmv::Matrix<T> a7x = a7;
    tmv::Matrix<std::complex<T> > c4x = c4;
    tmv::Matrix<std::complex<T> > c5x = c5;
    tmv::Matrix<std::complex<T> > c6x = c6;
    tmv::Matrix<std::complex<T> > c7x = c7;
    TestMatrixDivArith1<T>(dt,a4x,c4x,a1.View(),a4.View(),c1.View(),c4.View(),
                           "NonSquare/NonSquare");
    TestMatrixDivArith1<T>(dt,a5x,c5x,a1.View(),a5.View(),c1.View(),c5.View(),
                           "NonSquare/NonSquare");
    TestMatrixDivArith1<T>(dt,a6x,c6x,a1.View(),a6.View(),c1.View(),c6.View(),
                           "NonSquare/NonSquare");
    TestMatrixDivArith1<T>(dt,a7x,c7x,a1.View(),a7.View(),c1.View(),c7.View(),
                           "NonSquare/NonSquare");

#ifdef XTEST
    tmv::Matrix<T,stor> a8(4,0);
    tmv::Matrix<T,stor> a9(0,4);
    tmv::Matrix<T,stor> a10(6,0);
    tmv::Matrix<T,stor> a11(0,6);
    tmv::Matrix<std::complex<T>,stor> c8 = a8;
    tmv::Matrix<std::complex<T>,stor> c9 = a9;
    tmv::Matrix<std::complex<T>,stor> c10 = a10;
    tmv::Matrix<std::complex<T>,stor> c11 = a11;

    tmv::Matrix<T> a8x = a8;
    tmv::Matrix<T> a9x = a9;
    tmv::Matrix<T> a10x = a10;
    tmv::Matrix<T> a11x = a11;
    tmv::Matrix<std::complex<T> > c8x = c8;
    tmv::Matrix<std::complex<T> > c9x = c9;
    tmv::Matrix<std::complex<T> > c10x = c10;
    tmv::Matrix<std::complex<T> > c11x = c11;
    TestMatrixDivArith1<T>(dt,a8,c8,a1.View(),a8.View(),c1.View(),c8.View(),
                           "NonSquare/Degenerate");
    TestMatrixDivArith1<T>(dt,a9,c9,a1.View(),a9.View(),c1.View(),c9.View(),
                           "NonSquare/Degenerate");
    TestMatrixDivArith1<T>(dt,a10,c10,a1.View(),a10.View(),c1.View(),c10.View(),
                           "NonSquare/Degenerate");
    TestMatrixDivArith1<T>(dt,a11,c11,a1.View(),a11.View(),c1.View(),c11.View(),
                           "NonSquare/Degenerate");
#endif

    if (stor == tmv::ColMajor) {
#ifdef XTEST
        TestNonSquareDiv<T,tmv::RowMajor>(dt);
    } else {
#endif
        std::cout<<"NonSquare Matrix<"<<tmv::TMV_Text(T())<<"> Division using ";
        std::cout<<tmv::TMV_Text(dt)<<" passed all tests\n";
    }
}

template <class T, tmv::StorageType stor> 
static void TestSingularDiv(tmv::DivType dt)
{
    tmv::Matrix<T,stor> m(4,4);
    for(int i=0;i<4;++i) for(int j=0;j<4;++j) m(i,j) = T(2+4*i-5*j);
    m(2,2) = 30;

    tmv::Vector<T> x(4);
    x(0) = 2;
    x(1) = -10;
    x(2) = 5;
    x(3) = -5;

    m.DivideUsing(dt);
    m.SaveDiv();
    m.SetDiv();
    std::ostream* dbgout = showdiv ? &std::cout : 0;
    Assert(m.CheckDecomp(dbgout),"CheckDecomp");

    T eps = EPS * Norm(m) * Norm(m.Inverse());

    tmv::Vector<T> b = m * x;
    tmv::Vector<T> x2 = b/m;
    tmv::Vector<T> b2 = m*x2;
    Assert(Norm(b2-b) < eps*Norm(b),"Singular exact b/m");

    b = x * m;
    x2 = b%m;
    b2 = x2*m;
    Assert(Norm(b2-b) < eps*Norm(b),"Singular exact b%m");

    b(0) += 10;
    x = b/m;
    b2 = m*x;
    T refnorm = Norm(b2-b);
    T dxval = T(std::sqrt(eps))*Norm(x);
    tmv::Vector<T> dx = dxval*tmv::BasisVector<T>(4,0);
    b2 = m*(x+dx);
    Assert(Norm(b2-b) >= refnorm,"Singular Least Squares b/m (1)");
    b2 = m*(x-dx);
    Assert(Norm(b2-b) >= refnorm,"Singular Least Squares b/m (2)");
    dx = dxval*tmv::BasisVector<T>(4,1);
    b2 = m*(x+dx);
    Assert(Norm(b2-b) >= refnorm,"Singular Least Squares b/m (3)");
    b2 = m*(x-dx);
    Assert(Norm(b2-b) >= refnorm,"Singular Least Squares b/m (4)");
    dx = dxval*tmv::BasisVector<T>(4,2);
    b2 = m*(x+dx);
    Assert(Norm(b2-b) >= refnorm,"Singular Least Squares b/m (5)");
    b2 = m*(x-dx);
    Assert(Norm(b2-b) >= refnorm,"Singular Least Squares b/m (6)");
    dx = dxval*tmv::BasisVector<T>(4,3);
    b2 = m*(x+dx);
    Assert(Norm(b2-b) >= refnorm,"Singular Least Squares b/m (7)");
    b2 = m*(x-dx);
    Assert(Norm(b2-b) >= refnorm,"Singular Least Squares b/m (8)");

    tmv::Matrix<T> minv = m.Inverse();

    if (showacc) {
        std::cout<<"m = "<<m<<std::endl;
        std::cout<<"minv = "<<minv<<std::endl;
        std::cout<<"minv*m = "<<minv*m<<std::endl;
        std::cout<<"m*minv = "<<m*minv<<std::endl;
        std::cout<<"m*minv*m = "<<m*minv*m<<std::endl;
        std::cout<<"Norm(m*minv*m-m) = "<<Norm(m*minv*m-m)<<std::endl;
        std::cout<<"minv*m*minv = "<<minv*m*minv<<std::endl;
        std::cout<<"Norm(minv*m*minv-minv) = "<<
            Norm(minv*m*minv-minv)<<std::endl;
        std::cout<<"m*minv-(m*minv)T = "<<m*minv-Transpose(m*minv)<<std::endl;
        std::cout<<"Norm(m*minv-(m*minv)T) = "<<
            Norm(m*minv-(m*minv).Transpose())<<std::endl;
        std::cout<<"minv*m-(minv*m)T = "<<minv*m-Transpose(minv*m)<<std::endl;
        std::cout<<"Norm(minv*m-(minv*m)T) = "<<
            Norm(minv*m-(minv*m).Transpose())<<std::endl;
    }

    Assert(Norm(m*minv*m - m) < eps*Norm(m),"Singular Inverse M*X*M != M");
    Assert(Norm(minv*m*minv - minv) < eps*Norm(minv),
           "Singular Inverse X*M*X != X");
    Assert(Norm((m*minv)-(m*minv).Transpose()) < eps,
           "Singular Inverse M*X != (M*X)T");
    if (dt != tmv::QRP) { // QRP doesn't get this right.
        Assert(Norm((minv*m)-(minv*m).Transpose()) < eps,
               "Singular Inverse X*M != (X*M)T");
    }

    // Try big one with many singular values.
    tmv::Matrix<T,stor> mm(30,30);
    mm.DivideUsing(dt);
    mm.SaveDiv();
    for(int i=0;i<30;++i) for(int j=0;j<30;++j) mm(i,j) = T(4-17*i+23*j);
    mm(20,20) += 200;
    mm(12,12) += 500;
    mm(7,7) += 300;
    mm(28,28) += 700;
    mm(24,24) += 400;

    eps = EPS * Norm(mm) * Norm(mm.Inverse());

    tmv::Vector<T> xx(30);
    for(int i=0;i<30;++i) xx(i) = T(10+i);
    tmv::Vector<T> bb = mm*xx;
    tmv::Vector<T> xx2 = bb/mm;
    tmv::Vector<T> bb2 = mm*xx2;
    //Assert(mm.CheckDecomp(),"CheckDecomp");
    Assert(mm.CheckDecomp(dbgout),"CheckDecomp");
    if (showacc) {
        std::cout<<"Norm(bb2-bb) = "<<Norm(bb2-bb);
        std::cout<<", EPS*Norm(bb) = "<<eps*Norm(bb)<<std::endl;
    }
    Assert(Norm(bb2-bb) < eps*Norm(bb),"Singular exact bb/mm");

    bb = xx * mm;
    xx2 = bb%mm;
    bb2 = xx2*mm;
    if (showacc) {
        std::cout<<"Norm(bb2-bb) = "<<Norm(bb2-bb);
        std::cout<<", EPS*Norm(mm)*Norm(xx) = "<<eps*Norm(bb)<<std::endl;
    }
    Assert(Norm(bb2-bb) < eps*Norm(bb),"Singular exact bb%mm");

    // Similar, but complex
    tmv::Matrix<std::complex<T>,stor> cc(30,30);
    cc.DivideUsing(dt);
    cc.SaveDiv();
    for(int i=0;i<30;++i) for(int j=0;j<30;++j) cc(i,j) = T(4-17*i+23*j);
    cc(20,20) += std::complex<T>(200,-999);
    cc(12,12) += std::complex<T>(500,-104);
    cc(7,7) += std::complex<T>(300,123);
    cc(28,28) += std::complex<T>(700,231);
    cc(24,24) += std::complex<T>(400,-120);

    eps = EPS * Norm(cc) * Norm(cc.Inverse());

    tmv::Vector<std::complex<T> > cxx(30);
    for(int i=0;i<30;++i) cxx(i) = T(10+i);
    tmv::Vector<std::complex<T> > cbb = cc*cxx;
    tmv::Vector<std::complex<T> > cxx2 = cbb/cc;
    tmv::Vector<std::complex<T> > cbb2 = cc*cxx2;
    //Assert(cc.CheckDecomp(),"CheckDecomp");
    Assert(cc.CheckDecomp(dbgout),"CheckDecomp");
    if (showacc) {
        std::cout<<"Norm(cbb2-cbb) = "<<Norm(cbb2-cbb);
        std::cout<<", EPS*Norm(cbb) = "<<eps*Norm(cbb)<<std::endl;
    }
    Assert(Norm(cbb2-cbb) < eps*Norm(cbb),"Singular exact bb/cc");

    cbb = cxx * cc;
    cxx2 = cbb%cc;
    cbb2 = cxx2*cc;
    if (showacc) {
        std::cout<<"Norm(cbb2-cbb) = "<<Norm(cbb2-cbb);
        std::cout<<", EPS*Norm(cc)*Norm(cxx) = "<<eps*Norm(cbb)<<std::endl;
    }
    Assert(Norm(cbb2-cbb) < eps*Norm(cbb),"Singular exact bb%cc");

    if (stor == tmv::ColMajor) {
#ifdef XTEST
        TestSingularDiv<T,tmv::RowMajor>(dt);
    } else {
#endif
        std::cout<<"Singular Matrix<"<<tmv::TMV_Text(T())<<"> Division using ";
        std::cout<<tmv::TMV_Text(dt)<<" passed all tests\n";
    }
}

template <class T> 
void TestAllMatrixDiv()
{
    TestMatrixDecomp<T,tmv::ColMajor>();
    TestMatrixDecomp<T,tmv::RowMajor>();
    std::cout<<"Matrix<"<<tmv::TMV_Text(T())<<"> passed all ";
    std::cout<<"decomposition tests.\n";
    TestSquareDiv<T,tmv::ColMajor>(tmv::LU);
    TestSquareDiv<T,tmv::ColMajor>(tmv::QR);
    TestSquareDiv<T,tmv::ColMajor>(tmv::QRP);
    TestSquareDiv<T,tmv::ColMajor>(tmv::SV);
    TestNonSquareDiv<T,tmv::ColMajor>(tmv::QR);
    TestNonSquareDiv<T,tmv::ColMajor>(tmv::QRP);
    TestNonSquareDiv<T,tmv::ColMajor>(tmv::SV);
    TestSingularDiv<T,tmv::ColMajor>(tmv::QRP);
    TestSingularDiv<T,tmv::ColMajor>(tmv::SV);
}

#ifdef INST_DOUBLE
template void TestAllMatrixDiv<double>();
#endif
#ifdef INST_FLOAT
template void TestAllMatrixDiv<float>();
#endif
#ifdef INST_LONGDOUBLE
template void TestAllMatrixDiv<long double>();
#endif
