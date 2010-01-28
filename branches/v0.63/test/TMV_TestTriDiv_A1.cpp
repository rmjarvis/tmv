#include "TMV.h"
#include "TMV_Test.h"
#include "TMV_Test1.h"

template <class T1, class T2> 
inline bool CanLDivEq(
    const tmv::UpperTriMatrixView<T1>& a, const tmv::UpperTriMatrixView<T2>& b)
{ return a.size() == b.size() && !a.isunit(); }

template <class T1, class T2> 
inline bool CanRDivEq(
    const tmv::UpperTriMatrixView<T1>& a, const tmv::UpperTriMatrixView<T2>& b)
{ return a.size() == b.size() && !a.isunit(); }

template <class T1, class T2> 
inline bool CanLDivEq(
    const tmv::LowerTriMatrixView<T1>& a, const tmv::LowerTriMatrixView<T2>& b)
{ return a.size() == b.size() && !a.isunit(); }

template <class T1, class T2> 
inline bool CanRDivEq(
    const tmv::LowerTriMatrixView<T1>& a, const tmv::LowerTriMatrixView<T2>& b)
{ return a.size() == b.size() && !a.isunit(); }

#include "TMV_TestMatrixDivArith.h"

template <class T, tmv::DiagType D> 
void TestTriDiv() 
{
    const int N = 10;

    tmv::Matrix<T> m(N,N);
    for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) 
        m(i,j) = T(0.4+0.02*i-0.05*j);
    m.diag().addToAll(5);
    m.diag(1).addToAll(T(0.32));
    m.diag(-1).addToAll(T(0.91));

    tmv::UpperTriMatrix<T,D> a(m);
    m = a;
    a.saveDiv();
    m.saveDiv();

    tmv::Vector<T> b(N);
    for (int i=0;i<N;++i) b(i) = T(i+7);

    a.setDiv();
    m.divideUsing(tmv::LU);
    m.setDiv();

    if (showacc) {
        std::cout<<"b = "<<b<<std::endl;
        std::cout<<"a = "<<a<<std::endl;
        std::cout<<"m = "<<m<<std::endl;
    }

    T eps = EPS * Norm(m) * Norm(m.inverse());

    tmv::Vector<T> x1 = b/a;
    tmv::Vector<T> x2 = b/m;
    if (showacc) {
        std::cout<<"x1 = "<<x1<<std::endl;
        std::cout<<"x2 = "<<x2<<std::endl;
        std::cout<<"a*x1-b = "<<a*x1-b<<std::endl;
        std::cout<<"m*x2-b = "<<m*x2-b<<std::endl;
        std::cout<<"Norm(x1-x2) = "<<Norm(x1-x2)<<std::endl;
        std::cout<<"EPS*Norm(x1) = "<<eps*Norm(x1)<<std::endl;
    }
    Assert(Norm(x1-x2) < eps*Norm(x1),"Tri b/a");

    x1 = b%a;
    x2 = b%m;
    if (showacc) {
        std::cout<<"x1 = "<<x1<<std::endl;
        std::cout<<"x2 = "<<x2<<std::endl;
        std::cout<<"x1*a-b = "<<x1*a-b<<std::endl;
        std::cout<<"x2*m-b = "<<x2*m-b<<std::endl;
        std::cout<<"Norm(x1-x2) = "<<Norm(x1-x2)<<std::endl;
        std::cout<<"EPS*Norm(b) = "<<eps*Norm(x1)<<std::endl;
    }
    Assert(Norm(x1-x2) < eps*Norm(x1),"Tri b%a");

    tmv::UpperTriMatrix<T,D> ainv = a.inverse();
    tmv::Matrix<T> minv = m.inverse();
    if (showacc) {
        std::cout<<"ainv = "<<ainv<<std::endl;
        std::cout<<"minv = "<<minv<<std::endl;
        std::cout<<"Norm(ainv-minv) = "<<Norm(ainv-minv)<<std::endl;
        std::cout<<"EPS*Norm(ainv) = "<<eps*Norm(ainv)<<std::endl;
    }
    Assert(Norm(ainv-minv) < eps*Norm(ainv),"Tri Inverse");

    if (showacc) {
        std::cout<<"Det(a) = "<<Det(a)<<", Det(m) = "<<Det(m)<<std::endl;
        std::cout<<"abs(adet-mdet) = "<<std::abs(Det(a)-Det(m));
        std::cout<<"  EPS*abs(mdet) = "<<eps*std::abs(Det(m))<<std::endl;
        std::cout<<"LogDet(a) = "<<LogDet(a);
        std::cout<<", LogDet(m) = "<<LogDet(m)<<std::endl;
    }
    Assert(std::abs(Det(m)-Det(a)) < eps*std::abs(Det(m)),"Tri Det");
    T asign, msign;
    Assert(std::abs(m.logDet(&msign)-a.logDet(&asign)) < N*eps,"Tri LogDet");
    Assert(std::abs(asign-msign) < N*eps,"Tri LogDet - sign");
    Assert(std::abs(Det(a) - asign*std::exp(LogDet(a))) < 
           eps*std::abs(Det(m)),"Tri Det--LogDet");

    tmv::Matrix<std::complex<T> > cm(m);
    cm += std::complex<T>(10,2);
    cm.diag(1) *= std::complex<T>(T(-0.5),T(-0.8));
    cm.diag(-1) *= std::complex<T>(T(-0.7),T(0.1));

    tmv::UpperTriMatrix<std::complex<T>,D> ca(cm);
    cm = ca;
    ca.saveDiv();
    cm.saveDiv();

    cm.divideUsing(tmv::LU);
    cm.setDiv();
    ca.setDiv();

    T ceps = EPS * Norm(cm) * Norm(cm.inverse());

    if (showacc) {
        std::cout<<"Det(ca) = "<<Det(ca)<<", Det(cm) = "<<Det(cm)<<std::endl;
        std::cout<<"abs(cadet-cmdet) = "<<std::abs(Det(ca)-Det(cm));
        std::cout<<"  EPS*abs(cmdet) = "<<ceps*std::abs(Det(cm))<<std::endl;
    }
    Assert(std::abs(Det(ca)-Det(cm)) < ceps*std::abs(Det(cm)),"Tri CDet");
    std::complex<T> casign, cmsign;
    Assert(std::abs(cm.logDet(&cmsign)-ca.logDet(&casign)) < N*eps,
           "Tri CLogDet");
    Assert(std::abs(casign-cmsign) < N*eps,"Tri CLogDet - sign");
    Assert(std::abs(Det(ca) - casign*std::exp(LogDet(ca))) < 
           eps*std::abs(Det(cm)),"Tri CDet--LogDet");

    tmv::Vector<std::complex<T> > e(b);
    e(1) += std::complex<T>(-1,5);
    e(2) -= std::complex<T>(-1,5);

    tmv::Vector<std::complex<T> > y1 = b/ca;
    tmv::Vector<std::complex<T> > y2 = b/cm;
    if (showacc) {
        std::cout<<"y1 = "<<y1<<std::endl;
        std::cout<<"y2 = "<<y2<<std::endl;
        std::cout<<"Norm(diff) = "<<Norm(y1-y2)<<std::endl;
        std::cout<<"EPS*Norm(y1) = "<<ceps*Norm(y1)<<std::endl;
    }
    Assert(Norm(y1-y2) < ceps*Norm(y1),"Tri b/c");
    y1 = b%ca;
    y2 = b%cm;
    if (showacc) {
        std::cout<<"y1 = "<<y1<<std::endl;
        std::cout<<"y2 = "<<y1<<std::endl;
        std::cout<<"Norm(diff) = "<<Norm(y1-y2)<<std::endl;
        std::cout<<"EPS*Norm(y1) = "<<ceps*Norm(y1)<<std::endl;
    }
    Assert(Norm(y1-y2) < ceps*Norm(y1),"Tri b%c");

    // test complex / real
    y1 = e/a;
    y2 = e/m;
    Assert(Norm(y1-y2) < eps*Norm(y1),"Tri e/m");
    y1 = e%a;
    y2 = e%m;
    Assert(Norm(y1-y2) < eps*Norm(y1),"Tri e%m");

    // test complex / complex
    y1 = e/ca;
    y2 = e/cm;
    Assert(Norm(y1-y2) < ceps*Norm(y1),"Tri e/c");
    y1 = e%ca;
    y2 = e%cm;
    Assert(Norm(y1-y2) < ceps*Norm(y1),"Tri e%c");
}

template <class T> 
void TestTriDiv_A1() 
{
    const int N = 10;

    tmv::Matrix<T> m(N,N);
    for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) 
        m(i,j) = T(0.4+0.02*i-0.05*j);
    m.diag().addToAll(5);
    m.diag(1).addToAll(T(0.32));
    m.diag(-1).addToAll(T(0.91));

    tmv::Matrix<std::complex<T> > cm(m);
    cm += std::complex<T>(10,2);
    cm.diag(1) *= std::complex<T>(T(-0.5),T(-0.8));
    cm.diag(-1) *= std::complex<T>(T(-0.7),T(0.1));

    tmv::UpperTriMatrix<T,tmv::NonUnitDiag> a1(m);
    tmv::UpperTriMatrix<std::complex<T>,tmv::NonUnitDiag> ca1(cm);
    tmv::UpperTriMatrix<T,tmv::UnitDiag> a2(m);
    tmv::UpperTriMatrix<std::complex<T>,tmv::UnitDiag> ca2(cm);
    a1.saveDiv();
    a2.saveDiv();
    ca1.saveDiv();
    ca2.saveDiv();
    a1.setDiv();
    a2.setDiv();
    ca1.setDiv();
    ca2.setDiv();

    tmv::UpperTriMatrix<T,tmv::NonUnitDiag> a1x(m);
    tmv::UpperTriMatrix<std::complex<T>,tmv::NonUnitDiag> ca1x(cm);
    tmv::UpperTriMatrix<T,tmv::UnitDiag> a2x(m);
    tmv::UpperTriMatrix<std::complex<T>,tmv::UnitDiag> ca2x(cm);
    tmv::LowerTriMatrix<T,tmv::NonUnitDiag> b1x(m);
    tmv::LowerTriMatrix<std::complex<T>,tmv::NonUnitDiag> cb1x(cm);
    tmv::LowerTriMatrix<T,tmv::UnitDiag> b2x(m);
    tmv::LowerTriMatrix<std::complex<T>,tmv::UnitDiag> cb2x(cm);

    TestMatrixDivArith2<T>(tmv::LU,a2x,ca2x,a1.view(),a2.view(),
                           ca1.view(),ca2.view(),"U/U 1");
    TestMatrixDivArith2<T>(tmv::LU,b2x,cb2x,a1.transpose(),a2.transpose(),
                           ca1.transpose(),ca2.transpose(),"L/L 1");
    TestMatrixDivArith2<T>(tmv::LU,a1x,ca1x,a2.view(),a1.view(),
                           ca2.view(),ca1.view(),"U/U 2");
    TestMatrixDivArith2<T>(tmv::LU,b1x,cb1x,a2.transpose(),a1.transpose(),
                           ca2.transpose(),ca1.transpose(),"L/L 2");

#ifdef XTEST
    tmv::UpperTriMatrix<T,tmv::NonUnitDiag> a1b(m);
    tmv::UpperTriMatrix<std::complex<T>,tmv::NonUnitDiag> ca1b(cm);
    tmv::UpperTriMatrix<T,tmv::UnitDiag> a2b(m);
    tmv::UpperTriMatrix<std::complex<T>,tmv::UnitDiag> ca2b(cm);

    TestMatrixDivArith1<T>(tmv::LU,a1x,ca1x,a1.view(),a1b.view(),
                           ca1.view(),ca1b.view(),"U/U 3");
    TestMatrixDivArith1<T>(tmv::LU,b1x,cb1x,a1.transpose(),a1b.transpose(),
                           ca1.transpose(),ca1b.transpose(),"L/L 3");
    TestMatrixDivArith1<T>(tmv::LU,a2x,ca2x,a2.view(),a2b.view(),
                           ca2.view(),ca2b.view(),"U/U 4");
    TestMatrixDivArith1<T>(tmv::LU,b2x,cb2x,a2.transpose(),a2b.transpose(),
                           ca2.transpose(),ca2b.transpose(),"L/L 4");
#endif
}

template <class T> 
void TestAllTriDiv()
{
    TestTriDiv<T,tmv::NonUnitDiag>();
    TestTriDiv<T,tmv::UnitDiag>();
    TestTriDiv_A1<T>();
    TestTriDiv_A2<T>();
    TestTriDiv_B1<T>();
    TestTriDiv_B2<T>();
    TestTriDiv_C1<T>();
    TestTriDiv_C2<T>();
    std::cout<<"TriMatrix<"<<tmv::TMV_Text(T())<<"> Division passed all tests\n";
}

#ifdef INST_DOUBLE
template void TestAllTriDiv<double>();
#endif
#ifdef INST_FLOAT
template void TestAllTriDiv<float>();
#endif
#ifdef INST_LONGDOUBLE
template void TestAllTriDiv<long double>();
#endif
