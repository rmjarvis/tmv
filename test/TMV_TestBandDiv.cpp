
#define START 0

#include "TMV.h"
#include "TMV_Band.h"
#include "TMV_Test.h"
#include "TMV_Test_2.h"
#include "TMV_TestBandArith.h"

template <class T> 
void TestBandDiv(tmv::DivType dt)
{
    const int N = 10;

    std::vector<tmv::BandMatrixView<T> > b;
    std::vector<tmv::BandMatrixView<std::complex<T> > > cb;
    MakeBandList(b,cb);

    tmv::Matrix<T> a1(N,N);
    for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) a1(i,j) = T(3+i-2*j);
    a1.diag().addToAll(T(10)*N);
    a1 /= T(10);

    tmv::Matrix<std::complex<T> > ca1 = a1 * std::complex<T>(3,-4);

    tmv::Vector<T> v1(N);
    tmv::Vector<T> v2(N-1);
    for (int i=0; i<N; ++i) v1(i) = T(16-3*i);
    for (int i=0; i<N-1; ++i) v2(i) = T(-6+i); 
    tmv::Vector<std::complex<T> > cv1(N);
    tmv::Vector<std::complex<T> > cv2(N-1);
    for (int i=0; i<N; ++i) cv1(i) = std::complex<T>(16-3*i,i+4); 
    for (int i=0; i<N-1; ++i) cv2(i) = std::complex<T>(2*i-3,-6+i); 

    tmv::Matrix<T> a3 = a1.colRange(0,N/2);
    tmv::Matrix<std::complex<T> > ca3 = ca1.colRange(0,N/2);
    tmv::Matrix<T> a4 = a1.rowRange(0,N/2);
    tmv::Matrix<std::complex<T> > ca4 = ca1.rowRange(0,N/2);
    tmv::Matrix<T> a5 = a1.colRange(0,0);
    tmv::Matrix<std::complex<T> > ca5 = ca1.colRange(0,0);
    tmv::Matrix<T> a6 = a1.rowRange(0,0);
    tmv::Matrix<std::complex<T> > ca6 = ca1.rowRange(0,0);
    tmv::Matrix<T> a7 = a1;
    tmv::Matrix<std::complex<T> > ca7 = ca1;
    a7.diag().addToAll(T(10)*N);
    ca7.diag().addToAll(T(10)*N);

    for(size_t i=START;i<b.size();i++) {
        if (showstartdone) 
            std::cout<<"Start loop: i = "<<i<<"\nbi = "<<tmv::TMV_Text(b[i])<<
                "  "<<b[i]<<std::endl;
        tmv::BandMatrixView<T> bi = b[i];
        tmv::BandMatrixView<std::complex<T> > cbi = cb[i];
        if (dt == tmv::LU && !bi.isSquare()) continue;

        bi.saveDiv();
        cbi.saveDiv();

        tmv::Matrix<T> m(bi);
        m.saveDiv();
        bi.divideUsing(dt);
        bi.setDiv();
        m.divideUsing(dt);
        m.setDiv();

        std::ostream* divout = showdiv ? &std::cout : 0;
        Assert(bi.checkDecomp(divout),"CheckDecomp");
        T eps = m.rowsize()*EPS*Norm(m)*Norm(m.inverse());

        if (bi.colsize() == N) {
            tmv::Vector<T> x1 = v1/bi;
            tmv::Vector<T> x2 = v1/m;
            if (showacc) {
                std::cout<<"v/b: Norm(x1-x2) = "<<Norm(x1-x2)<<
                    "  "<<eps*Norm(x1)<<std::endl;
            }
            Assert(Norm(x1-x2) < eps*Norm(x1),"Band v/b");
        }

        if (bi.rowsize() == N) {
            tmv::Vector<T> x1 = v1%bi;
            tmv::Vector<T> x2 = v1%m;
            if (showacc) {
                std::cout<<"v%b: Norm(x1-x2) = "<<Norm(x1-x2)<<
                    "  "<<eps*Norm(x1)<<std::endl;
            }
            Assert(Norm(x1-x2) < eps*Norm(x1),"Band v%b");
        }

        tmv::Matrix<T,tmv::ColMajor> binv = bi.inverse();
        tmv::Matrix<T,tmv::ColMajor> minv = m.inverse();
        if (showacc) {
            std::cout<<"minv = "<<minv<<std::endl;
            std::cout<<"binv = "<<binv<<std::endl;
            std::cout<<"Norm(minv-binv) = "<<Norm(minv-binv)<<
                "  "<<eps*Norm(binv)<<std::endl;
        }
        Assert(Norm(binv-minv) < eps*Norm(binv),"Band Inverse");

        if (m.isSquare()) {
            if (showacc) {
                std::cout<<"Det(b) = "<<Det(bi)<<
                    ", Det(m) = "<<Det(m)<<std::endl;
                std::cout<<"abs(bdet-mdet) = "<<std::abs(Det(bi)-Det(m));
                std::cout<<"  EPS*abs(mdet) = "<<
                    eps*std::abs(Det(m))<<std::endl;
                std::cout<<"abs(abs(bdet)-abs(mdet)) = "<<
                    std::abs(std::abs(Det(bi))-std::abs(Det(m)));
                std::cout<<"  EPS*abs(mdet) = "<<
                    eps*std::abs(Det(m))<<std::endl;
            }
            Assert(std::abs(Det(m)-Det(bi)) < eps*std::abs(Det(m)+Norm(m)),
                   "Band Det");
            T msign, bsign;
            Assert(std::abs(m.logDet(&msign)-bi.logDet(&bsign)) < N*eps,
                   "Band LogDet");
            Assert(std::abs(msign-bsign) < N*eps,"Band LogDet - sign");
        }

        cbi.divideUsing(dt);
        cbi.setDiv();
        Assert(cbi.checkDecomp(divout),"CheckDecomp");

        tmv::Matrix<std::complex<T> > cm(cbi);
        cm.saveDiv();
        cm.divideUsing(dt);
        cm.setDiv();
        T ceps = EPS*Norm(cm)*Norm(cm.inverse());

        if (cm.isSquare()) {
            if (showacc) {
                std::cout<<"Det(cbi) = "<<Det(cbi)<<", Det(cm) = "<<
                    Det(cm)<<std::endl;
                std::cout<<"abs(cbidet-cmdet) = "<<std::abs(Det(cbi)-Det(cm));
                std::cout<<"  cbidet/cmdet = "<<Det(cbi)/Det(cm);
                std::cout<<"  EPS*abs(cmdet) = "<<
                    ceps*std::abs(Det(cm))<<std::endl;
                std::cout<<"abs(abs(bdet)-abs(mdet)) = "<<
                    std::abs(std::abs(Det(bi))-std::abs(Det(m)));
                std::cout<<"  EPS*abs(mdet) = "<<
                    ceps*std::abs(Det(m))<<std::endl;
            }
            Assert(std::abs(Det(cbi)-Det(cm)) < 
                   ceps*std::abs(Det(cm)+Norm(cm)),
                   "Band CDet");
            std::complex<T> cmsign, cbsign;
            Assert(std::abs(cm.logDet(&cmsign)-cbi.logDet(&cbsign)) < N*eps,
                   "Band CLogDet");
            Assert(std::abs(cmsign-cbsign) < N*eps,"Band CLogDet - sign");
        }

        tmv::Vector<std::complex<T> > cv(v1 * std::complex<T>(1,1));
        cv(1) += std::complex<T>(-1,5);
        cv(2) -= std::complex<T>(-1,5);

        if (m.colsize() == N) {
            // test real / complex
            tmv::Vector<std::complex<T> > y1 = v1/cbi;
            tmv::Vector<std::complex<T> > y2 = v1/cm;
            if (showacc) {
                std::cout<<"v/cb: Norm(y1-y2) = "<<Norm(y1-y2)<<
                    "  "<<ceps*Norm(y1)<<std::endl;
            }
            Assert(Norm(y1-y2) < ceps*Norm(y1),"Band v/cb");

            // test complex / real
            y1 = cv/bi;
            y2 = cv/m;
            if (showacc) {
                std::cout<<"cv/b: Norm(y1-y2) = "<<Norm(y1-y2)<<
                    "  "<<eps*Norm(y1)<<std::endl;
            }
            Assert(Norm(y1-y2) < eps*Norm(y1),"Band cv/b");

            // test complex / complex
            y1 = cv/cbi;
            y2 = cv/cm;
            if (showacc) {
                std::cout<<"cv/cb: Norm(y1-y2) = "<<Norm(y1-y2)<<
                    "  "<<ceps*Norm(y1)<<std::endl;
            }
            Assert(Norm(y1-y2) < ceps*Norm(y1),"Band cv/cb");
        }

        if (bi.rowsize() == N) {
            tmv::Vector<std::complex<T> > y1 = v1%cbi;
            tmv::Vector<std::complex<T> > y2 = v1%cm;
            if (showacc) {
                std::cout<<"v%cb: Norm(y1-y2) = "<<Norm(y1-y2)<<
                    "  "<<ceps*Norm(y1)<<std::endl;
            }
            Assert(Norm(y1-y2) < ceps*Norm(y1),"Band v%cb");

            y1 = cv%bi;
            y2 = cv%m;
            if (showacc) {
                std::cout<<"cv%b: Norm(y1-y2) = "<<Norm(y1-y2)<<
                    "  "<<eps*Norm(y1)<<std::endl;
            }
            Assert(Norm(y1-y2) < eps*Norm(y1),"Band cv%b");
            y1 = cv%cbi;
            y2 = cv%cm;
            if (showacc) {
                std::cout<<"cv%cb: Norm(y1-y2) = "<<Norm(y1-y2)<<
                    "  "<<ceps*Norm(y1)<<std::endl;
            }
            Assert(Norm(y1-y2) < ceps*Norm(y1),"Band cv%cb");
        }

    }

    TestBandDiv_A<T>(dt);
    TestBandDiv_B1<T>(dt);
    TestBandDiv_B2<T>(dt);
    TestBandDiv_C1<T>(dt);
    if (dt == tmv::LU) TestBandDiv_C2<T>(dt);
    TestBandDiv_D1<T>(dt);
    if (dt == tmv::LU) TestBandDiv_D2<T>(dt);

    std::cout<<"BandMatrix<"<<tmv::TMV_Text(T())<<"> Division using ";
    std::cout<<tmv::TMV_Text(dt)<<" passed all tests\n";
}

template <class T> void TestAllBandDiv()
{
    TestBandDecomp<T,tmv::ColMajor>();
    TestBandDecomp<T,tmv::RowMajor>();
    TestBandDecomp<T,tmv::DiagMajor>();
    std::cout<<"BandMatrix<"<<tmv::TMV_Text(T())<<"> passed all ";
    std::cout<<"decomposition tests.\n";
    TestBandDiv<T>(tmv::LU);
    TestBandDiv<T>(tmv::QR);
    TestBandDiv<T>(tmv::SV);
}

#ifdef TEST_DOUBLE
template void TestAllBandDiv<double>();
#endif
#ifdef TEST_FLOAT
template void TestAllBandDiv<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestAllBandDiv<long double>();
#endif

