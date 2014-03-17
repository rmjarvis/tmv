
#define START 0

#include "TMV.h"
#include "TMV_SymBand.h"
#include "TMV_Test.h"
#include "TMV_Test_2.h"
#include "TMV_TestSymBandArith.h"

#ifdef TMV_MEM_DEBUG
// See the discussion of this in TMV_TestTri.cpp.  But basically, there seems to be something
// in the std library exception class that doesn't interact well with the mmgr-style memory
// debugging.  So we skip those tests if we are doing MEM_DEBUG
#define NOTHROW
#endif

template <class T> 
static bool IsPosDef(const tmv::GenSymBandMatrix<T>& m)
{
#ifdef NOTHROW
    for(int i=1;i<=m.size();i++) {
        T d = Det(tmv::Matrix<T>(m.subSymBandMatrix(0,i)));
        if (tmv::TMV_REAL(d) <= 0) return false;
    }
    return true;
#else
    try {
        tmv::HermBandMatrix<T> m2 = m;
        CH_Decompose(m2.view());
    } catch (tmv::NonPosDef) {
#if (XTEST & 16)
        if (showacc)
            std::cout<<"caught\n";
        for(int i=1;i<=m.size();i++) {
            T d = Det(tmv::Matrix<T>(m.subSymBandMatrix(0,i)));
            if (showacc) 
                std::cout<<"Det(0.."<<i<<") = "<<d<<std::endl;
            if (!(tmv::TMV_REAL(d) > 0)) {
                return false;
            }
        }
        std::cout<<"m = "<<TMV_Text(m)<<"  "<<m<<std::endl;
        for(int i=1;i<=m.size();i++) {
            T d = Det(tmv::Matrix<T>(m.subSymBandMatrix(0,i)));
            std::cout<<"Det(0.."<<i<<") = "<<d<<std::endl;
        }
        Assert(tmv::TMV_FALSE,
               "Caught NonPosDef, but no determinants of sub-blocks are "
               "negative");
#endif
        return false;
    }
#if (XTEST & 16)
    if (showacc)
        std::cout<<"not caught\n";
    for(int i=1;i<=m.size();i++) {
        T d = Det(tmv::Matrix<T>(m.subSymBandMatrix(0,i)));
        if (showacc) 
            std::cout<<"Det(0.."<<i<<") = "<<d<<std::endl;
        if (tmv::TMV_REAL(d) < 0) {
            std::cout<<"m = "<<TMV_Text(m)<<"  "<<m<<std::endl;
            std::cout<<"Det(0.."<<i<<") = "<<d<<std::endl;
            Assert(tmv::TMV_FALSE,
                   "Didn't catch NonPosDef, but determinant of sub-block "
                   "is negative");
        }
    }
#endif
    return true;
#endif
}

template <class T> 
void TestSymBandDiv(tmv::DivType dt, PosDefCode pdc)
{
    const int N = 10;

    std::vector<tmv::SymBandMatrixView<T> > sb;
    std::vector<tmv::SymBandMatrixView<std::complex<T> > > csb;
    MakeSymBandList(sb,csb,pdc);

    tmv::Vector<T> v1(N);
    tmv::Vector<T> v2(N-1);
    for (int i=0; i<N; ++i) v1(i) = T(16-3*i); 
    for (int i=0; i<N-1; ++i) v2(i) = T(-7+2*i); 

    std::ostream* checkout = showdiv ? &std::cout : 0;

    for(size_t i=START;i<sb.size();i++) {
        if (i>START) break;
        tmv::SymBandMatrixView<T> si = sb[i];
        tmv::SymBandMatrixView<std::complex<T> > csi = csb[i];
        if (dt == tmv::CH && csi.issym()) continue;
        si.saveDiv();
        csi.saveDiv();
        if (showstartdone)
            std::cout<<"Start loop: i = "<<i<<", si = "<<tmv::TMV_Text(si)<<
                "  "<<si<<std::endl;

        Assert(IsPosDef(si) == (pdc==PosDef),"IsPosDef");
        if (csi.isherm()) {
            Assert(IsPosDef(csi) == (pdc==PosDef),"IsPosDef");
        }

        tmv::Matrix<T> m(si);
        m.saveDiv();
        if (dt == tmv::CH) m.divideUsing(tmv::LU);
        else m.divideUsing(dt);
        m.setDiv();
        Assert(m.checkDecomp(checkout),"CheckDecomp m"); 
        T eps = EPS;
        if (pdc == Sing) eps *= 1000;
        else eps *= Norm(m)*Norm(m.inverse());
        si.divideUsing(dt);
        si.setDiv();
        if (si.isherm()) {
            tmv::HermMatrix<T> six = si;
            six.divideUsing(dt);
            six.setDiv();
            Assert(six.checkDecomp(checkout),"CheckDecomp six(herm)"); 
        } else {
            tmv::SymBandMatrix<T> six = si;
            six.divideUsing(dt);
            six.setDiv();
            Assert(six.checkDecomp(checkout),"CheckDecomp six(sym)"); 
        }
        Assert(si.checkDecomp(checkout),"CheckDecomp si"); 

        tmv::Vector<T> x1 = v1/si;
        tmv::Vector<T> x2 = v1/m;
        if (showacc) {
            std::cout<<"v1 = "<<v1<<std::endl;
            std::cout<<"v1/si = "<<x1<<std::endl;
            std::cout<<"v1/m = "<<x2<<std::endl;
            std::cout<<"si*x1 = "<<si*x1<<std::endl;
            std::cout<<"m*x2 = "<<m*x2<<std::endl;
            std::cout<<"v/b: Norm(x1-x2) = "<<Norm(x1-x2)<<
                "  "<<eps*Norm(x1)<<std::endl;
        }
        Assert(Norm(x1-x2) < eps*Norm(x1),"SymBand v/b");

        x1 = v1%si;
        x2 = v1%m;
        if (showacc)
            std::cout<<"v%b: Norm(x1-x2) = "<<Norm(x1-x2)<<
                "  "<<eps*Norm(x1)<<std::endl;
        Assert(Norm(x1-x2) < eps*Norm(x1),"SymBand v%b");

        tmv::Matrix<T,tmv::ColMajor> sinv = si.inverse();
        tmv::Matrix<T,tmv::ColMajor> minv = m.inverse();
        if (showacc) {
            std::cout<<"sinv = "<<sinv<<std::endl;
            std::cout<<"minv = "<<minv<<std::endl;
            std::cout<<"Norm(minv-sinv) = "<<Norm(minv-sinv)<<
                "  "<<eps*Norm(sinv)<<std::endl;
        }
        Assert(Norm(sinv-minv) < eps*Norm(sinv),"SymBand Inverse");

        if (showacc) {
            std::cout<<"Det(si) = "<<Det(si)<<", Det(m) = "<<Det(m)<<std::endl;
            std::cout<<"abs(sdet-mdet) = "<<std::abs(Det(si)-Det(m));
            std::cout<<"  EPS*abs(mdet) = "<<eps*std::abs(Det(m))<<std::endl;
            std::cout<<"abs(abs(sdet)-abs(mdet)) = "<<
                std::abs(std::abs(Det(si))-std::abs(Det(m)));
            std::cout<<"  EPS*abs(mdet) = "<<eps*std::abs(Det(m))<<std::endl;
        }
        if (pdc != Sing) {
            Assert(std::abs(Det(m)-Det(si)) < eps*std::abs(Det(m)),
                   "SymBand Det");
            T msign,ssign;
            Assert(std::abs(m.logDet(&msign)-si.logDet(&ssign)) < 10*N*eps,
                   "SymBand LogDet");
            Assert(std::abs(msign-ssign) < 10*N*eps, "SymBand LogDet - sign");
        }

        tmv::Matrix<std::complex<T> > cm(csi);
        cm.saveDiv();
        if (dt == tmv::CH) cm.divideUsing(tmv::LU);
        else cm.divideUsing(dt);
        cm.setDiv();
        Assert(cm.checkDecomp(checkout),"CheckDecomp cm"); 
        csi.divideUsing(dt);
        csi.setDiv();
        if (csi.isherm()) {
            tmv::HermMatrix<std::complex<T> > csix = csi;
            csix.divideUsing(dt);
            csix.setDiv();
            Assert(csix.checkDecomp(checkout),"CheckDecomp csix(herm)"); 
        } else {
            tmv::SymBandMatrix<std::complex<T> > csix = csi;
            csix.divideUsing(dt);
            csix.setDiv();
            Assert(csix.checkDecomp(checkout),"CheckDecomp csix(sym)"); 
        }
        Assert(csi.checkDecomp(checkout),"CheckDecomp csi"); 

        T ceps = EPS;
        if (pdc == Sing) ceps *= 1000;
        else ceps *= Norm(cm)*Norm(cm.inverse());

        if (showacc) {
            std::cout<<"Det(csi) = "<<Det(csi)<<
                ", Det(cm) = "<<Det(cm)<<std::endl;
            std::cout<<"abs(csidet-cmdet) = "<<std::abs(Det(csi)-Det(cm));
            std::cout<<"  csidet/cmdet = "<<Det(csi)/Det(cm);
            std::cout<<"  EPS*abs(cmdet) = "<<
                ceps*std::abs(Det(cm))<<std::endl;
            std::cout<<"abs(abs(csdet)-abs(cmdet)) = "<<
                std::abs(std::abs(Det(csi))-std::abs(Det(cm)));
            std::cout<<"  EPS*abs(cmdet) = "<<
                ceps*std::abs(Det(cm))<<std::endl;
        }
        if (pdc != Sing) {
            Assert(std::abs(Det(csi)-Det(cm)) < ceps*std::abs(Det(cm)),
                   "SymBand CDet");
            std::complex<T> cmsign,cssign;
            Assert(std::abs(cm.logDet(&cmsign)-csi.logDet(&cssign)) < 10*N*eps,
                   "SymBand CLogDet");
            Assert(std::abs(cmsign-cssign) < 10*N*eps, 
                   "SymBand CLogDet - sign");
        }

        tmv::Vector<std::complex<T> > cv = v1 * std::complex<T>(1,1);
        cv(1) += std::complex<T>(-1,5);
        cv(2) -= std::complex<T>(-1,5);

        // test real / complex
        tmv::Vector<std::complex<T> > y1 = v1/csi;
        tmv::Vector<std::complex<T> > y2 = v1/cm;
        if (showacc) 
            std::cout<<"v/cs: Norm(y1-y2) = "<<Norm(y1-y2)<<"  "<<
                ceps*Norm(y1)<<std::endl;
        Assert(Norm(y1-y2) < ceps*Norm(y1),"SymBand v/cs");

        // test complex / real
        y1 = cv/si;
        y2 = cv/m;
        if (showacc) 
            std::cout<<"cv/b: Norm(y1-y2) = "<<Norm(y1-y2)<<"  "<<
                eps*Norm(y1)<<std::endl;
        Assert(Norm(y1-y2) < eps*Norm(y1),"SymBand cv/b");

        // test complex / complex
        y1 = cv/csi;
        y2 = cv/cm;
        if (showacc) 
            std::cout<<"cv/cs: Norm(y1-y2) = "<<Norm(y1-y2)<<"  "<<
                ceps*Norm(y1)<<std::endl;
        Assert(Norm(y1-y2) < ceps*Norm(y1),"SymBand cv/cs");

        y1 = v1%csi;
        y2 = v1%cm;
        if (showacc) 
            std::cout<<"v%cs: Norm(y1-y2) = "<<Norm(y1-y2)<<"  "<<
                ceps*Norm(y1)<<std::endl;
        Assert(Norm(y1-y2) < ceps*Norm(y1),"SymBand v%cs");

        y1 = cv%si;
        y2 = cv%m;
        if (showacc) 
            std::cout<<"cv%b: Norm(y1-y2) = "<<Norm(y1-y2)<<"  "<<
                eps*Norm(y1)<<std::endl;
        Assert(Norm(y1-y2) < eps*Norm(y1),"SymBand cv%b");
        y1 = cv%csi;
        y2 = cv%cm;
        if (showacc) 
            std::cout<<"cv%cs: Norm(y1-y2) = "<<Norm(y1-y2)<<"  "<<
                ceps*Norm(y1)<<std::endl;
        Assert(Norm(y1-y2) < ceps*Norm(y1),"SymBand cv%cs");
    }

    if (pdc != Sing) {
        TestSymBandDiv_A<T>(dt,pdc);
        TestSymBandDiv_B1<T>(dt,pdc);
        TestSymBandDiv_C1<T>(dt,pdc);
        TestSymBandDiv_D1<T>(dt,pdc);
        TestSymBandDiv_E1<T>(dt,pdc);
        TestSymBandDiv_F1<T>(dt,pdc);
    }
    if (pdc == PosDef) {
        if (dt != tmv::CH) TestSymBandDiv_B2<T>(dt,pdc);
        if (dt == tmv::LU) TestSymBandDiv_C2<T>(dt,pdc);
        if (dt == tmv::LU) TestSymBandDiv_D2<T>(dt,pdc);
        if (dt != tmv::CH) TestSymBandDiv_E2<T>(dt,pdc);
        TestSymBandDiv_F2<T>(dt,pdc);
    }

    std::cout<<PDLabel(pdc)<<" SymBandMatrix<"<<tmv::TMV_Text(T())<<
        "> Division using ";
    std::cout<<tmv::TMV_Text(dt)<<" passed all tests\n";
}

template <class T> 
void TestAllSymBandDiv()
{
    TestHermBandDecomp<T,tmv::Upper,tmv::ColMajor>();
    TestHermBandDecomp<T,tmv::Upper,tmv::RowMajor>();
    TestHermBandDecomp<T,tmv::Lower,tmv::ColMajor>();
    TestHermBandDecomp<T,tmv::Lower,tmv::RowMajor>();
    TestSymBandDecomp<T,tmv::Upper,tmv::ColMajor>();
    TestSymBandDecomp<T,tmv::Upper,tmv::RowMajor>();
    TestSymBandDecomp<T,tmv::Lower,tmv::ColMajor>();
    TestSymBandDecomp<T,tmv::Lower,tmv::RowMajor>();
    std::cout<<"SymBandMatrix<"<<tmv::TMV_Text(T())<<"> passed all ";
    std::cout<<"decomposition tests.\n";
    TestSymBandDiv<T>(tmv::CH,PosDef);
    TestSymBandDiv<T>(tmv::LU,PosDef);
    TestSymBandDiv<T>(tmv::LU,InDef);
    TestSymBandDiv<T>(tmv::SV,PosDef);
    TestSymBandDiv<T>(tmv::SV,InDef);
    TestSymBandDiv<T>(tmv::SV,Sing);
}

#ifdef TEST_DOUBLE
template void TestAllSymBandDiv<double>();
#endif
#ifdef TEST_FLOAT
template void TestAllSymBandDiv<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestAllSymBandDiv<long double>();
#endif
