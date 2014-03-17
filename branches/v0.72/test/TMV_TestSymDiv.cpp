
#define START 0

#include "TMV.h"
#include "TMV_Sym.h"
#include "TMV_Test.h"
#include "TMV_Test_2.h"
#include "TMV_TestSymArith.h"

#ifdef TMV_MEM_DEBUG
// See the discussion of this in TMV_TestTri.cpp.  But basically, there seems to be something
// in the std library exception class that doesn't interact well with the mmgr-style memory
// debugging.  So we skip those tests if we are doing MEM_DEBUG
#define NOTHROW
#endif

template <class T> 
static bool IsPosDef(const tmv::GenSymMatrix<T>& m)
{
#ifdef NOTHROW
    for(int i=1;i<=m.size();i++) {
        T d = Det(m.subSymMatrix(0,i));
        if (tmv::TMV_REAL(d) <= 0) return false;
    }
    return true;
#else
    try {
        tmv::HermMatrix<T> m2 = m;
        CH_Decompose(m2.view());
    }
    catch (tmv::NonPosDef) {
#if (XTEST & 16)
        if (showacc) {
            std::cout<<"caught nonposdef\n";
            std::cout<<"m.size ="<<m.size()<<std::endl;
        }
        for(int i=1;i<=m.size();i++) {
            if (showacc) {
                std::cout<<"i = "<<i<<std::endl;
            }
            T d = Det(tmv::Matrix<T>(m.subSymMatrix(0,i)));
            if (showacc) std::cout<<"Det(0.."<<i<<") = "<<d<<std::endl;
            if (!(tmv::TMV_REAL(d) > 0)) {
                if (showacc) {
                    std::cout<<"!>0 "<<tmv::TMV_REAL(d)<<std::endl;
                }
                return false;
            } else {
                if (showacc) {
                    std::cout<<">0 "<<tmv::TMV_REAL(d)<<std::endl;
                }
            }
        }
        std::cout<<"m = "<<TMV_Text(m)<<"  "<<m<<std::endl;
        for(int i=1;i<=m.size();i++) {
            T d = Det(tmv::Matrix<T>(m.subSymMatrix(0,i)));
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
        T d = Det(tmv::Matrix<T>(m.subSymMatrix(0,i)));
        if (showacc) 
            std::cout<<"Det(0.."<<i<<") = "<<d<<std::endl;
        if (tmv::TMV_REAL(d) < 0) {
            std::cout<<"m = "<<TMV_Text(m)<<"  "<<m<<std::endl;
            std::cout<<"Det(0.."<<i<<") = "<<d<<std::endl;
            Assert(tmv::TMV_FALSE,
                   "Didn't catch NonPosDef, but determinant of sub-block is "
                   "negative");
        }
    }
#endif
    return true;
#endif
}

template <class T> 
void TestSymDiv(tmv::DivType dt, PosDefCode pdc)
{
    const int N = 10;

    std::vector<tmv::SymMatrixView<T> > s;
    std::vector<tmv::SymMatrixView<std::complex<T> > > cs;
    MakeSymList(s,cs,pdc);

    tmv::Vector<T> v1(N);
    tmv::Vector<T> v2(N-1);
    for (int i=0; i<N; ++i) v1(i) = T(16-3*i); 
    for (int i=0; i<N-1; ++i) v2(i) = T(-7+2*i); 

    std::ostream* checkout = showdiv ? &std::cout : 0;

    for(size_t i=START;i<s.size();i++) {
        tmv::SymMatrixView<T> si = s[i];
        tmv::SymMatrixView<std::complex<T> > csi = cs[i];
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
            tmv::SymMatrix<T> six = si;
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
        Assert(Norm(x1-x2) < eps*Norm(x1),"Sym v/b");

        x1 = v1%si;
        x2 = v1%m;
        if (showacc)
            std::cout<<"v%b: Norm(x1-x2) = "<<Norm(x1-x2)<<
                "  "<<eps*Norm(x1)<<std::endl;
        Assert(Norm(x1-x2) < eps*Norm(x1),"Sym v%b");

        tmv::Matrix<T,tmv::ColMajor> sinv = si.inverse();
        tmv::Matrix<T,tmv::ColMajor> minv = m.inverse();
        if (showacc) {
            std::cout<<"sinv = "<<sinv<<std::endl;
            std::cout<<"minv = "<<minv<<std::endl;
            std::cout<<"Norm(minv-sinv) = "<<Norm(minv-sinv)<<
                "  "<<eps*Norm(sinv)<<std::endl;
        }
        Assert(Norm(sinv-minv) < eps*Norm(sinv),"Sym Inverse");

        if (pdc != Sing) {
            tmv::Matrix<T> m2 = m;
            if (showacc) {
                std::cout<<"Det(si) = "<<Det(si)<<
                    ", Det(m) = "<<Det(m2)<<std::endl;
                std::cout<<"abs(sdet-mdet) = "<<std::abs(Det(si)-Det(m2));
                std::cout<<"  EPS*abs(mdet) = "<<
                    eps*std::abs(Det(m2))<<std::endl;
                std::cout<<"abs(abs(sdet)-abs(mdet)) = "<<
                    std::abs(std::abs(Det(si))-std::abs(Det(m2)));
                std::cout<<"  EPS*abs(mdet) = "<<
                    eps*std::abs(Det(m2))<<std::endl;
                std::cout<<"LogDet(m) = "<<LogDet(m2)<<std::endl;
                std::cout<<"LogDet(si) = "<<LogDet(si)<<std::endl;
                std::cout<<"abs(diff) = "<<
                    std::abs(LogDet(m2)-LogDet(si))<<"   ";
                std::cout<<"eps = "<<eps<<std::endl;
                T msign,ssign;
                m.logDet(&msign);
                si.logDet(&ssign);
                std::cout<<"m sign = "<<msign<<std::endl;
                std::cout<<"si sign = "<<ssign<<std::endl;
                std::cout<<"abs(diff) = "<<std::abs(msign-ssign)<<"   ";
                std::cout<<"eps = "<<eps<<std::endl;
            }
            Assert(std::abs(Det(m2)-Det(si)) < eps*std::abs(Det(m2)),
                   "Sym Det");
            T msign,ssign;
            Assert(std::abs(m2.logDet(&msign)-si.logDet(&ssign)) < 10*N*eps,
                   "Sym LogDet");
            Assert(std::abs(msign-ssign) < 10*N*eps,"Sym LogDet - sign");
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
            tmv::SymMatrix<std::complex<T> > csix = csi;
            csix.divideUsing(dt);
            csix.setDiv();
            Assert(csix.checkDecomp(checkout),"CheckDecomp csix(sym)"); 
        }
        Assert(csi.checkDecomp(checkout),"CheckDecomp csi"); 

        T ceps = EPS;
        if (pdc == Sing) ceps *= 1000;
        else ceps *= Norm(cm)*Norm(cm.inverse());

        if (pdc != Sing) {
            tmv::Matrix<std::complex<T> > cm2 = cm;
            //std::cout<<"cm2 = "<<cm2<<std::endl;
            if (showacc) {
                std::cout<<"Det(csi) = "<<Det(csi)<<
                    ", Det(cm) = "<<Det(cm2)<<std::endl;
                std::cout<<"abs(csidet-cmdet) = "<<
                    std::abs(Det(csi)-Det(cm2));
                std::cout<<"  csidet/cmdet = "<<Det(csi)/Det(cm2);
                std::cout<<"  EPS*abs(cmdet) = "<<
                    ceps*std::abs(Det(cm2))<<std::endl;
                std::cout<<"abs(abs(csdet)-abs(cmdet)) = "<<
                    std::abs(std::abs(Det(csi))-std::abs(Det(cm2)));
                std::cout<<"  EPS*abs(cmdet) = "<<
                    ceps*std::abs(Det(cm2))<<std::endl;
                std::cout<<"LogDet(cm) = "<<LogDet(cm2)<<std::endl;
                std::cout<<"LogDet(csi) = "<<LogDet(csi)<<std::endl;
                std::cout<<"abs(diff) = "<<
                    std::abs(LogDet(cm2)-LogDet(csi))<<"   ";
                std::cout<<"eps = "<<eps<<std::endl;
                std::complex<T> cmsign,cssign;
                cm2.logDet(&cmsign);
                csi.logDet(&cssign);
                std::cout<<"cm sign = "<<cmsign<<std::endl;
                std::cout<<"csi sign = "<<cssign<<std::endl;
                std::cout<<"abs(diff) = "<<std::abs(cmsign-cssign)<<"   ";
                std::cout<<"eps = "<<eps<<std::endl;
            }
            Assert(std::abs(Det(csi)-Det(cm2)) < 
                   ceps*std::abs(Det(cm2)),"Sym CDet");
            std::complex<T> cmsign,cssign;
            Assert(std::abs(cm2.logDet(&cmsign)-csi.logDet(&cssign)) < 10*N*eps,
                   "Sym CLogDet");
            Assert(std::abs(cmsign-cssign) < 10*N*eps,"Sym CLogDet - sign");
        }

        tmv::Vector<std::complex<T> > cv = v1 * std::complex<T>(1,1);
        cv(1) += std::complex<T>(-1,5);
        cv(2) -= std::complex<T>(-1,5);

        // test real / complex
        tmv::Vector<std::complex<T> > y1 = v1/csi;
        tmv::Vector<std::complex<T> > y2 = v1/cm;
        if (showacc) 
            std::cout<<"v/cs: Norm(y1-y2) = "<<Norm(y1-y2)<<
                "  "<<ceps*Norm(y1)<<std::endl;
        Assert(Norm(y1-y2) < ceps*Norm(y1),"Sym v/cs");

        // test complex / real
        y1 = cv/si;
        y2 = cv/m;
        if (showacc) 
            std::cout<<"cv/b: Norm(y1-y2) = "<<Norm(y1-y2)<<
                "  "<<eps*Norm(y1)<<std::endl;
        Assert(Norm(y1-y2) < eps*Norm(y1),"Sym cv/b");

        // test complex / complex
        y1 = cv/csi;
        y2 = cv/cm;
        if (showacc) 
            std::cout<<"cv/cs: Norm(y1-y2) = "<<Norm(y1-y2)<<
                "  "<<ceps*Norm(y1)<<std::endl;
        Assert(Norm(y1-y2) < ceps*Norm(y1),"Sym cv/cs");

        y1 = v1%csi;
        y2 = v1%cm;
        if (showacc) 
            std::cout<<"v%cs: Norm(y1-y2) = "<<Norm(y1-y2)<<
                "  "<<ceps*Norm(y1)<<std::endl;
        Assert(Norm(y1-y2) < ceps*Norm(y1),"Sym v%cs");

        y1 = cv%si;
        y2 = cv%m;
        if (showacc) 
            std::cout<<"cv%b: Norm(y1-y2) = "<<Norm(y1-y2)<<
                "  "<<eps*Norm(y1)<<std::endl;
        Assert(Norm(y1-y2) < eps*Norm(y1),"Sym cv%b");
        y1 = cv%csi;
        y2 = cv%cm;
        if (showacc) 
            std::cout<<"cv%cs: Norm(y1-y2) = "<<Norm(y1-y2)<<
                "  "<<ceps*Norm(y1)<<std::endl;
        Assert(Norm(y1-y2) < ceps*Norm(y1),"Sym cv%cs");
    }

    if (pdc != Sing) {
        TestSymDiv_A<T>(dt,pdc);
        TestSymDiv_B1<T>(dt,pdc);
        TestSymDiv_C1<T>(dt,pdc);
        TestSymDiv_D1<T>(dt,pdc);
        TestSymDiv_E1<T>(dt,pdc);
    }
    if (pdc == PosDef) {
        if (dt != tmv::CH) TestSymDiv_B2<T>(dt,pdc);
        if (dt == tmv::LU) TestSymDiv_C2<T>(dt,pdc);
        if (dt == tmv::LU) TestSymDiv_D2<T>(dt,pdc);
        if (dt != tmv::CH) TestSymDiv_E2<T>(dt,pdc);
    }

    std::cout<<PDLabel(pdc)<<" SymMatrix<"<<tmv::TMV_Text(T())<<
        "> Division using ";
    std::cout<<tmv::TMV_Text(dt)<<" passed all tests\n";
}

template <class T> 
void TestAllSymDiv()
{
    TestHermDecomp<T,tmv::Lower,tmv::ColMajor>();
    TestHermDecomp<T,tmv::Lower,tmv::RowMajor>();
    TestHermDecomp<T,tmv::Upper,tmv::ColMajor>();
    TestHermDecomp<T,tmv::Upper,tmv::RowMajor>();
    TestSymDecomp<T,tmv::Upper,tmv::ColMajor>();
    TestSymDecomp<T,tmv::Upper,tmv::RowMajor>();
    TestSymDecomp<T,tmv::Lower,tmv::ColMajor>();
    TestSymDecomp<T,tmv::Lower,tmv::RowMajor>();
    TestPolar<T,tmv::RowMajor>();
    TestPolar<T,tmv::ColMajor>();
    std::cout<<"SymMatrix<"<<tmv::TMV_Text(T())<<"> passed all ";
    std::cout<<"decomposition tests.\n";
    TestSymDiv<T>(tmv::CH,PosDef);
    TestSymDiv<T>(tmv::LU,PosDef);
    TestSymDiv<T>(tmv::LU,InDef);
    TestSymDiv<T>(tmv::SV,PosDef);
    TestSymDiv<T>(tmv::SV,InDef);
    TestSymDiv<T>(tmv::SV,Sing);
}

#ifdef TEST_DOUBLE
template void TestAllSymDiv<double>();
#endif
#ifdef TEST_FLOAT
template void TestAllSymDiv<float>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestAllSymDiv<long double>();
#endif
