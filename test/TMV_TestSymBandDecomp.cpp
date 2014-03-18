
#define START 0

#include "TMV.h"
#include "TMV_SymBand.h"
#include "TMV_Test.h"
#include "TMV_Test_2.h"

#ifdef TMV_MEM_DEBUG
// See the discussion of this in TMV_TestTri.cpp.  But basically, there seems to be something
// in the std library exception class that doesn't interact well with the mmgr-style memory
// debugging.  So we skip those tests if we are doing MEM_DEBUG
#define NOTHROW
#endif

template <class T, tmv::UpLoType uplo, tmv::StorageType stor> 
void TestHermBandDecomp()
{
    typedef std::complex<T> CT;
    for (int mattype = START; mattype <= 17; mattype++) {
#if !(XTEST & 64) || defined(LAP)
        if (mattype % 6 >= 3) continue;
#endif
        if (showstartdone) {
            std::cout<<"HermBand: mattype = "<<mattype<<std::endl;
            std::cout<<"uplo, stor = "<<TMV_Text(uplo)<<
                "  "<<TMV_Text(stor)<<std::endl;
        }
        // mattype = 0  is PosDef, nlo = 8
        // mattype = 1  is Indef, nlo = 8
        // mattype = 2  is Singular, nlo = 8
        // mattype = 3  is Singular, nlo = 8, seriously bad defects
        // mattype = 4  is Singular, nlo = 8, nearly zero
        // mattype = 5  is Singular, nlo = 8, nearly overflow
        // mattype = 6  is PosDef, nlo = 1
        // mattype = 7  is Indef, nlo = 1
        // mattype = 8  is Singular, nlo = 1
        // mattype = 9  is Singular, nlo = 1, seriously bad defects
        // mattype = 10  is Singular, nlo = 1, nearly zero
        // mattype = 11  is Singular, nlo = 1, nearly overflow
        // mattype = 12  is PosDef, nlo = 0
        // mattype = 13  is Indef, nlo = 0
        // mattype = 14  is Singular, nlo = 0
        // mattype = 15  is Singular, nlo = 0, seriously bad defects
        // mattype = 16  is Singular, nlo = 0, nearly zero
        // mattype = 17  is Singular, nlo = 0, nearly overflow

        const int N = 200;
        const bool posdef = mattype % 6 == 0;
        const bool singular = mattype % 6 >= 2;
        const bool baddefect = mattype % 6 == 3;
        const bool nearunderflow = mattype % 6 == 4;
        const bool nearoverflow = mattype % 6 == 5;
        const int nlo = mattype / 6 == 0 ? 8 : mattype / 6 == 1 ? 1 : 0;
        const T Teps = std::numeric_limits<T>::epsilon();
        const T Tmin = std::numeric_limits<T>::min();
        const T Tmax = std::numeric_limits<T>::max();
        if (showstartdone) {
            std::cout<<"posdef = "<<posdef<<", singular = "<<singular<<
                ", nlo = "<<nlo<<std::endl;
        }

        tmv::HermBandMatrix<T,uplo|stor> m(N,nlo);
        //std::cout<<"m = "<<TMV_Text(m)<<std::endl;
        for(int i=0;i<N;++i) for(int j=0;j<N;++j) 
            if ((uplo == tmv::Upper && j>=i && j<=i+nlo) || 
                (uplo == tmv::Lower && i>=j && i<=j+nlo)) 
                m(i,j) = T(2.5+4*i-5*j);
        m /= T(10*N);
        if (posdef) {
            // pos def
            m += T(10);
            m.diag(0,N/4,N/2) *= T(10);
            m.subSymBandMatrix(N/2,3*N/4,nlo/2,2) *= T(100);
        } else if (singular) {
            // singular
            if (nlo > 0) {
                if (nlo > 2) m.row(2,0,2).setZero();
                else m.row(2,2-nlo,2).setZero();
                m.row(2,2,2+nlo+1).setZero();
                m.row(14,15-nlo,14) = m.row(15,15-nlo,14);
                m.row(14,15,14+nlo+1) = m.row(15,15,14+nlo+1);
                if (14-nlo >= 0) m(14,14-nlo) = T(0);
                if (15+nlo < N) m(15,15+nlo) = T(0);
                m(14,15) = m(14,14) = m(15,15);
            }
            m.subSymBandMatrix(N/2,3*N/4).setZero();
        } else {
            // indef
            if (nlo > 0) {
                m(2,1) = T(5);
                m.diag(1,4,7).addToAll(T(1));
                m.row(17,17-nlo,17).addToAll(T(20));
                m.diag(0,N/4,N/2) *= T(10);
                m.subSymBandMatrix(N/2,3*N/4,nlo/2,2) /= T(100);
            }
        }
        //std::cout<<"m = "<<m<<std::endl;

        tmv::HermBandMatrix<CT,uplo|stor> c(N,nlo);
        //std::cout<<"c = "<<TMV_Text(c)<<std::endl;
        for(int i=0;i<N;++i) for(int j=0;j<N;++j) 
            if ((uplo == tmv::Upper && j>=i && j<=i+nlo) || 
                (uplo == tmv::Lower && i>=j && i<=j+nlo)) 
                c(i,j) = CT(2.5+4*i-5*j,3-i);
        c.diag().imagPart().setZero();
        c /= T(10*N);
        if (posdef) {
            c += T(10);
            c.diag(0,N/4,N/2) *= T(10);
            c.subSymBandMatrix(N/2,3*N/4,nlo/2,2) *= T(100);
        } else if (singular) {
            if (nlo > 0) {
                if (nlo > 2) c.row(2,0,2).setZero();
                else c.row(2,2-nlo,2).setZero();
                c.row(2,2,2+nlo+1).setZero();
                c.row(14,15-nlo,14) = c.row(15,15-nlo,14);
                c.row(14,15,14+nlo+1) = c.row(15,15,14+nlo+1);
                if (14-nlo >= 0) c(14,14-nlo) = T(0);
                if (15+nlo < N) c(15,15+nlo) = T(0);
                c(14,15) = c(14,14) = c(15,15);
            }
            c.subSymBandMatrix(N/2,3*N/4).setZero();
        } else {
            // indef
            if (nlo > 0) {
                c(2,1) = T(5);
                c.diag(1,4,7).addToAll(T(1));
                c.row(17,17-nlo,17).addToAll(T(20));
                c.diag(0,N/4,N/2) *= T(10);
                c.subSymBandMatrix(N/2,3*N/4,nlo/2,2) /= T(100);
            }
        }
        //std::cout<<"c = "<<c<<std::endl;

        if (baddefect) {
            for(int i=10;i<N;i+=10) {
                m.subSymBandMatrix(i,N) *= T(-1.e-10);
                c.subSymBandMatrix(i,N) *= T(-1.e-10);
            }
        }

        if (nearunderflow) {
            T x = Tmin;
            m.setAllTo(x);
            if (nlo > 0)
                c.upperBandOff().setAllTo(CT(x,2*x));
            c.diag().setAllTo(x);
            for(int i=1;i<N;++i) {
                int start = i-nlo; if (start < 0) start = 0;
                m.col(i,start,i+1) /= T(i);
                c.col(i,start,i+1) /= T(i);
            }
            m(N/2,N/2) = x/Teps;
            c(N/2,N/2) = x/Teps;
            m.diag(0,N/5,2*N/5).setZero();
            c.diag(0,N/5,2*N/5).setZero();
            m.diag(0,N/2+1,N) *= T(-1);
            c.diag(0,N/2+1,N) *= T(-1);
            m.subSymBandMatrix(N-4,N).setAllTo(-x/Teps/T(8));
            c.subSymBandMatrix(N-4,N).setAllTo(-x/Teps/T(8));
        }

        if (nearoverflow) {
            T x = Tmax;
            x /= N;
            m.setAllTo(x);
            if (nlo > 0) 
                c.upperBandOff().setAllTo(CT(x,x/2));
            c.diag().setAllTo(x);
            for(int i=1;i<N;++i) {
                int start = i-nlo;  if (start < 0) start = 0;
                m.col(i,start,i+1) /= T(i);
                c.col(i,start,i+1) /= T(i);
            }
            m.diag(0,N/5,2*N/5).setZero();
            c.diag(0,N/5,2*N/5).setZero();
            m.diag(0,N/2,N) *= T(-1);
            c.diag(0,N/2,N) *= T(-1);
        }


        T eps = EPS;
        T ceps = EPS;
        T normm = Norm(m);
        T normc = Norm(c);
        if (showacc) std::cout<<"eps = "<<eps<<"  "<<ceps<<std::endl;
        if (!singular) {
            T kappa = normm * Norm(m.inverse());
            T ckappa = normc * Norm(c.inverse());
            if (showacc) std::cout<<"kappa = "<<kappa<<"  "<<ckappa<<std::endl;
            eps *= kappa;
            ceps *= ckappa;
        } else {
            if (showacc) std::cout<<"eps *= "<<T(10*N)<<std::endl;
            eps *= T(10*N);
            ceps *= T(10*N);
        }
        if (showacc) std::cout<<"eps => "<<eps<<"  "<<ceps<<std::endl;

        // CH Decomposition
#ifdef NOTHROW
        if (posdef) {
#else
            try {
#endif
                if (showstartdone) std::cout<<"CH"<<std::endl;
                tmv::BandMatrix<T> L(N,N,nlo,0);
                tmv::BandMatrix<T> LLt(N,N,nlo,nlo);
                if (nlo > 1) {
                    L = m.chd().getL();
                    LLt = L*L.adjoint();
                    if (showacc) {
                        std::cout<<"Norm(m-LLt) = "<<Norm(m-LLt)<<std::endl;
                        std::cout<<"cf eps*Norm(m) = "<<eps*normm<<std::endl;
                    }
                    Assert(Equal(m,LLt,eps*normm),"HermBand CH");
                }

                tmv::BandMatrix<CT> cL(N,N,nlo,0);
                tmv::BandMatrix<CT> cLLt(N,N,nlo,nlo);
                if (nlo > 1) {
                    cL = c.chd().getL();
                    cLLt = cL*cL.adjoint();
                    if (showacc) {
                        std::cout<<"Norm(c-cLcLt) = "<<Norm(c-cLLt)<<std::endl;
                        std::cout<<"cf eps*Norm(c) = "<<eps*normc<<std::endl;
                    }
                    Assert(Equal(c,cLLt,ceps*normc),"HermBand C CH");
                }

#if (XTEST & 16)
                tmv::HermBandMatrix<T,uplo|stor> m2 = m;
                CH_Decompose(m2);
                L = m2.lowerBand();
                LLt = L*L.adjoint();
                Assert(Equal(m,LLt,eps*normm),"HermBand CH2");

                tmv::HermBandMatrix<CT,uplo|stor> c2 = c;
                CH_Decompose(c2);
                cL = c2.lowerBand();
                cLLt = cL*cL.adjoint();
                Assert(Equal(c,cLLt,ceps*normc),"HermBand C CH2");

                c2.conjugate() = c;
                CH_Decompose(c2.conjugate());
                cL = c2.conjugate().lowerBand();
                cLLt = cL*cL.adjoint();
                Assert(Equal(c,cLLt,ceps*normc),"HermBand C CH2");
#endif
                std::cout<<"."; std::cout.flush();
#ifndef NOTHROW
            } catch(tmv::NonPosDef) {
                Assert(!posdef,"Caught NonPosDef, but mattype == pos def");
            }
#else
        }
#endif

        // LDL Decomposition
        if (nlo == 1) {
#ifdef NOTHROW
            if (posdef) {
#else
                try {
#endif
                    if (showstartdone) std::cout<<"LDL"<<std::endl;
                    tmv::BandMatrix<T> L(N,N,1,0);
                    tmv::DiagMatrix<T> D(N);
                    tmv::BandMatrix<T> LDL(N,N,1,1);
                    if (posdef) {
                        L = m.chd().getL();
                        D = m.chd().getD();
                        LDL = L*D*L.adjoint();
                        Assert(Equal(m,LDL,eps*normm),"HermBand LDL");
                    }

                    tmv::BandMatrix<CT> cL(N,N,nlo,0);
                    tmv::DiagMatrix<CT> cD(N);
                    tmv::BandMatrix<CT> cLDL(N,N,nlo,nlo);
                    if (posdef) {
                        cL = c.chd().getL();
                        cD = c.chd().getD();
                        cLDL = cL*cD*cL.adjoint();
                        Assert(Equal(c,cLDL,ceps*normc),"HermBand C LDL");
                    }

#if (XTEST & 16)
                    tmv::HermBandMatrix<T,uplo|stor> m2 = m;
                    LDL_Decompose(m2);
                    (L = m2.lowerBand()).diag().setAllTo(T(1));
                    D = DiagMatrixViewOf(m2.diag());
                    LDL = L*D*L.adjoint();
                    Assert(Equal(m,LDL,eps*normm),"HermBand LDL2");

                    // Alt calculation:
                    // (I + L) D (I + Lt)
                    // D + LD + DLt = LDLt
                    tmv::DiagMatrix<T> E(m2.diag(-1));
                    LDL.setZero();
                    LDL.diag(0,1,N) = (E*D.subDiagMatrix(0,N-1)*E).diag();
                    LDL.diag(-1) = (E*D.subDiagMatrix(0,N-1)).diag();
                    LDL += D;
                    LDL.diag(1) = LDL.diag(-1);
                    Assert(Equal(m,LDL,eps*normm),"HermBand LDL2 Alt");

                    tmv::HermBandMatrix<CT,uplo|stor> c2 = c;
                    LDL_Decompose(c2);
                    (cL = c2.lowerBand()).diag().setAllTo(T(1));
                    cD = DiagMatrixViewOf(c2.diag());
                    cLDL = cL*cD*cL.adjoint();
                    Assert(Equal(c,cLDL,ceps*normc),"HermBand C LDL2");
                    Assert(Norm(cD.imagPart()) <= ceps*normc,
                           "HermBand C LDL2 - D");

                    // Alt calculation:
                    tmv::DiagMatrix<CT> cE(c2.diag(-1));
                    cLDL.setZero();
                    cLDL.diag(0,1,N) =
                        (cE*cD.subDiagMatrix(0,N-1)*cE.conjugate()).diag();
                    cLDL.diag(-1) = (cE*cD.subDiagMatrix(0,N-1)).diag();
                    cLDL += cD;
                    cLDL.diag(1) = cLDL.diag(-1).conjugate();
                    Assert(Equal(c,cLDL,ceps*normc),"HermBand C LDL2 Alt");

                    c2.conjugate() = c;
                    LDL_Decompose(c2.conjugate());
                    (cL = c2.conjugate().lowerBand()).diag().setAllTo(T(1));
                    cD = DiagMatrixViewOf(c2.conjugate().diag());
                    cLDL = cL*cD*cL.adjoint();
                    Assert(Equal(c,cLDL,ceps*normc),"HermBand C LDL3");
                    Assert(Norm(cD.imagPart()) <= ceps*normc,
                           "HermBand C LDL2 - D");
#endif
                    std::cout<<"."; std::cout.flush();
#ifndef NOTHROW
                } catch (tmv::NonPosDef) { 
                    // The Lapack version throws whenever mattype is not posdef, 
                    // but native algorithm succeeds when mattype is indefinite,
                    // or even singular.  
                    Assert(!posdef,"caught NonPosDef, but mattype == posdef");
                }
#else
            }
#endif
        }

        // SV Decomposition
        {
            if (showstartdone) std::cout<<"SV"<<std::endl;
            tmv::Matrix<T> U = m.svd().getU();
            tmv::DiagMatrix<T> S = m.svd().getS();
            tmv::Matrix<T> Vt = m.svd().getVt();
            if (showacc) {
                std::cout<<"Norm(m-U*S*Vt) = "<<Norm(m-U*S*Vt)<<std::endl;
                std::cout<<"eps * normm = "<<eps*normm<<std::endl;
            }
            Assert(Equal(m,U*S*Vt,eps*normm),"HermBand SV");

            tmv::Matrix<CT> cU = c.svd().getU();
            tmv::DiagMatrix<T> cS = c.svd().getS();
            tmv::Matrix<CT> cVt = c.svd().getVt();
            if (showacc) {
                std::cout<<"Norm(c-cU*cS*cVt) = "<<Norm(c-cU*cS*cVt)<<std::endl;
                std::cout<<"ceps * normc = "<<ceps*normc<<std::endl;
                std::cout<<"Norm(cUt*cU-1) = "<<Norm(cU.adjoint()*cU-T(1))<<std::endl;
                std::cout<<"Norm(cU*cUt-1) = "<<Norm(cU*cU.adjoint()-T(1))<<std::endl;
                std::cout<<"Norm(cV*cVt-1) = "<<Norm(cVt.adjoint()*cVt-T(1))<<std::endl;
                std::cout<<"Norm(cVt*cV-1) = "<<Norm(cVt*cVt.adjoint()-T(1))<<std::endl;
            }
            Assert(Equal(c,cU*cS*cVt,ceps*normc),"HermBand C SV");

#if (XTEST & 16)
            const T x = 
                nearoverflow ? Tmax*Teps :
                nearunderflow ? Tmin/Teps : T(1);

            tmv::Matrix<T> U2(N,N);
            tmv::DiagMatrix<T> S2(N);
            tmv::Matrix<T> Vt2(N,N);
            SV_Decompose(m,U2,S2,Vt2);
            if (showacc) {
                std::cout<<"Norm(m-U*S*Vt) = "<<Norm(m-U2*S2*Vt2)<<std::endl;
                std::cout<<"eps * normm = "<<eps*normm<<std::endl;
            }
            Assert(Equal(m,U2*S2*Vt2,eps*normm),"HermBand SV2");

            SV_Decompose(m,S2);
            Assert(Equal(S2,S,eps*normm),"HermBand SV3");
            SV_Decompose(m,U2,S2);
            Assert(Equal(S2,S,eps*normm),"HermBand SV4 S");
            Assert(Equal(
                    (m/x)*tmv::Matrix<T>(m.adjoint()/x),
                    U2*tmv::DiagMatrix<T>(S2/x)*
                    tmv::DiagMatrix<T>(S2/x)*U2.adjoint(),
                    eps*(normm/x)*(normm/x)),"HermBand SV4 U");
            SV_Decompose(m,S2,Vt2);
            Assert(Equal(S2,S,eps*normm),"HermBand SV5 S");
            Assert(Equal(
                    (m.adjoint()/x)*tmv::Matrix<T>(m/x),
                    Vt2.adjoint()*tmv::DiagMatrix<T>(S2/x)*
                    tmv::DiagMatrix<T>(S2/x)*Vt2,
                    eps*(normm/x)*(normm/x)),"HermBand SV5 Vt");

            tmv::Matrix<CT> cU2(N,N);
            tmv::DiagMatrix<T> cS2(N);
            tmv::Matrix<CT> cVt2(N,N);
            SV_Decompose(c,cU2,cS2,cVt2);
            Assert(Equal(c,cU2*cS2*cVt2,ceps*normc),"HermBand C SV2");

            SV_Decompose(c,cS2);
            Assert(Equal(cS2,cS,ceps*normc),"HermBand C SV3");
            SV_Decompose(c,cU2,cS2);
            Assert(Equal(cS2,cS,ceps*normc),"HermBand C SV4 S");
            Assert(Equal(
                    tmv::Matrix<CT>(c/x)*tmv::Matrix<CT>(c.adjoint()/x),
                    cU2*tmv::DiagMatrix<CT>(cS2/x)*
                    tmv::DiagMatrix<CT>(cS2/x)*cU2.adjoint(),
                    ceps*(normc/x)*(normc/x)),"HermBand C SV4 U");
            SV_Decompose(c,cS2,cVt2);
            Assert(Equal(cS2,cS,ceps*normc),"HermBand C SV5 S");
            Assert(Equal(
                    tmv::Matrix<CT>(c.adjoint()/x)*tmv::Matrix<CT>(c/x),
                    cVt2.adjoint()*tmv::DiagMatrix<CT>(cS2/x)*
                    tmv::DiagMatrix<CT>(cS2/x)*cVt2,
                    ceps*(normc/x)*(normc/x)),"Herm C SV5 Vt");

            SV_Decompose(c.conjugate(),cS2);
            Assert(Equal(cS2,cS,ceps*normc),"HermBand C SV6");
            SV_Decompose(c.conjugate(),cU2,cS2);
            Assert(Equal(cS2,cS,ceps*normc),"HermBand C SV7 S");
            Assert(Equal(
                    tmv::Matrix<CT>(c.conjugate()/x)*tmv::Matrix<CT>(c.transpose()/x),
                    cU2*tmv::DiagMatrix<CT>(cS2/x)*
                    tmv::DiagMatrix<CT>(cS2/x)*cU2.adjoint(),
                    ceps*(normc/x)*(normc/x)),"HermBand C SV7 U");
            SV_Decompose(c.conjugate(),cU2.conjugate(),cS2);
            Assert(Equal(cS2,cS,ceps*normc),"HermBand C SV8 S");
            Assert(Equal(
                    tmv::Matrix<CT>(c.conjugate()/x)*tmv::Matrix<CT>(c.transpose()/x),
                    cU2.conjugate()*tmv::DiagMatrix<CT>(cS2/x)*
                    tmv::DiagMatrix<CT>(cS2/x)*cU2.transpose(),
                    ceps*(normc/x)*(normc/x)),"HermBand C SV8 U");
#endif
            std::cout<<"."; std::cout.flush();
        }

        // Eigen
        {
            if (showstartdone) std::cout<<"Eigen"<<std::endl;
            tmv::Matrix<T> V(N,N);
            tmv::Vector<T> L(N);
            Eigen(m,V,L);
            Assert(Equal(m*V,V*DiagMatrixViewOf(L),eps*normm),
                   "HermBand Eigen");

            tmv::Matrix<CT> cV(N,N);
            tmv::Vector<T> cL(N);
            Eigen(c,cV,cL);
            Assert(Equal(c*cV,cV*DiagMatrixViewOf(cL),ceps*normc),
                   "HermBand C Eigen");

#if (XTEST & 16)
            tmv::Vector<T> L2(N);
            Eigen(m,L2);
            Assert(Equal(L2,L,eps*normm),"HermBand Eigen2");

            tmv::Vector<T> cL2(N);
            Eigen(c,cL2);
            Assert(Equal(cL2,cL,ceps*normc),"HermBand C Eigen2");

            Eigen(c,cV.conjugate(),cL);
            Assert(Equal(
                    c*cV.conjugate(),cV.conjugate()*DiagMatrixViewOf(cL), 
                    ceps*normc),"HermBand C Eigen3");

            Eigen(c.conjugate(),cL2);
            Assert(Equal(cL2,cL,ceps*normc),"HermBand C Eigen4");
#endif
            std::cout<<"."; std::cout.flush();
        }

        // Square Root
#ifdef NOTHROW
        if (posdef) {
#else
            try {
#endif
                if (showstartdone) std::cout<<"Square Root"<<std::endl;
                tmv::HermMatrix<T,uplo|stor> S(N);
                SquareRoot(m,S);
                Assert(Equal(m,S*S,eps*normm),"HermBand Square Root");

                tmv::HermMatrix<CT,uplo|stor> cS(N);
                SquareRoot(c,cS);
                Assert(Equal(c,cS*cS,eps*normc),"HermBand C Square Root");

#if (XTEST & 16)
                SquareRoot(c,cS.conjugate());
                Assert(Equal(c,cS.conjugate()*cS.conjugate(),ceps*normc),
                       "HermBand C Square Root2");

                SquareRoot(c,cS.transpose());
                Assert(Equal(c,cS.transpose()*cS.transpose(),ceps*normc),
                       "HermBand C Square Root3");

                SquareRoot(c,cS.adjoint());
                Assert(Equal(c,cS.adjoint()*cS.adjoint(),ceps*normc),
                       "HermBand C Square Root4");

                SquareRoot(c.conjugate(),cS);
                Assert(Equal(c.conjugate(),cS*cS,ceps*normc),
                       "HermBand C Square Root5");

                SquareRoot(c.conjugate(),cS.conjugate());
                Assert(Equal(
                        c.conjugate(),cS.conjugate()*cS.conjugate(),
                        ceps*normc), "HermBand C Square Root6");

                SquareRoot(c.conjugate(),cS.transpose());
                Assert(Equal(
                        c.conjugate(),cS.transpose()*cS.transpose(),
                        ceps*normc), "HermBand C Square Root7");

                SquareRoot(c.conjugate(),cS.adjoint());
                Assert(Equal(
                        c.conjugate(),cS.adjoint()*cS.adjoint(),
                        ceps*normc), "HermBand C Square Root8");
#endif
                Assert(posdef, "didn't throw NonPosDef, and mattype != pos def"); 
                std::cout<<"."; std::cout.flush();
#ifndef NOTHROW
            } catch (tmv::NonPosDef) { 
                Assert(!posdef,"caught NonPosDef, but mattype == posdef"); 
            }
#else
        }
#endif

    }
}

template <class T, tmv::UpLoType uplo, tmv::StorageType stor> 
void TestSymBandDecomp()
{
    typedef std::complex<T> CT;
    for (int mattype = START; mattype <= 17; mattype++) {
#if !(XTEST & 64) || defined(LAP)
        if (mattype % 6 >= 3) continue;
#endif
        if (showstartdone) {
            std::cout<<"SymBand: mattype = "<<mattype<<std::endl;
            std::cout<<"uplo, stor = "<<TMV_Text(uplo)<<"  "<<
                TMV_Text(stor)<<std::endl;
        }
        // mattype = 0  is PosDef, nlo = 8
        // mattype = 1  is Indef, nlo = 8
        // mattype = 2  is Singular, nlo = 8
        // mattype = 3  is Singular, nlo = 8, seriously bad defects
        // mattype = 4  is Singular, nlo = 8, nearly zero
        // mattype = 5  is Singular, nlo = 8, nearly overflow
        // mattype = 6  is PosDef, nlo = 1
        // mattype = 7  is Indef, nlo = 1
        // mattype = 8  is Singular, nlo = 1
        // mattype = 9  is Singular, nlo = 1, seriously bad defects
        // mattype = 10  is Singular, nlo = 1, nearly zero
        // mattype = 11  is Singular, nlo = 1, nearly overflow
        // mattype = 12  is PosDef, nlo = 0
        // mattype = 13  is Indef, nlo = 0
        // mattype = 14  is Singular, nlo = 0
        // mattype = 15  is Singular, nlo = 0, seriously bad defects
        // mattype = 16  is Singular, nlo = 0, nearly zero
        // mattype = 17  is Singular, nlo = 0, nearly overflow

        const int N = 200;
        const bool posdef = mattype % 6 == 0;
        const bool singular = mattype % 6 >= 2;
        const bool baddefect = mattype % 6 == 3;
        const bool nearunderflow = mattype % 6 == 4;
        const bool nearoverflow = mattype % 6 == 5;
        const int nlo = mattype / 6 == 0 ? 8 : mattype / 6 == 1 ? 1 : 0;
        const T Teps = std::numeric_limits<T>::epsilon();
        const T Tmin = std::numeric_limits<T>::min();
        const T Tmax = std::numeric_limits<T>::max();
        if (showstartdone) {
            std::cout<<"posdef = "<<posdef<<", singular = "<<singular<<
                ", nlo = "<<nlo<<std::endl;
        }

        tmv::SymBandMatrix<T,uplo|stor> m(N,nlo);
        for(int i=0;i<N;++i) for(int j=0;j<N;++j) 
            if ((uplo == tmv::Upper && j>=i && j<=i+nlo) || 
                (uplo == tmv::Lower && i>=j && i<=j+nlo)) 
                m(i,j) = T(2.5+4*i-5*j);
        m /= T(10*N);
        if (posdef) {
            m += T(10);
            m.diag(0,N/4,N/2) *= T(10);
            m.subSymBandMatrix(N/2,3*N/4,nlo/2,2) *= T(100);
        } else if (singular) {
            if (nlo > 0) {
                if (nlo > 2) m.row(2,0,2).setZero();
                else m.row(2,2-nlo,2).setZero();
                m.row(2,2,2+nlo+1).setZero();
                m.row(14,15-nlo,14) = m.row(15,15-nlo,14);
                m.row(14,15,14+nlo+1) = m.row(15,15,14+nlo+1);
                if (14-nlo >= 0) m(14,14-nlo) = T(0);
                if (15+nlo < N) m(15,15+nlo) = T(0);
                m(14,15) = m(14,14) = m(15,15);
            }
            m.subSymBandMatrix(N/2,3*N/4).setZero();
        } else {
            // indef
            if (nlo > 0) {
                m(2,1) = T(5);
                m.diag(1,4,7).addToAll(T(1));
                m.row(17,17-nlo,17).addToAll(T(20));
                m.diag(0,N/4,N/2) *= T(10);
                m.subSymBandMatrix(N/2,3*N/4,nlo/2,2) /= T(100);
            }
        }
        //std::cout<<"m = "<<m<<std::endl;

        tmv::SymBandMatrix<CT,uplo|stor> c(N,nlo);
        for(int i=0;i<N;++i) for(int j=0;j<N;++j) 
            if ((uplo == tmv::Upper && j>=i && j<=i+nlo) || 
                (uplo == tmv::Lower && i>=j && i<=j+nlo)) 
                c(i,j) = CT(2+4*i-5*j,3-i);
        c /= T(10*N);
        if (posdef) {
            c += T(10);
            c.diag(0,N/4,N/2) *= T(10);
            c.subSymBandMatrix(N/2,3*N/4,nlo/2,2) *= T(100);
        } else if (singular) {
            if (nlo > 0) {
                if (nlo > 2) c.row(2,0,2).setZero();
                else c.row(2,2-nlo,2).setZero();
                c.row(2,2,2+nlo+1).setZero();
                c.row(14,15-nlo,14) = c.row(15,15-nlo,14);
                c.row(14,15,14+nlo+1) = c.row(15,15,14+nlo+1);
                if (14-nlo >= 0) c(14,14-nlo) = T(0);
                if (15+nlo < N) c(15,15+nlo) = T(0);
                c(14,15) = c(14,14) = c(15,15);
            }
            c.subSymBandMatrix(N/2,3*N/4).setZero();
        } else {
            // indef
            if (nlo > 0) {
                c(2,1) = T(5);
                c.diag(1,4,7).addToAll(T(1));
                c.row(17,17-nlo,17).addToAll(T(20));
                c.diag(0,N/4,N/2) *= T(10);
                c.subSymBandMatrix(N/2,3*N/4,nlo/2,2) /= T(100);
            }
        }
        //std::cout<<"c = "<<c<<std::endl;
        
        if (baddefect) {
            for(int i=10;i<N;i+=10) {
                m.subSymBandMatrix(i,N) *= T(-1.e-10);
                c.subSymBandMatrix(i,N) *= T(-1.e-10);
            }
        }

        if (nearunderflow) {
            T x = Tmin;
            m.setAllTo(x);
            c.upperBand().setAllTo(CT(x,2*x));
            for(int i=1;i<N;++i) {
                int start = i-nlo; if (start < 0) start = 0;
                m.col(i,start,i+1) /= T(i);
                c.col(i,start,i+1) /= T(i);
            }
            m(N/2,N/2) = x/Teps;
            c(N/2,N/2) = x/Teps;
            m.diag(0,N/5,2*N/5).setZero();
            c.diag(0,N/5,2*N/5).setZero();
            m.diag(0,N/2+1,N) *= T(-1);
            c.diag(0,N/2+1,N) *= T(-1);
            m.subSymBandMatrix(N-4,N).setAllTo(-x/Teps/T(8));
            c.subSymBandMatrix(N-4,N).setAllTo(-x/Teps/T(8));
        }

        if (nearoverflow) {
            T x = Tmax;
            x /= N;
            m.setAllTo(x);
            c.upperBand().setAllTo(CT(x,x/2));
            for(int i=1;i<N;++i) {
                int start = i-nlo;  if (start < 0) start = 0;
                m.col(i,start,i+1) /= T(i);
                c.col(i,start,i+1) /= T(i);
            }
            m.diag(0,N/5,2*N/5).setZero();
            c.diag(0,N/5,2*N/5).setZero();
            m.diag(0,N/2,N) *= T(-1);
            c.diag(0,N/2,N) *= T(-1);
        }

        T eps = EPS;
        T ceps = EPS;
        T normm = Norm(m);
        T normc = Norm(c);
        if (showacc) std::cout<<"eps = "<<eps<<"  "<<ceps<<std::endl;
        if (!singular) {
            T kappa = normm * Norm(m.inverse());
            T ckappa = normc * Norm(c.inverse());
            if (showacc) std::cout<<"kappa = "<<kappa<<"  "<<ckappa<<std::endl;
            eps *= kappa;
            ceps *= ckappa;
        } else {
            if (showacc) std::cout<<"eps *= "<<T(10*N)<<std::endl;
            eps *= T(10*N);
            ceps *= T(10*N);
        }
        if (showacc) std::cout<<"eps => "<<eps<<"  "<<ceps<<std::endl;

        // LDL Decomposition
        if (nlo == 1) {
#ifdef NOTHROW
            if (posdef) {
#else
                try {
#endif
                    if (showstartdone) std::cout<<"LDL"<<std::endl;
                    tmv::BandMatrix<T> L(N,N,1,0);
                    tmv::DiagMatrix<T> D(N);
                    tmv::BandMatrix<T> LDL(N,N,1,1);

                    tmv::SymBandMatrix<T,uplo|stor> m2 = m;
                    LDL_Decompose(m2);
                    (L = m2.lowerBand()).diag().setAllTo(T(1));
                    D = DiagMatrixViewOf(m2.diag());
                    LDL = L*D*L.transpose();
                    Assert(Equal(m,LDL,eps*normm),"SymBand LDL2");

#if (XTEST & 16)
                    // Alt calculation:
                    // (I + L) D (I + Lt)
                    // D + LD + DLt = LDLt
                    tmv::DiagMatrix<T> E(m2.diag(-1));
                    LDL.setZero();
                    LDL.diag(0,1,N) = (E*D.subDiagMatrix(0,N-1)*E).diag();
                    LDL.diag(-1) = (E*D.subDiagMatrix(0,N-1)).diag();
                    LDL += D;
                    LDL.diag(1) = LDL.diag(-1);
                    Assert(Equal(m,LDL,eps*normm),"SymBand LDL2 Alt");
#endif

                    tmv::BandMatrix<CT> cL(N,N,nlo,0);
                    tmv::DiagMatrix<CT> cD(N);
                    tmv::BandMatrix<CT> cLDL(N,N,nlo,nlo);

                    tmv::SymBandMatrix<CT,uplo|stor> c2 = c;
                    LDL_Decompose(c2);
                    (cL = c2.lowerBand()).diag().setAllTo(T(1));
                    cD = DiagMatrixViewOf(c2.diag());
                    cLDL = cL*cD*cL.transpose();
                    Assert(Equal(c,cLDL,ceps*normc),"SymBand C LDL2");

#if (XTEST & 16)
                    // Alt calculation:
                    tmv::DiagMatrix<CT> cE(c2.diag(-1));
                    cLDL.setZero();
                    cLDL.diag(0,1,N) = (cE*cD.subDiagMatrix(0,N-1)*cE).diag();
                    cLDL.diag(-1) = (cE*cD.subDiagMatrix(0,N-1)).diag();
                    cLDL += cD;
                    cLDL.diag(1) = cLDL.diag(-1);
                    Assert(Equal(c,cLDL,eps*normm),"SymBand C LDL2 Alt");

                    c2.conjugate() = c;
                    LDL_Decompose(c2.conjugate());
                    (cL = c2.conjugate().lowerBand()).diag().setAllTo(T(1));
                    cD = DiagMatrixViewOf(c2.conjugate().diag());
                    cLDL = cL*cD*cL.transpose();
                    Assert(Equal(c,cLDL,ceps*normc),"SymBand C LDL3");
#endif
                    std::cout<<"."; std::cout.flush();
#ifndef NOTHROW
                } catch (tmv::NonPosDef) {
                    Assert(!posdef,"caught NonPosDef, but mattype == posdef"); 
                }
#else
            }
#endif
        }

        // SV Decomposition
        {
            if (showstartdone) std::cout<<"SV"<<std::endl;
            tmv::Matrix<T> U = m.svd().getU();
            tmv::DiagMatrix<T> S = m.svd().getS();
            tmv::Matrix<T> Vt = m.svd().getVt();
            Assert(Equal(m,U*S*Vt,eps*normm),"SymBand SV");

            tmv::Matrix<CT> cU = c.symsvd().getU();
            tmv::DiagMatrix<T> cS = c.symsvd().getS();
            tmv::Matrix<CT> cVt = c.symsvd().getVt();
            Assert(Equal(c,cU*cS*cVt,ceps*normc),"SymBand C SV");

#if (XTEST & 16)
            const T x = 
                nearoverflow ? Tmax*Teps :
                nearunderflow ? Tmin/Teps : T(1);

            tmv::Matrix<T> U2(N,N);
            tmv::DiagMatrix<T> S2(N);
            tmv::Matrix<T> Vt2(N,N);
            SV_Decompose(m,U2,S2,Vt2);
            Assert(Equal(m,U2*S2*Vt2,eps*normm),"SymBand SV2");

            SV_Decompose(m,S2);
            Assert(Equal(S2,S,eps*normm),"SymBand SV3");
            SV_Decompose(m,U2,S2);
            Assert(Equal(S2,S,eps*normm),"SymBand SV4 S");
            Assert(Equal(
                    (m/x)*tmv::Matrix<T>(m.transpose()/x),
                    U2*tmv::DiagMatrix<T>(S2/x)*
                    tmv::DiagMatrix<T>(S2/x)*U2.transpose(),
                    eps*(normm/x)*(normm/x)),"SymBand SV4 U");
            SV_Decompose(m,S2,Vt2);
            Assert(Equal(S2,S,eps*normm),"SymBand SV5 S");
            Assert(Equal(
                    (m.adjoint()/x)*tmv::Matrix<T>(m/x),
                    Vt2.adjoint()*(S2/x)*tmv::DiagMatrix<T>(S2/x)*Vt2,
                    eps*(normm/x)*(normm/x)),"SymBand SV5 Vt");

            tmv::Matrix<CT> cU2(N,N);
            tmv::DiagMatrix<T> cS2(N);
            tmv::Matrix<CT> cVt2(N,N);
            SV_Decompose(c,cU2,cS2,cVt2);
            Assert(Equal(c,cU2*cS2*cVt2,ceps*normc),"SymBand C SV2");

            SV_Decompose(c,cS2);
            Assert(Equal(cS2,cS,ceps*normc),"SymBand C SV3");
            SV_Decompose(c,cU2,cS2);
            Assert(Equal(cS2,cS,ceps*normc),"SymBand C SV4 S");
            Assert(Equal(
                    tmv::Matrix<CT>(c/x)*tmv::Matrix<CT>(c.adjoint()/x),
                    cU2*tmv::DiagMatrix<CT>(cS2/x)*
                    tmv::DiagMatrix<CT>(cS2/x)*cU2.adjoint(),
                    ceps*(normc/x)*(normc/x)),"SymBand C SV4 U");
            SV_Decompose(c,cS2,cVt2);
            Assert(Equal(cS2,cS,ceps*normc),"SymBand C SV5 S");
            Assert(Equal(
                    tmv::Matrix<CT>(c.adjoint()/x)*tmv::Matrix<CT>(c/x),
                    cVt2.adjoint()*tmv::DiagMatrix<CT>(cS2/x)*
                    tmv::DiagMatrix<CT>(cS2/x)*cVt2,
                    ceps*(normc/x)*(normc/x)),"SymBand C SV5 Vt");

            SV_Decompose(c.conjugate(),cS2);
            Assert(Equal(cS2,cS,ceps*normc),"SymBand C SV6");
            SV_Decompose(c.conjugate(),cU2,cS2);
            Assert(Equal(cS2,cS,ceps*normc),"SymBand C SV7 S");
            Assert(Equal(
                    tmv::Matrix<CT>(c.conjugate()/x)*tmv::Matrix<CT>(c.transpose()/x),
                    cU2*tmv::DiagMatrix<CT>(cS2/x)*
                    tmv::DiagMatrix<CT>(cS2/x)*cU2.adjoint(),
                    ceps*(normc/x)*(normc/x)),"SymBand C SV7 U");
            SV_Decompose(c.conjugate(),cU2.conjugate(),cS2);
            Assert(Equal(cS2,cS,ceps*normc),"SymBand C SV8 S");
            Assert(Equal(
                    tmv::Matrix<CT>(c.conjugate()/x)*tmv::Matrix<CT>(c.transpose()/x),
                    cU2.conjugate()*tmv::DiagMatrix<CT>(cS2/x)*
                    tmv::DiagMatrix<CT>(cS2/x)*cU2.transpose(),
                    ceps*(normc/x)*(normc/x)),"SymBand C SV8 U");
#endif
            std::cout<<"."; std::cout.flush();
        }
    }
}

#ifdef TEST_DOUBLE
template void TestHermBandDecomp<double,tmv::Upper,tmv::ColMajor>();
template void TestHermBandDecomp<double,tmv::Upper,tmv::RowMajor>();
template void TestHermBandDecomp<double,tmv::Lower,tmv::ColMajor>();
template void TestHermBandDecomp<double,tmv::Lower,tmv::RowMajor>();
template void TestSymBandDecomp<double,tmv::Upper,tmv::ColMajor>();
template void TestSymBandDecomp<double,tmv::Upper,tmv::RowMajor>();
template void TestSymBandDecomp<double,tmv::Lower,tmv::ColMajor>();
template void TestSymBandDecomp<double,tmv::Lower,tmv::RowMajor>();
#endif
#ifdef TEST_FLOAT
template void TestHermBandDecomp<float,tmv::Upper,tmv::ColMajor>();
template void TestHermBandDecomp<float,tmv::Upper,tmv::RowMajor>();
template void TestHermBandDecomp<float,tmv::Lower,tmv::ColMajor>();
template void TestHermBandDecomp<float,tmv::Lower,tmv::RowMajor>();
template void TestSymBandDecomp<float,tmv::Upper,tmv::ColMajor>();
template void TestSymBandDecomp<float,tmv::Upper,tmv::RowMajor>();
template void TestSymBandDecomp<float,tmv::Lower,tmv::ColMajor>();
template void TestSymBandDecomp<float,tmv::Lower,tmv::RowMajor>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestHermBandDecomp<long double,tmv::Upper,tmv::ColMajor>();
template void TestHermBandDecomp<long double,tmv::Upper,tmv::RowMajor>();
template void TestHermBandDecomp<long double,tmv::Lower,tmv::ColMajor>();
template void TestHermBandDecomp<long double,tmv::Lower,tmv::RowMajor>();
template void TestSymBandDecomp<long double,tmv::Upper,tmv::ColMajor>();
template void TestSymBandDecomp<long double,tmv::Upper,tmv::RowMajor>();
template void TestSymBandDecomp<long double,tmv::Lower,tmv::ColMajor>();
template void TestSymBandDecomp<long double,tmv::Lower,tmv::RowMajor>();
#endif

