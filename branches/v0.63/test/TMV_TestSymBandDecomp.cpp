
#include "TMV_Test.h"
#include "TMV_Test2.h"
#include "TMV.h"
#include "TMV_SymBand.h"
#include "TMV_TestSymBandArith.h"

#include "TMV_TestMatrixDivArith.h"

template <class T, tmv::UpLoType uplo, tmv::StorageType stor> 
void TestHermBandDecomp()
{
    for (int mattype = 0; mattype < 9; mattype++) {
        if (showstartdone) {
            std::cout<<"mattype = "<<mattype<<std::endl;
            std::cout<<"uplo, stor = "<<TMV_Text(uplo)<<
                "  "<<TMV_Text(stor)<<std::endl;
        }
        // mattype = 0  is PosDef, nlo = 8
        // mattype = 1  is Indef, nlo = 8
        // mattype = 2  is Singular, nlo = 8
        // mattype = 3  is PosDef, nlo = 1
        // mattype = 4  is Indef, nlo = 1
        // mattype = 5  is Singular, nlo = 1
        // mattype = 6  is PosDef, nlo = 0
        // mattype = 7  is Indef, nlo = 0
        // mattype = 8  is Singular, nlo = 0

        const int N = 200;
        int nlo = 8;
        if (mattype / 3 == 1) nlo = 1;
        else if (mattype / 3 == 2) nlo = 0;

        tmv::HermBandMatrix<T,uplo,stor> m(N,nlo);
        //std::cout<<"m = "<<TMV_Text(m)<<std::endl;
        for(int i=0;i<N;++i) for(int j=0;j<N;++j) 
            if ((uplo == tmv::Upper && j>=i && j<=i+nlo) || 
                (uplo == tmv::Lower && i>=j && i<=j+nlo)) 
                m(i,j) = T(2.5+4*i-5*j);
        m /= T(10*N);
        if (mattype%3 == 0) {
            m += T(10);
            m.diag(0,N/4,N/2) *= T(10);
            m.SubSymBandMatrix(N/2,3*N/4,nlo/2,2) *= T(100);
        } else if (mattype%3 == 1) {
            if (nlo > 0) {
                m(2,1) = T(5);
                m.diag(1,4,7).AddToAll(T(1));
                m.row(17,17-nlo,17).AddToAll(T(20));
                m.diag(0,N/4,N/2) *= T(10);
                m.SubSymBandMatrix(N/2,3*N/4,nlo/2,2) /= T(100);
            }
        } else {
            if (nlo > 0) {
                if (nlo > 2) m.row(2,0,2).Zero();
                else m.row(2,2-nlo,2).Zero();
                m.row(2,2,2+nlo+1).Zero();
                m.row(14,15-nlo,14) = m.row(15,15-nlo,14);
                m.row(14,15,14+nlo+1) = m.row(15,15,14+nlo+1);
                if (14-nlo >= 0) m(14,14-nlo) = T(0);
                if (15+nlo < N) m(15,15+nlo) = T(0);
                m(14,15) = m(14,14) = m(15,15);
            }
            m.SubSymBandMatrix(N/2,3*N/4).Zero();
        }
        //std::cout<<"m = "<<m<<std::endl;

        tmv::HermBandMatrix<std::complex<T>,uplo,stor> c(N,nlo);
        //std::cout<<"c = "<<TMV_Text(c)<<std::endl;
        for(int i=0;i<N;++i) for(int j=0;j<N;++j) 
            if ((uplo == tmv::Upper && j>=i && j<=i+nlo) || 
                (uplo == tmv::Lower && i>=j && i<=j+nlo)) 
                c(i,j) = std::complex<T>(2.5+4*i-5*j,3-i);
        c.diag().Imag().Zero();
        c /= T(10*N);
        if (mattype%3 == 0) {
            c += T(10);
            c.diag(0,N/4,N/2) *= T(10);
            c.SubSymBandMatrix(N/2,3*N/4,nlo/2,2) *= T(100);
        } else if (mattype%3 == 1) {
            if (nlo > 0) {
                c(2,1) = T(5);
                c.diag(1,4,7).AddToAll(T(1));
                c.row(17,17-nlo,17).AddToAll(T(20));
                c.diag(0,N/4,N/2) *= T(10);
                c.SubSymBandMatrix(N/2,3*N/4,nlo/2,2) /= T(100);
            }
        } else {
            if (nlo > 0) {
                if (nlo > 2) c.row(2,0,2).Zero();
                else c.row(2,2-nlo,2).Zero();
                c.row(2,2,2+nlo+1).Zero();
                c.row(14,15-nlo,14) = c.row(15,15-nlo,14);
                c.row(14,15,14+nlo+1) = c.row(15,15,14+nlo+1);
                if (14-nlo >= 0) c(14,14-nlo) = T(0);
                if (15+nlo < N) c(15,15+nlo) = T(0);
                c(14,15) = c(14,14) = c(15,15);
            }
            c.SubSymBandMatrix(N/2,3*N/4).Zero();
        }
        //std::cout<<"c = "<<c<<std::endl;

        T eps = EPS;
        T ceps = EPS;
        if (showacc) std::cout<<"eps = "<<eps<<"  "<<ceps<<std::endl;
        if (mattype % 3 != 2) {
            T kappa = Norm(m) * Norm(m.Inverse());
            if (showacc) std::cout<<"kappa = "<<kappa<<std::endl;
            eps *= kappa;
            ceps *= kappa;
        } else {
            if (showacc) std::cout<<"eps *= "<<T(10*N)<<std::endl;
            eps *= T(10*N);
            ceps *= T(10*N);
        }
        if (showacc) std::cout<<"eps => "<<eps<<"  "<<ceps<<std::endl;


        // CH Decomposition
        if (mattype % 3 == 0) {
            if (showstartdone) std::cout<<"CH"<<std::endl;
            tmv::BandMatrix<T> L(N,N,nlo,0);
            tmv::BandMatrix<T> LLt(N,N,nlo,nlo);
            if (mattype == 0) {
                L = m.CHD().GetL();
                LLt = L*L.Adjoint();
                Assert(Norm(m-LLt) < eps*Norm(m),"HermBand CH");
            }

            tmv::BandMatrix<std::complex<T> > cL(N,N,nlo,0);
            tmv::BandMatrix<std::complex<T> > cLLt(N,N,nlo,nlo);
            if (mattype == 0) {
                cL = c.CHD().GetL();
                cLLt = cL*cL.Adjoint();
                Assert(Norm(c-cLLt) < ceps*Norm(c),"HermBand C CH");
            }

#ifdef XTEST
            tmv::HermBandMatrix<T,uplo,stor> m2 = m;
            CH_Decompose(m2.View());
            L = m2.LowerBand();
            LLt = L*L.Adjoint();
            Assert(Norm(m-LLt) < eps*Norm(m),"HermBand CH2");

            tmv::HermBandMatrix<std::complex<T>,uplo,stor> c2 = c;
            CH_Decompose(c2.View());
            cL = c2.LowerBand();
            cLLt = cL*cL.Adjoint();
            Assert(Norm(c-cLLt) < ceps*Norm(c),"HermBand C CH2");

            c2.Conjugate() = c;
            CH_Decompose(c2.Conjugate());
            cL = c2.Conjugate().LowerBand();
            cLLt = cL*cL.Adjoint();
            Assert(Norm(c-cLLt) < ceps*Norm(c),"HermBand C CH2");
#endif
            std::cout<<"."; std::cout.flush();
        }

        // LDL Decomposition
#ifdef NOTHROW
        if (mattype == 3) {
#else
            if (mattype /3 == 1) try {
#endif
                if (showstartdone) std::cout<<"LDL"<<std::endl;
                tmv::BandMatrix<T> L(N,N,1,0);
                tmv::DiagMatrix<T> D(N);
                tmv::BandMatrix<T> LDL(N,N,1,1);
                if (mattype == 3) {
                    L = m.CHD().GetL();
                    D = m.CHD().GetD();
                    LDL = L*D*L.Adjoint();
                    Assert(Norm(m-LDL) < eps*Norm(m),"HermBand LDL");
                }

                tmv::BandMatrix<std::complex<T> > cL(N,N,nlo,0);
                tmv::DiagMatrix<std::complex<T> > cD(N);
                tmv::BandMatrix<std::complex<T> > cLDL(N,N,nlo,nlo);
                if (mattype == 3) {
                    cL = c.CHD().GetL();
                    cD = c.CHD().GetD();
                    cLDL = cL*cD*cL.Adjoint();
                    Assert(Norm(c-cLDL) < ceps*Norm(c),"HermBand C LDL");
                }

#ifdef XTEST
                tmv::HermBandMatrix<T,uplo,stor> m2 = m;
                LDL_Decompose(m2.View());
                (L = m2.LowerBand()).diag().SetAllTo(T(1));
                D = DiagMatrixViewOf(m2.diag());
                LDL = L*D*L.Adjoint();
                Assert(Norm(m-LDL) < eps*Norm(m),"HermBand LDL2");

                // Alt calculation:
                // (I + L) D (I + Lt)
                // D + LD + DLt = LDLt
                tmv::DiagMatrix<T> E(m2.diag(-1));
                LDL.Zero();
                LDL.diag(0,1,N) = (E*D.SubDiagMatrix(0,N-1)*E).diag();
                LDL.diag(-1) = (E*D.SubDiagMatrix(0,N-1)).diag();
                LDL += D;
                LDL.diag(1) = LDL.diag(-1);
                Assert(Norm(m-LDL) < eps*Norm(m),"HermBand LDL2 Alt");

                tmv::HermBandMatrix<std::complex<T>,uplo,stor> c2 = c;
                LDL_Decompose(c2.View());
                (cL = c2.LowerBand()).diag().SetAllTo(T(1));
                cD = DiagMatrixViewOf(c2.diag());
                cLDL = cL*cD*cL.Adjoint();
                Assert(Norm(c-cLDL) < ceps*Norm(c),"HermBand C LDL2");
                Assert(Norm(cD.Imag()) < ceps*Norm(c),"HermBand C LDL2 - D");

                // Alt calculation:
                tmv::DiagMatrix<std::complex<T> > cE(c2.diag(-1));
                cLDL.Zero();
                cLDL.diag(0,1,N) =
                    (cE*cD.SubDiagMatrix(0,N-1)*cE.Conjugate()).diag();
                cLDL.diag(-1) = (cE*cD.SubDiagMatrix(0,N-1)).diag();
                cLDL += cD;
                cLDL.diag(1) = cLDL.diag(-1).Conjugate();
                Assert(Norm(c-cLDL) < ceps*Norm(c),"HermBand C LDL2 Alt");

                c2.Conjugate() = c;
                LDL_Decompose(c2.Conjugate());
                (cL = c2.Conjugate().LowerBand()).diag().SetAllTo(T(1));
                cD = DiagMatrixViewOf(c2.Conjugate().diag());
                cLDL = cL*cD*cL.Adjoint();
                Assert(Norm(c-cLDL) < ceps*Norm(c),"HermBand C LDL3");
                Assert(Norm(cD.Imag()) < ceps*Norm(c),"HermBand C LDL2 - D");
                // The Lapack version throws whenever mattype is not posdef, 
                // but native algorithm succeeds when mattype is indefinite.
                Assert(mattype%3 != 2,
                       "didn't throw NonPosDef, mattype != singular"); 
#endif
                std::cout<<"."; std::cout.flush();
#ifndef NOTHROW
            } catch (tmv::NonPosDef) { 
                Assert(mattype%3 != 0,"caught NonPosDef, mattype != posdef");
            }
#else
        }
#endif

        // SV Decomposition
        {
            if (showstartdone) std::cout<<"SV"<<std::endl;
            tmv::Matrix<T> U = m.SVD().GetU();
            tmv::DiagMatrix<T> S = m.SVD().GetS();
            tmv::Matrix<T> V = m.SVD().GetV();
            if (showacc) {
                std::cout<<"Norm(U) = "<<Norm(U)<<std::endl;
                std::cout<<"Norm(m-U*S*V) = "<<Norm(m-U*S*V)<<std::endl;
                std::cout<<"Norm(m) = "<<Norm(m)<<std::endl;
                std::cout<<"eps * Norm(m) = "<<eps*Norm(m)<<std::endl;
            }
            Assert(Norm(m-U*S*V) < eps*Norm(m),"HermBand SV");

            tmv::Matrix<std::complex<T> > cU = c.SVD().GetU();
            tmv::DiagMatrix<T> cS = c.SVD().GetS();
            tmv::Matrix<std::complex<T> > cV = c.SVD().GetV();
            if (showacc) {
                std::cout<<"Norm(cU) = "<<Norm(cU)<<std::endl;
                std::cout<<"Norm(c-cU*cS*cV) = "<<Norm(c-cU*cS*cV)<<std::endl;
                std::cout<<"Norm(c) = "<<Norm(c)<<std::endl;
                std::cout<<"ceps * Norm(c) = "<<ceps*Norm(c)<<std::endl;
            }
            Assert(Norm(c-cU*cS*cV) < ceps*Norm(c),"HermBand C SV");

#ifdef XTEST
            tmv::Matrix<T> U2(N,N);
            tmv::DiagMatrix<T> S2(N);
            tmv::Matrix<T> V2(N,N);
            SV_Decompose(m,U2.View(),S2.View(),V2.View());
            if (showacc) {
                std::cout<<"Norm(U) = "<<Norm(U2)<<std::endl;
                std::cout<<"Norm(m-U*S*V) = "<<Norm(m-U2*S2*V2)<<std::endl;
                std::cout<<"Norm(m) = "<<Norm(m)<<std::endl;
                std::cout<<"eps * Norm(m) = "<<eps*Norm(m)<<std::endl;
            }
            Assert(Norm(m-U2*S2*V2) < eps*Norm(m),"HermBand SV2");

            SV_Decompose(m,S2.View());
            Assert(Norm(S2-S) < eps*Norm(S),"HermBand SV3");
            SV_Decompose(m,U2.View(),S2.View());
            Assert(Norm(S2-S) < eps*Norm(S),"HermBand SV4 S");
            if (showacc) {
                //std::cout<<"U = "<<U2<<std::endl;
                //std::cout<<"S = "<<S2<<std::endl;
                std::cout<<"Norm(U*S*S*Ut-m*mt) = "<<
                    Norm(U2*S2*S2*U2.Adjoint()-m*m.Adjoint())<<std::endl;
                std::cout<<"Norm(m) = "<<Norm(m)<<std::endl;
                std::cout<<"eps * Norm(m)^2 = "<<eps*Norm(m)*Norm(m)<<std::endl;
                std::cout<<"eps = "<<eps<<std::endl;
                std::cout<<"Norm(U) = "<<Norm(U2)<<std::endl;
                std::cout<<"Norm(S) = "<<Norm(S2)<<std::endl;
                std::cout<<"eps * Norm(U)^2 * Norm(S)^2 = "<<
                    eps*Norm(U2)*Norm(U2)*Norm(S2)*Norm(S2)<<std::endl;
            }
            Assert(Norm(U2*S2*S2*U2.Adjoint()-m*m.Adjoint()) < 
                   eps*Norm(m)*Norm(m),"HermBand SV4 U");
            SV_Decompose(m,S2.View(),V2.View());
            Assert(Norm(S2-S) < eps*Norm(S),"HermBand SV5 S");
            Assert(Norm(m.Adjoint()*m-V2.Adjoint()*S2*S2*V2) < 
                   eps*Norm(m)*Norm(m),"HermBand SV5 V");

            tmv::Matrix<std::complex<T> > cU2(N,N);
            tmv::DiagMatrix<T> cS2(N);
            tmv::Matrix<std::complex<T> > cV2(N,N);
            SV_Decompose(c,cU2.View(),cS2.View(),cV2.View());
            Assert(Norm(c-cU2*cS2*cV2) < ceps*Norm(c),"HermBand C SV2");

            SV_Decompose(c,cS2.View());
            Assert(Norm(cS2-cS) < ceps*Norm(cS),"HermBand C SV3");
            SV_Decompose(c,cU2.View(),cS2.View());
            Assert(Norm(cS2-cS) < ceps*Norm(cS),"HermBand C SV4 S");
            Assert(Norm(cU2*cS2*cS2*cU2.Adjoint()-c*c.Adjoint()) < 
                   ceps*Norm(c)*Norm(c),"HermBand C SV4 U");
            SV_Decompose(c,cS2.View(),cV2.View());
            Assert(Norm(cS2-cS) < ceps*Norm(cS),"HermBand C SV5 S");
            Assert(Norm(c.Adjoint()*c-cV2.Adjoint()*cS2*cS2*cV2) < 
                   ceps*Norm(c)*Norm(c),"Herm C SV5 V");

            SV_Decompose(c.Conjugate(),cS2.View());
            Assert(Norm(cS2-cS) < ceps*Norm(cS),"HermBand C SV6");
            SV_Decompose(c.Conjugate(),cU2.View(),cS2.View());
            Assert(Norm(cS2-cS) < ceps*Norm(cS),"HermBand C SV7 S");
            Assert(Norm(cU2*cS2*cS2*cU2.Adjoint()-c.Conjugate()*c.Transpose()) <
                   ceps*Norm(c)*Norm(c),"HermBand C SV7 U");
            SV_Decompose(c.Conjugate(),cU2.Conjugate(),cS2.View());
            Assert(Norm(cS2-cS) < ceps*Norm(cS),"HermBand C SV8 S");
            Assert(Norm(cU2.Conjugate()*cS2*cS2*cU2.Transpose()-
                        c.Conjugate()*c.Transpose()) < ceps*Norm(c)*Norm(c),
                   "HermBand C SV8 U");
#endif
            std::cout<<"."; std::cout.flush();
        }

        // Eigen
        {
            if (showstartdone) std::cout<<"Eigen"<<std::endl;
            tmv::Matrix<T> V(N,N);
            tmv::Vector<T> L(N);
            Eigen(m,V.View(),L.View());
            Assert(Norm(m*V-V*DiagMatrixViewOf(L)) < eps*Norm(m),
                   "HermBand Eigen");

            tmv::Matrix<std::complex<T> > cV(N,N);
            tmv::Vector<T> cL(N);
            Eigen(c,cV.View(),cL.View());
            Assert(Norm(c*cV-cV*DiagMatrixViewOf(cL)) < ceps*Norm(c),
                   "HermBand C Eigen");

#ifdef XTEST
            tmv::Vector<T> L2(N);
            Eigen(m,L2.View());
            Assert(Norm(L2-L) < eps*Norm(L),"HermBand Eigen2");

            tmv::Vector<T> cL2(N);
            Eigen(c,cL2.View());
            Assert(Norm(cL2-cL) < ceps*Norm(cL),"HermBand C Eigen2");

            Eigen(c,cV.Conjugate(),cL.View());
            Assert(Norm(c*cV.Conjugate()-cV.Conjugate()*DiagMatrixViewOf(cL)) < 
                   ceps*Norm(c),"HermBand C Eigen3");

            Eigen(c.Conjugate(),cL2.View());
            Assert(Norm(cL2-cL) < ceps*Norm(cL),"HermBand C Eigen4");
#endif
            std::cout<<"."; std::cout.flush();
        }

        // Square Root
#ifdef NOTHROW
        if (mattype%3 == 0) {
#else
            try {
#endif
                if (showstartdone) std::cout<<"Square Root"<<std::endl;
                tmv::HermMatrix<T,uplo,stor> S(N);
                SquareRoot(m,S.View());
                //std::cout<<"SquareRoot of "<<m<<" is "<<S<<std::endl;
                Assert(Norm(m-S*S) < eps*Norm(m),"HermBand Square Root");

                tmv::HermMatrix<std::complex<T>,uplo,stor> cS(N);
                SquareRoot(c,cS.View());
                //std::cout<<"SquareRoot of "<<c<<" is "<<cS<<std::endl;
                Assert(Norm(c-cS*cS) < eps*Norm(c),"HermBand C Square Root");

#ifdef XTEST
                SquareRoot(c,cS.Conjugate());
                Assert(Norm(c-cS.Conjugate()*cS.Conjugate()) < ceps*Norm(c),
                       "HermBand C Square Root2");

                SquareRoot(c,cS.Transpose());
                Assert(Norm(c-cS.Transpose()*cS.Transpose()) < ceps*Norm(c),
                       "HermBand C Square Root3");

                SquareRoot(c,cS.Adjoint());
                Assert(Norm(c-cS.Adjoint()*cS.Adjoint()) < ceps*Norm(c),
                       "HermBand C Square Root4");

                SquareRoot(c.Conjugate(),cS.View());
                Assert(Norm(c.Conjugate()-cS*cS) < ceps*Norm(c),
                       "HermBand C Square Root5");

                SquareRoot(c.Conjugate(),cS.Conjugate());
                Assert(Norm(c.Conjugate()-cS.Conjugate()*cS.Conjugate()) < 
                       ceps*Norm(c),
                       "HermBand C Square Root6");

                SquareRoot(c.Conjugate(),cS.Transpose());
                Assert(Norm(c.Conjugate()-cS.Transpose()*cS.Transpose()) < 
                       ceps*Norm(c),
                       "HermBand C Square Root7");

                SquareRoot(c.Conjugate(),cS.Adjoint());
                Assert(Norm(c.Conjugate()-cS.Adjoint()*cS.Adjoint()) < 
                       ceps*Norm(c),
                       "HermBand C Square Root8");
#endif
                Assert(mattype%3 == 0,
                       "didn't throw NonPosDef, mattype == pos def"); 
                std::cout<<"."; std::cout.flush();
#ifndef NOTHROW
            } catch (tmv::NonPosDef)
            { 
                Assert(mattype%3 != 0,"caught NonPosDef, mattype != posdef"); 
            }
#else
        }
#endif
    }
}

template <class T, tmv::UpLoType uplo, tmv::StorageType stor> 
void TestSymBandDecomp()
{
    for (int mattype = 0; mattype < 9; mattype++) {
        if (showstartdone) {
            std::cout<<"mattype = "<<mattype<<std::endl;
            std::cout<<"uplo, stor = "<<TMV_Text(uplo)<<"  "<<TMV_Text(stor)<<std::endl;
        }
        // mattype = 0  is PosDef, nlo = 8
        // mattype = 1  is Indef, nlo = 8
        // mattype = 2  is Singular, nlo = 8
        // mattype = 3  is PosDef, nlo = 1
        // mattype = 4  is Indef, nlo = 1
        // mattype = 5  is Singular, nlo = 1
        // mattype = 6  is PosDef, nlo = 0
        // mattype = 7  is Indef, nlo = 0
        // mattype = 8  is Singular, nlo = 0

        const int N = 200;
        int nlo = 8;
        if (mattype / 3 == 1) nlo = 1;
        else if (mattype / 3 == 2) nlo = 0;

        tmv::SymBandMatrix<T,uplo,stor> m(N,nlo);
        for(int i=0;i<N;++i) for(int j=0;j<N;++j) 
            if ((uplo == tmv::Upper && j>=i && j<=i+nlo) || 
                (uplo == tmv::Lower && i>=j && i<=j+nlo)) 
                m(i,j) = T(2.5+4*i-5*j);
        m /= T(10*N);
        if (mattype%3 == 0) {
            m += T(10);
            m.diag(0,N/4,N/2) *= T(10);
            m.SubSymBandMatrix(N/2,3*N/4,nlo/2,2) *= T(100);
        } else if (mattype%3 == 1) {
            if (nlo > 0) {
                m(2,1) = T(5);
                m.diag(1,4,7).AddToAll(T(1));
                m.row(17,17-nlo,17).AddToAll(T(20));
                m.diag(0,N/4,N/2) *= T(10);
                m.SubSymBandMatrix(N/2,3*N/4,nlo/2,2) /= T(100);
            }
        } else {
            if (nlo > 0) {
                if (nlo > 2) m.row(2,0,2).Zero();
                else m.row(2,2-nlo,2).Zero();
                m.row(2,2,2+nlo+1).Zero();
                m.row(14,15-nlo,14) = m.row(15,15-nlo,14);
                m.row(14,15,14+nlo+1) = m.row(15,15,14+nlo+1);
                if (14-nlo >= 0) m(14,14-nlo) = T(0);
                if (15+nlo < N) m(15,15+nlo) = T(0);
                m(14,15) = m(14,14) = m(15,15);
            }
            m.SubSymBandMatrix(N/2,3*N/4).Zero();
        }
        //std::cout<<"m = "<<m<<std::endl;

        tmv::SymBandMatrix<std::complex<T>,uplo,stor> c(N,nlo);
        for(int i=0;i<N;++i) for(int j=0;j<N;++j) 
            if ((uplo == tmv::Upper && j>=i && j<=i+nlo) || 
                (uplo == tmv::Lower && i>=j && i<=j+nlo)) 
                c(i,j) = std::complex<T>(2+4*i-5*j,3-i);
        c /= T(10*N);
        if (mattype%3 == 0) {
            c += T(10);
            c.diag(0,N/4,N/2) *= T(10);
            c.SubSymBandMatrix(N/2,3*N/4,nlo/2,2) *= T(100);
        } else if (mattype%3 == 1) {
            if (nlo > 0) {
                c(2,1) = T(5);
                c.diag(1,4,7).AddToAll(T(1));
                c.row(17,17-nlo,17).AddToAll(T(20));
                c.diag(0,N/4,N/2) *= T(10);
                c.SubSymBandMatrix(N/2,3*N/4,nlo/2,2) /= T(100);
            }
        } else {
            if (nlo > 0) {
                if (nlo > 2) c.row(2,0,2).Zero();
                else c.row(2,2-nlo,2).Zero();
                c.row(2,2,2+nlo+1).Zero();
                c.row(14,15-nlo,14) = c.row(15,15-nlo,14);
                c.row(14,15,14+nlo+1) = c.row(15,15,14+nlo+1);
                if (14-nlo >= 0) c(14,14-nlo) = T(0);
                if (15+nlo < N) c(15,15+nlo) = T(0);
                c(14,15) = c(14,14) = c(15,15);
            }
            c.SubSymBandMatrix(N/2,3*N/4).Zero();
        }
        //std::cout<<"c = "<<c<<std::endl;

        T eps = EPS;
        T ceps = EPS;
        if (mattype %3 != 2) {
            eps *= Norm(m) * Norm(m.Inverse());
            ceps *= Norm(c) * Norm(c.Inverse());
        } else {
            eps *= T(10*N);
            ceps *= T(10*N);
        }

#ifdef XTEST
        // LDL Decomposition
#ifdef NOTHROW
        if (mattype == 3)
#else
            if (mattype /3 == 1) try 
#endif
            {
                if (showstartdone) std::cout<<"LDL"<<std::endl;
                tmv::BandMatrix<T> L(N,N,1,0);
                tmv::DiagMatrix<T> D(N);
                tmv::BandMatrix<T> LDL(N,N,1,1);

                tmv::SymBandMatrix<T,uplo,stor> m2 = m;
                LDL_Decompose(m2.View());
                (L = m2.LowerBand()).diag().SetAllTo(T(1));
                D = DiagMatrixViewOf(m2.diag());
                LDL = L*D*L.Transpose();
                Assert(Norm(m-LDL) < eps*Norm(m),"SymBand LDL2");

                // Alt calculation:
                // (I + L) D (I + Lt)
                // D + LD + DLt = LDLt
                tmv::DiagMatrix<T> E(m2.diag(-1));
                LDL.Zero();
                LDL.diag(0,1,N) = (E*D.SubDiagMatrix(0,N-1)*E).diag();
                LDL.diag(-1) = (E*D.SubDiagMatrix(0,N-1)).diag();
                LDL += D;
                LDL.diag(1) = LDL.diag(-1);
                Assert(Norm(m-LDL) < eps*Norm(m),"SymBand LDL2 Alt");

                tmv::BandMatrix<std::complex<T> > cL(N,N,nlo,0);
                tmv::DiagMatrix<std::complex<T> > cD(N);
                tmv::BandMatrix<std::complex<T> > cLDL(N,N,nlo,nlo);

                tmv::SymBandMatrix<std::complex<T>,uplo,stor> c2 = c;
                LDL_Decompose(c2.View());
                (cL = c2.LowerBand()).diag().SetAllTo(T(1));
                cD = DiagMatrixViewOf(c2.diag());
                cLDL = cL*cD*cL.Transpose();
                Assert(Norm(c-cLDL) < ceps*Norm(c),"SymBand C LDL2");

                // Alt calculation:
                tmv::DiagMatrix<std::complex<T> > cE(c2.diag(-1));
                cLDL.Zero();
                cLDL.diag(0,1,N) = (cE*cD.SubDiagMatrix(0,N-1)*cE).diag();
                cLDL.diag(-1) = (cE*cD.SubDiagMatrix(0,N-1)).diag();
                cLDL += cD;
                cLDL.diag(1) = cLDL.diag(-1);
                Assert(Norm(c-cLDL) < eps*Norm(m),"SymBand C LDL2 Alt");

                c2.Conjugate() = c;
                LDL_Decompose(c2.Conjugate());
                (cL = c2.Conjugate().LowerBand()).diag().SetAllTo(T(1));
                cD = DiagMatrixViewOf(c2.Conjugate().diag());
                cLDL = cL*cD*cL.Transpose();
                Assert(Norm(c-cLDL) < ceps*Norm(c),"SymBand C LDL3");
                Assert(mattype%3 != 2,"didn't throw NonPosDef, mattype != singular"); 
                std::cout<<"."; std::cout.flush();
            }
#ifndef NOTHROW
        catch (tmv::NonPosDef)
        { Assert(mattype%3 != 0,"caught NonPosDef, mattype != posdef"); }
#endif
#endif

        // SV Decomposition
        {
            if (showstartdone) std::cout<<"SV"<<std::endl;
            tmv::Matrix<T> U = m.SVD().GetU();
            tmv::DiagMatrix<T> S = m.SVD().GetS();
            tmv::Matrix<T> V = m.SVD().GetV();
            Assert(Norm(m-U*S*V) < eps*Norm(m),"SymBand SV");

            tmv::Matrix<std::complex<T> > cU = c.SymSVD().GetU();
            tmv::DiagMatrix<T> cS = c.SymSVD().GetS();
            tmv::Matrix<std::complex<T> > cV = c.SymSVD().GetV();
            Assert(Norm(c-cU*cS*cV) < ceps*Norm(c),"SymBand C SV");

#ifdef XTEST
            tmv::Matrix<T> U2(N,N);
            tmv::DiagMatrix<T> S2(N);
            tmv::Matrix<T> V2(N,N);
            SV_Decompose(m,U2.View(),S2.View(),V2.View());
            Assert(Norm(m-U2*S2*V2) < eps*Norm(m),"SymBand SV2");

            SV_Decompose(m,S2.View());
            Assert(Norm(S2-S) < eps*Norm(S),"SymBand SV3");
            SV_Decompose(m,U2.View(),S2.View());
            Assert(Norm(S2-S) < eps*Norm(S),"SymBand SV4 S");
            if (showacc) {
                std::cout<<"Norm(U*S*S*UT-m*mT) = "<<
                    Norm(U2*S2*S2*U2.Transpose()-m*m.Transpose())<<std::endl;
                std::cout<<"Norm(m) = "<<Norm(m)<<std::endl;
                std::cout<<"eps * Norm(m)^2 = "<<eps*Norm(m)*Norm(m)<<std::endl;
                std::cout<<"eps = "<<eps<<std::endl;
                std::cout<<"Norm(U) = "<<Norm(U2)<<std::endl;
                std::cout<<"Norm(S) = "<<Norm(S2)<<std::endl;
                std::cout<<"eps * Norm(U)^2 * Norm(S)^2 = "<<
                    eps*Norm(U2)*Norm(U2)*Norm(S2)*Norm(S2)<<std::endl;
            }
            Assert(Norm(U2*S2*S2*U2.Transpose()-m*m.Transpose()) <
                   eps*Norm(m)*Norm(m),"SymBand SV4 U");
            SV_Decompose(m,S2.View(),V2.View());
            Assert(Norm(S2-S) < eps*Norm(S),"SymBand SV5 S");
            Assert(Norm(m.Adjoint()*m-V2.Adjoint()*S2*S2*V2) < eps*Norm(m)*Norm(m),
                   "SymBand SV5 V");

            tmv::Matrix<std::complex<T> > cU2(N,N);
            tmv::DiagMatrix<T> cS2(N);
            tmv::Matrix<std::complex<T> > cV2(N,N);
            SV_Decompose(c,cU2.View(),cS2.View(),cV2.View());
            Assert(Norm(c-cU2*cS2*cV2) < ceps*Norm(c),"SymBand C SV2");

            SV_Decompose(c,cS2.View());
            Assert(Norm(cS2-cS) < ceps*Norm(cS),"SymBand C SV3");
            SV_Decompose(c,cU2.View(),cS2.View());
            Assert(Norm(cS2-cS) < ceps*Norm(cS),"SymBand C SV4 S");
            Assert(Norm(cU2*cS2*cS2*cU2.Adjoint()-c*c.Adjoint()) <
                   ceps*Norm(c*c.Adjoint()),"SymBand C SV4 U");
            SV_Decompose(c,cS2.View(),cV2.View());
            Assert(Norm(cS2-cS) < ceps*Norm(cS),"SymBand C SV5 S");
            Assert(Norm(c.Adjoint()*c-cV2.Adjoint()*cS2*cS2*cV2) < ceps*Norm(c)*Norm(c),
                   "SymBand C SV5 V");

            SV_Decompose(c.Conjugate(),cS2.View());
            Assert(Norm(cS2-cS) < ceps*Norm(cS),"SymBand C SV6");
            SV_Decompose(c.Conjugate(),cU2.View(),cS2.View());
            Assert(Norm(cS2-cS) < ceps*Norm(cS),"SymBand C SV7 S");
            Assert(Norm(cU2*cS2*cS2*cU2.Adjoint()-c.Conjugate()*c.Transpose()) <
                   ceps*Norm(c.Conjugate()*c.Transpose()),"SymBand C SV7 U");
            SV_Decompose(c.Conjugate(),cU2.Conjugate(),cS2.View());
            Assert(Norm(cS2-cS) < ceps*Norm(cS),"SymBand C SV8 S");
            Assert(Norm(cU2.Conjugate()*cS2*cS2*cU2.Transpose()-
                        c.Conjugate()*c.Transpose()) <
                   ceps*Norm(c.Conjugate()*c.Transpose()),"SymBand C SV8 U");
#endif
            std::cout<<"."; std::cout.flush();
        }
    }
}

#ifdef INST_DOUBLE
template void TestHermBandDecomp<double,tmv::Upper,tmv::ColMajor>();
template void TestHermBandDecomp<double,tmv::Upper,tmv::RowMajor>();
template void TestHermBandDecomp<double,tmv::Lower,tmv::ColMajor>();
template void TestHermBandDecomp<double,tmv::Lower,tmv::RowMajor>();
template void TestSymBandDecomp<double,tmv::Upper,tmv::ColMajor>();
template void TestSymBandDecomp<double,tmv::Upper,tmv::RowMajor>();
template void TestSymBandDecomp<double,tmv::Lower,tmv::ColMajor>();
template void TestSymBandDecomp<double,tmv::Lower,tmv::RowMajor>();
#endif
#ifdef INST_FLOAT
template void TestHermBandDecomp<float,tmv::Upper,tmv::ColMajor>();
template void TestHermBandDecomp<float,tmv::Upper,tmv::RowMajor>();
template void TestHermBandDecomp<float,tmv::Lower,tmv::ColMajor>();
template void TestHermBandDecomp<float,tmv::Lower,tmv::RowMajor>();
template void TestSymBandDecomp<float,tmv::Upper,tmv::ColMajor>();
template void TestSymBandDecomp<float,tmv::Upper,tmv::RowMajor>();
template void TestSymBandDecomp<float,tmv::Lower,tmv::ColMajor>();
template void TestSymBandDecomp<float,tmv::Lower,tmv::RowMajor>();
#endif
#ifdef INST_LONGDOUBLE
template void TestHermBandDecomp<long double,tmv::Upper,tmv::ColMajor>();
template void TestHermBandDecomp<long double,tmv::Upper,tmv::RowMajor>();
template void TestHermBandDecomp<long double,tmv::Lower,tmv::ColMajor>();
template void TestHermBandDecomp<long double,tmv::Lower,tmv::RowMajor>();
template void TestSymBandDecomp<long double,tmv::Upper,tmv::ColMajor>();
template void TestSymBandDecomp<long double,tmv::Upper,tmv::RowMajor>();
template void TestSymBandDecomp<long double,tmv::Lower,tmv::ColMajor>();
template void TestSymBandDecomp<long double,tmv::Lower,tmv::RowMajor>();
#endif

