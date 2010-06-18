
#define START 0

#include "../src/TMV_Blas.h"
#include "TMV.h"
#include "TMV_Sym.h"
#include "TMV_Band.h"
#include "TMV_SymBand.h"
#include "TMV_Test.h"
#include "TMV_Test2.h"
#include "TMV_TestSymArith.h"

template <class T, tmv::UpLoType uplo, tmv::StorageType stor> 
void TestHermDecomp()
{
    for (int mattype = START; mattype <= 5; mattype++) {
#if !(XTEST & 64)
        //if (mattype >= 3) break;
#endif
#ifdef LAP
        if (mattype >= 4) break;
#endif
        if (showstartdone) {
            std::cout<<"Herm: mattype = "<<mattype<<std::endl;
            std::cout<<"uplo, stor = "<<TMV_Text(uplo)<<
                "  "<<TMV_Text(stor)<<std::endl;
        }
        // mattype = 0  is PosDef
        // mattype = 1  is Indef
        // mattype = 2  is Singular
        // mattype = 3  is Singular with seriously bad defects
        // mattype = 4  is Singular and nearly zero
        // mattype = 5  is Singular and nearly overflow

        const int N = 200; // Must be multiple of 8
        const bool posdef = mattype == 0;
        const bool singular = mattype >= 2;
        const bool baddefect = mattype == 3;
        const bool nearunderflow = mattype == 4;
        const bool nearoverflow = mattype == 5;
        if (showstartdone) {
            std::cout<<"posdef = "<<posdef<<
                ", singular = "<<singular<<std::endl;
        }

        tmv::HermMatrix<T,uplo,stor> m(N);
        for(int i=0;i<N;++i) for(int j=0;j<N;++j) 
            if ((uplo == tmv::Upper && j>=i) || (uplo == tmv::Lower && i>=j)) 
                m(i,j) = T(2+4*i-5*j);
        if (posdef) {
            m /= T(10*N);
            m.diag().addToAll(T(1));
            for(int i=0;i<N;i++) m.diag()(i) += T(i);
            m.diag(0,N/4,N/2) *= T(10);
            m.subSymMatrix(N/2,3*N/4,2) *= T(100);
        } else if (singular) {
            m.row(2,0,2).setZero();
            m.row(2,2,N).setZero();
            m.row(4,0,4) = m.row(5,0,4);
            m.row(4,5,N) = m.row(5,5,N);
            m(4,5) = m(4,4) = m(5,5);
        } else {
            // indef
            m /= T(N);
            m.diag(0,N/10,N/5).setZero();
            m.col(0,N/3,3*N/7).addToAll(T(50));
            m.diag(1,4,7).addToAll(T(10));
            m.row(7,0,5).addToAll(T(10));
            m.diag(0,N/4,N/2) *= T(10);
            m.subSymMatrix(N/2,3*N/4,2) /= T(100);
        }

        tmv::HermMatrix<std::complex<T>,uplo,stor> c(N);
        for(int i=0;i<N;++i) for(int j=0;j<N;++j) 
            if ((uplo == tmv::Upper && j>=i) || (uplo == tmv::Lower && i>=j)) 
                c(i,j) = std::complex<T>(2+4*i-5*j,3-i);
        c.diag().imagPart().setZero();
        if (posdef) {
            c /= T(10*N);
            c.diag().addToAll(T(1));
            for(int i=0;i<N;i++) c.diag()(i) += T(i);
            c.diag(0,N/4,N/2) *= T(10);
            c.subSymMatrix(N/2,3*N/4,2) *= T(100);
        } else if (singular) {
            c.row(2,0,2).setZero();
            c.row(2,2,N).setZero();
            c.row(4,0,4) = c.row(5,0,4);
            c.row(4,5,N) = c.row(5,5,N);
            c(4,5) = c(4,4) = c(5,5);
        } else {
            // indef
            c /= T(N);
            c.diag(0,N/10,N/5).setZero();
            c.col(0,N/3,3*N/7).addToAll(T(50));
            c.diag(1,4,7).addToAll(T(10));
            c.row(7,0,5).addToAll(T(10));
            c.diag(0,N/4,N/2) *= T(10);
            c.subSymMatrix(N/2,3*N/4,2) /= T(100);
        }

        if (baddefect) {
            for(int i=10;i<N;i+=10) {
                m.subMatrix(0,i,i,N) *= T(1.e-10);
                m.subSymMatrix(i,N) *= T(-1.e-10);
                c.subMatrix(0,i,i,N) *= T(1.e-10);
                c.subSymMatrix(i,N) *= T(-1.e-10);
            }
        }

        if (nearunderflow) {
            T x = std::numeric_limits<T>::min();
            m.setAllTo(x);
            c.upperTri().offDiag().setAllTo(std::complex<T>(x,2*x));
            c.diag().setAllTo(x);
            for(int i=1;i<N;++i) {
                m.col(i,0,i+1) /= T(i);
                c.col(i,0,i+1) /= T(i);
            }
            m(N/2,N/2) = x/std::numeric_limits<T>::epsilon();
            c(N/2,N/2) = x/std::numeric_limits<T>::epsilon();
            m.diag(0,N/5,2*N/5).setZero();
            c.diag(0,N/5,2*N/5).setZero();
            m.diag(0,N/2+1,N) *= T(-1);
            c.diag(0,N/2+1,N) *= T(-1);
            m.subSymMatrix(N-4,N).setAllTo(
                -x/std::numeric_limits<T>::epsilon()/T(8));
            c.subSymMatrix(N-4,N).setAllTo(
                -x/std::numeric_limits<T>::epsilon()/T(8));
        }

        if (nearoverflow) {
            T x = std::numeric_limits<T>::max();
            x /= N;
            m.setAllTo(x);
            c.upperTri().offDiag().setAllTo(std::complex<T>(x,x/2));
            c.diag().setAllTo(x);
            for(int i=1;i<N;++i) {
                m.col(i,0,i+1) /= T(i);
                c.col(i,0,i+1) /= T(i);
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
        if (mattype < 2) {
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
            try  {
#endif
                if (showstartdone) std::cout<<"CH"<<std::endl;
                tmv::LowerTriMatrix<T> L = m.chd().getL();
                tmv::Matrix<T> LLt = L*L.adjoint();
                Assert(Norm(m-LLt) <= eps*normm,"Herm CH");

                tmv::LowerTriMatrix<std::complex<T> > cL = c.chd().getL();
                tmv::Matrix<std::complex<T> > cLLt = cL*cL.adjoint();
                Assert(Norm(c-cLLt) <= ceps*normc,"Herm C CH");

#if (XTEST & 16)
                tmv::HermMatrix<T,uplo,stor> m2 = m;
                CH_Decompose(m2.view());
                L = m2.lowerTri();
                LLt = L*L.adjoint();
                Assert(Norm(m-LLt) <= eps*normm,"Herm CH2");

                tmv::HermMatrix<std::complex<T>,uplo,stor> c2 = c;
                CH_Decompose(c2.view());
                cL = c2.lowerTri();
                cLLt = cL*cL.adjoint();
                Assert(Norm(c-cLLt) <= ceps*normc,"Herm C CH2");

                c2.conjugate() = c;
                CH_Decompose(c2.conjugate());
                cL = c2.conjugate().lowerTri();
                cLLt = cL*cL.adjoint();
                Assert(Norm(c-cLLt) <= ceps*normc,"Herm C CH2");
                Assert(posdef, "Didn't throw NonPosDef, and mattype != pos def");
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
        try {
            do {
#ifdef LAP
                if (baddefect || nearunderflow || nearoverflow) break;
#endif
                if (showstartdone) std::cout<<"LDL"<<std::endl;
                tmv::LowerTriMatrix<T,tmv::UnitDiag> L = m.lud().getL();
                tmv::BandMatrix<T> D = m.lud().getD();
                tmv::Permutation P = m.lud().getP();
                tmv::Matrix<T> LDL = P*L*D*L.adjoint()*P.transpose();
                Assert(Norm(m-LDL) <= eps*normm,"Herm LDL");

                tmv::LowerTriMatrix<std::complex<T>,tmv::UnitDiag> cL = c.lud().getL();
                tmv::BandMatrix<std::complex<T> > cD = c.lud().getD();
                tmv::Permutation cP = c.lud().getP();
                tmv::Matrix<std::complex<T> > cLDL = 
                    cP*cL*cD*cL.adjoint()*cP.transpose();
                Assert(Norm(c-cLDL) <= ceps*normc,"Herm C LDL");

#if (XTEST & 16)
                tmv::HermMatrix<T,uplo,stor> m2 = m;
                tmv::HermBandMatrix<T,uplo,stor> D2(N,1);
                tmv::Permutation P2(N);
                LDL_Decompose(m2.view(),D2.view(),P2);
                L = m2.lowerTri(tmv::UnitDiag);
                LDL = P2*L*D2*L.adjoint()*P2.transpose();
                Assert(Norm(m-LDL) <= eps*normm,"Herm LDL2");

                tmv::HermMatrix<std::complex<T>,uplo,stor> c2 = c;
                tmv::HermBandMatrix<std::complex<T>,uplo,stor> cD2(N,1);
                LDL_Decompose(c2.view(),cD2.view(),P2);
                cL = c2.lowerTri(tmv::UnitDiag);
                cLDL = P2*cL*cD2*cL.adjoint()*P2.transpose();
                Assert(Norm(c-cLDL) <= ceps*normc,"Herm C LDL2");

                c2.conjugate() = c;
                LDL_Decompose(c2.conjugate(),cD2.view(),P2);
                cL = c2.conjugate().lowerTri(tmv::UnitDiag);
                cLDL = P2*cL*cD2*cL.adjoint()*P2.transpose();
                Assert(Norm(c-cLDL) <= ceps*normc,"Herm C LDL3");

                c2 = c;
                LDL_Decompose(c2.view(),cD2.conjugate(),P2);
                cL = c2.lowerTri(tmv::UnitDiag);
                cLDL = P2*cL*cD2.conjugate()*cL.adjoint()*P2.transpose();
                Assert(Norm(c-cLDL) <= ceps*normc,"Herm C LDL4");

                c2.conjugate() = c;
                LDL_Decompose(c2.conjugate(),cD2.conjugate(),P2);
                cL = c2.conjugate().lowerTri(tmv::UnitDiag);
                cLDL = P2*cL*cD2.conjugate()*cL.adjoint()*P2.transpose();
                Assert(Norm(c-cLDL) <= ceps*normc,"Herm C LDL5");
#endif
                std::cout<<"."; std::cout.flush();
            } while (false);
        } catch (tmv::NonPosDef) {
            // The Lapack version throws whenever mattype is not posdef, 
            // but native algorithm succeeds when mattype is indefinite,
            // or even singular.  
            Assert(!posdef,"caught NonPosDef but mattype == posdef");
        }

        // SV Decomposition
        {
            if (showstartdone) std::cout<<"SV"<<std::endl;
            tmv::Matrix<T> U = m.svd().getU();
            tmv::DiagMatrix<T> S = m.svd().getS();
            tmv::Matrix<T> V = m.svd().getV();
            if (showacc) {
                std::cout<<"Norm(m-U*S*V) = "<<Norm(m-U*S*V)<<std::endl;
                std::cout<<"eps*Norm(m) = "<<eps*normm<<std::endl;
            }
            Assert(Norm(m-U*S*V) <= eps*normm,"Herm SV");

            tmv::Matrix<std::complex<T> > cU = c.svd().getU();
            tmv::DiagMatrix<T> cS = c.svd().getS();
            tmv::Matrix<std::complex<T> > cV = c.svd().getV();
            Assert(Norm(c-cU*cS*cV) <= ceps*normc,"Herm C SV");

#if (XTEST & 16)
            tmv::Matrix<T> U2(N,N);
            tmv::DiagMatrix<T> S2(N);
            tmv::Matrix<T> V2(N,N);
            SV_Decompose(m,U2.view(),S2.view(),V2.view());
            Assert(Norm(m-U2*S2*V2) <= eps*normm,"Herm SV2");

            tmv::HermMatrix<T,uplo,stor> m2 = m;
            SV_Decompose(m2.view(),S2.view());
            Assert(Norm(S2-S) <= eps*normm,"Herm SV3");
            SV_Decompose(m,U2.view(),S2.view());
            Assert(Norm(S2-S) <= eps*normm,"Herm SV4 S");
            if (!nearoverflow)
                Assert(Norm(m*m.transpose()-U2*S2*S2*U2.adjoint()) <= 
                       eps*normm*normm,"Herm SV4 U");
            SV_Decompose(m,S2.view(),V2.view());
            Assert(Norm(S2-S) <= eps*normm,"Herm SV5 S");
            if (!nearoverflow)
                Assert(Norm(m.adjoint()*m-V2.adjoint()*S2*S2*V2) <= 
                       eps*normm*normm,"Herm SV5 V");

            tmv::Matrix<std::complex<T> > cU2(N,N);
            tmv::DiagMatrix<T> cS2(N);
            tmv::Matrix<std::complex<T> > cV2(N,N);
            SV_Decompose(c,cU2.view(),cS2.view(),cV2.view());
            Assert(Norm(c-cU2*cS2*cV2) <= ceps*normc,"Herm C SV2");

            tmv::HermMatrix<std::complex<T>,uplo,stor> c2 = c;
            SV_Decompose(c2.view(),cS2.view());
            Assert(Norm(cS2-cS) <= ceps*normc,"Herm C SV3");
            SV_Decompose(c,cU2.view(),cS2.view());
            Assert(Norm(cS2-cS) <= ceps*normc,"Herm C SV4 S");
            if (!nearoverflow)
                Assert(Norm(c*c.adjoint()-cU2*cS2*cS2*cU2.adjoint()) <= 
                       ceps*normc*normc,"Herm C SV4 U");
            SV_Decompose(c,cS2.view(),cV2.view());
            Assert(Norm(cS2-cS) <= ceps*normc,"Herm C SV5 S");
            if (!nearoverflow)
                Assert(Norm(c.adjoint()*c-cV2.adjoint()*cS2*cS2*cV2) <= 
                       ceps*normc*normc,"Herm C SV5 V");

            c2 = c.conjugate();
            SV_Decompose(c2.conjugate(),cS2.view());
            Assert(Norm(cS2-cS) <= ceps*normc,"Herm C SV6");
            SV_Decompose(c,cU2.conjugate(),cS2.view());
            Assert(Norm(cS2-cS) <= ceps*normc,"Herm C SV4 S");
            if (!nearoverflow)
                Assert(Norm(c*c.adjoint()-
                            cU2.conjugate()*cS2*cS2*cU2.transpose()) <=
                       ceps*normc*normc, "Herm C SV4 U");
#endif
            std::cout<<"."; std::cout.flush();
        }

        // Eigen
        {
            if (showstartdone) std::cout<<"Eigen"<<std::endl;
            tmv::Matrix<T> V(N,N);
            tmv::Vector<T> L(N);
            Eigen(m,V.view(),L.view());
            Assert(Norm(m*V-V*DiagMatrixViewOf(L)) <= eps*normm,"Herm Eigen");

            tmv::Matrix<std::complex<T> > cV(N,N);
            tmv::Vector<T> cL(N);
            Eigen(c,cV.view(),cL.view());
            Assert(Norm(c*cV-cV*DiagMatrixViewOf(cL)) <= eps*normc,
                   "Herm C Eigen");

#if (XTEST & 16)
            tmv::Vector<T> L2(N);
            tmv::HermMatrix<T,uplo,stor> m2 = m;
            Eigen(m2.view(),L2.view());
            Assert(Norm(L2-L) <= eps*normm,"Herm Eigen2");

            tmv::Vector<T> cL2(N);
            tmv::HermMatrix<std::complex<T>,uplo,stor> c2 = c;
            Eigen(c2.view(),cL2.view());
            Assert(Norm(cL2-cL) <= eps*normc,"Herm C Eigen2");

            Eigen(c,cV.conjugate(),cL.view());
            Assert(Norm(c*cV.conjugate()-cV.conjugate()*DiagMatrixViewOf(cL)) <= 
                   eps*normc,"Herm C Eigen3");

            Eigen(c.conjugate(),cV.view(),cL.view());
            Assert(Norm(c.conjugate()*cV-cV*DiagMatrixViewOf(cL)) <= eps*normc,
                   "Herm C Eigen4");

            Eigen(c.conjugate(),cV.conjugate(),cL.view());
            Assert(Norm(c*cV-cV*DiagMatrixViewOf(cL)) <= eps*normc,
                   "Herm C Eigen5");

            c2.conjugate() = c;
            Eigen(c2.conjugate(),cL2.view());
            Assert(Norm(cL2-cL) <= eps*normc,"Herm C Eigen6");
#endif
            std::cout<<"."; std::cout.flush();
        }

        // Square Root
#ifdef NOTHROW
        if (posdef) {
#else
            try {
#endif
                if (showstartdone) std::cout<<"SquareRoot"<<std::endl;
                tmv::HermMatrix<T,uplo,stor> S = m;
                SquareRoot(S.view());
                Assert(Norm(m-S*S) <= eps*normm,"Herm Square Root");

                tmv::HermMatrix<std::complex<T>,uplo,stor> cS = c;
                SquareRoot(cS.view());

#if (XTEST & 16)
                cS.conjugate() = c;
                SquareRoot(cS.conjugate());
                Assert(Norm(c-cS.conjugate()*cS.conjugate()) <= eps*normc,
                       "Herm C Square Root 2");
#endif
                Assert(posdef,"Didn't throw NonPosDef, and mattype != pos def");
                std::cout<<"."; std::cout.flush();
#ifndef NOTHROW
            } catch(tmv::NonPosDef) { 
                Assert(!posdef,"Caught NonPosDef, but mattype == pos def"); 
            }
#else
        }
#endif
    }
}

template <class T, tmv::UpLoType uplo, tmv::StorageType stor> 
void TestSymDecomp()
{
    for (int mattype = START; mattype <= 5; mattype++) {
#if !(XTEST & 64)
        //if (mattype >= 3) break;
#endif
#ifdef LAP
        if (mattype >= 4) break;
#endif
        if (showstartdone) {
            std::cout<<"Sym: mattype = "<<mattype<<std::endl;
            std::cout<<"uplo, stor = "<<TMV_Text(uplo)<<
                "  "<<TMV_Text(stor)<<std::endl;
        }
        // mattype = 0  is Normal (Posdef if real)
        // mattype = 1  is Indef
        // mattype = 2  is Singular
        // mattype = 3  is Singular with seriously bad defects
        // mattype = 4  is Singular and nearly zero
        // mattype = 5  is Singular and nearly overflow

        const int N = 200;
        const bool posdef = mattype == 0;
        const bool baddefect = mattype == 3;
        const bool nearunderflow = mattype == 4;
        const bool nearoverflow = mattype == 5;
        // Note: posdef isn't really positive definite for complex matrices.
        const bool singular = mattype >= 2;
        if (showstartdone) {
            std::cout<<"posdef = "<<posdef<<
                ", singular = "<<singular<<std::endl;
        }

        tmv::SymMatrix<T,uplo,stor> m(N);
        for(int i=0;i<N;++i) for(int j=0;j<N;++j) 
            if ((uplo == tmv::Upper && j>=i) || (uplo == tmv::Lower && i>=j)) 
                m(i,j) = T(2+4*i-5*j);
        if (posdef) {
            m /= T(10*N);
            m.diag().addToAll(T(1));
            for(int i=0;i<N;i++) m.diag()(i) += T(i);
            m.diag(0,N/4,N/2) *= T(10);
            m.subSymMatrix(N/2,3*N/4,2) *= T(100);
        } else if (singular) {
            m.row(2,0,2).setZero();
            m.row(2,2,N).setZero();
            m.row(4,0,4) = m.row(5,0,4);
            m.row(4,5,N) = m.row(5,5,N);
            m(4,5) = m(4,4) = m(5,5);
        } else {
            // indef
            m /= T(N);
            m.diag(0,N/10,N/5).setZero();
            m.col(0,N/3,3*N/7).addToAll(T(50));
            m.diag(1,4,7).addToAll(T(10));
            m.row(7,0,5).addToAll(T(10));
            m.diag(0,N/4,N/2) *= T(10);
            m.subSymMatrix(N/2,3*N/4,2) /= T(100);
        } 

        tmv::SymMatrix<std::complex<T>,uplo,stor> c(N);
        for(int i=0;i<N;++i) for(int j=0;j<N;++j) 
            c(i,j) = std::complex<T>(2+4*i-5*j,3-i);
        if (posdef) {
            c /= T(10*N);
            c.diag().addToAll(T(1));
            for(int i=0;i<N;i++) c.diag()(i) += T(i);
            c.diag(0,N/4,N/2) *= T(10);
            c.subSymMatrix(N/2,3*N/4,2) *= T(100);
        } else if (singular) {
            c.row(2,0,2).setZero();
            c.row(2,2,N).setZero();
            c.row(4,0,4) = c.row(5,0,4);
            c.row(4,5,N) = c.row(5,5,N);
            c(4,5) = c(4,4) = c(5,5);
        } else {
            // indef
            c /= T(N);
            c.diag(0,N/10,N/5).setZero();
            c.col(0,N/3,3*N/7).addToAll(T(50));
            c.diag(1,4,7).addToAll(T(10));
            c.row(7,0,5).addToAll(T(10));
            c.diag(0,N/4,N/2) *= T(10);
            c.subSymMatrix(N/2,3*N/4,2) /= T(100);
        }

        if (baddefect) {
            for(int i=10;i<N;i+=10) {
                m.subMatrix(0,i,i,N) *= T(1.e-10);
                m.subSymMatrix(i,N) *= T(-1.e-10);
                c.subMatrix(0,i,i,N) *= T(1.e-10);
                c.subSymMatrix(i,N) *= T(-1.e-10);
            }
        }

        if (nearunderflow) {
            T x = std::numeric_limits<T>::min();
            m.setAllTo(x);
            c.upperTri().setAllTo(std::complex<T>(x,2*x));
            for(int i=1;i<N;++i) {
                m.col(i,0,i+1) /= T(i);
                c.col(i,0,i+1) /= T(i);
            }
            m(N/2,N/2) = x/std::numeric_limits<T>::epsilon();
            c(N/2,N/2) = x/std::numeric_limits<T>::epsilon();
            m.diag(0,N/5,2*N/5).setZero();
            c.diag(0,N/5,2*N/5).setZero();
            m.diag(0,N/2+1,N) *= T(-1);
            c.diag(0,N/2+1,N) *= T(-1);
            m.subSymMatrix(N-4,N).setAllTo(
                -x/std::numeric_limits<T>::epsilon()/T(8));
            c.subSymMatrix(N-4,N).setAllTo(
                -x/std::numeric_limits<T>::epsilon()/T(8));
        }

        if (nearoverflow) {
            T x = std::numeric_limits<T>::max();
            x /= N;
            m.setAllTo(x);
            c.upperTri().setAllTo(std::complex<T>(x,x/2));
            for(int i=1;i<N;++i) {
                m.col(i,0,i+1) /= T(i);
                c.col(i,0,i+1) /= T(i);
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
            T ckappa = normc* Norm(c.inverse());
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
        do {
#ifdef LAP
            if (baddefect || nearunderflow || nearoverflow) break;
#endif
            if (showstartdone) std::cout<<"LDL"<<std::endl;
            tmv::LowerTriMatrix<T,tmv::UnitDiag> L = m.lud().getL();
            tmv::BandMatrix<T> D = m.lud().getD();
            tmv::Permutation P = m.lud().getP();
            tmv::Matrix<T> LDL = P*L*D*L.transpose()*P.transpose();
            if (showacc) {
                std::cout<<"Norm(m-LDL) = "<<Norm(m-LDL)<<std::endl;
                std::cout<<"eps*Norm(m) = "<<eps*normm<<std::endl;
            }
            Assert(Norm(m-LDL) <= eps*normm,"Sym LDL");

            tmv::LowerTriMatrix<std::complex<T>,tmv::UnitDiag> cL = 
                c.lud().getL();
            tmv::BandMatrix<std::complex<T> > cD = c.lud().getD();
            tmv::Permutation cP = c.lud().getP();
            tmv::Matrix<std::complex<T> > cLDL = 
                cP*cL*cD*cL.transpose()*cP.transpose();
            if (showacc) {
                std::cout<<"Norm(c-cLDL) = "<<Norm(c-cLDL)<<std::endl;
                std::cout<<"ceps*Norm(c) = "<<ceps*normc<<std::endl;
            }
            Assert(Norm(c-cLDL) <= ceps*normc,"Sym C LDL");

#if (XTEST & 16)
            tmv::SymMatrix<T,uplo,stor> m2 = m;
            tmv::SymBandMatrix<T,uplo,stor> D2(N,1);
            tmv::Permutation P2(N);
            LDL_Decompose(m2.view(),D2.view(),P2);
            L = m2.lowerTri(tmv::UnitDiag);
            LDL = P2*L*D2*L.transpose()*P2.transpose();
            Assert(Norm(m-LDL) <= eps*normm,"Sym LDL2");

            tmv::SymMatrix<std::complex<T>,uplo,stor> c2 = c;
            tmv::SymBandMatrix<std::complex<T>,uplo,stor> cD2(N,1);
            LDL_Decompose(c2.view(),cD2.view(),P2);
            cL = c2.lowerTri(tmv::UnitDiag);
            cLDL = P2*cL*cD2*cL.transpose()*P2.transpose();
            Assert(Norm(c-cLDL) <= ceps*normc,"Sym C LDL2");

            c2.conjugate() = c;
            LDL_Decompose(c2.conjugate(),cD2.view(),P2);
            cL = c2.conjugate().lowerTri(tmv::UnitDiag);
            cLDL = P2*cL*cD2*cL.transpose()*P2.transpose();
            Assert(Norm(c-cLDL) <= ceps*normc,"Sym C LDL3");

            c2= c;
            LDL_Decompose(c2.view(),cD2.conjugate(),P2);
            cL = c2.lowerTri(tmv::UnitDiag);
            cLDL = P2*cL*cD2.conjugate()*cL.transpose()*P2.transpose();
            Assert(Norm(c-cLDL) <= ceps*normc,"Sym C LDL4");

            c2.conjugate() = c;
            LDL_Decompose(c2.conjugate(),cD2.conjugate(),P2);
            cL = c2.conjugate().lowerTri(tmv::UnitDiag);
            cLDL = P2*cL*cD2.conjugate()*cL.transpose()*P2.transpose();
            Assert(Norm(c-cLDL) <= ceps*normc,"Sym C LDL5");
#endif
            std::cout<<"."; std::cout.flush();
        } while (false);

        // SV Decomposition
        {
            if (showstartdone) std::cout<<"SV"<<std::endl;
            tmv::Matrix<T> U = m.svd().getU();
            tmv::DiagMatrix<T> S = m.svd().getS();
            tmv::Matrix<T> V = m.svd().getV();
            if (showacc) {
                std::cout<<"Norm(m-U*S*V) = "<<Norm(m-U*S*V)<<std::endl;
                std::cout<<"eps*Norm(m) = "<<eps*normm<<std::endl;
            }
            Assert(Norm(m-U*S*V) <= eps*normm,"Sym SV");

            tmv::Matrix<std::complex<T> > cU = c.symsvd().getU();
            tmv::DiagMatrix<T> cS = c.symsvd().getS();
            tmv::Matrix<std::complex<T> > cV = c.symsvd().getV();
            if (showacc) {
                std::cout<<"Norm(c-U*S*V) = "<<Norm(c-cU*cS*cV)<<std::endl;
                std::cout<<"eps*Norm(c) = "<<eps*normc<<std::endl;
            }
            Assert(Norm(c-cU*cS*cV) <= ceps*normc,"Sym C SV");

#if (XTEST & 16)
            tmv::Matrix<T> U2(N,N);
            tmv::DiagMatrix<T> S2(N);
            tmv::Matrix<T> V2(N,N);
            SV_Decompose(m,U2.view(),S2.view(),V2.view());
            Assert(Norm(m-U2*S2*V2) <= eps*normm,"Sym SV2");

            tmv::SymMatrix<T,uplo,stor> m2 = m;
            SV_Decompose(m2.view(),S2.view());
            Assert(Norm(S2-S) <= eps*normm,"Sym SV3");
            SV_Decompose(m,U2.view(),S2.view());
            Assert(Norm(S2-S) <= eps*normm,"Sym SV4 S");
            if (!nearoverflow)
                Assert(Norm(m*m.adjoint()-U2*S2*S2*U2.adjoint()) <=
                       eps*normm*normm, "Sym SV4 U");
            SV_Decompose(m,S2.view(),V2.view());
            Assert(Norm(S2-S) <= eps*normm,"Sym SV5 S");
            if (!nearoverflow)
                Assert(Norm(m.adjoint()*m-V2.adjoint()*S2*S2*V2) <=
                       eps*normm*normm, "Sym SV5 V");

            tmv::Matrix<std::complex<T> > cU2(N,N);
            tmv::DiagMatrix<T> cS2(N);
            tmv::Matrix<std::complex<T> > cV2(N,N);
            SV_Decompose(c,cU2.view(),cS2.view(),cV2.view());
            Assert(Norm(c-cU2*cS2*cV2) <= ceps*normc,"Sym C SV2");

            tmv::SymMatrix<std::complex<T>,uplo,stor> c2 = c;
            SV_Decompose(c2.view(),cS2.view());
            Assert(Norm(cS2-cS) <= ceps*normc,"Sym C SV3");
            SV_Decompose(c,cU2.view(),cS2.view());
            Assert(Norm(cS2-cS) <= ceps*normc,"Sym C SV4 S");
            if (!nearoverflow)
                Assert(Norm(c*c.adjoint()-cU2*cS2*cS2*cU2.adjoint()) <= 
                       ceps*normc*normc,"Sym C SV4 U");
            SV_Decompose(c,cS2.view(),cV2.view());
            Assert(Norm(cS2-cS) <= ceps*normc,"Sym C SV5 S");
            if (!nearoverflow)
                Assert(Norm(c.adjoint()*c-cV2.adjoint()*cS2*cS2*cV2) <= 
                       ceps*normc*normc,"Sym C SV5 V");

            c2 = c.conjugate();
            SV_Decompose(c2.conjugate(),cS2.view());
            Assert(Norm(cS2-cS) <= ceps*normc,"Sym C SV6");
            SV_Decompose(c,cU2.conjugate(),cS2.view());
            Assert(Norm(cS2-cS) <= ceps*normc,"Sym C SV7 S");
            if (!nearoverflow)
                Assert(Norm(c*c.adjoint()-
                            cU2.conjugate()*cS2*cS2*cU2.transpose()) <= 
                       ceps*normc*normc,"Sym C SV7 U");
#endif
            std::cout<<"."; std::cout.flush();
        }
    }
}

template <class T, tmv::StorageType stor> 
void TestPolar()
{
    if (showstartdone) std::cout<<"PolarDecomp "<<TMV_Text(stor)<<std::endl;

    for (int mattype = START; mattype <= 6; mattype++) {
#if !(XTEST & 64) || defined(LAP)
        // Some LAPACK packages don't manage the bad defect case here
        // so skip that as well, not just over/underflow.
        //if (mattype >= 4) break;
#endif
        if (showstartdone) {
            std::cout<<"Polar: mattype = "<<mattype<<std::endl;
            std::cout<<"stor = "<<TMV_Text(stor)<<std::endl;
        }
        // mattype = 0  is Square
        // mattype = 1  is NonSquare slightly tall
        // mattype = 2  is NonSquare very tall
        // mattype = 3  is Singular
        // mattype = 4  is Singular with seriously bad defects
        // mattype = 5  is Singular and nearly zero
        // mattype = 6  is Singular and nearly overflow

        const int N = 200;
        int M = N;
        if (mattype == 1) M = 211;
        else if (mattype == 2) M = 545;
        const bool singular = mattype >= 3;
        const bool baddefect = mattype == 4;
        const bool nearunderflow = mattype == 5;
        const bool nearoverflow = mattype == 6;

        tmv::Matrix<T,stor> m(M,N);
        for(int i=0;i<M;++i) for(int j=0;j<N;++j) m(i,j) = T(2+4*i-5*j);
        m(0,0) = T(14);
        m(1,0) = T(-2);
        m(2,0) = T(7);
        m(3,4) = T(-10);
        if (!singular) m.diag() *= T(30);
        else { m.col(1).setZero(); m.row(7) = m.row(6); }

        tmv::Matrix<std::complex<T>,stor> c(M,N);
        for(int i=0;i<M;++i) for(int j=0;j<N;++j) 
            c(i,j) = std::complex<T>(2+4*i-5*j,3-i);
        c(0,0) *= T(14);
        c(1,0) *= T(-2);
        c(2,0) *= T(7);
        c(7,6) *= T(-10);
        if (!singular) c.diag() *= T(30);
        else { c.col(1).setZero(); c.row(7) = c.row(6); }

        if (baddefect) {
            for(int i=10;i<N;i+=10) {
                m.colRange(i,N) *= T(1.e-10);
                c.colRange(i,N) *= T(1.e-10);
                m.rowRange(i+5,N) *= T(1.e-10);
                c.rowRange(i+5,N) *= T(1.e-10);
            }
        }

        if (nearunderflow) {
            T x = std::numeric_limits<T>::min();
            m.setAllTo(x);
            c.setAllTo(std::complex<T>(x,2*x));
            for(int i=1;i<N;++i) {
                m.col(i) /= T(i);
                c.col(i) /= T(i);
            }
            m(N/2,N/2) = x/std::numeric_limits<T>::epsilon();
            c(N/2,N/2) = x/std::numeric_limits<T>::epsilon();
            m.diag(0,N/5,2*N/5).setZero();
            c.diag(0,N/5,2*N/5).setZero();
            m.diag(0,N/2+1,N) *= T(-1);
            c.diag(0,N/2+1,N) *= T(-1);
            m.subMatrix(N-4,N,N-4,N).setAllTo(
                -x/std::numeric_limits<T>::epsilon()/T(8));
            c.subMatrix(N-4,N,N-4,N).setAllTo(
                -x/std::numeric_limits<T>::epsilon()/T(8));
        }

        if (nearoverflow) {
            T x = std::numeric_limits<T>::max();
            x /= N;
            m.setAllTo(x);
            c.setAllTo(std::complex<T>(x,x/2));
            for(int i=1;i<N;++i) {
                m.col(i) /= T(i);
                c.col(i) /= T(i);
            }
        }

        T eps = EPS;
        T ceps = EPS;
        T normm = Norm(m);
        T normc = Norm(c);
        if (showacc) {
            std::cout<<"eps = "<<eps<<"  "<<ceps<<std::endl;
            std::cout<<"norm = "<<normm<<"  "<<normc<<std::endl;
        }
        if (!singular) {
            T kappa = normm * Norm(m.inverse());
            T ckappa = normc* Norm(c.inverse());
            if (showacc) std::cout<<"kappa = "<<kappa<<"  "<<ckappa<<std::endl;
            eps *= kappa;
            ceps *= ckappa;
        } else {
            if (showacc) std::cout<<"eps *= "<<T(10*N)<<std::endl;
            eps *= T(10*N);
            ceps *= T(10*N);
        }
        if (showacc) std::cout<<"eps => "<<eps<<"  "<<ceps<<std::endl;

        // Matrix Polar Decomposition
        do {
            if (showstartdone) std::cout<<"Matrix Polar"<<std::endl;
            tmv::HermMatrix<T,tmv::Lower,tmv::ColMajor> P(N);
            tmv::Matrix<T,stor> U = m;
            PolarDecompose(U.view(),P.view());
            if (showacc) {
                std::cout<<"Norm(m-UP) = "<<Norm(m-U*P)<<std::endl;
                std::cout<<"eps*Norm(m) = "<<eps*normm<<std::endl;
                std::cout<<"Norm(UtU-1) = "<<Norm(U.adjoint()*U-T(1))<<std::endl;
                std::cout<<"eps = "<<eps<<std::endl;
            }
            Assert(Norm(m-U*P) <= eps*normm,"Polar");
            Assert(Norm(U.adjoint()*U-T(1)) <= eps,"Polar UtU");

            tmv::HermMatrix<std::complex<T>,tmv::Lower,tmv::ColMajor> cP(N);
            tmv::Matrix<std::complex<T>,stor> cU = c;
            PolarDecompose(cU.view(),cP.view());
            Assert(Norm(c-cU*cP) <= ceps*normc,"C Polar");
            Assert(Norm(cU.adjoint()*cU-T(1)) <= ceps,"C Polar UtU");

#if (XTEST & 16)
            U = m;
            PolarDecompose(U.view(),P.transpose());
            Assert(Norm(m-U*P) <= eps*normm,"Polar2");
            Assert(Norm(U.adjoint()*U-T(1)) <= eps,"Polar2 UtU");

            cU = c;
            PolarDecompose(cU.view(),cP.adjoint());
            Assert(Norm(c-cU*cP) <= ceps*normc,"C Polar2");
            Assert(Norm(cU.adjoint()*cU-T(1)) <= ceps,"C Polar2 UtU");

            cU = c;
            PolarDecompose(cU.view(),cP.transpose());
            Assert(Norm(c-cU*cP.transpose()) <= ceps*normc,"C Polar3");
            Assert(Norm(cU.adjoint()*cU-T(1)) <= ceps,"C Polar3 UtU");

            cU = c;
            PolarDecompose(cU.view(),cP.conjugate());
            Assert(Norm(c-cU*cP.transpose()) <= ceps*normc,"C Polar4");
            Assert(Norm(cU.adjoint()*cU-T(1)) <= ceps,"C Polar4 UtU");

            cU.conjugate() = c;
            PolarDecompose(cU.conjugate(),cP.view());
            Assert(Norm(c-cU.conjugate()*cP) <= ceps*normc,"C Polar5");
            Assert(Norm(cU.adjoint()*cU-T(1)) <= ceps,"C Polar5 UtU");

            cU.conjugate() = c;
            PolarDecompose(cU.conjugate(),cP.adjoint());
            Assert(Norm(c-cU.conjugate()*cP) <= ceps*normc,"C Polar6");
            Assert(Norm(cU.adjoint()*cU-T(1)) <= ceps,"C Polar6 UtU");

            cU.conjugate() = c;
            PolarDecompose(cU.conjugate(),cP.transpose());
            Assert(Norm(c-cU.conjugate()*cP.transpose()) <= ceps*normc,
                   "C Polar7");
            Assert(Norm(cU.adjoint()*cU-T(1)) <= ceps,"C Polar7 UtU");

            cU.conjugate() = c;
            PolarDecompose(cU.conjugate(),cP.conjugate());
            Assert(Norm(c-cU.conjugate()*cP.transpose()) <= ceps*normc,
                   "C Polar8");
            Assert(Norm(cU.adjoint()*cU-T(1)) <= ceps,"C Polar8 UtU");
#endif
            std::cout<<"."; std::cout.flush();
        } while (false);

        // BandMatrix Polar Decomposition
        do {
            if (showstartdone) std::cout<<"Band Polar"<<std::endl;
            tmv::BandMatrixView<T> b(m.view(),5,11);
            tmv::BandMatrixView<std::complex<T> > cb(c.view(),5,11);
            T normb = Norm(b);
            T normcb = Norm(cb);

            tmv::HermMatrix<T,tmv::Lower,tmv::ColMajor> P(N);
            tmv::Matrix<T,stor> U(b.colsize(),N);
            PolarDecompose(b,U.view(),P.view());
            if (showacc) {
                std::cout<<"Norm(b-UP) = "<<Norm(b-U*P)<<"  "<<eps*normb<<std::endl;
                std::cout<<"Norm(UtU-1) = "<<Norm(U.adjoint()*U-T(1))<<"  "<<eps<<std::endl;
            }
            Assert(Norm(b-U*P) <= eps*normb,"Band Polar");
            Assert(Norm(U.adjoint()*U-T(1)) <= eps,"Band Polar UtU");

            tmv::HermMatrix<std::complex<T>,tmv::Lower,tmv::ColMajor> cP(N);
            tmv::Matrix<std::complex<T>,stor> cU(cb.colsize(),N);
            PolarDecompose(cb,cU.view(),cP.view());
            if (showacc) {
                std::cout<<"Norm(cb-UP) = "<<Norm(cb-cU*cP)<<"  "<<ceps*normcb<<std::endl;
                std::cout<<"Norm(UtU-1) = "<<Norm(cU.adjoint()*cU-T(1))<<"  "<<ceps<<std::endl;
            }
            Assert(Norm(cb-cU*cP) <= ceps*normcb,"C Band Polar");
            Assert(Norm(cU.adjoint()*cU-T(1)) <= ceps,"C Band Polar UtU");

#if (XTEST & 16)
            PolarDecompose(b,U.view(),P.transpose());
            Assert(Norm(b-U*P) <= eps*normb,"Band Polar2");
            Assert(Norm(U.adjoint()*U-T(1)) <= eps,"Band Polar2 UtU");

            PolarDecompose(cb,cU.view(),cP.adjoint());
            Assert(Norm(cb-cU*cP) <= ceps*normcb,"C Band Polar2");
            Assert(Norm(cU.adjoint()*cU-T(1)) <= ceps,"C Band Polar2 UtU");

            PolarDecompose(cb,cU.view(),cP.transpose());
            Assert(Norm(cb-cU*cP.transpose()) <= ceps*normcb,"C Band Polar3");
            Assert(Norm(cU.adjoint()*cU-T(1)) <= ceps,"C Band Polar3 UtU");

            PolarDecompose(cb,cU.view(),cP.conjugate());
            Assert(Norm(cb-cU*cP.transpose()) <= ceps*normcb,"C Band Polar4");
            Assert(Norm(cU.adjoint()*cU-T(1)) <= ceps,"C Band Polar4 UtU");

            PolarDecompose(cb,cU.conjugate(),cP.view());
            Assert(Norm(cb-cU.conjugate()*cP) <= ceps*normcb,"C Band Polar5");
            Assert(Norm(cU.adjoint()*cU-T(1)) <= ceps,"C Band Polar5 UtU");

            PolarDecompose(cb,cU.conjugate(),cP.adjoint());
            Assert(Norm(cb-cU.conjugate()*cP) <= ceps*normcb,"C Band Polar6");
            Assert(Norm(cU.adjoint()*cU-T(1)) <= ceps,"C Band Polar6 UtU");

            PolarDecompose(cb,cU.conjugate(),cP.transpose());
            Assert(Norm(cb-cU.conjugate()*cP.transpose()) <= ceps*normcb,
                   "C Band Polar7");
            Assert(Norm(cU.adjoint()*cU-T(1)) <= ceps,"C Band Polar7 UtU");

            PolarDecompose(cb,cU.conjugate(),cP.conjugate());
            Assert(Norm(cb-cU.conjugate()*cP.transpose()) <= ceps*normcb,
                   "C Band Polar8");
            Assert(Norm(cU.adjoint()*cU-T(1)) <= ceps,"C Band Polar8 UtU");
#endif
            std::cout<<"."; std::cout.flush();
        } while (false);
    }
}

#ifdef TEST_DOUBLE
template void TestHermDecomp<double,tmv::Upper,tmv::ColMajor>();
template void TestHermDecomp<double,tmv::Upper,tmv::RowMajor>();
template void TestHermDecomp<double,tmv::Lower,tmv::ColMajor>();
template void TestHermDecomp<double,tmv::Lower,tmv::RowMajor>();
template void TestSymDecomp<double,tmv::Upper,tmv::ColMajor>();
template void TestSymDecomp<double,tmv::Upper,tmv::RowMajor>();
template void TestSymDecomp<double,tmv::Lower,tmv::ColMajor>();
template void TestSymDecomp<double,tmv::Lower,tmv::RowMajor>();
template void TestPolar<double,tmv::ColMajor>();
template void TestPolar<double,tmv::RowMajor>();
#endif
#ifdef TEST_FLOAT
template void TestHermDecomp<float,tmv::Upper,tmv::ColMajor>();
template void TestHermDecomp<float,tmv::Upper,tmv::RowMajor>();
template void TestHermDecomp<float,tmv::Lower,tmv::ColMajor>();
template void TestHermDecomp<float,tmv::Lower,tmv::RowMajor>();
template void TestSymDecomp<float,tmv::Upper,tmv::ColMajor>();
template void TestSymDecomp<float,tmv::Upper,tmv::RowMajor>();
template void TestSymDecomp<float,tmv::Lower,tmv::ColMajor>();
template void TestSymDecomp<float,tmv::Lower,tmv::RowMajor>();
template void TestPolar<float,tmv::ColMajor>();
template void TestPolar<float,tmv::RowMajor>();
#endif
#ifdef TEST_LONGDOUBLE
template void TestHermDecomp<long double,tmv::Upper,tmv::ColMajor>();
template void TestHermDecomp<long double,tmv::Upper,tmv::RowMajor>();
template void TestHermDecomp<long double,tmv::Lower,tmv::ColMajor>();
template void TestHermDecomp<long double,tmv::Lower,tmv::RowMajor>();
template void TestSymDecomp<long double,tmv::Upper,tmv::ColMajor>();
template void TestSymDecomp<long double,tmv::Upper,tmv::RowMajor>();
template void TestSymDecomp<long double,tmv::Lower,tmv::ColMajor>();
template void TestSymDecomp<long double,tmv::Lower,tmv::RowMajor>();
template void TestPolar<long double,tmv::ColMajor>();
template void TestPolar<long double,tmv::RowMajor>();
#endif



