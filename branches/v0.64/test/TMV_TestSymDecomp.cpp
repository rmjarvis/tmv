
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
    for (int mattype = 0; mattype < 4; mattype++) {
        if (showstartdone) {
            std::cout<<"Herm: mattype = "<<mattype<<std::endl;
            std::cout<<"uplo, stor = "<<TMV_Text(uplo)<<
                "  "<<TMV_Text(stor)<<std::endl;
        }
        // mattype = 0  is PosDef
        // mattype = 1  is Indef
        // mattype = 2  is Singular
        // mattype = 3  is Singular with seriously bad defects

        const int N = 200; // Must be multiple of 8
        const bool posdef = mattype == 0;
        const bool singular = mattype == 2 || mattype == 3;
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

        if (mattype == 3) {
            for(int i=10;i<N;i+=10) {
                m.subMatrix(0,i,i,N) *= T(1.e-10);
                m.subSymMatrix(i,N) *= T(1.e-10);
                c.subMatrix(0,i,i,N) *= T(1.e-10);
                c.subSymMatrix(i,N) *= T(1.e-10);
            }
        }

        T eps = EPS;
        T ceps = EPS;
        if (showacc) std::cout<<"eps = "<<eps<<"  "<<ceps<<std::endl;
        if (mattype < 2) {
            T kappa = Norm(m) * Norm(m.inverse());
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
#ifdef NOTHROW
        if (posdef) {
#else
            try  {
#endif
                tmv::LowerTriMatrix<T> L = m.chd().getL();
                tmv::Matrix<T> LLt = L*L.adjoint();
                Assert(Norm(m-LLt) < eps*Norm(m),"Herm CH");

                tmv::LowerTriMatrix<std::complex<T> > cL = c.chd().getL();
                tmv::Matrix<std::complex<T> > cLLt = cL*cL.adjoint();
                Assert(Norm(c-cLLt) < ceps*Norm(c),"Herm C CH");

#ifdef XTEST
                tmv::HermMatrix<T,uplo,stor> m2 = m;
                CH_Decompose(m2.view());
                L = m2.lowerTri();
                LLt = L*L.adjoint();
                Assert(Norm(m-LLt) < eps*Norm(m),"Herm CH2");

                tmv::HermMatrix<std::complex<T>,uplo,stor> c2 = c;
                CH_Decompose(c2.view());
                cL = c2.lowerTri();
                cLLt = cL*cL.adjoint();
                Assert(Norm(c-cLLt) < ceps*Norm(c),"Herm C CH2");

                c2.conjugate() = c;
                CH_Decompose(c2.conjugate());
                cL = c2.conjugate().lowerTri();
                cLLt = cL*cL.adjoint();
                Assert(Norm(c-cLLt) < ceps*Norm(c),"Herm C CH2");
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
            tmv::LowerTriMatrix<T,tmv::UnitDiag> L = m.lud().getL();
            tmv::BandMatrix<T> D = m.lud().getD();
            const int* p = m.lud().getP();
            tmv::Matrix<T> LDL = L*D*L.adjoint();
            LDL.reversePermuteRows(p);
            LDL.reversePermuteCols(p);
            Assert(Norm(m-LDL) < eps*Norm(m),"Herm LDL");

            tmv::LowerTriMatrix<std::complex<T>,tmv::UnitDiag> cL = c.lud().getL();
            tmv::BandMatrix<std::complex<T> > cD = c.lud().getD();
            p = c.lud().getP();
            tmv::Matrix<std::complex<T> > cLDL = cL*cD*cL.adjoint();
            cLDL.reversePermuteRows(p);
            cLDL.reversePermuteCols(p);
            Assert(Norm(c-cLDL) < ceps*Norm(c),"Herm C LDL");

#ifdef XTEST
            tmv::HermMatrix<T,uplo,stor> m2 = m;
            tmv::HermBandMatrix<T,uplo,stor> D2(N,1);
            int p2[N];
            LDL_Decompose(m2.view(),D2.view(),p2);
            L = m2.lowerTri(tmv::UnitDiag);
            LDL = L*D2*L.adjoint();
            LDL.reversePermuteRows(p2);
            LDL.reversePermuteCols(p2);
            Assert(Norm(m-LDL) < eps*Norm(m),"Herm LDL2");

            tmv::HermMatrix<std::complex<T>,uplo,stor> c2 = c;
            tmv::HermBandMatrix<std::complex<T>,uplo,stor> cD2(N,1);
            LDL_Decompose(c2.view(),cD2.view(),p2);
            cL = c2.lowerTri(tmv::UnitDiag);
            cLDL = cL*cD2*cL.adjoint();
            cLDL.reversePermuteRows(p2);
            cLDL.reversePermuteCols(p2);
            Assert(Norm(c-cLDL) < ceps*Norm(c),"Herm C LDL2");

            c2.conjugate() = c;
            LDL_Decompose(c2.conjugate(),cD2.view(),p2);
            cL = c2.conjugate().lowerTri(tmv::UnitDiag);
            cLDL = cL*cD2*cL.adjoint();
            cLDL.reversePermuteRows(p2);
            cLDL.reversePermuteCols(p2);
            Assert(Norm(c-cLDL) < ceps*Norm(c),"Herm C LDL3");

            c2 = c;
            LDL_Decompose(c2.view(),cD2.conjugate(),p2);
            cL = c2.lowerTri(tmv::UnitDiag);
            cLDL = cL*cD2.conjugate()*cL.adjoint();
            cLDL.reversePermuteRows(p2);
            cLDL.reversePermuteCols(p2);
            Assert(Norm(c-cLDL) < ceps*Norm(c),"Herm C LDL4");

            c2.conjugate() = c;
            LDL_Decompose(c2.conjugate(),cD2.conjugate(),p2);
            cL = c2.conjugate().lowerTri(tmv::UnitDiag);
            cLDL = cL*cD2.conjugate()*cL.adjoint();
            cLDL.reversePermuteRows(p2);
            cLDL.reversePermuteCols(p2);
            Assert(Norm(c-cLDL) < ceps*Norm(c),"Herm C LDL5");
#endif
            std::cout<<"."; std::cout.flush();
        } catch (tmv::NonPosDef) {
            // The Lapack version throws whenever mattype is not posdef, 
            // but native algorithm succeeds when mattype is indefinite,
            // or even singular.  
            Assert(!posdef,"caught NonPosDef but mattype == posdef");
        }

        // SV Decomposition
        {
            tmv::Matrix<T> U = m.svd().getU();
            tmv::DiagMatrix<T> S = m.svd().getS();
            tmv::Matrix<T> V = m.svd().getV();
            Assert(Norm(m-U*S*V) < eps*Norm(m),"Herm SV");

            tmv::Matrix<std::complex<T> > cU = c.svd().getU();
            tmv::DiagMatrix<T> cS = c.svd().getS();
            tmv::Matrix<std::complex<T> > cV = c.svd().getV();
            Assert(Norm(c-cU*cS*cV) < ceps*Norm(c),"Herm C SV");

#ifdef XTEST
            tmv::Matrix<T> U2(N,N);
            tmv::DiagMatrix<T> S2(N);
            tmv::Matrix<T> V2(N,N);
            SV_Decompose(m,U2.view(),S2.view(),V2.view());
            Assert(Norm(m-U2*S2*V2) < eps*Norm(m),"Herm SV2");

            tmv::HermMatrix<T,uplo,stor> m2 = m;
            SV_Decompose(m2.view(),S2.view());
            Assert(Norm(S2-S) < eps*Norm(S),"Herm SV3");
            SV_Decompose(m,U2.view(),S2.view());
            Assert(Norm(S2-S) < eps*Norm(S),"Herm SV4 S");
            Assert(Norm(U2*S2*S2*U2.adjoint()-m*m.adjoint()) < 
                   eps*Norm(m)*Norm(m),"Herm SV4 U");
            SV_Decompose(m,S2.view(),V2.view());
            Assert(Norm(S2-S) < eps*Norm(S),"Herm SV5 S");
            Assert(Norm(m.adjoint()*m-V2.adjoint()*S2*S2*V2) < 
                   eps*Norm(m)*Norm(m),"Herm SV5 V");

            tmv::Matrix<std::complex<T> > cU2(N,N);
            tmv::DiagMatrix<T> cS2(N);
            tmv::Matrix<std::complex<T> > cV2(N,N);
            SV_Decompose(c,cU2.view(),cS2.view(),cV2.view());
            Assert(Norm(c-cU2*cS2*cV2) < ceps*Norm(c),"Herm C SV2");

            tmv::HermMatrix<std::complex<T>,uplo,stor> c2 = c;
            SV_Decompose(c2.view(),cS2.view());
            Assert(Norm(cS2-cS) < ceps*Norm(cS),"Herm C SV3");
            SV_Decompose(c,cU2.view(),cS2.view());
            Assert(Norm(cS2-cS) < ceps*Norm(cS),"Herm C SV4 S");
            Assert(Norm(cU2*cS2*cS2*cU2.adjoint()-c*c.adjoint()) < 
                   ceps*Norm(c)*Norm(c),"Herm C SV4 U");
            SV_Decompose(c,cS2.view(),cV2.view());
            Assert(Norm(cS2-cS) < ceps*Norm(cS),"Herm C SV5 S");
            Assert(Norm(c.adjoint()*c-cV2.adjoint()*cS2*cS2*cV2) < 
                   ceps*Norm(c)*Norm(c),"Herm C SV5 V");

            c2 = c.conjugate();
            SV_Decompose(c2.conjugate(),cS2.view());
            Assert(Norm(cS2-cS) < ceps*Norm(cS),"Herm C SV6");
            SV_Decompose(c,cU2.conjugate(),cS2.view());
            Assert(Norm(cS2-cS) < ceps*Norm(cS),"Herm C SV4 S");
            Assert(Norm(cU2.conjugate()*cS2*cS2*cU2.transpose() - 
                        c*c.adjoint()) < ceps*Norm(c)*Norm(c),
                   "Herm C SV4 U");
#endif
            std::cout<<"."; std::cout.flush();
        }

        // Eigen
        {
            tmv::Matrix<T> V(N,N);
            tmv::Vector<T> L(N);
            Eigen(m,V.view(),L.view());
            Assert(Norm(m*V-V*DiagMatrixViewOf(L)) < eps*Norm(m),"Herm Eigen");

            tmv::Matrix<std::complex<T> > cV(N,N);
            tmv::Vector<T> cL(N);
            Eigen(c,cV.view(),cL.view());
            Assert(Norm(c*cV-cV*DiagMatrixViewOf(cL)) < eps*Norm(c),
                   "Herm C Eigen");

#ifdef XTEST
            tmv::Vector<T> L2(N);
            tmv::HermMatrix<T,uplo,stor> m2 = m;
            Eigen(m2.view(),L2.view());
            Assert(Norm(L2-L) < eps*Norm(L),"Herm Eigen2");

            tmv::Vector<T> cL2(N);
            tmv::HermMatrix<std::complex<T>,uplo,stor> c2 = c;
            Eigen(c2.view(),cL2.view());
            Assert(Norm(cL2-cL) < eps*Norm(cL),"Herm C Eigen2");

            Eigen(c,cV.conjugate(),cL.view());
            Assert(Norm(c*cV.conjugate()-cV.conjugate()*DiagMatrixViewOf(cL)) < 
                   eps*Norm(c),"Herm C Eigen3");

            Eigen(c.conjugate(),cV.view(),cL.view());
            Assert(Norm(c.conjugate()*cV-cV*DiagMatrixViewOf(cL)) < eps*Norm(c),
                   "Herm C Eigen4");

            Eigen(c.conjugate(),cV.conjugate(),cL.view());
            Assert(Norm(c*cV-cV*DiagMatrixViewOf(cL)) < eps*Norm(c),
                   "Herm C Eigen5");

            c2.conjugate() = c;
            Eigen(c2.conjugate(),cL2.view());
            Assert(Norm(cL2-cL) < eps*Norm(cL),"Herm C Eigen6");
#endif
            std::cout<<"."; std::cout.flush();
        }

        // Square Root
#ifdef NOTHROW
        if (posdef) {
#else
            try {
#endif
                tmv::HermMatrix<T,uplo,stor> S = m;
                SquareRoot(S.view());
                Assert(Norm(m-S*S) < eps*Norm(m),"Herm Square Root");

                tmv::HermMatrix<std::complex<T>,uplo,stor> cS = c;
                SquareRoot(cS.view());

#ifdef XTEST
                cS.conjugate() = c;
                SquareRoot(cS.conjugate());
                Assert(Norm(c-cS.conjugate()*cS.conjugate()) < eps*Norm(c),
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
    for (int mattype = 0; mattype < 3; mattype++) {
        if (showstartdone) {
            std::cout<<"Symm: mattype = "<<mattype<<std::endl;
            std::cout<<"uplo, stor = "<<TMV_Text(uplo)<<
                "  "<<TMV_Text(stor)<<std::endl;
        }
        // mattype = 0  is Normal
        // mattype = 1  is Singular
        // mattype = 2  is Singular with seriously bad defects

        const int N = 200;
        const bool posdef = mattype == 0;
        // Note: posdef isn't really positive definite for complex matrices.
        const bool singular = mattype == 2 || mattype == 3;
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

        if (mattype == 2) {
            for(int i=10;i<N;i+=10) {
                m.subMatrix(0,i,i,N) *= T(1.e-10);
                m.subSymMatrix(i,N) *= T(1.e-10);
                c.subMatrix(0,i,i,N) *= T(1.e-10);
                c.subSymMatrix(i,N) *= T(1.e-10);
            }
        }

        T eps = EPS;
        T ceps = EPS;
        if (showacc) std::cout<<"eps = "<<eps<<"  "<<ceps<<std::endl;
        if (!singular) {
            T kappa = Norm(m) * Norm(m.inverse());
            if (showacc) std::cout<<"kappa = "<<kappa<<std::endl;
            eps *= kappa;
            ceps *= kappa;
        } else {
            if (showacc) std::cout<<"eps *= "<<T(10*N)<<std::endl;
            eps *= T(10*N);
            ceps *= T(10*N);
        }
        if (showacc) std::cout<<"eps => "<<eps<<"  "<<ceps<<std::endl;

        // LDL Decomposition
        {
            tmv::LowerTriMatrix<T,tmv::UnitDiag> L = m.lud().getL();
            tmv::BandMatrix<T> D = m.lud().getD();
            const int* p = m.lud().getP();
            tmv::Matrix<T> LDL = L*D*L.transpose();
            LDL.reversePermuteRows(p);
            LDL.reversePermuteCols(p);
            Assert(Norm(m-LDL) < eps*Norm(m),"Sym LDL");

            tmv::LowerTriMatrix<std::complex<T>,tmv::UnitDiag> cL = c.lud().getL();
            tmv::BandMatrix<std::complex<T> > cD = c.lud().getD();
            p = c.lud().getP();
            tmv::Matrix<std::complex<T> > cLDL = cL*cD*cL.transpose();
            cLDL.reversePermuteRows(p);
            cLDL.reversePermuteCols(p);
            Assert(Norm(c-cLDL) < ceps*Norm(c),"Sym C LDL");

#ifdef XTEST
            tmv::SymMatrix<T,uplo,stor> m2 = m;
            tmv::SymBandMatrix<T,uplo,stor> D2(N,1);
            int p2[N];
            LDL_Decompose(m2.view(),D2.view(),p2);
            L = m2.lowerTri(tmv::UnitDiag);
            LDL = L*D2*L.transpose();
            LDL.reversePermuteRows(p2);
            LDL.reversePermuteCols(p2);
            Assert(Norm(m-LDL) < eps*Norm(m),"Sym LDL2");

            tmv::SymMatrix<std::complex<T>,uplo,stor> c2 = c;
            tmv::SymBandMatrix<std::complex<T>,uplo,stor> cD2(N,1);
            LDL_Decompose(c2.view(),cD2.view(),p2);
            cL = c2.lowerTri(tmv::UnitDiag);
            cLDL = cL*cD2*cL.transpose();
            cLDL.reversePermuteRows(p2);
            cLDL.reversePermuteCols(p2);
            Assert(Norm(c-cLDL) < ceps*Norm(c),"Sym C LDL2");

            c2.conjugate() = c;
            LDL_Decompose(c2.conjugate(),cD2.view(),p2);
            cL = c2.conjugate().lowerTri(tmv::UnitDiag);
            cLDL = cL*cD2*cL.transpose();
            cLDL.reversePermuteRows(p2);
            cLDL.reversePermuteCols(p2);
            Assert(Norm(c-cLDL) < ceps*Norm(c),"Sym C LDL3");

            c2= c;
            LDL_Decompose(c2.view(),cD2.conjugate(),p2);
            cL = c2.lowerTri(tmv::UnitDiag);
            cLDL = cL*cD2.conjugate()*cL.transpose();
            cLDL.reversePermuteRows(p2);
            cLDL.reversePermuteCols(p2);
            Assert(Norm(c-cLDL) < ceps*Norm(c),"Sym C LDL4");

            c2.conjugate() = c;
            LDL_Decompose(c2.conjugate(),cD2.conjugate(),p2);
            cL = c2.conjugate().lowerTri(tmv::UnitDiag);
            cLDL = cL*cD2.conjugate()*cL.transpose();
            cLDL.reversePermuteRows(p2);
            cLDL.reversePermuteCols(p2);
            Assert(Norm(c-cLDL) < ceps*Norm(c),"Sym C LDL5");
#endif
            std::cout<<"."; std::cout.flush();
        }

        // SV Decomposition
        {
            tmv::Matrix<T> U = m.svd().getU();
            tmv::DiagMatrix<T> S = m.svd().getS();
            tmv::Matrix<T> V = m.svd().getV();
            Assert(Norm(m-U*S*V) < eps*Norm(m),"Sym SV");

            tmv::Matrix<std::complex<T> > cU = c.symsvd().getU();
            tmv::DiagMatrix<T> cS = c.symsvd().getS();
            tmv::Matrix<std::complex<T> > cV = c.symsvd().getV();
            Assert(Norm(c-cU*cS*cV) < ceps*Norm(c),"Sym C SV");

#ifdef XTEST
            tmv::Matrix<T> U2(N,N);
            tmv::DiagMatrix<T> S2(N);
            tmv::Matrix<T> V2(N,N);
            SV_Decompose(m,U2.view(),S2.view(),V2.view());
            Assert(Norm(m-U2*S2*V2) < eps*Norm(m),"Sym SV2");

            tmv::SymMatrix<T,uplo,stor> m2 = m;
            SV_Decompose(m2.view(),S2.view());
            Assert(Norm(S2-S) < eps*Norm(S),"Sym SV3");
            SV_Decompose(m,U2.view(),S2.view());
            Assert(Norm(S2-S) < eps*Norm(S),"Sym SV4 S");
            Assert(Norm(U2-U) < eps*Norm(U),"Sym SV4 U");
            SV_Decompose(m,S2.view(),V2.view());
            Assert(Norm(S2-S) < eps*Norm(S),"Sym SV5 S");
            Assert(Norm(m.adjoint()*m-V2.adjoint()*S2*S2*V2) <
                   eps*Norm(m)*Norm(m), "Sym SV5 V");

            tmv::Matrix<std::complex<T> > cU2(N,N);
            tmv::DiagMatrix<T> cS2(N);
            tmv::Matrix<std::complex<T> > cV2(N,N);
            SV_Decompose(c,cU2.view(),cS2.view(),cV2.view());
            Assert(Norm(c-cU2*cS2*cV2) < ceps*Norm(c),"Sym C SV2");

            tmv::SymMatrix<std::complex<T>,uplo,stor> c2 = c;
            SV_Decompose(c2.view(),cS2.view());
            Assert(Norm(cS2-cS) < ceps*Norm(cS),"Sym C SV3");
            SV_Decompose(c,cU2.view(),cS2.view());
            Assert(Norm(cS2-cS) < ceps*Norm(cS),"Sym C SV4 S");
            Assert(Norm(cU2-cU) < ceps*Norm(cU),"Sym C SV4 U");
            SV_Decompose(c,cS2.view(),cV2.view());
            Assert(Norm(cS2-cS) < ceps*Norm(cS),"Sym C SV5 S");
            Assert(Norm(c.adjoint()*c-cV2.adjoint()*cS2*cS2*cV2) < 
                   ceps*Norm(c)*Norm(c), "Sym C SV5 V");

            c2 = c.conjugate();
            SV_Decompose(c2.conjugate(),cS2.view());
            Assert(Norm(cS2-cS) < ceps*Norm(cS),"Sym C SV6");
            SV_Decompose(c,cU2.conjugate(),cS2.view());
            Assert(Norm(cS2-cS) < ceps*Norm(cS),"Sym C SV7 S");
            Assert(Norm(cU2.conjugate()-cU) < ceps*Norm(cU),"Sym C SV7 U");
#endif
            std::cout<<"."; std::cout.flush();
        }
    }
}

template <class T, tmv::StorageType stor> 
void TestPolar()
{
    if (showstartdone) std::cout<<"PolarDecomp "<<TMV_Text(stor)<<std::endl;

    for (int mattype = 0; mattype < 5; mattype++) {
        if (showstartdone) {
            std::cout<<"Polar: mattype = "<<mattype<<std::endl;
            std::cout<<"stor = "<<TMV_Text(stor)<<std::endl;
        }
        // mattype = 0  is Square
        // mattype = 1  is NonSquare slightly tall
        // mattype = 2  is NonSquare very tall
        // mattype = 3  is Singular
        // mattype = 4  is Singular with seriously bad defects

        const int N = 200;
        int M = N;
        if (mattype == 1) M = 211;
        else if (mattype == 2) M = 545;
        const bool singular = mattype == 3 || mattype == 4;

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

        if (mattype == 4) {
            for(int i=10;i<N;i+=10) {
                m.colRange(i,N) *= T(1.e-10);
                c.colRange(i,N) *= T(1.e-10);
                m.rowRange(i+5,N) *= T(1.e-10);
                c.rowRange(i+5,N) *= T(1.e-10);
            }
        }

        T eps = EPS;
        T ceps = EPS;
        if (showacc) std::cout<<"eps = "<<eps<<"  "<<ceps<<std::endl;
        if (!singular) {
            T kappa = Norm(m) * Norm(m.inverse());
            if (showacc) std::cout<<"kappa = "<<kappa<<std::endl;
            eps *= kappa;
            ceps *= kappa;
        } else {
            if (showacc) std::cout<<"eps *= "<<T(10*N)<<std::endl;
            eps *= T(10*N);
            ceps *= T(10*N);
        }
        if (showacc) std::cout<<"eps => "<<eps<<"  "<<ceps<<std::endl;

        // Matrix Polar Decomposition
        {
            tmv::HermMatrix<T,tmv::Lower,tmv::ColMajor> P(N);
            tmv::Matrix<T,stor> U = m;
            PolarDecompose(U.view(),P.view());
            Assert(Norm(m-U*P) < eps*Norm(m),"Polar");
            Assert(Norm(U.adjoint()*U-T(1)) < eps,"Polar UtU");

            tmv::HermMatrix<std::complex<T>,tmv::Lower,tmv::ColMajor> cP(N);
            tmv::Matrix<std::complex<T>,stor> cU = c;
            PolarDecompose(cU.view(),cP.view());
            Assert(Norm(c-cU*cP) < ceps*Norm(c),"C Polar");
            Assert(Norm(cU.adjoint()*cU-T(1)) < ceps,"C Polar UtU");

#ifdef XTEST
            U = m;
            PolarDecompose(U.view(),P.transpose());
            Assert(Norm(m-U*P) < eps*Norm(m),"Polar2");
            Assert(Norm(U.adjoint()*U-T(1)) < eps,"Polar2 UtU");

            cU = c;
            PolarDecompose(cU.view(),cP.adjoint());
            Assert(Norm(c-cU*cP) < ceps*Norm(c),"C Polar2");
            Assert(Norm(cU.adjoint()*cU-T(1)) < ceps,"C Polar2 UtU");

            cU = c;
            PolarDecompose(cU.view(),cP.transpose());
            Assert(Norm(c-cU*cP.transpose()) < ceps*Norm(c),"C Polar3");
            Assert(Norm(cU.adjoint()*cU-T(1)) < ceps,"C Polar3 UtU");

            cU = c;
            PolarDecompose(cU.view(),cP.conjugate());
            Assert(Norm(c-cU*cP.transpose()) < ceps*Norm(c),"C Polar4");
            Assert(Norm(cU.adjoint()*cU-T(1)) < ceps,"C Polar4 UtU");

            cU.conjugate() = c;
            PolarDecompose(cU.conjugate(),cP.view());
            Assert(Norm(c-cU.conjugate()*cP) < ceps*Norm(c),"C Polar5");
            Assert(Norm(cU.adjoint()*cU-T(1)) < ceps,"C Polar5 UtU");

            cU.conjugate() = c;
            PolarDecompose(cU.conjugate(),cP.adjoint());
            Assert(Norm(c-cU.conjugate()*cP) < ceps*Norm(c),"C Polar6");
            Assert(Norm(cU.adjoint()*cU-T(1)) < ceps,"C Polar6 UtU");

            cU.conjugate() = c;
            PolarDecompose(cU.conjugate(),cP.transpose());
            Assert(Norm(c-cU.conjugate()*cP.transpose()) < ceps*Norm(c),
                   "C Polar7");
            Assert(Norm(cU.adjoint()*cU-T(1)) < ceps,"C Polar7 UtU");

            cU.conjugate() = c;
            PolarDecompose(cU.conjugate(),cP.conjugate());
            Assert(Norm(c-cU.conjugate()*cP.transpose()) < ceps*Norm(c),
                   "C Polar8");
            Assert(Norm(cU.adjoint()*cU-T(1)) < ceps,"C Polar8 UtU");
#endif
            std::cout<<"."; std::cout.flush();
        }

        // BandMatrix Polar Decomposition
        {
            tmv::BandMatrixView<T> b(m.view(),5,11);
            tmv::BandMatrixView<std::complex<T> > cb(c.view(),5,11);

            tmv::HermMatrix<T,tmv::Lower,tmv::ColMajor> P(N);
            tmv::Matrix<T,stor> U(b.colsize(),N);
            PolarDecompose(b,U.view(),P.view());
            Assert(Norm(b-U*P) < eps*Norm(b),"Band Polar");
            Assert(Norm(U.adjoint()*U-T(1)) < eps,"Band Polar UtU");

            tmv::HermMatrix<std::complex<T>,tmv::Lower,tmv::ColMajor> cP(N);
            tmv::Matrix<std::complex<T>,stor> cU(cb.colsize(),N);
            PolarDecompose(cb,cU.view(),cP.view());
            Assert(Norm(cb-cU*cP) < ceps*Norm(cb),"C Band Polar");
            Assert(Norm(cU.adjoint()*cU-T(1)) < ceps,"C Band Polar UtU");

#ifdef XTEST
            PolarDecompose(b,U.view(),P.transpose());
            Assert(Norm(b-U*P) < eps*Norm(b),"Band Polar2");
            Assert(Norm(U.adjoint()*U-T(1)) < eps,"Band Polar2 UtU");

            PolarDecompose(cb,cU.view(),cP.adjoint());
            Assert(Norm(cb-cU*cP) < ceps*Norm(cb),"C Band Polar2");
            Assert(Norm(cU.adjoint()*cU-T(1)) < ceps,"C Band Polar2 UtU");

            PolarDecompose(cb,cU.view(),cP.transpose());
            Assert(Norm(cb-cU*cP.transpose()) < ceps*Norm(cb),"C Band Polar3");
            Assert(Norm(cU.adjoint()*cU-T(1)) < ceps,"C Band Polar3 UtU");

            PolarDecompose(cb,cU.view(),cP.conjugate());
            Assert(Norm(cb-cU*cP.transpose()) < ceps*Norm(cb),"C Band Polar4");
            Assert(Norm(cU.adjoint()*cU-T(1)) < ceps,"C Band Polar4 UtU");

            PolarDecompose(cb,cU.conjugate(),cP.view());
            Assert(Norm(cb-cU.conjugate()*cP) < ceps*Norm(cb),"C Band Polar5");
            Assert(Norm(cU.adjoint()*cU-T(1)) < ceps,"C Band Polar5 UtU");

            PolarDecompose(cb,cU.conjugate(),cP.adjoint());
            Assert(Norm(cb-cU.conjugate()*cP) < ceps*Norm(cb),"C Band Polar6");
            Assert(Norm(cU.adjoint()*cU-T(1)) < ceps,"C Band Polar6 UtU");

            PolarDecompose(cb,cU.conjugate(),cP.transpose());
            Assert(Norm(cb-cU.conjugate()*cP.transpose()) < ceps*Norm(cb),
                   "C Band Polar7");
            Assert(Norm(cU.adjoint()*cU-T(1)) < ceps,"C Band Polar7 UtU");

            PolarDecompose(cb,cU.conjugate(),cP.conjugate());
            Assert(Norm(cb-cU.conjugate()*cP.transpose()) < ceps*Norm(cb),"C Band Polar8");
            Assert(Norm(cU.adjoint()*cU-T(1)) < ceps,"C Band Polar8 UtU");
#endif
            std::cout<<"."; std::cout.flush();
        }
    }
}

#ifdef INST_DOUBLE
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
#ifdef INST_FLOAT
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
#ifdef INST_LONGDOUBLE
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


