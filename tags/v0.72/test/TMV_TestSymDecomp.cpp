
#define START 0

#include "TMV.h"
#include "TMV_Sym.h"
#include "TMV_Band.h"
#include "TMV_SymBand.h"
#include "TMV_Test.h"
#include "TMV_Test_2.h"
#include "TMV_TestSymArith.h"

#ifdef TMV_MEM_DEBUG
// See the discussion of this in TMV_TestTri.cpp.  But basically, there seems to be something
// in the std library exception class that doesn't interact well with the mmgr-style memory
// debugging.  So we skip those tests if we are doing MEM_DEBUG
#define NOTHROW
#endif

template <class T, tmv::UpLoType uplo, tmv::StorageType stor> 
void TestHermDecomp()
{
    typedef std::complex<T> CT;
    for (int mattype = START; mattype <= 5; mattype++) {
#if !(XTEST & 64) || defined(LAP)
        if (mattype >= 3) break;
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
        const T Teps = std::numeric_limits<T>::epsilon();
        const T Tmin = std::numeric_limits<T>::min();
        const T Tmax = std::numeric_limits<T>::max();
        if (showstartdone) {
            std::cout<<"posdef = "<<posdef<<
                ", singular = "<<singular<<std::endl;
        }

        tmv::HermMatrix<T,uplo|stor> m(N);
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

        tmv::HermMatrix<CT,uplo|stor> c(N);
        for(int i=0;i<N;++i) for(int j=0;j<N;++j) 
            if ((uplo == tmv::Upper && j>=i) || (uplo == tmv::Lower && i>=j)) 
                c(i,j) = CT(2+4*i-5*j,3-i);
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
            T x = Tmin;
            m.setAllTo(x);
            c.upperTri().offDiag().setAllTo(CT(x,2*x));
            c.diag().setAllTo(x);
            for(int i=1;i<N;++i) {
                m.col(i,0,i+1) /= T(i);
                c.col(i,0,i+1) /= T(i);
            }
            m(N/2,N/2) = x/Teps;
            c(N/2,N/2) = x/Teps;
            m.diag(0,N/5,2*N/5).setZero();
            c.diag(0,N/5,2*N/5).setZero();
            m.diag(0,N/2+1,N) *= T(-1);
            c.diag(0,N/2+1,N) *= T(-1);
            m.subSymMatrix(N-4,N).setAllTo(-x/Teps/T(8));
            c.subSymMatrix(N-4,N).setAllTo(-x/Teps/T(8));
        }

        if (nearoverflow) {
            T x = Tmax;
            x /= N;
            m.setAllTo(x);
            c.upperTri().offDiag().setAllTo(CT(x,x/2));
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
            try {
#endif
                if (showstartdone) std::cout<<"CH"<<std::endl;
                tmv::LowerTriMatrix<T> L = m.chd().getL();
                tmv::Matrix<T> LLt = L*L.adjoint();
                Assert(Equal(m,LLt,eps*normm),"Herm CH");

                tmv::LowerTriMatrix<CT> cL = c.chd().getL();
                tmv::Matrix<CT> cLLt = cL*cL.adjoint();
                Assert(Equal(c,cLLt,ceps*normc),"Herm C CH");

#if (XTEST & 16)
                tmv::HermMatrix<T,uplo|stor> m2 = m;
                CH_Decompose(m2);
                L = m2.lowerTri();
                LLt = L*L.adjoint();
                Assert(Equal(m,LLt,eps*normm),"Herm CH2");

                tmv::HermMatrix<CT,uplo|stor> c2 = c;
                CH_Decompose(c2);
                cL = c2.lowerTri();
                cLLt = cL*cL.adjoint();
                Assert(Equal(c,cLLt,ceps*normc),"Herm C CH2");

                c2.conjugate() = c;
                CH_Decompose(c2.conjugate());
                cL = c2.conjugate().lowerTri();
                cLLt = cL*cL.adjoint();
                Assert(Equal(c,cLLt,ceps*normc),"Herm C CH2");
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
#ifdef NOTHROW
        if (!singular) {
#else
            try  {
#endif
                if (showstartdone) std::cout<<"LDL"<<std::endl;
                tmv::LowerTriMatrix<T,tmv::UnitDiag> L = m.lud().getL();
                tmv::BandMatrix<T> D = m.lud().getD();
                tmv::Permutation P = m.lud().getP();
                tmv::Matrix<T> LDL = P*L*D*L.adjoint()*P.transpose();
                Assert(Equal(m,LDL,eps*normm),"Herm LDL");

                tmv::LowerTriMatrix<CT,tmv::UnitDiag> cL = c.lud().getL();
                tmv::BandMatrix<CT> cD = c.lud().getD();
                tmv::Permutation cP = c.lud().getP();
                tmv::Matrix<CT> cLDL = 
                    cP*cL*cD*cL.adjoint()*cP.transpose();
                Assert(Equal(c,cLDL,ceps*normc),"Herm C LDL");

#if (XTEST & 16)
                tmv::HermMatrix<T,uplo|stor> m2 = m;
                tmv::HermBandMatrix<T,uplo|stor> D2(N,1);
                tmv::Permutation P2(N);
                LDL_Decompose(m2,D2,P2);
                L = m2.unitLowerTri();
                LDL = P2*L*D2*L.adjoint()*P2.transpose();
                Assert(Equal(m,LDL,eps*normm),"Herm LDL2");

                tmv::HermMatrix<CT,uplo|stor> c2 = c;
                tmv::HermBandMatrix<CT,uplo|stor> cD2(N,1);
                LDL_Decompose(c2,cD2,P2);
                cL = c2.unitLowerTri();
                cLDL = P2*cL*cD2*cL.adjoint()*P2.transpose();
                Assert(Equal(c,cLDL,ceps*normc),"Herm C LDL2");

                c2.conjugate() = c;
                LDL_Decompose(c2.conjugate(),cD2,P2);
                cL = c2.conjugate().unitLowerTri();
                cLDL = P2*cL*cD2*cL.adjoint()*P2.transpose();
                Assert(Equal(c,cLDL,ceps*normc),"Herm C LDL3");

                c2 = c;
                LDL_Decompose(c2,cD2.conjugate(),P2);
                cL = c2.unitLowerTri();
                cLDL = P2*cL*cD2.conjugate()*cL.adjoint()*P2.transpose();
                Assert(Equal(c,cLDL,ceps*normc),"Herm C LDL4");

                c2.conjugate() = c;
                LDL_Decompose(c2.conjugate(),cD2.conjugate(),P2);
                cL = c2.conjugate().unitLowerTri();
                cLDL = P2*cL*cD2.conjugate()*cL.adjoint()*P2.transpose();
                Assert(Equal(c,cLDL,ceps*normc),"Herm C LDL5");
#endif
                std::cout<<"."; std::cout.flush();
#ifndef NOTHROW
            } catch (tmv::Singular) {
                // The Lapack version throws when matrix is exactly singular,
                // but native algorithm succeeds for all matrices.
                Assert(singular,"caught Singular but mattype != singular");
            }
#else
        }
#endif

        // SV Decomposition
        {
            if (showstartdone) std::cout<<"SV"<<std::endl;
            tmv::Matrix<T> U = m.svd().getU();
            tmv::DiagMatrix<T> S = m.svd().getS();
            tmv::Matrix<T> Vt = m.svd().getVt();
            if (showacc) {
                std::cout<<"Norm(m-U*S*Vt) = "<<Norm(m-U*S*Vt)<<std::endl;
                std::cout<<"eps*Norm(m) = "<<eps*normm<<std::endl;
            }
            Assert(Equal(m,U*S*Vt,eps*normm),"Herm SV");

            tmv::Matrix<CT> cU = c.svd().getU();
            tmv::DiagMatrix<T> cS = c.svd().getS();
            tmv::Matrix<CT> cVt = c.svd().getVt();
            Assert(Equal(c,cU*cS*cVt,ceps*normc),"Herm C SV");

#if (XTEST & 16)
            const T x = 
                nearoverflow ? Tmax*Teps :
                nearunderflow ? Tmin/Teps : T(1);

            tmv::Matrix<T> U2(N,N);
            tmv::DiagMatrix<T> S2(N);
            tmv::Matrix<T> Vt2(N,N);
            SV_Decompose(m,U2,S2,Vt2);
            Assert(Equal(m,U2*S2*Vt2,eps*normm),"Herm SV2");

            tmv::HermMatrix<T,uplo|stor> m2 = m;
            SV_Decompose(m2,S2);
            Assert(Equal(S2,S,eps*normm),"Herm SV3");
            SV_Decompose(m,U2,S2);
            Assert(Equal(S2,S,eps*normm),"Herm SV4 S");
            Assert(Equal(
                    tmv::Matrix<T>(m/x)*tmv::Matrix<T>(m.transpose()/x),
                    U2*tmv::DiagMatrix<T>(S2/x)*
                    tmv::DiagMatrix<T>(S2/x)*U2.adjoint(),
                    eps*(normm/x)*(normm/x)),"Herm SV4 U");
            SV_Decompose(m,S2,Vt2);
            Assert(Equal(S2,S,eps*normm),"Herm SV5 S");
            Assert(Equal(
                    tmv::Matrix<T>(m.adjoint()/x)*tmv::Matrix<T>(m/x),
                    Vt2.adjoint()*tmv::DiagMatrix<T>(S2/x)*
                    tmv::DiagMatrix<T>(S2/x)*Vt2,
                    eps*(normm/x)*(normm/x)),"Herm SV5 Vt");

            tmv::Matrix<CT> cU2(N,N);
            tmv::DiagMatrix<T> cS2(N);
            tmv::Matrix<CT> cVt2(N,N);
            SV_Decompose(c,cU2,cS2,cVt2);
            Assert(Equal(c,cU2*cS2*cVt2,ceps*normc),"Herm C SV2");

            tmv::HermMatrix<CT,uplo|stor> c2 = c;
            SV_Decompose(c2,cS2);
            Assert(Equal(cS2,cS,ceps*normc),"Herm C SV3");
            SV_Decompose(c,cU2,cS2);
            Assert(Equal(cS2,cS,ceps*normc),"Herm C SV4 S");
            Assert(Equal(
                    tmv::Matrix<CT>(c/x)*tmv::Matrix<CT>(c.adjoint()/x),
                    cU2*tmv::DiagMatrix<CT>(cS2/x)*
                    tmv::DiagMatrix<CT>(cS2/x)*cU2.adjoint(),
                    ceps*(normc/x)*(normc/x)),"Herm C SV4 U");
            SV_Decompose(c,cS2,cVt2);
            Assert(Equal(cS2,cS,ceps*normc),"Herm C SV5 S");
            Assert(Equal(
                    tmv::Matrix<CT>(c.adjoint()/x)*tmv::Matrix<CT>(c/x),
                    cVt2.adjoint()*tmv::DiagMatrix<CT>(cS2/x)*
                    tmv::DiagMatrix<CT>(cS2/x)*cVt2,
                    ceps*(normc/x)*(normc/x)),"Herm C SV5 Vt");

            c2 = c.conjugate();
            SV_Decompose(c2.conjugate(),cS2);
            Assert(Equal(cS2,cS,ceps*normc),"Herm C SV6");
            SV_Decompose(c,cU2.conjugate(),cS2);
            Assert(Equal(cS2,cS,ceps*normc),"Herm C SV4 S");
            Assert(Equal(
                    tmv::Matrix<CT>(c/x)*tmv::Matrix<CT>(c.adjoint()/x),
                    cU2.conjugate()*tmv::DiagMatrix<CT>(cS2/x)*
                    tmv::DiagMatrix<CT>(cS2/x)*cU2.transpose(),
                    ceps*(normc/x)*(normc/x)), "Herm C SV4 U");
#endif
            std::cout<<"."; std::cout.flush();
        }

        // Eigen
        {
            if (showstartdone) std::cout<<"Eigen"<<std::endl;
            tmv::Matrix<T> V(N,N);
            tmv::Vector<T> L(N);
            Eigen(m,V,L);
            Assert(Equal(m*V,V*DiagMatrixViewOf(L),eps*normm),"Herm Eigen");

            tmv::Matrix<CT> cV(N,N);
            tmv::Vector<T> cL(N);
            Eigen(c,cV,cL);
            Assert(Equal(c*cV,cV*DiagMatrixViewOf(cL),eps*normc),
                   "Herm C Eigen");

#if (XTEST & 16)
            tmv::Vector<T> L2(N);
            Eigen(m,L2);
            Assert(Equal(L2,L,eps*normm),"Herm Eigen2");

            tmv::Vector<T> cL2(N);
            Eigen(c,cL2);
            Assert(Equal(cL2,cL,eps*normc),"Herm C Eigen2");

            Eigen(c,cV.conjugate(),cL);
            Assert(Equal(
                    c*cV.conjugate(),cV.conjugate()*DiagMatrixViewOf(cL),
                    eps*normc),"Herm C Eigen3");

            Eigen(c.conjugate(),cV,cL);
            Assert(Equal(c.conjugate()*cV,cV*DiagMatrixViewOf(cL),eps*normc),
                   "Herm C Eigen4");

            Eigen(c.conjugate(),cV.conjugate(),cL);
            Assert(Equal(c*cV,cV*DiagMatrixViewOf(cL),eps*normc),
                   "Herm C Eigen5");

            Eigen(c.conjugate(),cL2);
            Assert(Equal(cL2,cL,eps*normc),"Herm C Eigen6");
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
                tmv::HermMatrix<T,uplo|stor> S = m;
                SquareRoot(S);
                Assert(Equal(m,S*S,eps*normm),"Herm Square Root");

                tmv::HermMatrix<CT,uplo|stor> cS = c;
                SquareRoot(cS);

#if (XTEST & 16)
                cS.conjugate() = c;
                SquareRoot(cS.conjugate());
                Assert(Equal(c,cS.conjugate()*cS.conjugate(),eps*normc),
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
    typedef std::complex<T> CT;
    for (int mattype = START; mattype <= 5; mattype++) {
#if !(XTEST & 64) || defined(LAP)
        if (mattype >= 3) break;
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
        const T Teps = std::numeric_limits<T>::epsilon();
        const T Tmin = std::numeric_limits<T>::min();
        const T Tmax = std::numeric_limits<T>::max();
        if (showstartdone) {
            std::cout<<"posdef = "<<posdef<<
                ", singular = "<<singular<<std::endl;
        }

        tmv::SymMatrix<T,uplo|stor> m(N);
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

        tmv::SymMatrix<CT,uplo|stor> c(N);
        for(int i=0;i<N;++i) for(int j=0;j<N;++j) 
            c(i,j) = CT(2+4*i-5*j,3-i);
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
            T x = Tmin;
            m.setAllTo(x);
            c.upperTri().setAllTo(CT(x,2*x));
            for(int i=1;i<N;++i) {
                m.col(i,0,i+1) /= T(i);
                c.col(i,0,i+1) /= T(i);
            }
            m(N/2,N/2) = x/Teps;
            c(N/2,N/2) = x/Teps;
            m.diag(0,N/5,2*N/5).setZero();
            c.diag(0,N/5,2*N/5).setZero();
            m.diag(0,N/2+1,N) *= T(-1);
            c.diag(0,N/2+1,N) *= T(-1);
            m.subSymMatrix(N-4,N).setAllTo(-x/Teps/T(8));
            c.subSymMatrix(N-4,N).setAllTo(-x/Teps/T(8));
        }

        if (nearoverflow) {
            T x = Tmax;
            x /= N;
            m.setAllTo(x);
            c.upperTri().setAllTo(CT(x,x/2));
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
#ifdef NOTHROW
        if (!singular) {
#else
            try  {
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
                Assert(Equal(m,LDL,eps*normm),"Sym LDL");

                tmv::LowerTriMatrix<CT,tmv::UnitDiag> cL = c.lud().getL();
                tmv::BandMatrix<CT> cD = c.lud().getD();
                tmv::Permutation cP = c.lud().getP();
                tmv::Matrix<CT> cLDL = 
                    cP*cL*cD*cL.transpose()*cP.transpose();
                if (showacc) {
                    std::cout<<"Norm(c-cLDL) = "<<Norm(c-cLDL)<<std::endl;
                    std::cout<<"ceps*Norm(c) = "<<ceps*normc<<std::endl;
                }
                Assert(Equal(c,cLDL,ceps*normc),"Sym C LDL");

#if (XTEST & 16)
                tmv::SymMatrix<T,uplo|stor> m2 = m;
                tmv::SymBandMatrix<T,uplo|stor> D2(N,1);
                tmv::Permutation P2(N);
                LDL_Decompose(m2,D2,P2);
                L = m2.unitLowerTri();
                LDL = P2*L*D2*L.transpose()*P2.transpose();
                Assert(Equal(m,LDL,eps*normm),"Sym LDL2");

                tmv::SymMatrix<CT,uplo|stor> c2 = c;
                tmv::SymBandMatrix<CT,uplo|stor> cD2(N,1);
                LDL_Decompose(c2,cD2,P2);
                cL = c2.unitLowerTri();
                cLDL = P2*cL*cD2*cL.transpose()*P2.transpose();
                Assert(Equal(c,cLDL,ceps*normc),"Sym C LDL2");

                c2.conjugate() = c;
                LDL_Decompose(c2.conjugate(),cD2,P2);
                cL = c2.conjugate().unitLowerTri();
                cLDL = P2*cL*cD2*cL.transpose()*P2.transpose();
                Assert(Equal(c,cLDL,ceps*normc),"Sym C LDL3");

                c2= c;
                LDL_Decompose(c2,cD2.conjugate(),P2);
                cL = c2.unitLowerTri();
                cLDL = P2*cL*cD2.conjugate()*cL.transpose()*P2.transpose();
                Assert(Equal(c,cLDL,ceps*normc),"Sym C LDL4");

                c2.conjugate() = c;
                LDL_Decompose(c2.conjugate(),cD2.conjugate(),P2);
                cL = c2.conjugate().unitLowerTri();
                cLDL = P2*cL*cD2.conjugate()*cL.transpose()*P2.transpose();
                Assert(Equal(c,cLDL,ceps*normc),"Sym C LDL5");
#endif
                std::cout<<"."; std::cout.flush();
#ifndef NOTHROW
            } catch (tmv::Singular) {
                // The Lapack version throws when matrix is exactly singular,
                // but native algorithm succeeds for all matrices.
                Assert(singular,"caught Singular but mattype != singular");
            }
#else
        }
#endif

        // SV Decomposition
        {
            if (showstartdone) std::cout<<"SV"<<std::endl;
            tmv::Matrix<T> U = m.svd().getU();
            tmv::DiagMatrix<T> S = m.svd().getS();
            tmv::Matrix<T> Vt = m.svd().getVt();
            if (showacc) {
                std::cout<<"Norm(m-U*S*Vt) = "<<Norm(m-U*S*Vt)<<std::endl;
                std::cout<<"eps*Norm(m) = "<<eps*normm<<std::endl;
            }
            Assert(Equal(m,U*S*Vt,eps*normm),"Sym SV");

            tmv::Matrix<CT> cU = c.symsvd().getU();
            tmv::DiagMatrix<T> cS = c.symsvd().getS();
            tmv::Matrix<CT> cVt = c.symsvd().getVt();
            if (showacc) {
                std::cout<<"Norm(c-U*S*Vt) = "<<Norm(c-cU*cS*cVt)<<std::endl;
                std::cout<<"eps*Norm(c) = "<<eps*normc<<std::endl;
            }
            Assert(Equal(c,cU*cS*cVt,ceps*normc),"Sym C SV");

#if (XTEST & 16)
            const T x = 
                nearoverflow ? Tmax*Teps :
                nearunderflow ? Tmin/Teps : T(1);

            tmv::Matrix<T> U2(N,N);
            tmv::DiagMatrix<T> S2(N);
            tmv::Matrix<T> Vt2(N,N);
            SV_Decompose(m,U2,S2,Vt2);
            Assert(Equal(m,U2*S2*Vt2,eps*normm),"Sym SV2");

            tmv::SymMatrix<T,uplo|stor> m2 = m;
            SV_Decompose(m2,S2);
            Assert(Equal(S2,S,eps*normm),"Sym SV3");
            SV_Decompose(m,U2,S2);
            Assert(Equal(S2,S,eps*normm),"Sym SV4 S");
            Assert(Equal(
                    tmv::Matrix<T>(m/x)*tmv::Matrix<T>(m.adjoint()/x),
                    U2*tmv::DiagMatrix<T>(S2/x)*
                    tmv::DiagMatrix<T>(S2/x)*U2.adjoint(),
                    eps*(normm/x)*(normm/x)), "Sym SV4 U");
            SV_Decompose(m,S2,Vt2);
            Assert(Equal(S2,S,eps*normm),"Sym SV5 S");
            Assert(Equal(
                    tmv::Matrix<T>(m.adjoint()/x)*tmv::Matrix<T>(m/x),
                    Vt2.adjoint()*tmv::DiagMatrix<T>(S2/x)*
                    tmv::DiagMatrix<T>(S2/x)*Vt2,
                    eps*(normm/x)*(normm/x)), "Sym SV5 Vt");

            tmv::Matrix<CT> cU2(N,N);
            tmv::DiagMatrix<T> cS2(N);
            tmv::Matrix<CT> cVt2(N,N);
            SV_Decompose(c,cU2,cS2,cVt2);
            Assert(Equal(c,cU2*cS2*cVt2,ceps*normc),"Sym C SV2");

            tmv::SymMatrix<CT,uplo|stor> c2 = c;
            SV_Decompose(c2,cS2);
            Assert(Equal(cS2,cS,ceps*normc),"Sym C SV3");
            SV_Decompose(c,cU2,cS2);
            Assert(Equal(cS2,cS,ceps*normc),"Sym C SV4 S");
            Assert(Equal(
                    tmv::Matrix<CT>(c/x)*tmv::Matrix<CT>(c.adjoint()/x),
                    cU2*tmv::DiagMatrix<CT>(cS2/x)*
                    tmv::DiagMatrix<CT>(cS2/x)*cU2.adjoint(),
                    ceps*(normc/x)*(normc/x)),"Sym C SV4 U");
            SV_Decompose(c,cS2,cVt2);
            Assert(Equal(cS2,cS,ceps*normc),"Sym C SV5 S");
            Assert(Equal(tmv::Matrix<CT>(c.adjoint()/x)*tmv::Matrix<CT>(c/x),
                        cVt2.adjoint()*tmv::DiagMatrix<CT>(cS2/x)*
                        tmv::DiagMatrix<CT>(cS2/x)*cVt2,
                        ceps*(normc/x)*(normc/x)),"Sym C SV5 Vt");

            c2 = c.conjugate();
            SV_Decompose(c2.conjugate(),cS2);
            Assert(Equal(cS2,cS,ceps*normc),"Sym C SV6");
            SV_Decompose(c,cU2.conjugate(),cS2);
            Assert(Equal(cS2,cS,ceps*normc),"Sym C SV7 S");
            Assert(Equal(
                    tmv::Matrix<CT>(c/x)*tmv::Matrix<CT>(c.adjoint()/x),
                    cU2.conjugate()*tmv::DiagMatrix<CT>(cS2/x)*
                    tmv::DiagMatrix<CT>(cS2/x)*cU2.transpose(),
                    ceps*(normc/x)*(normc/x)),"Sym C SV7 U");
#endif
            std::cout<<"."; std::cout.flush();
        }
    }
}

template <class T, tmv::StorageType stor> 
void TestPolar()
{
    typedef std::complex<T> CT;
    if (showstartdone) std::cout<<"PolarDecomp "<<TMV_Text(stor)<<std::endl;

    for (int mattype = START; mattype <= 6; mattype++) {
#if !(XTEST & 64) || defined(LAP)
        if (mattype >= 4) break;
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
        const T Teps = std::numeric_limits<T>::epsilon();
        const T Tmin = std::numeric_limits<T>::min();
        const T Tmax = std::numeric_limits<T>::max();

        tmv::Matrix<T,stor> m(M,N);
        for(int i=0;i<M;++i) for(int j=0;j<N;++j) m(i,j) = T(2+4*i-5*j);
        m(0,0) = T(14);
        m(1,0) = T(-2);
        m(2,0) = T(7);
        m(3,4) = T(-10);
        if (!singular) m.diag() *= T(30);
        else { m.col(1).setZero(); m.row(7) = m.row(6); }

        tmv::Matrix<CT,stor> c(M,N);
        for(int i=0;i<M;++i) for(int j=0;j<N;++j) 
            c(i,j) = CT(2+4*i-5*j,3-i);
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
            T x = Tmin;
            m.setAllTo(x);
            c.setAllTo(CT(x,2*x));
            for(int i=1;i<N;++i) {
                m.col(i) /= T(i);
                c.col(i) /= T(i);
            }
            m(N/2,N/2) = x/Teps;
            c(N/2,N/2) = x/Teps;
            m.diag(0,N/5,2*N/5).setZero();
            c.diag(0,N/5,2*N/5).setZero();
            m.diag(0,N/2+1,N) *= T(-1);
            c.diag(0,N/2+1,N) *= T(-1);
            m.subMatrix(N-4,N,N-4,N).setAllTo(-x/Teps/T(8));
            c.subMatrix(N-4,N,N-4,N).setAllTo(-x/Teps/T(8));
        }

        if (nearoverflow) {
            T x = Tmax;
            x /= N;
            m.setAllTo(x);
            c.setAllTo(CT(x,x/2));
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
            tmv::HermMatrix<T,tmv::Lower|tmv::ColMajor> P(N);
            tmv::Matrix<T,stor> U = m;
            PolarDecompose(U,P);
            if (showacc) {
                std::cout<<"Norm(m-UP) = "<<Norm(m-U*P)<<std::endl;
                std::cout<<"eps*Norm(m) = "<<eps*normm<<std::endl;
                std::cout<<"Norm(UtU-1) = "<<Norm(U.adjoint()*U-T(1))<<std::endl;
                std::cout<<"eps = "<<eps<<std::endl;
            }
            Assert(Equal(m,U*P,eps*normm),"Polar");
            Assert(Equal(U.adjoint()*U,T(1),eps),"Polar UtU");

            tmv::HermMatrix<CT,tmv::Lower|tmv::ColMajor> cP(N);
            tmv::Matrix<CT,stor> cU = c;
            PolarDecompose(cU,cP);
            Assert(Equal(c,cU*cP,ceps*normc),"C Polar");
            Assert(Equal(cU.adjoint()*cU,T(1),ceps),"C Polar UtU");

#if (XTEST & 16)
            U = m;
            PolarDecompose(U,P.transpose());
            Assert(Equal(m,U*P,eps*normm),"Polar2");
            Assert(Equal(U.adjoint()*U,T(1),eps),"Polar2 UtU");

            cU = c;
            PolarDecompose(cU,cP.adjoint());
            Assert(Equal(c,cU*cP,ceps*normc),"C Polar2");
            Assert(Equal(cU.adjoint()*cU,T(1),ceps),"C Polar2 UtU");

            cU = c;
            PolarDecompose(cU,cP.transpose());
            Assert(Equal(c,cU*cP.transpose(),ceps*normc),"C Polar3");
            Assert(Equal(cU.adjoint()*cU,T(1),ceps),"C Polar3 UtU");

            cU = c;
            PolarDecompose(cU,cP.conjugate());
            Assert(Equal(c,cU*cP.transpose(),ceps*normc),"C Polar4");
            Assert(Equal(cU.adjoint()*cU,T(1),ceps),"C Polar4 UtU");

            cU.conjugate() = c;
            PolarDecompose(cU.conjugate(),cP);
            Assert(Equal(c,cU.conjugate()*cP,ceps*normc),"C Polar5");
            Assert(Equal(cU.adjoint()*cU,T(1),ceps),"C Polar5 UtU");

            cU.conjugate() = c;
            PolarDecompose(cU.conjugate(),cP.adjoint());
            Assert(Equal(c,cU.conjugate()*cP,ceps*normc),"C Polar6");
            Assert(Equal(cU.adjoint()*cU,T(1),ceps),"C Polar6 UtU");

            cU.conjugate() = c;
            PolarDecompose(cU.conjugate(),cP.transpose());
            Assert(Equal(c,cU.conjugate()*cP.transpose(),ceps*normc),
                   "C Polar7");
            Assert(Equal(cU.adjoint()*cU,T(1),ceps),"C Polar7 UtU");

            cU.conjugate() = c;
            PolarDecompose(cU.conjugate(),cP.conjugate());
            Assert(Equal(c,cU.conjugate()*cP.transpose(),ceps*normc),
                   "C Polar8");
            Assert(Equal(cU.adjoint()*cU,T(1),ceps),"C Polar8 UtU");
#endif
            std::cout<<"."; std::cout.flush();
        } while (false);

        // BandMatrix Polar Decomposition
        do {
            if (showstartdone) std::cout<<"Band Polar"<<std::endl;
            tmv::BandMatrixView<T> b(m.view(),5,11);
            tmv::BandMatrixView<CT> cb(c.view(),5,11);
            T normb = Norm(b);
            T normcb = Norm(cb);

            tmv::HermMatrix<T,tmv::Lower|tmv::ColMajor> P(N);
            tmv::Matrix<T,stor> U(b.colsize(),N);
            PolarDecompose(b,U,P);
            if (showacc) {
                std::cout<<"Norm(b-UP) = "<<Norm(b-U*P)<<"  "<<eps*normb<<std::endl;
                std::cout<<"Norm(UtU-1) = "<<Norm(U.adjoint()*U-T(1))<<"  "<<eps<<std::endl;
            }
            Assert(Equal(b,U*P,eps*normb),"Band Polar");
            Assert(Equal(U.adjoint()*U,T(1),eps),"Band Polar UtU");

            tmv::HermMatrix<CT,tmv::Lower|tmv::ColMajor> cP(N);
            tmv::Matrix<CT,stor> cU(cb.colsize(),N);
            PolarDecompose(cb,cU,cP);
            if (showacc) {
                std::cout<<"Norm(cb-UP) = "<<Norm(cb-cU*cP)<<
                    "  "<<ceps*normcb<<std::endl;
                std::cout<<"Norm(UtU-1) = "<<Norm(cU.adjoint()*cU-T(1))<<
                    "  "<<ceps<<std::endl;
            }
            Assert(Equal(cb,cU*cP,ceps*normcb),"C Band Polar");
            Assert(Equal(cU.adjoint()*cU,T(1),ceps),"C Band Polar UtU");

#if (XTEST & 16)
            PolarDecompose(b,U,P.transpose());
            Assert(Equal(b,U*P,eps*normb),"Band Polar2");
            Assert(Equal(U.adjoint()*U,T(1),eps),"Band Polar2 UtU");

            PolarDecompose(cb,cU,cP.adjoint());
            Assert(Equal(cb,cU*cP,ceps*normcb),"C Band Polar2");
            Assert(Equal(cU.adjoint()*cU,T(1),ceps),"C Band Polar2 UtU");

            PolarDecompose(cb,cU,cP.transpose());
            Assert(Equal(cb,cU*cP.transpose(),ceps*normcb),"C Band Polar3");
            Assert(Equal(cU.adjoint()*cU,T(1),ceps),"C Band Polar3 UtU");

            PolarDecompose(cb,cU,cP.conjugate());
            Assert(Equal(cb,cU*cP.transpose(),ceps*normcb),"C Band Polar4");
            Assert(Equal(cU.adjoint()*cU,T(1),ceps),"C Band Polar4 UtU");

            PolarDecompose(cb,cU.conjugate(),cP);
            Assert(Equal(cb,cU.conjugate()*cP,ceps*normcb),"C Band Polar5");
            Assert(Equal(cU.adjoint()*cU,T(1),ceps),"C Band Polar5 UtU");

            PolarDecompose(cb,cU.conjugate(),cP.adjoint());
            Assert(Equal(cb,cU.conjugate()*cP,ceps*normcb),"C Band Polar6");
            Assert(Equal(cU.adjoint()*cU,T(1),ceps),"C Band Polar6 UtU");

            PolarDecompose(cb,cU.conjugate(),cP.transpose());
            Assert(Equal(cb,cU.conjugate()*cP.transpose(),ceps*normcb),
                   "C Band Polar7");
            Assert(Equal(cU.adjoint()*cU,T(1),ceps),"C Band Polar7 UtU");

            PolarDecompose(cb,cU.conjugate(),cP.conjugate());
            Assert(Equal(cb,cU.conjugate()*cP.transpose(),ceps*normcb),
                   "C Band Polar8");
            Assert(Equal(cU.adjoint()*cU,T(1),ceps),"C Band Polar8 UtU");
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



