// vim:et:ts=2:sw=2:ci:cino=f0,g0,t0,+0:

#define START 0

#include "TMV_Test.h"
#include "TMV_Test2.h"
#include "TMV.h"
#include "TMV_SymBand.h"
#include "TMV_TestSymArith.h"

template <class T, tmv::UpLoType uplo, tmv::StorageType stor> 
void TestHermDecomp()
{
  for (int mattype = 0; mattype < 3; mattype++) {
    //std::cout<<"mattype = "<<mattype<<std::endl;
    // mattype = 0  is PosDef
    // mattype = 1  is Indef
    // mattype = 2  is Singular

    const int N = 200;

    tmv::HermMatrix<T,uplo,stor> m(N);
    for(int i=0;i<N;++i) for(int j=0;j<N;++j) 
      if ((uplo == tmv::Upper && j>=i) || (uplo == tmv::Lower && i>=j)) 
        m(i,j) = T(2+4*i-5*j);
    if (mattype == 0) {
      m /= T(10*N);
      m.diag().AddToAll(T(1));
      for(int i=0;i<N;i++) m.diag()(i) += T(i);
      m.diag(0,N/4,N/2) *= T(10);
      m.SubSymMatrix(N/2,3*N/4,2) *= T(100);
    } else if (mattype == 1) {
      m /= T(N);
      m.diag(0,N/10,N/5).Zero();
      m.col(0,N/3,3*N/7).AddToAll(T(50));
      m.diag(1,4,7).AddToAll(T(10));
      m.row(7,0,5).AddToAll(T(10));
      m.diag(0,N/4,N/2) *= T(10);
      m.SubSymMatrix(N/2,3*N/4,2) /= T(100);
    } else {
      m.row(2,0,2).Zero();
      m.row(2,2,8).Zero();
      m.row(4,0,4) = m.row(5,0,4);
      m.row(4,5,8) = m.row(5,5,8);
      m(4,5) = m(4,4) = m(5,5);
    }

    tmv::HermMatrix<std::complex<T>,uplo,stor> c(N);
    for(int i=0;i<N;++i) for(int j=0;j<N;++j) 
      if ((uplo == tmv::Upper && j>=i) || (uplo == tmv::Lower && i>=j)) 
        c(i,j) = std::complex<T>(2+4*i-5*j,3-i);
    c.diag().Imag().Zero();
    if (mattype == 0) {
      c /= T(10*N);
      c.diag().AddToAll(T(1));
      for(int i=0;i<N;i++) c.diag()(i) += T(i);
      c.diag(0,N/4,N/2) *= T(10);
      c.SubSymMatrix(N/2,3*N/4,2) *= T(100);
    } else if (mattype == 1) {
      c /= T(N);
      c.diag(0,N/10,N/5).Zero();
      c.col(0,N/3,3*N/7).AddToAll(T(50));
      c.diag(1,4,7).AddToAll(T(10));
      c.row(7,0,5).AddToAll(T(10));
      c.diag(0,N/4,N/2) *= T(10);
      c.SubSymMatrix(N/2,3*N/4,2) /= T(100);
    } else {
      c.row(2,0,2).Zero();
      c.row(2,2,8).Zero();
      c.row(4,0,4) = c.row(5,0,4);
      c.row(4,5,8) = c.row(5,5,8);
      c(4,5) = c(4,4) = c(5,5);
    }


    T eps = EPS;
    T ceps = EPS;
    if (showacc) std::cout<<"eps = "<<eps<<"  "<<ceps<<std::endl;
    if (mattype != 2) {
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

    //std::cout<<"m = "<<m<<std::endl;
    //std::cout<<"c = "<<c<<std::endl;

    // CH Decomposition
#ifdef NOTHROW
    if (mattype == 0)
#else
      try 
#endif
      {
        tmv::LowerTriMatrix<T> L = m.CHD().GetL();
        tmv::Matrix<T> LLt = L*L.Adjoint();
        Assert(Norm(m-LLt) < eps*Norm(m),"Herm CH");

        tmv::LowerTriMatrix<std::complex<T> > cL = c.CHD().GetL();
        tmv::Matrix<std::complex<T> > cLLt = cL*cL.Adjoint();
        Assert(Norm(c-cLLt) < ceps*Norm(c),"Herm C CH");

#ifdef XTEST
        tmv::HermMatrix<T,uplo,stor> m2 = m;
        CH_Decompose(m2.View());
        L = m2.LowerTri();
        LLt = L*L.Adjoint();
        Assert(Norm(m-LLt) < eps*Norm(m),"Herm CH2");

        tmv::HermMatrix<std::complex<T>,uplo,stor> c2 = c;
        CH_Decompose(c2.View());
        cL = c2.LowerTri();
        cLLt = cL*cL.Adjoint();
        Assert(Norm(c-cLLt) < ceps*Norm(c),"Herm C CH2");

        c2.Conjugate() = c;
        CH_Decompose(c2.Conjugate());
        cL = c2.Conjugate().LowerTri();
        cLLt = cL*cL.Adjoint();
        Assert(Norm(c-cLLt) < ceps*Norm(c),"Herm C CH2");
        Assert(mattype == 0,"Didn't throw NonPosDef, mattype != pos def");
#endif
        std::cout<<"."; std::cout.flush();
      }
#ifndef NOTHROW
    catch(tmv::NonPosDef) 
    { Assert(mattype != 0,"Caught NonPosDef, mattype == pos def"); }
#endif

    // LDL Decomposition
    {
      tmv::LowerTriMatrix<T,tmv::UnitDiag> L = m.LUD().GetL();
      tmv::BandMatrix<T> D = m.LUD().GetD();
      const int* p = m.LUD().GetP();
      tmv::Matrix<T> LDL = L*D*L.Adjoint();
      LDL.ReversePermuteRows(p);
      LDL.ReversePermuteCols(p);
      Assert(Norm(m-LDL) < eps*Norm(m),"Herm LDL");

      tmv::LowerTriMatrix<std::complex<T>,tmv::UnitDiag> cL = c.LUD().GetL();
      tmv::BandMatrix<std::complex<T> > cD = c.LUD().GetD();
      p = c.LUD().GetP();
      tmv::Matrix<std::complex<T> > cLDL = cL*cD*cL.Adjoint();
      cLDL.ReversePermuteRows(p);
      cLDL.ReversePermuteCols(p);
      Assert(Norm(c-cLDL) < ceps*Norm(c),"Herm C LDL");

#ifdef XTEST
      tmv::HermMatrix<T,uplo,stor> m2 = m;
      tmv::HermBandMatrix<T,uplo,stor> D2(N,1);
      int p2[N];
      LDL_Decompose(m2.View(),D2.View(),p2);
      L = m2.LowerTri(tmv::UnitDiag);
      LDL = L*D2*L.Adjoint();
      LDL.ReversePermuteRows(p2);
      LDL.ReversePermuteCols(p2);
      Assert(Norm(m-LDL) < eps*Norm(m),"Herm LDL2");

      tmv::HermMatrix<std::complex<T>,uplo,stor> c2 = c;
      tmv::HermBandMatrix<std::complex<T>,uplo,stor> cD2(N,1);
      LDL_Decompose(c2.View(),cD2.View(),p2);
      cL = c2.LowerTri(tmv::UnitDiag);
      cLDL = cL*cD2*cL.Adjoint();
      cLDL.ReversePermuteRows(p2);
      cLDL.ReversePermuteCols(p2);
      Assert(Norm(c-cLDL) < ceps*Norm(c),"Herm C LDL2");

      c2.Conjugate() = c;
      LDL_Decompose(c2.Conjugate(),cD2.View(),p2);
      cL = c2.Conjugate().LowerTri(tmv::UnitDiag);
      cLDL = cL*cD2*cL.Adjoint();
      cLDL.ReversePermuteRows(p2);
      cLDL.ReversePermuteCols(p2);
      Assert(Norm(c-cLDL) < ceps*Norm(c),"Herm C LDL3");

      c2 = c;
      LDL_Decompose(c2.View(),cD2.Conjugate(),p2);
      cL = c2.LowerTri(tmv::UnitDiag);
      cLDL = cL*cD2.Conjugate()*cL.Adjoint();
      cLDL.ReversePermuteRows(p2);
      cLDL.ReversePermuteCols(p2);
      Assert(Norm(c-cLDL) < ceps*Norm(c),"Herm C LDL4");

      c2.Conjugate() = c;
      LDL_Decompose(c2.Conjugate(),cD2.Conjugate(),p2);
      cL = c2.Conjugate().LowerTri(tmv::UnitDiag);
      cLDL = cL*cD2.Conjugate()*cL.Adjoint();
      cLDL.ReversePermuteRows(p2);
      cLDL.ReversePermuteCols(p2);
      Assert(Norm(c-cLDL) < ceps*Norm(c),"Herm C LDL5");
#endif
      std::cout<<"."; std::cout.flush();
    }

    // SV Decomposition
    {
      tmv::Matrix<T> U = m.SVD().GetU();
      tmv::DiagMatrix<T> S = m.SVD().GetS();
      tmv::Matrix<T> V = m.SVD().GetV();
      Assert(Norm(m-U*S*V) < eps*Norm(m),"Herm SV");

      tmv::Matrix<std::complex<T> > cU = c.SVD().GetU();
      tmv::DiagMatrix<T> cS = c.SVD().GetS();
      tmv::Matrix<std::complex<T> > cV = c.SVD().GetV();
      Assert(Norm(c-cU*cS*cV) < ceps*Norm(c),"Herm C SV");

#ifdef XTEST
      tmv::Matrix<T> U2(N,N);
      tmv::DiagMatrix<T> S2(N);
      tmv::Matrix<T> V2(N,N);
      SV_Decompose(m,U2.View(),S2.View(),V2.View());
      Assert(Norm(m-U2*S2*V2) < eps*Norm(m),"Herm SV2");

      tmv::HermMatrix<T,uplo,stor> m2 = m;
      SV_Decompose(m2.View(),S2.View());
      Assert(Norm(S2-S) < eps*Norm(S),"Herm SV3");
      SV_Decompose(m,U2.View(),S2.View());
      Assert(Norm(S2-S) < eps*Norm(S),"Herm SV4 S");
      Assert(Norm(U2*S2*S2*U2.Adjoint()-m*m.Adjoint()) < 
          eps*Norm(m)*Norm(m),"Herm SV4 U");
      SV_Decompose(m,S2.View(),V2.View());
      Assert(Norm(S2-S) < eps*Norm(S),"Herm SV5 S");
      Assert(Norm(m.Adjoint()*m-V2.Adjoint()*S2*S2*V2) < 
          eps*Norm(m)*Norm(m),"Herm SV5 V");

      tmv::Matrix<std::complex<T> > cU2(N,N);
      tmv::DiagMatrix<T> cS2(N);
      tmv::Matrix<std::complex<T> > cV2(N,N);
      SV_Decompose(c,cU2.View(),cS2.View(),cV2.View());
      Assert(Norm(c-cU2*cS2*cV2) < ceps*Norm(c),"Herm C SV2");

      tmv::HermMatrix<std::complex<T>,uplo,stor> c2 = c;
      SV_Decompose(c2.View(),cS2.View());
      Assert(Norm(cS2-cS) < ceps*Norm(cS),"Herm C SV3");
      SV_Decompose(c,cU2.View(),cS2.View());
      Assert(Norm(cS2-cS) < ceps*Norm(cS),"Herm C SV4 S");
      Assert(Norm(cU2*cS2*cS2*cU2.Adjoint()-c*c.Adjoint()) < 
          ceps*Norm(c)*Norm(c),"Herm C SV4 U");
      SV_Decompose(c,cS2.View(),cV2.View());
      Assert(Norm(cS2-cS) < ceps*Norm(cS),"Herm C SV5 S");
      Assert(Norm(c.Adjoint()*c-cV2.Adjoint()*cS2*cS2*cV2) < 
          ceps*Norm(c)*Norm(c),"Herm C SV5 V");

      c2 = c.Conjugate();
      SV_Decompose(c2.Conjugate(),cS2.View());
      Assert(Norm(cS2-cS) < ceps*Norm(cS),"Herm C SV6");
      SV_Decompose(c,cU2.Conjugate(),cS2.View());
      Assert(Norm(cS2-cS) < ceps*Norm(cS),"Herm C SV4 S");
      Assert(Norm(cU2.Conjugate()*cS2*cS2*cU2.Transpose()-c*c.Adjoint()) < 
          ceps*Norm(c)*Norm(c),"Herm C SV4 U");
#endif
      std::cout<<"."; std::cout.flush();
    }

    // Eigen
    {
      tmv::Matrix<T> V(N,N);
      tmv::Vector<T> L(N);
      Eigen(m,V.View(),L.View());
      Assert(Norm(m*V-V*DiagMatrixViewOf(L)) < eps*Norm(m),"Herm Eigen");

      tmv::Matrix<std::complex<T> > cV(N,N);
      tmv::Vector<T> cL(N);
      Eigen(c,cV.View(),cL.View());
      Assert(Norm(c*cV-cV*DiagMatrixViewOf(cL)) < eps*Norm(c),"Herm C Eigen");

#ifdef XTEST
      tmv::Vector<T> L2(N);
      tmv::HermMatrix<T,uplo,stor> m2 = m;
      Eigen(m2.View(),L2.View());
      Assert(Norm(L2-L) < eps*Norm(L),"Herm Eigen2");

      tmv::Vector<T> cL2(N);
      tmv::HermMatrix<std::complex<T>,uplo,stor> c2 = c;
      Eigen(c2.View(),cL2.View());
      Assert(Norm(cL2-cL) < eps*Norm(cL),"Herm C Eigen2");

      Eigen(c,cV.Conjugate(),cL.View());
      Assert(Norm(c*cV.Conjugate()-cV.Conjugate()*DiagMatrixViewOf(cL)) < 
          eps*Norm(c),"Herm C Eigen3");

      Eigen(c.Conjugate(),cV.View(),cL.View());
      Assert(Norm(c.Conjugate()*cV-cV*DiagMatrixViewOf(cL)) < eps*Norm(c),
          "Herm C Eigen4");

      Eigen(c.Conjugate(),cV.Conjugate(),cL.View());
      Assert(Norm(c*cV-cV*DiagMatrixViewOf(cL)) < eps*Norm(c),"Herm C Eigen5");

      c2.Conjugate() = c;
      Eigen(c2.Conjugate(),cL2.View());
      Assert(Norm(cL2-cL) < eps*Norm(cL),"Herm C Eigen6");
#endif
      std::cout<<"."; std::cout.flush();
    }

    // Square Root
#ifdef NOTHROW
    if (mattype == 0)
#else
      try 
#endif
      {
        tmv::HermMatrix<T,uplo,stor> S = m;
        SquareRoot(S.View());
        Assert(Norm(m-S*S) < eps*Norm(m),"Herm Square Root");

        tmv::HermMatrix<std::complex<T>,uplo,stor> cS = c;
        SquareRoot(cS.View());

#ifdef XTEST
        cS.Conjugate() = c;
        SquareRoot(cS.Conjugate());
        Assert(Norm(c-cS.Conjugate()*cS.Conjugate()) < eps*Norm(c),
            "Herm C Square Root 2");
#endif
        Assert(mattype == 0,"Didn't throw NonPosDef, mattype != pos def");
        std::cout<<"."; std::cout.flush();
      }
#ifndef NOTHROW
    catch(tmv::NonPosDef) 
    { Assert(mattype != 0,"Caught NonPosDef, mattype == pos def"); }
#endif
  }
}

template <class T, tmv::UpLoType uplo, tmv::StorageType stor> 
void TestSymDecomp()
{
  for (int mattype = 0; mattype < 2; mattype++) {
    // mattype = 0  is Normal
    // mattype = 1  is Singular

    const int N = 200;

    tmv::SymMatrix<T,uplo,stor> m(N);
    for(int i=0;i<N;++i) for(int j=0;j<N;++j) 
      if ((uplo == tmv::Upper && j>=i) || (uplo == tmv::Lower && i>=j)) 
        m(i,j) = T(2+4*i-5*j);
    if (mattype == 0) {
      m /= T(10*N);
      m.diag().AddToAll(T(1));
      for(int i=0;i<N;i++) m.diag()(i) += T(i);
      m.diag(0,N/4,N/2) *= T(10);
      m.SubSymMatrix(N/2,3*N/4,2) *= T(100);
    } else if (mattype == 1) {
      m /= T(N);
      m.diag(0,N/10,N/5).Zero();
      m.col(0,N/3,3*N/7).AddToAll(T(50));
      m.diag(1,4,7).AddToAll(T(10));
      m.row(7,0,5).AddToAll(T(10));
      m.diag(0,N/4,N/2) *= T(10);
      m.SubSymMatrix(N/2,3*N/4,2) /= T(100);
    } else {
      m.row(2,0,2).Zero();
      m.row(2,2,8).Zero();
      m.row(4,0,4) = m.row(5,0,4);
      m.row(4,5,8) = m.row(5,5,8);
      m(4,5) = m(4,4) = m(5,5);
    }


    tmv::SymMatrix<std::complex<T>,uplo,stor> c(N);
    for(int i=0;i<N;++i) for(int j=0;j<N;++j) 
      c(i,j) = std::complex<T>(2+4*i-5*j,3-i);
    if (mattype == 0) {
      c /= T(10*N);
      c.diag().AddToAll(T(1));
      for(int i=0;i<N;i++) c.diag()(i) += T(i);
      c.diag(0,N/4,N/2) *= T(10);
      c.SubSymMatrix(N/2,3*N/4,2) *= T(100);
    } else if (mattype == 1) {
      c /= T(N);
      c.diag(0,N/10,N/5).Zero();
      c.col(0,N/3,3*N/7).AddToAll(T(50));
      c.diag(1,4,7).AddToAll(T(10));
      c.row(7,0,5).AddToAll(T(10));
      c.diag(0,N/4,N/2) *= T(10);
      c.SubSymMatrix(N/2,3*N/4,2) /= T(100);
    } else {
      c.row(2,0,2).Zero();
      c.row(2,2,8).Zero();
      c.row(4,0,4) = c.row(5,0,4);
      c.row(4,5,8) = c.row(5,5,8);
      c(4,5) = c(4,4) = c(5,5);
    }


    T eps = EPS;
    T ceps = EPS;
    if (mattype != 1) {
      eps *= Norm(m) * Norm(m.Inverse());
      ceps *= Norm(c) * Norm(c.Inverse());
    } else {
      eps *= T(10*N);
      ceps *= T(10*N);
    }

    // LDL Decomposition
    {
      tmv::LowerTriMatrix<T,tmv::UnitDiag> L = m.LUD().GetL();
      tmv::BandMatrix<T> D = m.LUD().GetD();
      const int* p = m.LUD().GetP();
      tmv::Matrix<T> LDL = L*D*L.Transpose();
      LDL.ReversePermuteRows(p);
      LDL.ReversePermuteCols(p);
      Assert(Norm(m-LDL) < eps*Norm(m),"Sym LDL");

      tmv::LowerTriMatrix<std::complex<T>,tmv::UnitDiag> cL = c.LUD().GetL();
      tmv::BandMatrix<std::complex<T> > cD = c.LUD().GetD();
      p = c.LUD().GetP();
      tmv::Matrix<std::complex<T> > cLDL = cL*cD*cL.Transpose();
      cLDL.ReversePermuteRows(p);
      cLDL.ReversePermuteCols(p);
      Assert(Norm(c-cLDL) < ceps*Norm(c),"Sym C LDL");

#ifdef XTEST
      tmv::SymMatrix<T,uplo,stor> m2 = m;
      tmv::SymBandMatrix<T,uplo,stor> D2(N,1);
      int p2[N];
      LDL_Decompose(m2.View(),D2.View(),p2);
      L = m2.LowerTri(tmv::UnitDiag);
      LDL = L*D2*L.Transpose();
      LDL.ReversePermuteRows(p2);
      LDL.ReversePermuteCols(p2);
      Assert(Norm(m-LDL) < eps*Norm(m),"Sym LDL2");

      tmv::SymMatrix<std::complex<T>,uplo,stor> c2 = c;
      tmv::SymBandMatrix<std::complex<T>,uplo,stor> cD2(N,1);
      LDL_Decompose(c2.View(),cD2.View(),p2);
      cL = c2.LowerTri(tmv::UnitDiag);
      cLDL = cL*cD2*cL.Transpose();
      cLDL.ReversePermuteRows(p2);
      cLDL.ReversePermuteCols(p2);
      Assert(Norm(c-cLDL) < ceps*Norm(c),"Sym C LDL2");

      c2.Conjugate() = c;
      LDL_Decompose(c2.Conjugate(),cD2.View(),p2);
      cL = c2.Conjugate().LowerTri(tmv::UnitDiag);
      cLDL = cL*cD2*cL.Transpose();
      cLDL.ReversePermuteRows(p2);
      cLDL.ReversePermuteCols(p2);
      Assert(Norm(c-cLDL) < ceps*Norm(c),"Sym C LDL3");

      c2= c;
      LDL_Decompose(c2.View(),cD2.Conjugate(),p2);
      cL = c2.LowerTri(tmv::UnitDiag);
      cLDL = cL*cD2.Conjugate()*cL.Transpose();
      cLDL.ReversePermuteRows(p2);
      cLDL.ReversePermuteCols(p2);
      Assert(Norm(c-cLDL) < ceps*Norm(c),"Sym C LDL4");

      c2.Conjugate() = c;
      LDL_Decompose(c2.Conjugate(),cD2.Conjugate(),p2);
      cL = c2.Conjugate().LowerTri(tmv::UnitDiag);
      cLDL = cL*cD2.Conjugate()*cL.Transpose();
      cLDL.ReversePermuteRows(p2);
      cLDL.ReversePermuteCols(p2);
      Assert(Norm(c-cLDL) < ceps*Norm(c),"Sym C LDL5");
#endif
      std::cout<<"."; std::cout.flush();
    }

    // SV Decomposition
    {
      tmv::Matrix<T> U = m.SVD().GetU();
      tmv::DiagMatrix<T> S = m.SVD().GetS();
      tmv::Matrix<T> V = m.SVD().GetV();
      Assert(Norm(m-U*S*V) < eps*Norm(m),"Sym SV");

      tmv::Matrix<std::complex<T> > cU = c.SymSVD().GetU();
      tmv::DiagMatrix<T> cS = c.SymSVD().GetS();
      tmv::Matrix<std::complex<T> > cV = c.SymSVD().GetV();
      Assert(Norm(c-cU*cS*cV) < ceps*Norm(c),"Sym C SV");

#ifdef XTEST
      tmv::Matrix<T> U2(N,N);
      tmv::DiagMatrix<T> S2(N);
      tmv::Matrix<T> V2(N,N);
      SV_Decompose(m,U2.View(),S2.View(),V2.View());
      Assert(Norm(m-U2*S2*V2) < eps*Norm(m),"Sym SV2");

      tmv::SymMatrix<T,uplo,stor> m2 = m;
      SV_Decompose(m2.View(),S2.View());
      Assert(Norm(S2-S) < eps*Norm(S),"Sym SV3");
      SV_Decompose(m,U2.View(),S2.View());
      Assert(Norm(S2-S) < eps*Norm(S),"Sym SV4 S");
      Assert(Norm(U2-U) < eps*Norm(U),"Sym SV4 U");
      SV_Decompose(m,S2.View(),V2.View());
      Assert(Norm(S2-S) < eps*Norm(S),"Sym SV5 S");
      Assert(Norm(m.Adjoint()*m-V2.Adjoint()*S2*S2*V2) < eps*Norm(m)*Norm(m),
          "Sym SV5 V");

      tmv::Matrix<std::complex<T> > cU2(N,N);
      tmv::DiagMatrix<T> cS2(N);
      tmv::Matrix<std::complex<T> > cV2(N,N);
      SV_Decompose(c,cU2.View(),cS2.View(),cV2.View());
      Assert(Norm(c-cU2*cS2*cV2) < ceps*Norm(c),"Sym C SV2");

      tmv::SymMatrix<std::complex<T>,uplo,stor> c2 = c;
      SV_Decompose(c2.View(),cS2.View());
      Assert(Norm(cS2-cS) < ceps*Norm(cS),"Sym C SV3");
      SV_Decompose(c,cU2.View(),cS2.View());
      Assert(Norm(cS2-cS) < ceps*Norm(cS),"Sym C SV4 S");
      Assert(Norm(cU2-cU) < ceps*Norm(cU),"Sym C SV4 U");
      SV_Decompose(c,cS2.View(),cV2.View());
      Assert(Norm(cS2-cS) < ceps*Norm(cS),"Sym C SV5 S");
      Assert(Norm(c.Adjoint()*c-cV2.Adjoint()*cS2*cS2*cV2) < ceps*Norm(c)*Norm(c),
          "Sym C SV5 V");

      c2 = c.Conjugate();
      SV_Decompose(c2.Conjugate(),cS2.View());
      Assert(Norm(cS2-cS) < ceps*Norm(cS),"Sym C SV6");
      SV_Decompose(c,cU2.Conjugate(),cS2.View());
      Assert(Norm(cS2-cS) < ceps*Norm(cS),"Sym C SV7 S");
      Assert(Norm(cU2.Conjugate()-cU) < ceps*Norm(cU),"Sym C SV7 U");
#endif
      std::cout<<"."; std::cout.flush();
    }
  }
}

template <class T, tmv::StorageType stor> 
void TestPolar()
{
  for (int mattype = 0; mattype < 4; mattype++) {
    //std::cout<<"mattype = "<<mattype<<std::endl;
    // mattype = 0  is Square
    // mattype = 1  is NonSquare slightly tall
    // mattype = 2  is NonSquare very tall
    // mattype = 3  is Singular

    const int N = 200;
    int M = N;
    if (mattype == 1) M = 211;
    else if (mattype == 2) M = 545;

    tmv::Matrix<T,stor> m(M,N);
    for(int i=0;i<M;++i) for(int j=0;j<N;++j) m(i,j) = T(2+4*i-5*j);
    m(0,0) = T(14);
    m(1,0) = T(-2);
    m(2,0) = T(7);
    m(3,4) = T(-10);
    if (mattype != 3) m.diag() *= T(30);

    tmv::Matrix<std::complex<T>,stor> c(M,N);
    for(int i=0;i<M;++i) for(int j=0;j<N;++j) 
      c(i,j) = std::complex<T>(2+4*i-5*j,3-i);
    c(0,0) *= T(14);
    c(1,0) *= T(-2);
    c(2,0) *= T(7);
    c(7,6) *= T(-10);
    if (mattype != 3) c.diag() *= T(30);

    T eps = EPS;
    T ceps = EPS;
    if (mattype != 3) {
      eps *= Norm(m) * Norm(m.Inverse());
      ceps *= Norm(c) * Norm(c.Inverse());
    } else {
      eps *= T(100);
      ceps *= T(100);
    }

    // Matrix Polar Decomposition
    {
      tmv::HermMatrix<T,tmv::Lower,tmv::ColMajor> P(N);
      tmv::Matrix<T,stor> U = m;
      Polar_Decompose(U.View(),P.View());
      Assert(Norm(m-U*P) < eps*Norm(m),"Polar");
      Assert(Norm(U.Adjoint()*U-T(1)) < eps,"Polar UtU");

      tmv::HermMatrix<std::complex<T>,tmv::Lower,tmv::ColMajor> cP(N);
      tmv::Matrix<std::complex<T>,stor> cU = c;
      Polar_Decompose(cU.View(),cP.View());
      Assert(Norm(c-cU*cP) < ceps*Norm(c),"C Polar");
      Assert(Norm(cU.Adjoint()*cU-T(1)) < ceps,"C Polar UtU");

#ifdef XTEST
      U = m;
      Polar_Decompose(U.View(),P.Transpose());
      Assert(Norm(m-U*P) < eps*Norm(m),"Polar2");
      Assert(Norm(U.Adjoint()*U-T(1)) < eps,"Polar2 UtU");

      cU = c;
      Polar_Decompose(cU.View(),cP.Adjoint());
      Assert(Norm(c-cU*cP) < ceps*Norm(c),"C Polar2");
      Assert(Norm(cU.Adjoint()*cU-T(1)) < ceps,"C Polar2 UtU");

      cU = c;
      Polar_Decompose(cU.View(),cP.Transpose());
      Assert(Norm(c-cU*cP.Transpose()) < ceps*Norm(c),"C Polar3");
      Assert(Norm(cU.Adjoint()*cU-T(1)) < ceps,"C Polar3 UtU");

      cU = c;
      Polar_Decompose(cU.View(),cP.Conjugate());
      Assert(Norm(c-cU*cP.Transpose()) < ceps*Norm(c),"C Polar4");
      Assert(Norm(cU.Adjoint()*cU-T(1)) < ceps,"C Polar4 UtU");

      cU.Conjugate() = c;
      Polar_Decompose(cU.Conjugate(),cP.View());
      Assert(Norm(c-cU.Conjugate()*cP) < ceps*Norm(c),"C Polar5");
      Assert(Norm(cU.Adjoint()*cU-T(1)) < ceps,"C Polar5 UtU");

      cU.Conjugate() = c;
      Polar_Decompose(cU.Conjugate(),cP.Adjoint());
      Assert(Norm(c-cU.Conjugate()*cP) < ceps*Norm(c),"C Polar6");
      Assert(Norm(cU.Adjoint()*cU-T(1)) < ceps,"C Polar6 UtU");

      cU.Conjugate() = c;
      Polar_Decompose(cU.Conjugate(),cP.Transpose());
      Assert(Norm(c-cU.Conjugate()*cP.Transpose()) < ceps*Norm(c),"C Polar7");
      Assert(Norm(cU.Adjoint()*cU-T(1)) < ceps,"C Polar7 UtU");

      cU.Conjugate() = c;
      Polar_Decompose(cU.Conjugate(),cP.Conjugate());
      Assert(Norm(c-cU.Conjugate()*cP.Transpose()) < ceps*Norm(c),"C Polar8");
      Assert(Norm(cU.Adjoint()*cU-T(1)) < ceps,"C Polar8 UtU");
#endif
      std::cout<<"."; std::cout.flush();
    }

    // BandMatrix Polar Decomposition
    {
      tmv::BandMatrixView<T> b(m.View(),5,11);
      tmv::BandMatrixView<std::complex<T> > cb(c.View(),5,11);

      tmv::HermMatrix<T,tmv::Lower,tmv::ColMajor> P(N);
      tmv::Matrix<T,stor> U(b.colsize(),N);
      Polar_Decompose(b,U.View(),P.View());
      //std::cout<<"Norm(b-UP) = "<<Norm(b-U*P)<<std::endl;
      //std::cout<<"eps*Norm(b) = "<<eps*Norm(b)<<std::endl;
      Assert(Norm(b-U*P) < eps*Norm(b),"Band Polar");
      Assert(Norm(U.Adjoint()*U-T(1)) < eps,"Band Polar UtU");

      tmv::HermMatrix<std::complex<T>,tmv::Lower,tmv::ColMajor> cP(N);
      tmv::Matrix<std::complex<T>,stor> cU(cb.colsize(),N);
      Polar_Decompose(cb,cU.View(),cP.View());
      Assert(Norm(cb-cU*cP) < ceps*Norm(cb),"C Band Polar");
      Assert(Norm(cU.Adjoint()*cU-T(1)) < ceps,"C Band Polar UtU");

#ifdef XTEST
      Polar_Decompose(b,U.View(),P.Transpose());
      Assert(Norm(b-U*P) < eps*Norm(b),"Band Polar2");
      Assert(Norm(U.Adjoint()*U-T(1)) < eps,"Band Polar2 UtU");

      Polar_Decompose(cb,cU.View(),cP.Adjoint());
      Assert(Norm(cb-cU*cP) < ceps*Norm(cb),"C Band Polar2");
      Assert(Norm(cU.Adjoint()*cU-T(1)) < ceps,"C Band Polar2 UtU");

      Polar_Decompose(cb,cU.View(),cP.Transpose());
      Assert(Norm(cb-cU*cP.Transpose()) < ceps*Norm(cb),"C Band Polar3");
      Assert(Norm(cU.Adjoint()*cU-T(1)) < ceps,"C Band Polar3 UtU");

      Polar_Decompose(cb,cU.View(),cP.Conjugate());
      Assert(Norm(cb-cU*cP.Transpose()) < ceps*Norm(cb),"C Band Polar4");
      Assert(Norm(cU.Adjoint()*cU-T(1)) < ceps,"C Band Polar4 UtU");

      Polar_Decompose(cb,cU.Conjugate(),cP.View());
      Assert(Norm(cb-cU.Conjugate()*cP) < ceps*Norm(cb),"C Band Polar5");
      Assert(Norm(cU.Adjoint()*cU-T(1)) < ceps,"C Band Polar5 UtU");

      Polar_Decompose(cb,cU.Conjugate(),cP.Adjoint());
      Assert(Norm(cb-cU.Conjugate()*cP) < ceps*Norm(cb),"C Band Polar6");
      Assert(Norm(cU.Adjoint()*cU-T(1)) < ceps,"C Band Polar6 UtU");

      Polar_Decompose(cb,cU.Conjugate(),cP.Transpose());
      Assert(Norm(cb-cU.Conjugate()*cP.Transpose()) < ceps*Norm(cb),"C Band Polar7");
      Assert(Norm(cU.Adjoint()*cU-T(1)) < ceps,"C Band Polar7 UtU");

      Polar_Decompose(cb,cU.Conjugate(),cP.Conjugate());
      Assert(Norm(cb-cU.Conjugate()*cP.Transpose()) < ceps*Norm(cb),"C Band Polar8");
      Assert(Norm(cU.Adjoint()*cU-T(1)) < ceps,"C Band Polar8 UtU");
#endif
      std::cout<<"."; std::cout.flush();
    }
  }
}

template <class T> void TestAllSymDiv()
{
  TestHermDecomp<T,tmv::Upper,tmv::ColMajor>();
  TestHermDecomp<T,tmv::Upper,tmv::RowMajor>();
  TestHermDecomp<T,tmv::Lower,tmv::ColMajor>();
  TestHermDecomp<T,tmv::Lower,tmv::RowMajor>();
  TestSymDecomp<T,tmv::Upper,tmv::ColMajor>();
  TestSymDecomp<T,tmv::Upper,tmv::RowMajor>();
  TestSymDecomp<T,tmv::Lower,tmv::ColMajor>();
  TestSymDecomp<T,tmv::Lower,tmv::RowMajor>();
  TestPolar<T,tmv::RowMajor>();
  TestPolar<T,tmv::ColMajor>();
  std::cout<<"SymMatrix<"<<tmv::TypeText(T())<<"> passed all ";
  std::cout<<"decomposition tests.\n";
  TestSymDiv<T>(tmv::CH,PosDef);
  TestSymDiv<T>(tmv::LU,PosDef);
  TestSymDiv<T>(tmv::LU,InDef);
  TestSymDiv<T>(tmv::SV,PosDef);
  TestSymDiv<T>(tmv::SV,InDef);
  TestSymDiv<T>(tmv::SV,Sing);
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



