
#include "TMV.h"
#include "TMV_Householder.h"
#include "TMV_Givens.h"
#include "TMV_Diag.h"

namespace tmv {

  //
  // LDiv
  //

  template <class T1, class T2, class T3> void SV_LDiv(
      const GenMatrix<T1>& U, const GenVector<RealType(T1)>& S, 
      const GenMatrix<T1>& V, size_t kmax,
      const GenMatrix<T2>& m, const MatrixView<T3>& x)
  {
    // A x = m
    // U S V x = m
    // x = Vt S^-1 Ut m
    TMVAssert(m.colsize() == U.colsize());
    TMVAssert(x.colsize() == V.rowsize());
    TMVAssert(x.rowsize() == m.rowsize());
    TMVAssert(kmax <= V.rowsize());
    TMVAssert(kmax <= U.colsize());
    x.Rows(0,kmax) = U.QuickAdjoint().Rows(0,kmax) * m;
    x.Rows(kmax,x.colsize()).Zero();
    x.Rows(0,kmax) /= DiagMatrixViewOf(S.SubVector(0,kmax));
    x = V.QuickAdjoint().Cols(0,kmax) * x.Rows(0,kmax);
  }

  //
  // RDiv
  //

  template <class T1, class T2, class T3> void SV_RDiv(
      const GenMatrix<T1>& U, const GenVector<RealType(T1)>& S, 
      const GenMatrix<T1>& V, size_t kmax,
      const GenMatrix<T2>& m, const MatrixView<T3>& x) 
  {
    // x A = m
    // x U S V = m
    // x = m Vt S^-1 Ut
    TMVAssert(m.rowsize() == V.rowsize());
    TMVAssert(x.rowsize() == U.colsize());
    TMVAssert(x.colsize() == m.colsize());
    TMVAssert(kmax <= U.colsize());
    TMVAssert(kmax <= V.rowsize());
    x.Cols(0,kmax) = m * V.QuickAdjoint().Cols(0,kmax);
    x.Cols(kmax,x.rowsize()).Zero();
    x.Cols(0,kmax) %= DiagMatrixViewOf(S.SubVector(0,kmax));
    x = x.Cols(0,kmax) * U.QuickAdjoint().Rows(0,kmax);
  }

  template <class T> void NonLapBidiagonalize(
      const MatrixView<T>& U, const VectorView<RealType(T)>& D,
      const VectorView<RealType(T)>& E, const MatrixView<T>& V,
      T& detuv)
  {
    // Decompose A (input as U) into U B V
    // The Bidiagonal Matrix B is stored as two vectors: D, E
    // D is the diagonal, E is the super-diagonal
    // We use Householder reflections to reduce A to the bidiagonal form:
    const size_t M = U.colsize();
    const size_t N = U.rowsize();

    TMVAssert(N <= M);
    TMVAssert(V.colsize() == N);
    TMVAssert(V.rowsize() == N);
    TMVAssert(D.size() == N);
    TMVAssert(N > 0);
    TMVAssert(E.size() == N-1);

    Vector<T> Ubeta(N);
    Vector<T> Vbeta(N-1);
    for(size_t j=0;j<N-1;++j) {
      Ubeta(j) = Householder_Reflect(U.SubMatrix(j,M,j,N),detuv);
      Vbeta(j) = Householder_Reflect(U.QuickTranspose().SubMatrix(j+1,N,j,M),
	  detuv);
    }
    Ubeta(N-1) = Householder_Reflect(U.SubMatrix(N-1,M,N-1,N),detuv);

    // The bidiagonal of U is the bidiagonal we want, so copy it to D,E
    if (IsComplex(T())) {
      TMVAssert(NormInf(U.diag().Imag()) == RealType(T)(0));
      TMVAssert(NormInf(U.diag(1).Imag()) == RealType(T)(0));
    }
    D = U.diag().Real();
    E = U.diag(1).Real();
    U.diag().Zero();
    U.diag(1).Zero();

    // Now U stores Householder vectors for U in lower diagonal columns (HLi)
    // and Householder vectors for V in upper diagonal rows (HRi)
    // The Householder matrices for U are actually the adjoints of the 
    // matrices that bidiagonalize A, and for V are the transposes:
    // B = HLn-1t ... HL1t HL0t A HR0T HR1T ... HRn-2T
    // Using the fact that H Ht = I, we get A = U B V with:
    // U = HL0 ... HLn-1  and  V = HRn-2* ... HR0*
    // Or, VT = HR0t ... HRn-2t
    V.SetToIdentity();
    Householder_Unpack(U.SubMatrix(N-1,M,N-1,N),Ubeta(N-1));
    for (int j=N-2;j>=0;--j) {
      V.row(j+1,j+2,N) = U.row(j,j+2,N);
      Householder_Unpack(V.QuickTranspose().SubMatrix(j+1,N,j+1,N),Vbeta(j));
      U.row(j,j+2,N).Zero();
      Householder_Unpack(U.SubMatrix(j,M,j,N),Ubeta(j));
    }
  }

#ifdef LAP
  template <class T> void LapBidiagonalize(
      const MatrixView<T>& U, const VectorView<RealType(T)>& D,
      const VectorView<RealType(T)>& E, const MatrixView<T>& V, 
      T& detuv)
  { NonLapBidiagonalize(U,D,E,V,detuv); }
  template <> void LapBidiagonalize(
      const MatrixView<double>& U, const VectorView<double>& D,
      const VectorView<double>& E, const MatrixView<double>& V, 
      double& detuv)
  {
    TMVAssert(U.iscm());
    TMVAssert(V.iscm());
    TMVAssert(D.step() == 1);
    TMVAssert(E.step() == 1);
    TMVAssert(U.ct()==NonConj);
    TMVAssert(D.ct()==NonConj);
    TMVAssert(E.ct()==NonConj);
    TMVAssert(V.ct()==NonConj);
    int m = U.colsize();
    int n = U.rowsize();
    int ldu = U.stepj();
    Vector<double> tauq(n);
    Vector<double> taup(n);
    int lwork = (m+n)*LAP_BLOCKSIZE;
    double* work = LAP_DWork(lwork);
    int info;
    dgebrd(&m,&n,U.ptr(),&ldu,D.ptr(),E.ptr(),tauq.ptr(),taup.ptr(),
	work,&lwork,&info);
    LAP_Results(info,int(work[0]),m,n,lwork,"dgebrd");
    for(size_t i=0;i<tauq.size();++i) if (tauq(i) != double(0)) detuv = -detuv;
    for(size_t i=0;i<taup.size()-1;++i) if (taup(i) != double(0)) detuv = -detuv;
    for(size_t j=1;j<U.rowsize();++j) V.col(j,0,j) = U.col(j,0,j);
    char c = 'P';
    int ldv = V.stepj();
    dorgbr(&c,&n,&n,&n,V.ptr(),&ldv,taup.ptr(),work,&lwork,&info);
    LAP_Results(info,int(work[0]),m,n,lwork,"dorgbr");
    c = 'Q';
    dorgbr(&c,&m,&n,&n,U.ptr(),&ldu,tauq.ptr(),work,&lwork,&info);
    LAP_Results(info,int(work[0]),m,n,lwork,"dorgbr");
  }
  template <> void LapBidiagonalize(
      const MatrixView<complex<double> >& U, 
      const VectorView<double>& D, const VectorView<double>& E,
      const MatrixView<complex<double> >& V, 
      complex<double>& detuv)
  {
    TMVAssert(U.iscm());
    TMVAssert(V.iscm());
    TMVAssert(D.step() == 1);
    TMVAssert(E.step() == 1);
    TMVAssert(U.ct()==NonConj);
    TMVAssert(D.ct()==NonConj);
    TMVAssert(E.ct()==NonConj);
    TMVAssert(V.ct()==NonConj);
    int m = U.colsize();
    int n = U.rowsize();
    int ldu = U.stepj();
    Vector<complex<double> > tauq(n);
    Vector<complex<double> > taup(n);
    int lwork = (m+n)*LAP_BLOCKSIZE;
    complex<double>* work = LAP_ZWork(lwork);
    int info;
    zgebrd(&m,&n,LAP_Complex(U.ptr()),&ldu,D.ptr(),E.ptr(),
	LAP_Complex(tauq.ptr()),LAP_Complex(taup.ptr()),
	LAP_Complex(work),&lwork,&info);
    LAP_Results(info,int(real(work[0])),m,n,lwork,"zgebrd");
    for(size_t i=0;i<tauq.size();++i) if (tauq(i) != double(0)) {
      detuv *= (tauq(i)*tauq(i))/norm(tauq(i));
    }
    for(size_t i=0;i<taup.size()-1;++i) if (taup(i) != double(0)) {
      detuv *= conj(taup(i)*taup(i))/norm(taup(i));
    }
    for(size_t i=0;i<taup.size()-1;++i) if (taup(i) != double(0)) detuv = -detuv;
    for(size_t j=1;j<U.rowsize();++j) V.col(j,0,j) = U.col(j,0,j);
    char c = 'P';
    int ldv = V.stepj();
    zungbr(&c,&n,&n,&n,LAP_Complex(V.ptr()),&ldv,LAP_Complex(taup.ptr()),
	LAP_Complex(work),&lwork,&info);
    LAP_Results(info,int(real(work[0])),m,n,lwork,"zungbr");
    c = 'Q';
    zungbr(&c,&m,&n,&n,LAP_Complex(U.ptr()),&ldu,LAP_Complex(tauq.ptr()),
	LAP_Complex(work),&lwork,&info);
    LAP_Results(info,int(real(work[0])),m,n,lwork,"zungbr");
  }
#ifndef NOFLOAT
  template <> void LapBidiagonalize(
      const MatrixView<float>& U, const VectorView<float>& D,
      const VectorView<float>& E, const MatrixView<float>& V, 
      float& detuv)
  {
    TMVAssert(U.iscm());
    TMVAssert(V.iscm());
    TMVAssert(D.step() == 1);
    TMVAssert(E.step() == 1);
    TMVAssert(U.ct()==NonConj);
    TMVAssert(D.ct()==NonConj);
    TMVAssert(E.ct()==NonConj);
    TMVAssert(V.ct()==NonConj);
    int m = U.colsize();
    int n = U.rowsize();
    int ldu = U.stepj();
    Vector<float> tauq(n);
    Vector<float> taup(n);
    int lwork = (m+n)*LAP_BLOCKSIZE;
    float* work = LAP_SWork(lwork);
    int info;
    sgebrd(&m,&n,U.ptr(),&ldu,D.ptr(),E.ptr(),tauq.ptr(),taup.ptr(),
	work,&lwork,&info);
    LAP_Results(info,int(work[0]),m,n,lwork,"sgebrd");
    for(size_t i=0;i<tauq.size();++i) if (tauq(i) != float(0)) detuv = -detuv;
    for(size_t i=0;i<taup.size()-1;++i) if (taup(i) != float(0)) detuv = -detuv;
    for(size_t j=1;j<U.rowsize();++j) V.col(j,0,j) = U.col(j,0,j);
    char c = 'P';
    int ldv = V.stepj();
    sorgbr(&c,&n,&n,&n,V.ptr(),&ldv,taup.ptr(),work,&lwork,&info);
    LAP_Results(info,int(work[0]),m,n,lwork,"sorgbr");
    c = 'Q';
    sorgbr(&c,&m,&n,&n,U.ptr(),&ldu,tauq.ptr(),work,&lwork,&info);
    LAP_Results(info,int(work[0]),m,n,lwork,"sorgbr");
  }
  template <> void LapBidiagonalize(
      const MatrixView<complex<float> >& U, 
      const VectorView<float>& D, const VectorView<float>& E,
      const MatrixView<complex<float> >& V, complex<float>& detuv)
  {
    TMVAssert(U.iscm());
    TMVAssert(V.iscm());
    TMVAssert(D.step() == 1);
    TMVAssert(E.step() == 1);
    TMVAssert(U.ct()==NonConj);
    TMVAssert(D.ct()==NonConj);
    TMVAssert(E.ct()==NonConj);
    TMVAssert(V.ct()==NonConj);
    int m = U.colsize();
    int n = U.rowsize();
    int ldu = U.stepj();
    Vector<complex<float> > tauq(n);
    Vector<complex<float> > taup(n);
    int lwork = (m+n)*LAP_BLOCKSIZE;
    complex<float>* work = LAP_CWork(lwork);
    int info;
    cgebrd(&m,&n,LAP_Complex(U.ptr()),&ldu,D.ptr(),E.ptr(),
	LAP_Complex(tauq.ptr()),LAP_Complex(taup.ptr()),
	LAP_Complex(work),&lwork,&info);
    LAP_Results(info,int(real(work[0])),m,n,lwork,"cgebrd");
    for(size_t i=0;i<tauq.size();++i) if (tauq(i) != float(0)) {
      detuv *= (tauq(i)*tauq(i))/norm(tauq(i));
    }
    for(size_t i=0;i<taup.size()-1;++i) if (taup(i) != float(0)) {
      detuv *= conj(taup(i)*taup(i))/norm(taup(i));
    }
    for(size_t i=0;i<taup.size()-1;++i) if (taup(i) != float(0)) detuv = -detuv;
    for(size_t j=1;j<U.rowsize();++j) V.col(j,0,j) = U.col(j,0,j);
    char c = 'P';
    int ldv = V.stepj();
    cungbr(&c,&n,&n,&n,LAP_Complex(V.ptr()),&ldv,LAP_Complex(taup.ptr()),
	LAP_Complex(work),&lwork,&info);
    LAP_Results(info,int(real(work[0])),m,n,lwork,"cungbr");
    c = 'Q';
    cungbr(&c,&m,&n,&n,LAP_Complex(U.ptr()),&ldu,LAP_Complex(tauq.ptr()),
	LAP_Complex(work),&lwork,&info);
    LAP_Results(info,int(real(work[0])),m,n,lwork,"cungbr");
  }
#endif
#endif
  template <class T> void Bidiagonalize(
      const MatrixView<T>& U, const VectorView<RealType(T)>& D,
      const VectorView<RealType(T)>& E, const MatrixView<T>& V, 
      T& detuv)
  {
    TMVAssert(U.colsize() >= U.rowsize());
    TMVAssert(U.rowsize() == D.size());
    TMVAssert(U.rowsize() == V.colsize());
    TMVAssert(V.rowsize() == V.colsize());
    TMVAssert(U.isrm() || U.iscm());
    TMVAssert(V.isrm() || V.iscm());
    TMVAssert(U.ct()==NonConj);
    TMVAssert(D.ct()==NonConj);
    TMVAssert(E.ct()==NonConj);
    TMVAssert(V.ct()==NonConj);
    
    if (U.rowsize() > 0) {
      TMVAssert(E.size() == D.size()-1);
#ifdef LAP
      if (U.iscm() && V.iscm()) 
	LapBidiagonalize(U,D,E,V,detuv);
      else {
#ifdef TMVDEBUG
	cout<<"LAP Bidiag: U,V are wrong:\n";
	cout<<"U isrm = "<<U.isrm()<<", step = "<<U.stepi()<<endl;
	cout<<"V isrm = "<<V.isrm()<<", step = "<<V.stepi()<<endl;
	cout<<"D.step = "<<D.step()<<", E.step = "<<E.step()<<endl;
	cout<<"U_input = "<<U<<endl;
#endif
	NonLapBidiagonalize(U,D,E,V,detuv);
      }
#else
      NonLapBidiagonalize(U,D,E,V,detuv);
#endif
    }
  }

  template <class T> void BidiagonalChopSmallElements(
      const VectorView<T>& D, const VectorView<T>& E)
  {
    // This routines sets to 0 any elements in E or D which
    // are essentially 0, given the machine precision:
    // if |E(i)| < Epsilon() * (|D(i)| + |D(i+1)|) then E(i) <- 0
    // if |D(i)| < Epsilon() * |B| then D(i) <- 0
    TMVAssert(E.size() == D.size()-1);
    TMVAssert(D.step() == 1);
    TMVAssert(E.step() == 1);
    RealType(T) eps = Epsilon<T>();
    RealType(T) dthresh = eps * SQRT(NormSq(D) + NormSq(E));

    VIt<T,Unit,NonConj> Di = D.begin();
    VIt<T,Unit,NonConj> Ei = E.begin();
    const VIt<T,Unit,NonConj> E_end = E.end();
    for(;Ei!=E_end;++Di,++Ei) {
      if (abs(*Di) < dthresh) *Di = T(0);
      if (abs(*Ei) < eps * (abs(*Di)+abs(*(Di+1)))) *Ei = T(0);
    }
    if (abs(*Di) < dthresh) *Di = T(0);
  }

  template <class T, class TU> void BidiagonalZeroFirstRow(
      const MatrixView<TU>& U, const VectorView<T>& D, 
      const VectorView<T>& E)
  {
    // Input D,E form a bidiagonal matrix with the first element of D = 0:
    // (eg. for N = 5)
    //     [ 0 x 0 0 0 ]
    //     [ 0 x x 0 0 ]
    // B = [ 0 0 x x 0 ]
    //     [ 0 0 0 x x ]
    //     [ 0 0 0 0 x ]
    // Zero out the first row maintaining the constancy of U B
    // using Givens transformations.
    const size_t N = D.size();
    if (N <= 1) return; 
    TMVAssert(E.size() == N-1);
    TMVAssert(U.rowsize() == N);
    TMVAssert(D.step() == 1);
    TMVAssert(E.step() == 1);

    VIt<T,Unit,NonConj> Di = D.begin();
    VIt<T,Unit,NonConj> Ei = E.begin();
    TMVAssert(*Di == T(0));

    T x = *Ei;
    if (x != T(0)) {
      *Ei = T(0);
      ++Ei; ++Di;
      // Loop Invariant: x = B(0,i)
      for(size_t i=1; i<N; ++i,++Di,++Ei) {
	Givens<T> G = Givens_Rotate(*Di,x);
	// Make new B = G B
	if (i<N-1) G.Mult(*Ei,x);
	// Make new U = U Gt
	G.ConjMult(U.ColPair(i,0).QuickTranspose());
      }
    }
  }

  template <class T, class TV> void BidiagonalZeroLastCol(
      const VectorView<T>& D, const VectorView<T>& E,
      const MatrixView<TV>& V)
  {
    // Input D,E form a bidiagonal matrix with the last element of D = 0:
    // (eg. for N = 5)
    //     [ x x 0 0 0 ]
    //     [ 0 x x 0 0 ]
    // B = [ 0 0 x x 0 ]
    //     [ 0 0 0 x x ]
    //     [ 0 0 0 0 0 ]
    // Zero out the last col maintaining the constancy of B V
    // using Givens transformations.
    const size_t N = D.size();
    if (N <= 1) return; 
    TMVAssert(E.size() == N-1);
    TMVAssert(V.colsize() == N);
    TMVAssert(D.step() == 1);
    TMVAssert(E.step() == 1);
    TMVAssert(D(N-1) == T(0));

    VIt<T,Unit,NonConj> Di = D.begin()+N-2;
    VIt<T,Unit,NonConj> Ei = E.begin()+N-2;

    T x = *Ei;
    if (x != T(0)) {
      *Ei = T(0);
      // Loop Invariant: x = B(i,N-1)
      for(int i=N-2; i>=0; --i,--Di) {
	Givens<T> G = Givens_Rotate(*Di,x);
	// Make new B = B GT
	if (i>0) G.Mult(*(--Ei),x);
	// Make new V = G* V 
	G.ConjMult(V.RowPair(i,N-1));
      }
    }
  }

  template <class T> RealType(T) BidiagonalTrailingEigenValue(
      const VectorView<T>& D, const VectorView<T>& E)
  {
    // Return the Wilkinson choice for an eigenvalue of T = BtB, namely
    // the eigenvalue of the trailing 2x2 block of T which is closer
    // to the last diagonal element of T.
    //
    // Trailing 2x2 block =  (Use i = N-2, j = N-1)
    // [ a  b ] = [ |Di|^2 + |Ei-1|^2      Di Ei      ]
    // [ b* c ]   [     Di* Ei*       |Dj|^2 + |Ei|^2 ]
    // 
    // mu = c - d +- sqrt(d^2+|b|^2), where d = (c-a)/2
    // if d>0 we use +, if d<0 we use -.
    // 
    // For stability when |b| is small, we rearrange this to:
    // mu = c + |b|^2/(d +- sqrt(d^2+|b|^2))
    //    = c +- |b|^2/|d|(1 + sqrt(1+|b|^2/d^2))
    const size_t N = D.size();
    TMVAssert(E.size() == N-1);
    TMVAssert(N > 1);

    RealType(T) a = NORM(D(N-2)) + (N>2 ? NORM(E(N-3)) : RealType(T)(0));
    RealType(T) c = NORM(D(N-1)) + NORM(E(N-2));
    T b = D(N-2)*E(N-2);
    RealType(T) bsq = NORM(b);
    RealType(T) d = (c-a)/2;
    RealType(T) dsq = SQR(d);
    if (dsq < bsq) {
      RealType(T) x = SQRT(dsq+bsq);
      if (d > 0) return c - d + x;
      else return c - d - x;
    } else {
      RealType(T) x = bsq/abs(d)/(1 + SQRT(1+bsq/dsq));
      if (d > 0) return c + x;
      else return c - x;
    }
  }

  template <class TUV, class T> void ReduceUnredBidiagonal(
      const MatrixView<TUV>& U, const VectorView<T>& D,
      const VectorView<T>& E, const MatrixView<TUV>& V)
  {
    // Reduce the superdiagonal elements of unreduced Bidiagonal Matrix B 
    // (given by D,E) while maintaining U B V. 
    // Note: the input B must be unreduced - ie. all entries are non-zero.
    const size_t N = D.size();
    TMVAssert(N>0);
    TMVAssert(E.size() == N-1);
    TMVAssert(U.rowsize() == N);
    TMVAssert(V.colsize() == N);
    TMVAssert(D.step()==1);
    TMVAssert(E.step()==1);

    if (N == 1) return;

    // The reduction is based on the QR algorithm to diagonalize the
    // unreduced symmetric tridiagonal matrix T = BtB
    // The basic idea is as follows:
    // (see Golub and van Loan, chapter 8 for a derivation)
    //
    // if T is a symmetric tridiagonal matrix
    // and mu is (approximately) an eigenvalue of T
    // and the QR decomposition of T - mu I = V R
    // Then, T' = R V + mu I will be tridiagonal with the last 
    // subdiagonal element small.
    // (Note: T' = Vt (T-muI) V + muI = Vt T V.)
    //
    // Wilkinson (1968) suggested that a good choice for mu is
    // the eigenvalue of the trailing 2x2 block of T that is 
    // closer to the trailing diagonal element.
    //
    // Rather than explicitly forming T = BtB and doing this
    // procedure, Golub and van Load show that it can be done
    // in place.
    // If T' = Vt T V, 
    // then Bt'B' = Vt Bt B V
    // B' = U B V for some U
    //
    // So, start with the first Givens matrix in the QR algorithm for T:
    // G0 [ T00 - mu ] = [ x ]
    //    [   T10    ]   [ 0 ]
    // We apply this to the right of B which has the effect: (for N=5)
    //              [ x x 0 0 0 ]
    //              [ + x x 0 0 ]
    // B <- B G0T = [ 0 0 x x 0 ]
    //              [ 0 0 0 x x ]
    //              [ 0 0 0 0 x ]
    // The + is the element which screws up the bidiagonal structure.
    // The rest of the procedure simply involves chasing this + down
    // the diagonal using Givens rotations.
    // For each Givens rotation we use, we also multiply U or V by the
    // adjoint to maintain the constancy of U B V.
    //
    // At the end of this procedure, E(N-1) should be smaller than it was.
    // Note: This procedure works exactly if N=2.
    VIt<T,Unit,NonConj> Di = D.begin();
    VIt<T,Unit,NonConj> Ei = E.begin();

    T mu = BidiagonalTrailingEigenValue(D,E);
    T y = NORM(*Di) - mu;  // = T00 - mu
    T x = CONJ(*Di)*(*Ei);  // = T10
    Givens<T> G = Givens_Rotate(y,x);
    for(size_t i=1;i<N;++i) {
      G.Mult(*Di,*Ei);
      G.ConjMult(V.RowPair(i-1,i));
      TMVAssert(x==T(0));
      G.Mult(x,*(++Di)); // x = B(i,i-1)
      G = Givens_Rotate(*(Di-1),x);
      G.Mult(*Ei,*Di);
      G.ConjMult(U.ColPair(i-1,i).QuickTranspose());
      if (i < N-1) {
	TMVAssert(x==T(0));
	G.Mult(x,*(++Ei)); // x = B(i-1,i+1)
	G = Givens_Rotate(*(Ei-1),x);
      } 
    }
  }

  template <class TUV, class T> void ReduceBidiagonal(
      const MatrixView<TUV>& U, const VectorView<T>& D,
      const VectorView<T>& E, const MatrixView<TUV>& V)
  {
    // Reduce the superdiagonal elements of Bidiagonal Matrix B 
    // (given by D,E) while maintaining U B V. 
    const size_t N = D.size();
    TMVAssert(E.size() == N-1);
    TMVAssert(U.rowsize() == N);
    TMVAssert(V.colsize() == N);

    // The input E(i) are all assumed to be non-zero.
    // If there are any zeros in D, we can zero the corresponding
    // E's (above and right) directly, so look for these first.
    // Loop invariant: all D(i) with p<=i<q are non-zero.
    size_t p=0; 
    for(size_t q=0; q<N; ++q) {
      if (D(q) == T(0)) {
	if (p<q) {
          BidiagonalZeroLastCol(D.SubVector(p,q+1),E.SubVector(p,q),
	      V.SubMatrix(p,q+1,0,V.rowsize()));
          ReduceUnredBidiagonal(U.SubMatrix(0,U.colsize(),p,q), 
	      D.SubVector(p,q), E.SubVector(p,q-1), 
	      V.SubMatrix(p,q,0,V.rowsize()));
	}
	if (q<N-1) {
          BidiagonalZeroFirstRow(U.SubMatrix(0,U.colsize(),q,N),
	      D.SubVector(q,N), E.SubVector(q,N-1));
	}
	p=q+1;
      }
    }
    if (p<N) {
      ReduceUnredBidiagonal(U.SubMatrix(0,U.colsize(),p,N), D.SubVector(p,N),
 	  E.SubVector(p,N-1), V.SubMatrix(p,N,0,V.rowsize()));
    }
  }

  template <class T> void NonLapSV_Decompose_From_Bidiagonal(
      const MatrixView<T>& U, const VectorView<RealType(T)>& D, 
      const VectorView<RealType(T)>& E, const MatrixView<T>& V)
  {
    TMVAssert(U.colsize() >= U.rowsize());
    TMVAssert(U.rowsize() == D.size());
    TMVAssert(U.rowsize() == E.size()+1);
    TMVAssert(U.rowsize() == V.colsize());
    TMVAssert(V.rowsize() == V.colsize());

    const size_t N = U.rowsize();
    const size_t M = U.colsize();

    // We successively reduce the superdiagonal of B (E) to 0
    // using a sequence of Givens rotations. 
    // The reduction procedure tends to push the values up and left, so it 
    // makes sense to start at the lower right and work back up the matrix.
    // We also set to zero any very small values based on machine precision.
    // Loop invariant: all E(i) with i>=q are 0.
    // Initially q = N-1. (ie. All E(i) are potentially non-zero.)
    // When q = 0, we are done.
    BidiagonalChopSmallElements(D,E);
    for(size_t q = N-1; q>0; ) {
      if (E(q-1) == T(0)) --q;
      else {
	size_t p=q-1;
	while (p > 0 && (E(p-1) != T(0))) --p; 
	// Set p such that E(p-1) = 0 and all E(i) with p<=i<q are non-zero.
	ReduceBidiagonal(U.SubMatrix(0,M,p,q+1),D.SubVector(p,q+1),
	    E.SubVector(p,q),V.SubMatrix(p,q+1,0,N));
	BidiagonalChopSmallElements(D,E);
      }
    }

    // Make all of the singular values positive
    VIt<RealType(T),Unit,NonConj> Di = D.begin();
    for(size_t i=0;i<N;++i,++Di) if (*Di < 0) {
      *Di = -(*Di);
      V.row(i) = -V.row(i);
    }

    // Now A = U * S * V
    // Sort output singular values 
    Permutation sortp = SortPermutation(D,DESCEND);
    D = sortp * D;
    U = U * sortp.Transpose();
    V = sortp * V;
  }

#ifdef XLAP // Old Lap version
  template <class T> void LapSV_Decompose_From_Bidiagonal(
      const MatrixView<T>& U, const VectorView<RealType(T)>& D, 
      const VectorView<RealType(T)>& E, const MatrixView<T>& V)
  { NonLapSV_Decompose_From_Bidiagonal(U,D,E,V); }
  template <> void LapSV_Decompose_From_Bidiagonal(
      const MatrixView<double>& U, const VectorView<double>& D, 
      const VectorView<double>& E, const MatrixView<double>& V)
  {
    TMVAssert(U.colsize() >= U.rowsize());
    TMVAssert(U.rowsize() == D.size());
    TMVAssert(U.rowsize() == E.size()+1);
    TMVAssert(U.rowsize() == V.colsize());
    TMVAssert(V.rowsize() == V.colsize());
    TMVAssert(U.iscm());
    TMVAssert(V.iscm());
    TMVAssert(U.ct()==NonConj);
    TMVAssert(D.ct()==NonConj);
    TMVAssert(E.ct()==NonConj);
    TMVAssert(V.ct()==NonConj);
    char u = 'U';
    int n = D.size();
    int m = U.colsize();
    int o = 0;
    int ldv = V.stepj();
    int ldu = U.stepj();
    int lwork = 4*n;
    double* work = LAP_DWork(lwork);
    int info;
    dbdsqr(&u,&n,&n,&m,&o,D.ptr(),E.ptr(),V.ptr(),&ldv,U.ptr(),&ldu,0,&n,
	work,&info);
    LAP_Results(info,"dbdsqr");
  }
  template <> void LapSV_Decompose_From_Bidiagonal(
      const MatrixView<complex<double> >& U, 
      const VectorView<double>& D, const VectorView<double>& E,
      const MatrixView<complex<double> >& V)
  {
    TMVAssert(U.colsize() >= U.rowsize());
    TMVAssert(U.rowsize() == D.size());
    TMVAssert(U.rowsize() == E.size()+1);
    TMVAssert(U.rowsize() == V.colsize());
    TMVAssert(V.rowsize() == V.colsize());
    TMVAssert(U.iscm());
    TMVAssert(V.iscm());
    TMVAssert(U.ct()==NonConj);
    TMVAssert(D.ct()==NonConj);
    TMVAssert(E.ct()==NonConj);
    TMVAssert(V.ct()==NonConj);
    char u = 'U';
    int n = D.size();
    int m = U.colsize();
    int o = 0;
    int ldv = V.stepj();
    int ldu = U.stepj();
    int lwork = 4*n;
    double* work = LAP_DWork(lwork);
    int info;
    zbdsqr(&u,&n,&n,&m,&o,D.ptr(),E.ptr(),LAP_Complex(V.ptr()),&ldv,
	LAP_Complex(U.ptr()),&ldu,0,&n,work,&info);
    LAP_Results(info,"zbdsqr");
  }
#ifndef NOFLOAT
  template <> void LapSV_Decompose_From_Bidiagonal(
      const MatrixView<float>& U, const VectorView<float>& D, 
      const VectorView<float>& E, const MatrixView<float>& V)
  {
    TMVAssert(U.colsize() >= U.rowsize());
    TMVAssert(U.rowsize() == D.size());
    TMVAssert(U.rowsize() == E.size()+1);
    TMVAssert(U.rowsize() == V.colsize());
    TMVAssert(V.rowsize() == V.colsize());
    TMVAssert(U.iscm());
    TMVAssert(V.iscm());
    TMVAssert(U.ct()==NonConj);
    TMVAssert(D.ct()==NonConj);
    TMVAssert(E.ct()==NonConj);
    TMVAssert(V.ct()==NonConj);
    char u = 'U';
    int n = D.size();
    int m = U.colsize();
    int o = 0;
    int ldv = V.stepj();
    int ldu = U.stepj();
    int lwork = 4*n;
    float* work = LAP_SWork(lwork);
    int info;
    sbdsqr(&u,&n,&n,&m,&o,D.ptr(),E.ptr(),V.ptr(),&ldv,U.ptr(),&ldu,0,&n,
	work,&info);
    LAP_Results(info,"sbdsqr");
  }
  template <> void LapSV_Decompose_From_Bidiagonal(
      const MatrixView<complex<float> >& U, 
      const VectorView<float>& D, const VectorView<float>& E,
      const MatrixView<complex<float> >& V)
  {
    TMVAssert(U.colsize() >= U.rowsize());
    TMVAssert(U.rowsize() == D.size());
    TMVAssert(U.rowsize() == E.size()+1);
    TMVAssert(U.rowsize() == V.colsize());
    TMVAssert(V.rowsize() == V.colsize());
    TMVAssert(U.iscm());
    TMVAssert(V.iscm());
    TMVAssert(U.ct()==NonConj);
    TMVAssert(D.ct()==NonConj);
    TMVAssert(E.ct()==NonConj);
    TMVAssert(V.ct()==NonConj);
    char u = 'U';
    int n = D.size();
    int m = U.colsize();
    int o = 0;
    int ldv = V.stepj();
    int ldu = U.stepj();
    int lwork = 4*n;
    float* work = LAP_SWork(lwork);
    int info;
    cbdsqr(&u,&n,&n,&m,&o,D.ptr(),E.ptr(),LAP_Complex(V.ptr()),&ldv,
	LAP_Complex(U.ptr()),&ldu,0,&n,work,&info);
    LAP_Results(info,"cbdsqr");
  }
#endif
#endif // OLD
#ifdef LAP // New LAP SVD using divide and conquer technique:
  template <class T> void LapSV_Decompose_From_Bidiagonal(
      const MatrixView<T>& U, const VectorView<RealType(T)>& D, 
      const VectorView<RealType(T)>& E, const MatrixView<T>& V)
  { NonLapSV_Decompose_From_Bidiagonal(U,D,E,V); }
  template <> void LapSV_Decompose_From_Bidiagonal(
      const MatrixView<double>& U, const VectorView<double>& D, 
      const VectorView<double>& E, const MatrixView<double>& V)
  {
    TMVAssert(U.colsize() >= U.rowsize());
    TMVAssert(U.rowsize() == D.size());
    TMVAssert(U.rowsize() == E.size()+1);
    TMVAssert(U.rowsize() == V.colsize());
    TMVAssert(V.rowsize() == V.colsize());
    TMVAssert(U.iscm());
    TMVAssert(V.iscm());
    TMVAssert(U.ct()==NonConj);
    TMVAssert(D.ct()==NonConj);
    TMVAssert(E.ct()==NonConj);
    TMVAssert(V.ct()==NonConj);
    char u = 'U';
    char c = 'I';
    int n = D.size();
    int lwork = (3*n+4)*n;
    double* work = LAP_DWork(lwork);
    lwork = 8*n;
    int* iwork = LAP_IWork(lwork);
    Matrix<double,ColMajor> U1(n,n);
    Matrix<double,ColMajor> V1(n,n);
    int ldu = U1.stepj();
    int ldv = V1.stepj();
    int info;
    dbdsdc(&u,&c,&n,D.ptr(),E.ptr(),U1.ptr(),&ldu,V1.ptr(),&ldv,0,0,
	work,iwork,&info);
    LAP_Results(info,"dbdsdc");
    U = U*U1;
    V = V1*V;
  }
  template <> void LapSV_Decompose_From_Bidiagonal(
      const MatrixView<complex<double> >& U, 
      const VectorView<double>& D, const VectorView<double>& E,
      const MatrixView<complex<double> >& V)
  {
    TMVAssert(U.colsize() >= U.rowsize());
    TMVAssert(U.rowsize() == D.size());
    TMVAssert(U.rowsize() == E.size()+1);
    TMVAssert(U.rowsize() == V.colsize());
    TMVAssert(V.rowsize() == V.colsize());
    TMVAssert(U.iscm());
    TMVAssert(V.iscm());
    TMVAssert(U.ct()==NonConj);
    TMVAssert(D.ct()==NonConj);
    TMVAssert(E.ct()==NonConj);
    TMVAssert(V.ct()==NonConj);
    char u = 'U';
    char c = 'I';
    int n = D.size();
    int lwork = (3*n+4)*n;
    double* work = LAP_DWork(lwork);
    lwork = 8*n;
    int* iwork = LAP_IWork(lwork);
    Matrix<double,ColMajor> U1(n,n);
    Matrix<double,ColMajor> V1(n,n);
    int ldu = U1.stepj();
    int ldv = V1.stepj();
    int info;
    dbdsdc(&u,&c,&n,D.ptr(),E.ptr(),U1.ptr(),&ldu,V1.ptr(),&ldv,0,0,
	work,iwork,&info);
    LAP_Results(info,"dbdsdc");
    U = U*U1;
    V = V1*V;
  }
#ifndef NOFLOAT
  template <> void LapSV_Decompose_From_Bidiagonal(
      const MatrixView<float>& U, const VectorView<float>& D, 
      const VectorView<float>& E, const MatrixView<float>& V)
  {
    TMVAssert(U.colsize() >= U.rowsize());
    TMVAssert(U.rowsize() == D.size());
    TMVAssert(U.rowsize() == E.size()+1);
    TMVAssert(U.rowsize() == V.colsize());
    TMVAssert(V.rowsize() == V.colsize());
    TMVAssert(U.iscm());
    TMVAssert(V.iscm());
    TMVAssert(U.ct()==NonConj);
    TMVAssert(D.ct()==NonConj);
    TMVAssert(E.ct()==NonConj);
    TMVAssert(V.ct()==NonConj);
    char u = 'U';
    char c = 'I';
    int n = D.size();
    int lwork = (3*n+4)*n;
    float* work = LAP_SWork(lwork);
    lwork = 8*n;
    int* iwork = LAP_IWork(lwork);
    Matrix<float,ColMajor> U1(n,n);
    Matrix<float,ColMajor> V1(n,n);
    int ldu = U1.stepj();
    int ldv = V1.stepj();
    int info;
    sbdsdc(&u,&c,&n,D.ptr(),E.ptr(),U1.ptr(),&ldu,V1.ptr(),&ldv,0,0,
	work,iwork,&info);
    LAP_Results(info,"sbdsdc");
    U = U*U1;
    V = V1*V;
  }
  template <> void LapSV_Decompose_From_Bidiagonal(
      const MatrixView<complex<float> >& U, 
      const VectorView<float>& D, const VectorView<float>& E,
      const MatrixView<complex<float> >& V)
  {
    TMVAssert(U.colsize() >= U.rowsize());
    TMVAssert(U.rowsize() == D.size());
    TMVAssert(U.rowsize() == E.size()+1);
    TMVAssert(U.rowsize() == V.colsize());
    TMVAssert(V.rowsize() == V.colsize());
    TMVAssert(U.iscm());
    TMVAssert(V.iscm());
    TMVAssert(U.ct()==NonConj);
    TMVAssert(D.ct()==NonConj);
    TMVAssert(E.ct()==NonConj);
    TMVAssert(V.ct()==NonConj);
    char u = 'U';
    char c = 'I';
    int n = D.size();
    int lwork = (3*n+4)*n;
    float* work = LAP_SWork(lwork);
    lwork = 8*n;
    int* iwork = LAP_IWork(lwork);
    Matrix<float,ColMajor> U1(n,n);
    Matrix<float,ColMajor> V1(n,n);
    int ldu = U1.stepj();
    int ldv = V1.stepj();
    int info;
    sbdsdc(&u,&c,&n,D.ptr(),E.ptr(),U1.ptr(),&ldu,V1.ptr(),&ldv,0,0,
	work,iwork,&info);
    LAP_Results(info,"sbdsdc");
    U = U*U1;
    V = V1*V;
  }
#endif
#endif // LAP

  template <class T> void SV_Decompose_From_Bidiagonal(
      const MatrixView<T>& U, const VectorView<RealType(T)>& D, 
      const VectorView<RealType(T)>& E, const MatrixView<T>& V)
  {
    TMVAssert(U.colsize() >= U.rowsize());
    TMVAssert(U.rowsize() == D.size());
    TMVAssert(U.rowsize() == V.colsize());
    TMVAssert(V.rowsize() == V.colsize());
    TMVAssert(U.ct()==NonConj);
    TMVAssert(D.ct()==NonConj);
    TMVAssert(E.ct()==NonConj);
    TMVAssert(V.ct()==NonConj);
    
    if (U.rowsize() > 0) {
      TMVAssert(E.size() == D.size()-1);
#ifdef LAP
      if (U.iscm() && D.step() == 1 && E.step() == 1 && V.iscm()) {
	LapSV_Decompose_From_Bidiagonal(U,D,E,V);
      }
      else {
#ifdef TMVDEBUG
	cout<<"LAP SVDecomp_from_Bidiag: U,V are wrong:\n";
	cout<<"U isrm = "<<U.isrm()<<", step = "<<U.stepi()<<endl;
	cout<<"V isrm = "<<V.isrm()<<", step = "<<V.stepj()<<endl;
	cout<<"D.step = "<<D.step()<<", E.step = "<<E.step()<<endl;
	cout<<"A = "<<U<<endl;
#endif
	NonLapSV_Decompose_From_Bidiagonal(U,D,E,V);
      }
#else
      NonLapSV_Decompose_From_Bidiagonal(U,D,E,V);
#endif
    }
  }

  template <class T> void SV_Decompose(
      const MatrixView<T>& U, const VectorView<RealType(T)>& S, 
      const MatrixView<T>& V, T& det)
  {
    // Decompose A (input as U) into U S V
    // where S is a diagonal real matrix, and U,V are unitary matrices.
    // A,U are M x N (M >= N)
    // S,V are N x N
    // The determinant is returned in det.
    // (Technically, det is multiplied by the determinant, so det should
    // be set to 1 on entry.)
    const size_t M = U.colsize();
    const size_t N = U.rowsize();

    TMVAssert(N <= M);
    TMVAssert(V.colsize() == N);
    TMVAssert(V.rowsize() == N);
    TMVAssert(S.size() == N);

    // If M is much larger than N (technically M > 5/3 N), then it is quicker
    // to start by doing a QR decomposition and then do SVD on the square
    // R matrix.  Thus, the final U of the SVD is Q (from the QR decomp)
    // times U from R's SVD.
    if (M > 5*N/3) {
      Matrix<T,ColMajor> R(N,N);
      QR_Decompose(U,R.QuickView(),det);
      SV_Decompose(R.QuickView(),S,V,det);
      // Now R is a Unitary Matrix U'.  Need to multiply U by U'
      U = U*R;
      return;
    }

    // First we reduce A to bidiagonal form: A = U * B * V
    // using a series of Householder transformations.
    // The diagonal of the Bidiagonal Matrix B is stored in D.
    // The superdiagonal is stored in E.
    Vector<RealType(T)> E(N-1);
    Bidiagonalize(U,S,E.View(),V,det);

    // The determinant of B is just the product of the diagonal elements:
    det *= DiagMatrixViewOf(S).Det();

    // The rest of the procedure we separate out, since Decompositions of some
    // sparse matrices use the same procedure from here on out.
    SV_Decompose_From_Bidiagonal(U,S,E.View(),V);
  }

  template <class T> Matrix<T,ColMajor> SVDiv<T>::Inverse() const
  { 
    if (istrans) {
      // A^-1 = (Vt S^-1 Ut)T = U* S^-1 V*
      Matrix<T,ColMajor> SinvV = V.QuickConjugate().Rows(0,kmax) /
	DiagMatrixViewOf(S.SubVector(0,kmax));
      return U.QuickConjugate().Cols(0,kmax) * SinvV;
    } else {
      // A^-1 = Vt S^-1 Ut
      Matrix<T,ColMajor> SinvUt = U.QuickAdjoint().Rows(0,kmax) /
	DiagMatrixViewOf(S.SubVector(0,kmax));
      return V.QuickAdjoint().Cols(0,kmax) * SinvUt;
    }
  }

  template <class T> Matrix<T,ColMajor> SVDiv<T>::DoInverseATA() const
  {
    // A = U S V
    // At = Vt S Ut
    // AtA = Vt S^2 V
    // (AtA)^-1 = Vt S^-2 V
    //
    // if istrans:
    // AT = U S V
    // At = U* S V*
    // AAt = VT S^2 V*
    // (AAt)^-1 = VT S^-2 V*
    //
    Matrix<T,ColMajor> SinvV = V.Rows(0,kmax) /
      DiagMatrixViewOf(S.SubVector(0,kmax));
    if (istrans)
      return SinvV.QuickTranspose() * SinvV.QuickConjugate();
    else
      return SinvV.QuickAdjoint() * SinvV;
  }

  template <class T> void SVDiv<T>::Thresh(RealType(T) toler,
      ostream* debugout) const
  {
    TMVAssert(toler < RealType(T)(1));
    RealType(T) thresh = S(0)*toler;
    for(kmax=S.size(); kmax>0 && S(kmax-1)<=thresh; --kmax);
    if(debugout) {
      (*debugout)<<"S = "<<S<<endl;
      (*debugout)<<"Smax = "<<S(0)<<", thresh = "<<thresh<<std::endl;
      (*debugout)<<"kmax = "<<kmax<<" (S.size = "<<S.size()<<")"<<endl;
    }
  }

  template <class T> void SVDiv<T>::Top(size_t neigen,
      ostream* debugout) const
  {
    TMVAssert(neigen > 0);
    TMVAssert(neigen <= S.size());
    kmax = neigen;
    if(debugout) {
      (*debugout)<<"S = "<<S<<endl;
      (*debugout)<<"kmax = "<<kmax<<" (S.size = "<<S.size()<<")"<<endl;
    }
  }

#define InstFile "TMV_SVDiv.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


