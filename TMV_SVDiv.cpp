
#include "TMV.h"
#include "TMV_Householder.h"
#include "TMV_Givens.h"
#include "TMV_Diag.h"
#include "TMV_Tri.h"

namespace tmv {

  template <class T> const MatrixView<T>* ZMV() 
  { return (const MatrixView<T>*)(0); }

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
    TMVAssert(m.colsize() == U.colsize()); // = M
    TMVAssert(x.colsize() == V.rowsize()); // = N
    TMVAssert(x.rowsize() == m.rowsize()); // = R
    TMVAssert(kmax <= V.rowsize()); // = K
    TMVAssert(kmax <= U.colsize());
    Matrix<T3> m2 = U.QuickAdjoint().Rows(0,kmax) * m; // KxR
    m2 /= DiagMatrixViewOf(S.SubVector(0,kmax));
    x = V.QuickAdjoint().Cols(0,kmax) * m2; // NxR
  }

  template <class T1, class T2, class T3> void SV_LDiv(
      const GenMatrix<T1>& UV,
      const GenVector<T1>& Ubeta, const GenVector<T1>& Vbeta, 
      const GenMatrix<T1>* UV2, const GenVector<T1>* Qbeta,
      const GenMatrix<RealType(T1)>& U1, const GenMatrix<RealType(T1)>& V1,
      const GenVector<RealType(T1)>& S, size_t kmax,
      const GenMatrix<T2>& m, const MatrixView<T3>& x)
  {
    // A x = m
    // U U1 S V1 V x = m
    // x = Vt V1t S^-1 U1t Ut m
    TMVAssert(m.colsize() == UV.colsize()); // = M
    TMVAssert(x.colsize() == UV.rowsize()); // = N
    TMVAssert(x.rowsize() == m.rowsize()); // = R
    TMVAssert(kmax <= UV.rowsize());
    TMVAssert(UV.rowsize() <= UV.colsize());
    
    const size_t N = UV.rowsize();

    Matrix<T3,ColMajor> m2 = m; // MxR
    if (UV2) {
      TMVAssert(Qbeta);
      Q_LDivEq(UV,*Qbeta,m2.QuickView());
      Q_LDivEq(*UV2,Ubeta,m2.Rows(0,N));
    } else {
      Q_LDivEq(UV,Ubeta,m2.QuickView());
    }
    Matrix<T3,ColMajor> m3 = U1.QuickAdjoint().Rows(0,kmax) * m2.Rows(0,N); 
      // = KxR
    m3 /= DiagMatrixViewOf(S.SubVector(0,kmax));
    x = V1.QuickAdjoint().Cols(0,kmax) * m3; // NxR
    if (UV2) {
      Q_RDivEq(UV2->SubMatrix(0,N-1,1,N).QuickTranspose(),Vbeta,
	  x.Rows(1,N).QuickTranspose());
    } else {
      Q_RDivEq(UV.SubMatrix(0,N-1,1,N).QuickTranspose(),Vbeta,
	  x.Rows(1,N).QuickTranspose());
    }
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
    TMVAssert(m.rowsize() == V.rowsize()); // = N
    TMVAssert(x.rowsize() == U.colsize()); // = M
    TMVAssert(x.colsize() == m.colsize()); // = R
    TMVAssert(kmax <= U.colsize()); // = K
    TMVAssert(kmax <= V.rowsize());

    Matrix<T3,ColMajor> m2 = m * V.QuickAdjoint().Cols(0,kmax); // = RxK
    m2 %= DiagMatrixViewOf(S.SubVector(0,kmax));
    x = m2 * U.QuickAdjoint().Rows(0,kmax); // = RxM
  }

  template <class T1, class T2, class T3> void SV_RDiv(
      const GenMatrix<T1>& UV,
      const GenVector<T1>& Ubeta, const GenVector<T1>& Vbeta, 
      const GenMatrix<T1>* UV2, const GenVector<T1>* Qbeta,
      const GenMatrix<RealType(T1)>& U1, const GenMatrix<RealType(T1)>& V1,
      const GenVector<RealType(T1)>& S, size_t kmax,
      const GenMatrix<T2>& m, const MatrixView<T3>& x)
  {
    // x A = m
    // x U U1 S V1 V = m
    // x = m Vt V1t S^-1 U1t Ut
    TMVAssert(m.rowsize() == UV.rowsize()); // = N
    TMVAssert(x.rowsize() == UV.colsize()); // = M
    TMVAssert(x.colsize() == m.colsize()); // = R
    TMVAssert(kmax <= UV.rowsize()); // = K
    TMVAssert(UV.rowsize() <= UV.colsize());
    
    const size_t N = UV.rowsize();

    Matrix<T3,ColMajor> m2 = m; // = RxN
    if (UV2) {
      Q_LDivEq(UV2->SubMatrix(0,N-1,1,N).QuickTranspose(),Vbeta,
	  m2.Cols(1,N).QuickTranspose());
    } else {
      Q_LDivEq(UV.SubMatrix(0,N-1,1,N).QuickTranspose(),Vbeta,
	  m2.Cols(1,N).QuickTranspose());
    }
    Matrix<T3,ColMajor> m3 = m2 * V1.QuickAdjoint().Cols(0,kmax); // = RxK
    m3 %= DiagMatrixViewOf(S.SubVector(0,kmax));
    x.Cols(0,N) = m3 * U1.QuickAdjoint().Rows(0,kmax); // = RxN
    x.Cols(N,x.rowsize()).Zero(); // x = RxM
    if (UV2) {
      TMVAssert(Qbeta);
      Q_RDivEq(*UV2,Ubeta,x.Cols(0,N));
      Q_RDivEq(UV,*Qbeta,x);
    } else {
      Q_RDivEq(UV,Ubeta,x);
    }
  }

  // MJ: Write level 3 version of this? 
  template <class T> void NonLapBidiagonalize(
      const MatrixView<T>& UV, const VectorView<T>& Ubeta,
      const VectorView<T>& Vbeta, const VectorView<RealType(T)>& D,
      const VectorView<RealType(T)>& E, T& det)
  {
    // Decompose A (input as UV) into U B V
    // The Bidiagonal Matrix B is stored as two vectors: D, E
    // D is the diagonal, E is the super-diagonal
    // UV along with Ubeta and Vbeta hold the U and V matrices.
    const size_t M = UV.colsize();
    const size_t N = UV.rowsize();

    TMVAssert(N <= M);
    TMVAssert(N > 0);
    TMVAssert(Ubeta.size() == N);
    TMVAssert(Vbeta.size() == N-1);
    TMVAssert(D.size() == N);
    TMVAssert(E.size() == N-1);

    // We use Householder reflections to reduce A to the bidiagonal form:
    for(size_t j=0;j<N-1;++j) {
      Ubeta(j) = Householder_Reflect(UV.SubMatrix(j,M,j,N),det);
      Vbeta(j) = Householder_Reflect(UV.QuickTranspose().SubMatrix(j+1,N,j,M),
	  det);
    }
    Ubeta(N-1) = Householder_Reflect(UV.SubMatrix(N-1,M,N-1,N),det);

    // The bidiagonal of UV is the bidiagonal we want, so copy it to D,E
    if (IsComplex(T())) {
      TMVAssert(NormInf(UV.diag().Imag()) == RealType(T)(0));
      TMVAssert(NormInf(UV.diag(1).Imag()) == RealType(T)(0));
    }
    D = UV.diag().Real();
    E = UV.diag(1).Real();

    // The determinant of B is just the product of the diagonal elements:
    if (det!=T(0)) det *= DiagMatrixViewOf(D).Det();
  }
#ifdef LAP
  template <class T> void LapBidiagonalize(
      const MatrixView<T>& UV, const VectorView<T>& Ubeta,
      const VectorView<T>& Vbeta, const VectorView<RealType(T)>& D,
      const VectorView<RealType(T)>& E, T& det)
  { NonLapBidiagonalize(UV,Ubeta,Vbeta,D,E,det); }
  template <> void LapBidiagonalize(
      const MatrixView<double>& UV, const VectorView<double>& Ubeta,
      const VectorView<double>& Vbeta, const VectorView<double>& D,
      const VectorView<double>& E, double& det)
  {
    TMVAssert(UV.iscm());
    TMVAssert(D.step() == 1);
    TMVAssert(E.step() == 1);
    TMVAssert(UV.colsize() >= UV.rowsize());
    TMVAssert(Ubeta.size() == UV.rowsize());
    TMVAssert(Vbeta.size() == UV.rowsize()-1);
    TMVAssert(UV.ct()==NonConj);
    TMVAssert(D.ct()==NonConj);
    TMVAssert(E.ct()==NonConj);

    int m = UV.colsize();
    int n = UV.rowsize();
    int ldu = UV.stepj();
    Vector<double> Vbeta2(n);  
    // Stupid LAPack requires an extra element in the Vbeta vector
    // which it sets to 0 (!!!) rather than ignores.
    // So we need to create a temporary Vector which is size n.
    int lwork = (m+n)*LAP_BLOCKSIZE;
    double* work = LAP_DWork(lwork);
    int info;
    dgebrd(&m,&n,UV.ptr(),&ldu,D.ptr(),E.ptr(),Ubeta.ptr(),Vbeta2.ptr(),
	work,&lwork,&info);
    Vbeta = Vbeta2.SubVector(0,n-1);
    LAP_Results(info,int(work[0]),m,n,lwork,"dgebrd");
    if (det) {
      for(size_t i=0;i<Ubeta.size();++i) if (Ubeta(i) != double(0)) 
	det = -det;
      for(size_t i=0;i<Vbeta.size();++i) if (Vbeta(i) != double(0)) 
	det = -det;
      det *= DiagMatrixViewOf(D).Det();
    }
  }
  template <> void LapBidiagonalize(
      const MatrixView<complex<double> >& UV, 
      const VectorView<complex<double> >& Ubeta, 
      const VectorView<complex<double> >& Vbeta,
      const VectorView<double>& D, const VectorView<double>& E,
      complex<double>& det)
  {
    TMVAssert(UV.iscm());
    TMVAssert(D.step() == 1);
    TMVAssert(E.step() == 1);
    TMVAssert(UV.colsize() >= UV.rowsize());
    TMVAssert(Ubeta.size() == UV.rowsize());
    TMVAssert(Vbeta.size() == UV.rowsize()-1);
    TMVAssert(UV.ct()==NonConj);
    TMVAssert(D.ct()==NonConj);
    TMVAssert(E.ct()==NonConj);
    int m = UV.colsize();
    int n = UV.rowsize();
    int ldu = UV.stepj();
    Vector<complex<double> > Vbeta2(n);
    int lwork = (m+n)*LAP_BLOCKSIZE;
    complex<double>* work = LAP_ZWork(lwork);
    int info;

    zgebrd(&m,&n,LAP_Complex(UV.ptr()),&ldu,D.ptr(),E.ptr(),
	LAP_Complex(Ubeta.ptr()),LAP_Complex(Vbeta2.ptr()),
	LAP_Complex(work),&lwork,&info);
    Ubeta.ConjugateSelf();
    Vbeta = Vbeta2.SubVector(0,n-1);
    LAP_Results(info,int(real(work[0])),m,n,lwork,"zgebrd");
    if (det!=double(0)) {
      for(size_t i=0;i<Ubeta.size();++i) if (Ubeta(i) != double(0)) {
	det *= conj(Ubeta(i)*Ubeta(i))/norm(Ubeta(i));
      }
      for(size_t i=0;i<Vbeta.size();++i) if (Vbeta(i) != double(0)) {
	det *= -conj(Vbeta(i)*Vbeta(i))/norm(Vbeta(i));
      }
      det *= DiagMatrixViewOf(D).Det();
    }
  }
#ifndef NOFLOAT
  template <> void LapBidiagonalize(
      const MatrixView<float>& UV, const VectorView<float>& Ubeta,
      const VectorView<float>& Vbeta, const VectorView<float>& D,
      const VectorView<float>& E, float& det)
  {
    TMVAssert(UV.iscm());
    TMVAssert(D.step() == 1);
    TMVAssert(E.step() == 1);
    TMVAssert(UV.colsize() >= UV.rowsize());
    TMVAssert(Ubeta.size() == UV.rowsize());
    TMVAssert(Vbeta.size() == UV.rowsize()-1);
    TMVAssert(UV.ct()==NonConj);
    TMVAssert(D.ct()==NonConj);
    TMVAssert(E.ct()==NonConj);

    int m = UV.colsize();
    int n = UV.rowsize();
    int ldu = UV.stepj();
    Vector<float> Vbeta2(n);  
    int lwork = (m+n)*LAP_BLOCKSIZE;
    float* work = LAP_SWork(lwork);
    int info;
    sgebrd(&m,&n,UV.ptr(),&ldu,D.ptr(),E.ptr(),Ubeta.ptr(),Vbeta2.ptr(),
	work,&lwork,&info);
    Vbeta = Vbeta2.SubVector(0,n-1);
    LAP_Results(info,int(work[0]),m,n,lwork,"dgebrd");
    if (det) {
      for(size_t i=0;i<Ubeta.size();++i) if (Ubeta(i) != float(0)) 
	det = -det;
      for(size_t i=0;i<Vbeta.size();++i) if (Vbeta(i) != float(0)) 
	det = -det;
      det *= DiagMatrixViewOf(D).Det();
    }
  }
  template <> void LapBidiagonalize(
      const MatrixView<complex<float> >& UV, 
      const VectorView<complex<float> >& Ubeta, 
      const VectorView<complex<float> >& Vbeta,
      const VectorView<float>& D, const VectorView<float>& E,
      complex<float>& det)
  {
    TMVAssert(UV.iscm());
    TMVAssert(D.step() == 1);
    TMVAssert(E.step() == 1);
    TMVAssert(UV.colsize() >= UV.rowsize());
    TMVAssert(Ubeta.size() == UV.rowsize());
    TMVAssert(Vbeta.size() == UV.rowsize()-1);
    TMVAssert(UV.ct()==NonConj);
    TMVAssert(D.ct()==NonConj);
    TMVAssert(E.ct()==NonConj);
    int m = UV.colsize();
    int n = UV.rowsize();
    int ldu = UV.stepj();
    Vector<complex<float> > Vbeta2(n);
    int lwork = (m+n)*LAP_BLOCKSIZE;
    complex<float>* work = LAP_CWork(lwork);
    int info;

    cgebrd(&m,&n,LAP_Complex(UV.ptr()),&ldu,D.ptr(),E.ptr(),
	LAP_Complex(Ubeta.ptr()),LAP_Complex(Vbeta2.ptr()),
	LAP_Complex(work),&lwork,&info);
    Ubeta.ConjugateSelf();
    Vbeta = Vbeta2.SubVector(0,n-1);
    LAP_Results(info,int(real(work[0])),m,n,lwork,"zgebrd");
    if (det!=float(0)) {
      for(size_t i=0;i<Ubeta.size();++i) if (Ubeta(i) != float(0)) {
	det *= conj(Ubeta(i)*Ubeta(i))/norm(Ubeta(i));
      }
      for(size_t i=0;i<Vbeta.size();++i) if (Vbeta(i) != float(0)) {
	det *= -conj(Vbeta(i)*Vbeta(i))/norm(Vbeta(i));
      }
      det *= DiagMatrixViewOf(D).Det();
    }
  }
#endif // NOFLOAT
#endif // LAP

  template <class T> void Bidiagonalize(
      const MatrixView<T>& UV, const VectorView<T>& Ubeta,
      const VectorView<T>& Vbeta, const VectorView<RealType(T)>& D,
      const VectorView<RealType(T)>& E, T& det)
  {
    TMVAssert(UV.colsize() >= UV.rowsize());
    TMVAssert(UV.rowsize() == D.size());
    TMVAssert(Ubeta.size() == UV.rowsize());
    TMVAssert(Vbeta.size() == UV.rowsize()-1);
    TMVAssert(UV.isrm() || UV.iscm());
    TMVAssert(UV.ct()==NonConj);
    TMVAssert(D.ct()==NonConj);
    TMVAssert(E.ct()==NonConj);
    
    if (UV.rowsize() > 0) {
#ifdef LAP
      if (UV.iscm()) {
	LapBidiagonalize(UV,Ubeta,Vbeta,D,E,det);
      }
      else {
#ifdef TMVDEBUG
	cout<<"LAP Bidiag: UV is wrong:\n";
	cout<<"UV isrm = "<<UV.isrm()<<", step = "<<UV.stepi()<<endl;
	cout<<"D.step = "<<D.step()<<", E.step = "<<E.step()<<endl;
	cout<<"A_input = "<<UV<<endl;
#endif
	NonLapBidiagonalize(UV,Ubeta,Vbeta,D,E,det);
      }
#else
      NonLapBidiagonalize(UV,Ubeta,Vbeta,D,E,det);
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
      const MatrixView<TU>* U, const VectorView<T>& D, 
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
    if (U) TMVAssert(U->rowsize() == N);
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
	if (U) G.ConjMult(U->ColPair(i,0).QuickTranspose());
      }
    }
  }

  template <class T, class TV> void BidiagonalZeroLastCol(
      const VectorView<T>& D, const VectorView<T>& E,
      const MatrixView<TV>* V)
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
    if (V) { TMVAssert(V->colsize() == N); }
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
	if (V) G.ConjMult(V->RowPair(i,N-1));
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
      const MatrixView<TUV>* U, const VectorView<T>& D,
      const VectorView<T>& E, const MatrixView<TUV>* V)
  {
    // Reduce the superdiagonal elements of unreduced Bidiagonal Matrix B 
    // (given by D,E) while maintaining U B V. 
    // Note: the input B must be unreduced - ie. all entries are non-zero.
    const size_t N = D.size();
    TMVAssert(N>0);
    TMVAssert(E.size() == N-1);
    if (U) { TMVAssert(U->rowsize() == N); }
    if (V) { TMVAssert(V->colsize() == N); }
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
      if (V) G.ConjMult(V->RowPair(i-1,i));
      TMVAssert(x==T(0));
      G.Mult(x,*(++Di)); // x = B(i,i-1)
      G = Givens_Rotate(*(Di-1),x);
      G.Mult(*Ei,*Di);
      if (U) G.ConjMult(U->ColPair(i-1,i).QuickTranspose());
      if (i < N-1) {
	TMVAssert(x==T(0));
	G.Mult(x,*(++Ei)); // x = B(i-1,i+1)
	G = Givens_Rotate(*(Ei-1),x);
      } 
    }
  }

  template <class TUV, class T> void ReduceBidiagonal(
      const MatrixView<TUV>* U, const VectorView<T>& D,
      const VectorView<T>& E, const MatrixView<TUV>* V)
  {
    // Reduce the superdiagonal elements of Bidiagonal Matrix B 
    // (given by D,E) while maintaining U B V. 
    const size_t N = D.size();
    TMVAssert(E.size() == N-1);
    if (U) { TMVAssert(U->rowsize() == N); }
    if (V) { TMVAssert(V->colsize() == N); }

    // The input E(i) are all assumed to be non-zero.
    // If there are any zeros in D, we can zero the corresponding
    // E's (above and right) directly, so look for these first.
    // Loop invariant: all D(i) with p<=i<q are non-zero.
    size_t p=0; 
    for(size_t q=0; q<N; ++q) {
      if (D(q) == T(0)) {
	if (p<q) {
	  if (V) {
	    MatrixView<TUV> V1 = V->Rows(0,q+1);
	    BidiagonalZeroLastCol(D.SubVector(p,q+1), E.SubVector(p,q), &V1);
	    MatrixView<TUV> V2 = V->Rows(p,q);
	    if (U) {
	      MatrixView<TUV> U1 = U->Cols(p,q);
	      ReduceUnredBidiagonal(&U1, D.SubVector(p,q), E.SubVector(p,q-1), 
		  &V2);
	    } else {
	      ReduceUnredBidiagonal(ZMV<TUV>(), D.SubVector(p,q), 
		  E.SubVector(p,q-1), &V2);
	    }
	  } else {
	    BidiagonalZeroLastCol(D.SubVector(p,q+1), E.SubVector(p,q),
		ZMV<TUV>());
	    if (U) {
	      MatrixView<TUV> U1 = U->Cols(p,q);
	      ReduceUnredBidiagonal(&U1, D.SubVector(p,q), E.SubVector(p,q-1), 
		  ZMV<TUV>());
	    } else {
	      ReduceUnredBidiagonal(ZMV<TUV>(), D.SubVector(p,q),
		  E.SubVector(p,q-1), ZMV<TUV>());
	    }
	  }
	}
	if (q<N-1) {
	  if (U) {
	    MatrixView<TUV> U1 = U->Cols(q,N);
	    BidiagonalZeroFirstRow(&U1, D.SubVector(q,N), E.SubVector(q,N-1));
	  } else {
	    BidiagonalZeroFirstRow(ZMV<TUV>(), D.SubVector(q,N), 
		E.SubVector(q,N-1));
	  }
	}
	p=q+1;
      }
    }
    if (p<N) {
      if (U) {
	MatrixView<TUV> U1 = U->Cols(p,N);
	if (V) {
	  MatrixView<TUV> V1 = V->Rows(p,N);
	  ReduceUnredBidiagonal(&U1, D.SubVector(p,N), E.SubVector(p,N-1), &V1);
	} else {
	  ReduceUnredBidiagonal(&U1, D.SubVector(p,N), E.SubVector(p,N-1),
	      ZMV<TUV>());
	}
      } else {
	if (V) {
	  MatrixView<TUV> V1 = V->Rows(p,N);
	  ReduceUnredBidiagonal(ZMV<TUV>(), D.SubVector(p,N), 
	      E.SubVector(p,N-1), &V1);
	} else {
	  ReduceUnredBidiagonal(ZMV<TUV>(), D.SubVector(p,N), 
	      E.SubVector(p,N-1), ZMV<TUV>());
	}
      }
    }
  }

  template <class T> void NonLapSV_Decompose_From_Bidiagonal(
      const MatrixView<T>* U, const VectorView<RealType(T)>& D, 
      const VectorView<RealType(T)>& E, const MatrixView<T>* V,
      bool SetUV)
  {
    TMVAssert(D.size() == E.size()+1);
    if (U) {
      TMVAssert(U->colsize() >= U->rowsize());
      TMVAssert(U->rowsize() == D.size());
    }
    if (V) {
      TMVAssert(V->rowsize() == V->colsize());
      TMVAssert(V->rowsize() == D.size());
    }

    const size_t N = D.size();
    const size_t M = U ? U->colsize() : 0;

    if (SetUV) {
      TMVAssert(U && V);
      U->SetToIdentity();
      V->SetToIdentity();
    }

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
	if (U)
	  if (V) {
	    MatrixView<T> U1 = U->SubMatrix(0,M,p,q+1);
	    MatrixView<T> V1 = V->SubMatrix(p,q+1,0,N);
	    ReduceBidiagonal(&U1,D.SubVector(p,q+1),E.SubVector(p,q),&V1);
	  } else {
	    MatrixView<T> U1 = U->SubMatrix(0,M,p,q+1);
	    ReduceBidiagonal(&U1,D.SubVector(p,q+1),
		E.SubVector(p,q),(const MatrixView<T>*)(0));
	  }
	else
	  if (V) {
	    MatrixView<T> V1 = V->SubMatrix(p,q+1,0,N);
	    ReduceBidiagonal((const MatrixView<T>*)(0),D.SubVector(p,q+1),
		E.SubVector(p,q),&V1);
	  } else
	    ReduceBidiagonal((const MatrixView<T>*)(0),D.SubVector(p,q+1),
		E.SubVector(p,q),(const MatrixView<T>*)(0));
	BidiagonalChopSmallElements(D,E);
      }
    }

    // Make all of the singular values positive
    VIt<RealType(T),Unit,NonConj> Di = D.begin();
    for(size_t i=0;i<N;++i,++Di) if (*Di < 0) {
      *Di = -(*Di);
      if (V) V->row(i) = -V->row(i);
    }

    // Now A = U * S * V
    // Sort output singular values 
    Permutation sortp = SortPermutation(D,DESCEND);
    D = sortp * D;
    if (U) *U = *U * sortp.Transpose();
    if (V) *V = sortp * (*V);
  }

#ifdef LAP 
  template <class T> void LapSV_Decompose_From_Bidiagonal(
      const MatrixView<T>* U, const VectorView<RealType(T)>& D, 
      const VectorView<RealType(T)>& E, const MatrixView<T>* V,
      bool SetUV)
  { NonLapSV_Decompose_From_Bidiagonal(U,D,E,V,SetUV); }
  template <> void LapSV_Decompose_From_Bidiagonal(
      const MatrixView<double>* U, const VectorView<double>& D, 
      const VectorView<double>& E, const MatrixView<double>* V,
      bool SetUV)
  {
    TMVAssert(D.size() == E.size()+1);
    TMVAssert(D.ct()==NonConj);
    TMVAssert(E.ct()==NonConj);
    if (U) {
      TMVAssert(U->colsize() >= U->rowsize());
      TMVAssert(U->rowsize() == D.size());
      TMVAssert(U->iscm());
      TMVAssert(U->ct()==NonConj);
    }
    if (V) { 
      TMVAssert(V->rowsize() == V->colsize()); 
      TMVAssert(V->rowsize() == D.size()); 
      TMVAssert(V->iscm());
      TMVAssert(V->ct()==NonConj);
    }
    char u = 'U';
    int n = D.size();
    int lwork = (3*n+4)*n;
    double* work = LAP_DWork(lwork);
    lwork = 8*n;
    int* iwork = LAP_IWork(lwork);
    int info;
    if (SetUV) {
      char c = 'I';
      TMVAssert(U && V);
      int ldu = U->stepj();
      int ldv = V->stepj();
      dbdsdc(&u,&c,&n,D.ptr(),E.ptr(),U->ptr(),&ldu,V->ptr(),&ldv,0,0,
	  work,iwork,&info);
    } else if (U || V) {
      char c = 'I';
      Matrix<double,ColMajor> U1(n,n);
      Matrix<double,ColMajor> V1(n,n);
      int ldu = U1.stepj();
      int ldv = V1.stepj();
      dbdsdc(&u,&c,&n,D.ptr(),E.ptr(),U1.ptr(),&ldu,V1.ptr(),&ldv,0,0,
	  work,iwork,&info);
      if (U) *U = *U*U1;
      if (V) *V = V1*(*V);
    } else {
      int ldu = n;
      int ldv = n;
      char c = 'N';
      dbdsdc(&u,&c,&n,D.ptr(),E.ptr(),0,&ldu,0,&ldv,0,0,
	  work,iwork,&info);
    }
    LAP_Results(info,"dbdsdc");
  }
  template <> void LapSV_Decompose_From_Bidiagonal(
      const MatrixView<complex<double> >* U, const VectorView<double>& D, 
      const VectorView<double>& E, const MatrixView<complex<double> >* V,
      bool SetUV)
  {
    TMVAssert(D.size() == E.size()+1);
    TMVAssert(D.ct()==NonConj);
    TMVAssert(E.ct()==NonConj);
    if (U) {
      TMVAssert(U->colsize() >= U->rowsize());
      TMVAssert(U->rowsize() == D.size());
      TMVAssert(U->iscm());
      TMVAssert(U->ct()==NonConj);
    }
    if (V) { 
      TMVAssert(V->rowsize() == V->colsize()); 
      TMVAssert(V->rowsize() == D.size()); 
      TMVAssert(V->iscm());
      TMVAssert(V->ct()==NonConj);
    }
    TMVAssert(!SetUV);

    char u = 'U';
    int n = D.size();
    int lwork = (3*n+4)*n;
    double* work = LAP_DWork(lwork);
    lwork = 8*n;
    int* iwork = LAP_IWork(lwork);
    int info;
    if (U || V) {
      char c = 'I';
      Matrix<double,ColMajor> U1(n,n);
      Matrix<double,ColMajor> V1(n,n);
      int ldu = U1.stepj();
      int ldv = V1.stepj();
      dbdsdc(&u,&c,&n,D.ptr(),E.ptr(),U1.ptr(),&ldu,V1.ptr(),&ldv,0,0,
	  work,iwork,&info);
      if (U) *U = *U*U1;
      if (V) *V = V1*(*V);
    } else {
      int ldu = n;
      int ldv = n;
      char c = 'N';
      dbdsdc(&u,&c,&n,D.ptr(),E.ptr(),0,&ldu,0,&ldv,0,0,
	  work,iwork,&info);
    }
    LAP_Results(info,"dbdsdc");
  }
#ifndef NOFLOAT
  template <> void LapSV_Decompose_From_Bidiagonal(
      const MatrixView<float>* U, const VectorView<float>& D, 
      const VectorView<float>& E, const MatrixView<float>* V,
      bool SetUV)
  {
    TMVAssert(D.size() == E.size()+1);
    TMVAssert(D.ct()==NonConj);
    TMVAssert(E.ct()==NonConj);
    if (U) {
      TMVAssert(U->colsize() >= U->rowsize());
      TMVAssert(U->rowsize() == D.size());
      TMVAssert(U->iscm());
      TMVAssert(U->ct()==NonConj);
    }
    if (V) { 
      TMVAssert(V->rowsize() == V->colsize()); 
      TMVAssert(V->rowsize() == D.size()); 
      TMVAssert(V->iscm());
      TMVAssert(V->ct()==NonConj);
    }
    char u = 'U';
    int n = D.size();
    int lwork = (3*n+4)*n;
    float* work = LAP_SWork(lwork);
    lwork = 8*n;
    int* iwork = LAP_IWork(lwork);
    int info;
    if (SetUV) {
      char c = 'I';
      TMVAssert(U && V);
      int ldu = U->stepj();
      int ldv = V->stepj();
      sbdsdc(&u,&c,&n,D.ptr(),E.ptr(),U->ptr(),&ldu,V->ptr(),&ldv,0,0,
	  work,iwork,&info);
    } else if (U || V) {
      char c = 'I';
      Matrix<float,ColMajor> U1(n,n);
      Matrix<float,ColMajor> V1(n,n);
      int ldu = U1.stepj();
      int ldv = V1.stepj();
      sbdsdc(&u,&c,&n,D.ptr(),E.ptr(),U1.ptr(),&ldu,V1.ptr(),&ldv,0,0,
	  work,iwork,&info);
      if (U) *U = *U*U1;
      if (V) *V = V1*(*V);
    } else {
      char c = 'N';
      int ldu = n;
      int ldv = n;
      sbdsdc(&u,&c,&n,D.ptr(),E.ptr(),0,&ldu,0,&ldv,0,0,
	  work,iwork,&info);
    }
    LAP_Results(info,"sbdsdc");
  }
  template <> void LapSV_Decompose_From_Bidiagonal(
      const MatrixView<complex<float> >* U, const VectorView<float>& D, 
      const VectorView<float>& E, const MatrixView<complex<float> >* V,
      bool SetUV)
  {
    TMVAssert(D.size() == E.size()+1);
    TMVAssert(D.ct()==NonConj);
    TMVAssert(E.ct()==NonConj);
    if (U) {
      TMVAssert(U->colsize() >= U->rowsize());
      TMVAssert(U->rowsize() == D.size());
      TMVAssert(U->iscm());
      TMVAssert(U->ct()==NonConj);
    }
    if (V) { 
      TMVAssert(V->rowsize() == V->colsize()); 
      TMVAssert(V->rowsize() == D.size()); 
      TMVAssert(V->iscm());
      TMVAssert(V->ct()==NonConj);
    }
    TMVAssert(!SetUV);
    char u = 'U';
    int n = D.size();
    int lwork = (3*n+4)*n;
    float* work = LAP_SWork(lwork);
    lwork = 8*n;
    int* iwork = LAP_IWork(lwork);
    int info;
    if (U || V) {
      char c = 'I';
      Matrix<float,ColMajor> U1(n,n);
      Matrix<float,ColMajor> V1(n,n);
      int ldu = U1.stepj();
      int ldv = V1.stepj();
      sbdsdc(&u,&c,&n,D.ptr(),E.ptr(),U1.ptr(),&ldu,V1.ptr(),&ldv,0,0,
	  work,iwork,&info);
      if (U) *U = *U*U1;
      if (V) *V = V1*(*V);
    } else {
      char c = 'N';
      int ldu = n;
      int ldv = n;
      sbdsdc(&u,&c,&n,D.ptr(),E.ptr(),0,&ldu,0,&ldv,0,0,
	  work,iwork,&info);
    }
    LAP_Results(info,"sbdsdc");
  }
#endif
#endif // LAP

  template <class T> void SV_Decompose_From_Bidiagonal(
      const MatrixView<T>* U, const VectorView<RealType(T)>& D, 
      const VectorView<RealType(T)>& E, const MatrixView<T>* V,
      bool SetUV)
  {
    TMVAssert(D.size() == E.size()+1);
    TMVAssert(D.ct()==NonConj);
    TMVAssert(E.ct()==NonConj);
    if (U) {
      TMVAssert(U->colsize() >= U->rowsize());
      TMVAssert(U->rowsize() == D.size());
      TMVAssert(U->ct()==NonConj);
    }
    if (V) { 
      TMVAssert(V->rowsize() == V->colsize()); 
      TMVAssert(V->rowsize() == D.size()); 
      TMVAssert(V->ct()==NonConj);
    }
    
    if (D.size() > 0) {
#ifdef LAP
      // MJ: Write the row major versions (U = 'L' in Lap routines)
      if ((!U || U->iscm()) && D.step() == 1 && E.step() == 1 && 
	  (!V || V->iscm())) {
	LapSV_Decompose_From_Bidiagonal(U,D,E,V,SetUV);
      }
      else {
#ifdef TMVDEBUG
	cout<<"LAP SVDecomp_from_Bidiag: U,V are wrong:\n";
	if (U) cout<<"U iscm = "<<U->iscm()<<", step = "<<U->stepj()<<endl;
	if (V) cout<<"V iscm = "<<V->iscm()<<", step = "<<V->stepj()<<endl;
	cout<<"D.step = "<<D.step()<<", E.step = "<<E.step()<<endl;
#endif
	NonLapSV_Decompose_From_Bidiagonal(U,D,E,V,SetUV);
      }
#else
      NonLapSV_Decompose_From_Bidiagonal(U,D,E,V,SetUV);
#endif
    }
  }

  template <class T> void SV_Decompose(
      const MatrixView<T>& U, const VectorView<RealType(T)>& S, 
      const MatrixView<T>* V, T& det, bool StoreU)
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
    if (N == 0) return;

    TMVAssert(N <= M);
    if (V) {
      TMVAssert(V->colsize() == N);
      TMVAssert(V->rowsize() == N);
    }
    TMVAssert(S.size() == N);

    // If M is much larger than N (technically M > 5/3 N), then it is quicker
    // to start by doing a QR decomposition and then do SVD on the square
    // R matrix.  Thus, the final U of the SVD is Q (from the QR decomp)
    // times U from R's SVD.
    if (M > 5*N/3) {
      if (StoreU) {
	Matrix<T,ColMajor> R(N,N);
	QR_Decompose(U,R.QuickView(),det);
	SV_Decompose(R.QuickView(),S,V,det,StoreU);
	// Now R is a Unitary Matrix U'.  Need to multiply U by U'
	U = U*R;
      } else {
	Vector<T> Qbeta(N);
	QR_Decompose(U,Qbeta.View(),det);
	LowerTriMatrixViewOf(U.Rows(0,N)).OffDiag().Zero();
	SV_Decompose(U.Rows(0,N),S,V,det,StoreU);
      }
    } else {
      // First we reduce A to bidiagonal form: A = U * B * V
      // using a series of Householder transformations.
      // The diagonal of the Bidiagonal Matrix B is stored in D.
      // The superdiagonal is stored in E.
      Vector<RealType(T)> E(N-1);
      Vector<T> Ubeta(N);
      Vector<T> Vbeta(N-1);
      Bidiagonalize(U,Ubeta.View(),Vbeta.View(),S,E.View(),det);
      // Now UV stores Householder vectors for U in lower diagonal columns 
      // (HLi) and Householder vectors for V in upper diagonal rows (HRi)
      // The Householder matrices for U are actually the adjoints of the 
      // matrices that bidiagonalize A, and for V are the transposes:
      // U = HLn-1t ... HL1t HL0t A HR0T HR1T ... HRn-2T
      // Using the fact that H Ht = I, we get A = U B V with:
      // U = HL0 ... HLn-1 
      if (V) {
	V->row(0).MakeBasis(0);
	V->Rows(1,N) = U.Rows(0,N-1);
	V->col(0,1,N).Zero();
	GetQFromQR(V->SubMatrix(1,N,1,N).QuickTranspose(),Vbeta);
      }
      if (StoreU) GetQFromQR(U,Ubeta);

      if (StoreU) SV_Decompose_From_Bidiagonal(&U,S,E.View(),V);
      else SV_Decompose_From_Bidiagonal((const MatrixView<T>*)(0),S,E.View(),V);
    }
  }

  template <class T> void SV_Decompose(
      const MatrixView<T>& UV,
      const VectorView<T>& Ubeta, const VectorView<T>& Vbeta,
      Matrix<T,ColMajor>*& UV2, Vector<T>*& Qbeta,
      const MatrixView<RealType(T)>& U1, const MatrixView<RealType(T)>& V1,
      const VectorView<RealType(T)>& S, T& det)
  {
    // Decompose A (input as UV) into U U1 S V1 V
    // with U,V stored in UV, Ubeta, Vbeta
    // or Q U U1 S V1 V with Q stored in UV, Qbeta and U,V stored
    // in UV2, Ubeta, Vbeta
    
    const size_t M = UV.colsize();
    const size_t N = UV.rowsize();
    if (N == 0) return;

    TMVAssert(N <= M);
    TMVAssert(Ubeta.size() == N);
    TMVAssert(Vbeta.size() == N-1);
    TMVAssert(U1.colsize() == N);
    TMVAssert(U1.rowsize() == N);
    TMVAssert(V1.colsize() == N);
    TMVAssert(V1.rowsize() == N);
    TMVAssert(S.size() == N);

    if (M > 5*N/3) {
      UV2 = new Matrix<T,ColMajor>(N,N);
      Qbeta = new Vector<T>(N);

      QR_Decompose(UV,Qbeta->View(),det);
      *UV2 = UpperTriMatrixViewOf(UV);
      Matrix<T,ColMajor>* temp1=0;
      Vector<T>* temp2=0;
      SV_Decompose(UV2->QuickView(),Ubeta,Vbeta,temp1,temp2,U1,V1,S,det);
      TMVAssert(!temp1 && !temp2);
    } else {
      Vector<RealType(T)> E(N-1);
      Bidiagonalize(UV,Ubeta,Vbeta,S,E.View(),det);

      SV_Decompose_From_Bidiagonal(&U1,S,E.View(),&V1,true);
    }
  }

  template <class T> void SV_Decompose(
      const MatrixView<T>& U, const VectorView<RealType(T)>& S, 
      const MatrixView<T>& V, T& det, bool StoreU)
  { SV_Decompose(U,S,&V,det,StoreU); }

  template <class T> void SV_Decompose(
      const MatrixView<T>& U, const VectorView<RealType(T)>& S, 
      T& det, bool StoreU)
  { SV_Decompose(U,S,(const MatrixView<T>*)(0),det,StoreU); }

  template <class T> Matrix<T,ColMajor> SVFDiv<T>::Inverse() const
  { 
    TMVAssert(U && V);
    if (istrans) {
      // A^-1 = (Vt S^-1 Ut)T = U* S^-1 V*
      Matrix<T,ColMajor> SinvV = V->QuickConjugate().Rows(0,kmax) /
	DiagMatrixViewOf(S.SubVector(0,kmax));
      return U->QuickConjugate().Cols(0,kmax) * SinvV;
    } else {
      // A^-1 = Vt S^-1 Ut
      Matrix<T,ColMajor> SinvUt = U->QuickAdjoint().Rows(0,kmax) /
	DiagMatrixViewOf(S.SubVector(0,kmax));
      return V->QuickAdjoint().Cols(0,kmax) * SinvUt;
    }
  }

  template <class T> Matrix<T,ColMajor> SVFDiv<T>::DoInverseATA() const
  {
    TMVAssert(V);
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
    Matrix<T,ColMajor> SinvV = V->Rows(0,kmax) /
      DiagMatrixViewOf(S.SubVector(0,kmax));
    if (istrans)
      return SinvV.QuickTranspose() * SinvV.QuickConjugate();
    else
      return SinvV.QuickAdjoint() * SinvV;
  }

  template <class T> void SVFDiv<T>::Thresh(RealType(T) toler,
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

  template <class T> void SVFDiv<T>::Top(size_t neigen,
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

  template <class T> Matrix<T,ColMajor> SVDiv<T>::GetU() const
  {
    if (UV2) {
      TMVAssert(Qbeta);
      Matrix<T,ColMajor> U=UV;
      GetQFromQR(U.QuickView(),*Qbeta);
      Q_LDivEq(*UV2,Ubeta,U.Adjoint());
      return U*U1;
    } else {
      Matrix<T,ColMajor> U=UV;
      GetQFromQR(U.QuickView(),Ubeta);
      return U*U1;
    }
  }

  template <class T> Matrix<T,ColMajor> SVDiv<T>::GetV() const
  {
    const size_t N = UV.rowsize();
    Matrix<T,ColMajor> V(N,N);
    V.row(0).MakeBasis(0);
    if (UV2)
      V.Rows(1,N) = UV2->Rows(0,N-1);
    else
      V.Rows(1,N) = UV.Rows(0,N-1);
    V.col(0,1,N).Zero();
    GetQFromQR(V.SubMatrix(1,N,1,N).QuickTranspose(),Vbeta);
    return V1*V;
  }

  template <class T> Matrix<T,ColMajor> SVDiv<T>::Inverse() const
  { 
    if (istrans) {
      Matrix<T,ColMajor> inv(UV.colsize(),UV.rowsize());
      DoLDiv(Eye<T,ColMajor>(UV.rowsize()),inv.QuickView());
      return inv;
    } else {
      Matrix<T,ColMajor> inv(UV.rowsize(),UV.colsize());
      DoRDiv(Eye<T,ColMajor>(UV.rowsize()),inv.QuickView());
      return inv;
    }
  }

  template <class T> Matrix<T,ColMajor> SVDiv<T>::DoInverseATA() const
  {
    // A = U U1 S V1 V
    // At = Vt V1t S U1t Ut
    // AtA = Vt V1t S^2 V1 V
    // (AtA)^-1 = Vt V1t S^-2 V1 V
    //
    // if istrans:
    // AT = U S V
    // At = U* S V*
    // AAt = VT S^2 V*
    // (AAt)^-1 = VT S^-2 V*
    
    const size_t N = UV.rowsize();

    Matrix<T,ColMajor> VtSinv = V1.Rows(0,kmax).Adjoint() %
      DiagMatrixViewOf(S.SubVector(0,kmax));
    if (UV2) {
      Q_RDivEq(UV2->SubMatrix(0,N-1,1,N).QuickTranspose(),Vbeta,
	  VtSinv.Rows(1,N).QuickTranspose());
    } else {
      Q_RDivEq(UV.SubMatrix(0,N-1,1,N).QuickTranspose(),Vbeta,
	  VtSinv.Rows(1,N).QuickTranspose());
    }
    if (istrans)
      return VtSinv.QuickConjugate() * VtSinv.QuickTranspose();
    else
      return VtSinv * VtSinv.QuickAdjoint();
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


