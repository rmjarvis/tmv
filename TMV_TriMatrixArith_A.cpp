
#include "TMV_Tri.h"
#include "TMV_MatrixArith_A.h"

//#define XDEBUG

namespace tmv {

  //
  // Note: All routines here are designed to work even if x,y,A have 
  // the same storage.  So things like x += L*x will be efficient.
  // This also lets you do things like L *= L without a temporary.
  //

  //
  // MultMV
  //

  template <bool rm, bool cx, bool ca, class T, class Ta, class Tx>
    inline void DoRowAddMultMV(
	const GenUpperTriMatrix<Ta>& A, const GenVector<Tx>& x,
	const VectorView<T>& y)
    {
      //cerr<<"RowAddMultMV Upper\n";
      TMVAssert(A.size() == x.size());
      TMVAssert(A.size() == y.size());
      TMVAssert(y.size() > 0);
      TMVAssert(x.step()==1);
      TMVAssert(y.step()==1);
      TMVAssert(y.ct() == NonConj);
      TMVAssert(!A.isunit());
      TMVAssert(cx == x.isconj());
      TMVAssert(rm == A.isrm());
      TMVAssert(ca == A.isconj());

      const size_t N = x.size();

      const int sj = rm ? 1 : A.stepj();
      const int ds = A.stepi()+sj;
      const Tx* xi = x.cptr();
      T* yi = y.ptr();
      const Ta* Aii = A.cptr();
      size_t len = N;

      for(; len>0; --len,++xi,++yi,Aii+=ds) {
	// *yi += A.row(i,ii,N) * x.SubVector(ii,N);
	const Tx* xj = xi;
	const Ta* Aij = Aii;
	T temp(0);
	for(size_t j=len;j>0;--j,++xj,(rm?++Aij:Aij+=sj))
	  temp += (cx ? CONJ(*xj) : *xj) * (ca ? CONJ(*Aij) : *Aij);
	*yi += temp;
      }
    }

  template <bool rm, class T, class Ta, class Tx> inline void RowAddMultMV(
      const GenUpperTriMatrix<Ta>& A, const GenVector<Tx>& x,
      const VectorView<T>& y)
  { 
    if (x.isconj())
      if (A.isconj())
	DoRowAddMultMV<rm,true,true>(A,x,y);
      else
	DoRowAddMultMV<rm,true,false>(A,x,y);
    else
      if (A.isconj())
	DoRowAddMultMV<rm,false,true>(A,x,y);
      else
	DoRowAddMultMV<rm,false,false>(A,x,y);
  }

  template <bool cm, bool cx, bool ca, class T, class Ta, class Tx> 
    inline void DoColAddMultMV(
	const GenUpperTriMatrix<Ta>& A, const GenVector<Tx>& x,
	const VectorView<T>& y)
    {
      //cerr<<"ColAddMultMV Upper\n";
      TMVAssert(A.size() == x.size());
      TMVAssert(A.size() == y.size());
      TMVAssert(y.size() > 0);
      TMVAssert(x.step()==1);
      TMVAssert(y.step()==1);
      TMVAssert(y.ct() == NonConj);
      TMVAssert(!A.isunit());
      TMVAssert(cx == x.isconj());
      TMVAssert(cm == A.iscm());
      TMVAssert(ca == A.isconj());

      const size_t N = x.size();

      T* y0 = y.ptr();
      const Tx* xj = x.cptr();
      const int si = cm ? 1 : A.stepi();
      const Ta* A0j = A.cptr();

      for(size_t j=0; j<N; ++j,++xj,A0j+=A.stepj()) if (*xj != T(0)) {
	// y.SubVector(0,j+1) += *xj * A.col(j,0,j+1);
	const Ta* Aij = A0j;
	T* yi = y0;
	const Tx xxj = (cx ? CONJ(*xj) : *xj);
	for(size_t i=j+1;i>0;--i,++yi,(cm?++Aij:Aij+=si))
	  *yi += xxj * (ca ? CONJ(*Aij) : *Aij);
      }
    }

  template <bool cm, class T, class Ta, class Tx> inline void ColAddMultMV(
      const GenUpperTriMatrix<Ta>& A, const GenVector<Tx>& x,
      const VectorView<T>& y)
  { 
    if (x.isconj())
      if (A.isconj())
	DoColAddMultMV<cm,true,true>(A,x,y);
      else
	DoColAddMultMV<cm,true,false>(A,x,y);
    else
      if (A.isconj())
	DoColAddMultMV<cm,false,true>(A,x,y);
      else
	DoColAddMultMV<cm,false,false>(A,x,y);
  }

  template <bool rm, bool cx, bool ca, class T, class Ta, class Tx> 
    inline void DoRowAddMultMV(
	const GenLowerTriMatrix<Ta>& A, const GenVector<Tx>& x,
	const VectorView<T>& y)
    {
      //cerr<<"RowAddMultMV Lower\n";
      TMVAssert(A.size() == x.size());
      TMVAssert(A.size() == y.size());
      TMVAssert(y.size() > 0);
      TMVAssert(x.step()==1);
      TMVAssert(y.step()==1);
      TMVAssert(y.ct() == NonConj);
      TMVAssert(!A.isunit());
      TMVAssert(cx == x.isconj());
      TMVAssert(rm == A.isrm());
      TMVAssert(ca == A.isconj());

      const size_t N = x.size();
      const int sj = rm ? 1 : A.stepj();
      const int si = A.stepi();

      const Tx* x0 = x.cptr();
      T* yi = y.ptr() + N-1;
      const Ta* Ai0 = A.cptr()+(N-1)*A.stepi();

      for(size_t len=N; len>0; --len,--yi,Ai0-=si) {
	// i = N-1..0
	// *yi += A.row(i,0,ii) * x.SubVector(0,ii);
	const Ta* Aij = Ai0;
	const Tx* xj = x0;
	T temp(0);
	for(size_t j=len;j>0;--j,++xj,(rm?++Aij:Aij+=sj))
	  temp += (cx ? CONJ(*xj) : *xj) * (ca ? CONJ(*Aij) : *Aij);
	*yi += temp;
      }
    }

  template <bool rm, class T, class Ta, class Tx> inline void RowAddMultMV(
      const GenLowerTriMatrix<Ta>& A, const GenVector<Tx>& x,
      const VectorView<T>& y)
  { 
    if (x.isconj())
      if (A.isconj())
	DoRowAddMultMV<rm,true,true>(A,x,y);
      else
	DoRowAddMultMV<rm,true,false>(A,x,y);
    else
      if (A.isconj())
	DoRowAddMultMV<rm,false,true>(A,x,y);
      else
	DoRowAddMultMV<rm,false,false>(A,x,y);
  }

  template <bool cm, bool cx, bool ca, class T, class Ta, class Tx> 
    inline void DoColAddMultMV(
	const GenLowerTriMatrix<Ta>& A, const GenVector<Tx>& x,
	const VectorView<T>& y)
    {
      //cerr<<"ColAddMultMV Lower\n";
      TMVAssert(A.size() == x.size());
      TMVAssert(A.size() == y.size());
      TMVAssert(y.size() > 0);
      TMVAssert(x.step()==1);
      TMVAssert(y.step()==1);
      TMVAssert(y.ct() == NonConj);
      TMVAssert(!A.isunit());
      TMVAssert(cx == x.isconj());
      TMVAssert(cm == A.iscm());
      TMVAssert(ca == A.isconj());

      const size_t N = x.size();

      T* yj = y.ptr() + N-1;
      const Tx* xj = x.cptr() + N-1;
      const int si = cm ? 1 : A.stepi();
      const int ds = A.stepj()+si;
      const Ta* Ajj = A.cptr()+(N-1)*ds;

      for(size_t j=N,len=1;j>0;--j,++len,--yj,--xj,Ajj-=ds) if (*xj!=T(0)) {
	// y.SubVector(j,N) += *xj * A.col(j,j,N);
	T* yi = yj;
	const Ta* Aij = Ajj;
	const Tx xxj = (cx ? CONJ(*xj) : *xj);
	for (size_t i=len;i>0;--i,++yi,(cm?++Aij:Aij+=si))
	  *yi += xxj * (ca ? CONJ(*Aij) : *Aij);
      }
    }

  template <bool cm, class T, class Ta, class Tx> inline void ColAddMultMV(
      const GenLowerTriMatrix<Ta>& A, const GenVector<Tx>& x,
      const VectorView<T>& y)
  { 
    if (x.isconj())
      if (A.isconj())
	DoColAddMultMV<cm,true,true>(A,x,y);
      else
	DoColAddMultMV<cm,true,false>(A,x,y);
    else
      if (A.isconj())
	DoColAddMultMV<cm,false,true>(A,x,y);
      else
	DoColAddMultMV<cm,false,false>(A,x,y);
  }

  template <class T, class Ta, class Tx> inline void DoAddMultMV(
      const GenUpperTriMatrix<Ta>& A, const GenVector<Tx>& x,
      const VectorView<T>& y)
  {
    if (A.isunit()) {
      const size_t N = y.size();
      if (x.SameStorageAs(y) && N>1) {
	Vector<Tx> xx = x;
	AddMultMV(A.OffDiag(),xx.SubVector(1,N),y.SubVector(0,N-1));
	y += xx;
      } else {
	if (N > 1) 
	  AddMultMV(A.OffDiag(),x.SubVector(1,N),y.SubVector(0,N-1));
	y += x;
      }
    } else {
      if (A.isrm()) RowAddMultMV<true>(A,x,y);
      else if (A.iscm()) ColAddMultMV<true>(A,x,y);
      else RowAddMultMV<false>(A,x,y);
    }
  }

  template <class T, class Ta, class Tx> inline void AddMultMV(
      const GenUpperTriMatrix<Ta>& A, const GenVector<Tx>& x,
      const VectorView<T>& y)
    // y += A * x
  {
    TMVAssert(A.size() == x.size());
    TMVAssert(A.size() == y.size());
    TMVAssert(y.size() > 0);
    TMVAssert(x.step() == 1);
    TMVAssert(y.step() == 1);
    TMVAssert(y.ct() == NonConj);

    //      [ A11 A12 A13 ] [ 0  ]   [ A12 x2 ]
    // y += [  0  A22 A23 ] [ x2 ] = [ A22 x2 ]
    //      [  0   0  A33 ] [ 0  ]   [   0    ]

    const size_t N = x.size(); // = A.size()
    size_t j2 = N;
    for(const Tx* x2=x.cptr()+N-1; j2>0 && *x2==Tx(0); --j2,--x2);
    if (j2 == 0) return;
    size_t j1 = 0;
    for(const Tx* x1=x.cptr(); *x1==Tx(0); ++j1,++x1);
    if (j1 == 0 && j2 == N) {
      DoAddMultMV(A,x,y);
    } else {
      TMVAssert(j1 < j2);
      ConstUpperTriMatrixView<Ta> A22 = A.SubTriMatrix(j1,j2);
      ConstVectorView<Tx> x2 = x.SubVector(j1,j2);
      VectorView<T> y2 = y.SubVector(j1,j2);

      if (j1 != 0) {
	ConstMatrixView<Ta> A12 = A.SubMatrix(0,j1,j1,j2);
	VectorView<T> y1 = y.SubVector(0,j1);
	if (x.isconj())
	  UnitAMultMV1<false,true>(A12,x2,y1);
	else
	  UnitAMultMV1<false,false>(A12,x2,y1);
      }
      DoAddMultMV(A22,x2,y2);
    }
  }

  template <class T, class Ta, class Tx> inline void DoAddMultMV(
      const GenLowerTriMatrix<Ta>& A, const GenVector<Tx>& x,
      const VectorView<T>& y)
    // y += A * x
  {
    TMVAssert(A.size() == x.size());
    TMVAssert(A.size() == y.size());
    TMVAssert(y.size() > 0);
    TMVAssert(x.step() == 1);
    TMVAssert(y.step() == 1);

    if (A.isunit()) {
      const size_t N = y.size();
      if (x.SameStorageAs(y) && N>1) {
	Vector<Tx> xx = x;
	AddMultMV(A.OffDiag(),x.SubVector(0,N-1),y.SubVector(1,N));
	AddVV(T(1),xx,y);
      } else {
	if (N > 1) 
	  AddMultMV(A.OffDiag(),x.SubVector(0,N-1),y.SubVector(1,N));
	AddVV(T(1),x,y);
      }
    } else {
      if (A.isrm()) RowAddMultMV<true>(A,x,y);
      else if (A.iscm()) ColAddMultMV<true>(A,x,y);
      else RowAddMultMV<false>(A,x,y);
    }
  }

  template <class T, class Ta, class Tx> inline void AddMultMV(
      const GenLowerTriMatrix<Ta>& A, const GenVector<Tx>& x,
      const VectorView<T>& y)
    // y += A * x
  {
    TMVAssert(A.size() == x.size());
    TMVAssert(A.size() == y.size());
    TMVAssert(y.size() > 0);
    TMVAssert(x.step() == 1);
    TMVAssert(y.step() == 1);
    TMVAssert(y.ct() == NonConj);

    //      [ A11  0   0  ] [ 0  ]   [   0    ]
    // y += [ A21 A22  0  ] [ x2 ] = [ A22 x2 ]
    //      [ A31 A32 A33 ] [ 0  ]   [ A32 x2 ]

    const size_t N = x.size(); // = A.size()
    size_t j2 = N;
    for(const Tx* x2=x.cptr()+N-1; j2>0 && *x2==Tx(0); --j2,--x2);
    if (j2 == 0) return;
    size_t j1 = 0;
    for(const Tx* x1=x.cptr(); *x1==Tx(0); ++j1,++x1);
    if (j1 == 0 && j2 == N) {
      DoAddMultMV(A,x,y);
    } else {
      TMVAssert(j1 < j2);
      ConstLowerTriMatrixView<Ta> A22 = A.SubTriMatrix(j1,j2);
      ConstVectorView<Tx> x2 = x.SubVector(j1,j2);
      VectorView<T> y2 = y.SubVector(j1,j2);

      if (j2 != N) {
	ConstMatrixView<Ta> A32 = A.SubMatrix(j2,N,j1,j2);
	VectorView<T> y3 = y.SubVector(j2,N);
	if (x.isconj())
	  UnitAMultMV1<false,true>(A32,x2,y3);
	else
	  UnitAMultMV1<false,false>(A32,x2,y3);
      }
      DoAddMultMV(A22,x2,y2);
    }
  }

  // 
  // MultEqMV
  //

  template <bool rm, bool ca, bool ua, class T, class Ta> inline void DoRowMultEqMV(
      const GenUpperTriMatrix<Ta>& A, const VectorView<T>& x)
  {
    //cerr<<"RowMultEqMV Upper\n";
    TMVAssert(x.step()==1);
    TMVAssert(A.size() == x.size());
    TMVAssert(x.size() > 0);
    TMVAssert(x.ct() == NonConj);
    TMVAssert(rm == A.isrm());
    TMVAssert(ca == A.isconj());
    TMVAssert(ua == A.isunit());

    const size_t N = x.size();

    const int sj = rm ? 1 : A.stepj();
    const int ds = A.stepi()+sj;
    T* xi = x.ptr();
    const Ta* Aii = A.cptr();
    size_t len = N-1;

    for(; len>0; --len,++xi,Aii+=ds) {
      // i = 0..N-2
      // x(i) = A.row(i,i,N) * x.SubVector(i,N);
      if (!ua) *xi *= (ca ? CONJ(*Aii) : *Aii);
      const T* xj = xi+1;
      const Ta* Aij = Aii+sj;
      for(size_t j=len;j>0;--j,++xj,(rm?++Aij:Aij+=sj))
	*xi += (*xj) * (ca ? CONJ(*Aij) : *Aij);
    }
    if (!ua) *xi *= (ca ? CONJ(*Aii) : *Aii);
  }

  template <bool rm, class T, class Ta> inline void RowMultEqMV(
      const GenUpperTriMatrix<Ta>& A, const VectorView<T>& x)
  { 
    if (A.isconj())
      if (A.isunit())
	DoRowMultEqMV<rm,true,true>(A,x);
      else
	DoRowMultEqMV<rm,true,false>(A,x);
    else
      if (A.isunit())
	DoRowMultEqMV<rm,false,true>(A,x);
      else
	DoRowMultEqMV<rm,false,false>(A,x);
  }

  template <bool cm, bool ca, bool ua, class T, class Ta> inline void DoColMultEqMV(
      const GenUpperTriMatrix<Ta>& A, const VectorView<T>& x)
  {
    //cerr<<"ColMultEqMV Upper\n";
    TMVAssert(x.step()==1);
    TMVAssert(A.size() == x.size());
    TMVAssert(x.size() > 0);
    TMVAssert(x.ct() == NonConj);
    TMVAssert(cm == A.iscm());
    TMVAssert(ca == A.isconj());
    TMVAssert(ua == A.isunit());

    const size_t N = x.size();

    T* x0 = x.ptr();
    const T* xj = x0+1;
    const int si = cm ? 1 : A.stepi();
    const Ta* A0j = A.cptr();

    if (!ua) *x0 *= (ca ? CONJ(*A0j) : *A0j);
    A0j += A.stepj();

    for(size_t len=1; len<N; ++len,++xj,A0j+=A.stepj()) if (*xj != T(0)) {
      // j = 1..N-1
      // x.SubVector(0,j) += *xj * A.col(j,0,j);
      const Ta* Aij = A0j;
      T* xi = x0;
      for(size_t i=len;i>0;--i,++xi,(cm?++Aij:Aij+=si))
	*xi += *xj * (ca ? CONJ(*Aij) : *Aij);
      // Now Aij == Ajj, xi == xj
      // so this next statement is really *xj *= *Ajj
      if (!ua) *xi *= (ca ? CONJ(*Aij) : *Aij);
    }
  }

  template <bool cm, class T, class Ta> inline void ColMultEqMV(
      const GenUpperTriMatrix<Ta>& A, const VectorView<T>& x)
  { 
    if (A.isconj())
      if (A.isunit())
	DoColMultEqMV<cm,true,true>(A,x);
      else
	DoColMultEqMV<cm,true,false>(A,x);
    else
      if (A.isunit())
	DoColMultEqMV<cm,false,true>(A,x);
      else
	DoColMultEqMV<cm,false,false>(A,x);
  }

  template <bool rm, bool ca, bool ua, class T, class Ta> inline void DoRowMultEqMV(
      const GenLowerTriMatrix<Ta>& A, const VectorView<T>& x)
  {
    //cerr<<"RowMultEqMV Lower\n";
    TMVAssert(x.step()==1);
    TMVAssert(A.size() == x.size());
    TMVAssert(x.size() > 0);
    TMVAssert(x.ct() == NonConj);
    TMVAssert(rm == A.isrm());
    TMVAssert(ca == A.isconj());
    TMVAssert(ua == A.isunit());

    const size_t N = x.size();
    const int si = A.stepi();
    const int sj = rm ? 1 : A.stepj();
    const int ds = si+sj;

    const T* x0 = x.cptr();
    T* xi = x.ptr() + N-1;
    const Ta* Ai0 = A.cptr()+(N-1)*si;
    const Ta* Aii = Ai0 + (N-1)*sj;

    for(size_t len=N-1; len>0; --len,--xi,Ai0-=si,Aii-=ds) {
      // i = N-1..1
      // x(i) = A.row(i,0,i+1) * x.SubVector(0,i+1);
      if (!ua) *xi *= (ca ? CONJ(*Aii) : *Aii);
      const Ta* Aij = Ai0;
      const T* xj = x0;
      for(size_t j=len;j>0;--j,++xj,(rm?++Aij:Aij+=sj))
	*xi += *xj * (ca ? CONJ(*Aij) : *Aij);
    }
    if (!ua) *xi *= (ca ? CONJ(*Aii) : *Aii);
  }

  template <bool rm, class T, class Ta> inline void RowMultEqMV(
      const GenLowerTriMatrix<Ta>& A, const VectorView<T>& x)
  { 
    if (A.isconj())
      if (A.isunit())
	DoRowMultEqMV<rm,true,true>(A,x);
      else
	DoRowMultEqMV<rm,true,false>(A,x);
    else
      if (A.isunit())
	DoRowMultEqMV<rm,false,true>(A,x);
      else
	DoRowMultEqMV<rm,false,false>(A,x);
  }

  template <bool cm, bool ca, bool ua, class T, class Ta> inline void DoColMultEqMV(
      const GenLowerTriMatrix<Ta>& A, const VectorView<T>& x)
  {
    //cerr<<"ColMultEqMV Lower\n";
    TMVAssert(A.size() == x.size());
    TMVAssert(x.size() > 0);
    TMVAssert(x.ct() == NonConj);
    TMVAssert(x.step() == 1);
    TMVAssert(cm == A.iscm());
    TMVAssert(ca == A.isconj());
    TMVAssert(ua == A.isunit());

    const size_t N = x.size();

    T* xj = x.ptr() + N-2;
    const int si = cm ? 1 : A.stepi();
    const int ds = A.stepj()+si;
    const Ta* Ajj = A.cptr()+(N-2)*ds;

    if (!ua) *(xj+1) *= (ca ? CONJ(*(Ajj+ds)) : *(Ajj+ds));
    for(size_t j=N-1,len=1;j>0;--j,++len,--xj,Ajj-=ds) if (*xj!=T(0)) {
      // Actual j = N-2..0
      // x.SubVector(j+1,N) += *xj * A.col(j,j+1,N);
      T* xi = xj+1;
      const Ta* Aij = Ajj+si;
      for (size_t i=len;i>0;--i,++xi,(cm?++Aij:Aij+=si))
	*xi += *xj * (ca ? CONJ(*Aij) : *Aij);
      if (!ua) *xj *= (ca ? CONJ(*Ajj) : *Ajj);
    }
  }

  template <bool cm, class T, class Ta> inline void ColMultEqMV(
      const GenLowerTriMatrix<Ta>& A, const VectorView<T>& x)
  { 
    if (A.isconj())
      if (A.isunit())
	DoColMultEqMV<cm,true,true>(A,x);
      else
	DoColMultEqMV<cm,true,false>(A,x);
    else
      if (A.isunit())
	DoColMultEqMV<cm,false,true>(A,x);
      else
	DoColMultEqMV<cm,false,false>(A,x);
  }

  template <class T, class Ta> inline void DoMultEqMV(
      const GenUpperTriMatrix<Ta>& A, const VectorView<T>& x)
    // x = A * x
  {
    if (A.isrm()) RowMultEqMV<true>(A,x);
    else if (A.iscm()) ColMultEqMV<true>(A,x);
    else RowMultEqMV<false>(A,x);
  }

  template <class T, class Ta> inline void NonBlasMultEqMV(
      const GenUpperTriMatrix<Ta>& A, const VectorView<T>& x)
  {
    TMVAssert(A.size() == x.size());
    TMVAssert(x.size() > 0);
    TMVAssert(x.step() == 1);
    TMVAssert(x.ct() == NonConj);

    //     [ A11 A12 A13 ] [ 0  ]   [ A12 x2 ]
    // x = [  0  A22 A23 ] [ x2 ] = [ A22 x2 ]
    //     [  0   0  A33 ] [ 0  ]   [   0    ]

    const size_t N = x.size(); // = A.size()
    size_t j2 = N;
    for(const T* x2=x.cptr()+N-1; j2>0 && *x2==T(0); --j2,--x2);
    if (j2 == 0) return;
    size_t j1 = 0;
    for(const T* x1=x.cptr(); *x1==T(0); ++j1,++x1);
    if (j1 == 0 && j2 == N) {
      DoMultEqMV(A,x);
    } else {
      TMVAssert(j1 < j2);
      ConstUpperTriMatrixView<Ta> A22 = A.SubTriMatrix(j1,j2);
      VectorView<T> x2 = x.SubVector(j1,j2);

      if (j1 != 0) {
	ConstMatrixView<Ta> A12 = A.SubMatrix(0,j1,j1,j2);
	VectorView<T> x1 = x.SubVector(0,j1);
	UnitAMultMV1<false,false>(A12,x2,x1);
      }
      DoMultEqMV(A22,x2);
    }
  }

  template <class T, class Ta> inline void DoMultEqMV(
      const GenLowerTriMatrix<Ta>& A, const VectorView<T>& x)
  {
    if (A.isrm()) RowMultEqMV<true>(A,x);
    else if (A.iscm()) ColMultEqMV<true>(A,x);
    else RowMultEqMV<false>(A,x);
  }

  template <class T, class Ta> inline void NonBlasMultEqMV(
      const GenLowerTriMatrix<Ta>& A, const VectorView<T>& x)
    // x = A * x
  {
    TMVAssert(A.size() == x.size());
    TMVAssert(x.size() > 0);
    TMVAssert(x.step() == 1);
    TMVAssert(x.ct() == NonConj);

    //     [ A11  0   0  ] [ 0  ]   [   0    ]
    // x = [ A21 A22  0  ] [ x2 ] = [ A22 x2 ]
    //     [ A31 A32 A33 ] [ 0  ]   [ A32 x2 ]

    const size_t N = x.size(); // = A.size()
    size_t j2 = N;
    for(const T* x2=x.cptr()+N-1; j2>0 && *x2==T(0); --j2,--x2);
    if (j2 == 0) return;
    size_t j1 = 0;
    for(const T* x1=x.cptr(); *x1==T(0); ++j1,++x1);
    if (j1 == 0 && j2 == N) {
      DoMultEqMV(A,x);
    } else {
      TMVAssert(j1 < j2);
      ConstLowerTriMatrixView<Ta> A22 = A.SubTriMatrix(j1,j2);
      VectorView<T> x2 = x.SubVector(j1,j2);

      if (j2 != N) {
	ConstMatrixView<Ta> A32 = A.SubMatrix(j2,N,j1,j2);
	VectorView<T> x3 = x.SubVector(j2,N);
	UnitAMultMV1<false,false>(A32,x2,x3);
      }
      DoMultEqMV(A22,x2);
    }
  }

#ifdef BLAS
  template <class T, class Ta> inline void BlasMultEqMV(
      const GenUpperTriMatrix<Ta>& A, const VectorView<T>& x)
  { NonBlasMultEqMV(A,x); }
  template <class T, class Ta> inline void BlasMultEqMV(
      const GenLowerTriMatrix<Ta>& A, const VectorView<T>& x)
  { NonBlasMultEqMV(A,x); }
  template <> inline void BlasMultEqMV( 
      const GenUpperTriMatrix<double>& A, const VectorView<double>& x)
  {
    TMVAssert(A.size() == x.size());
    TMVAssert(x.size() > 0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(x.ct() == NonConj);

    int n=A.size();
    int lda=A.isrm()?A.stepi():A.stepj();
    int xs=x.step();
    BLASNAME(dtrmv) (BLASCM A.iscm()?BLASCH_UP:BLASCH_LO,
	A.iscm()?BLASCH_NT:BLASCH_T, A.isunit()?BLASCH_U:BLASCH_NU,
	BLASV(n),BLASP(A.cptr()),BLASV(lda),
	BLASP(x.ptr()),BLASV(xs) BLAS1 BLAS1 BLAS1);
  }
  template <> inline void BlasMultEqMV( 
      const GenLowerTriMatrix<double>& A, const VectorView<double>& x)
  {
    TMVAssert(A.size() == x.size());
    TMVAssert(x.size() > 0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(x.ct() == NonConj);

    int n=A.size();
    int lda=A.isrm()?A.stepi():A.stepj();
    int xs=x.step();
    BLASNAME(dtrmv) (BLASCM A.iscm()?BLASCH_LO:BLASCH_UP,
	A.iscm()?BLASCH_NT:BLASCH_T, A.isunit()?BLASCH_U:BLASCH_NU,
	BLASV(n),BLASP(A.cptr()),BLASV(lda),
	BLASP(x.ptr()),BLASV(xs) BLAS1 BLAS1 BLAS1);
  }
  template <> inline void BlasMultEqMV(
      const GenUpperTriMatrix<complex<double> >& A,
      const VectorView<complex<double> >& x)
  {
    TMVAssert(A.size() == x.size());
    TMVAssert(x.size() > 0);
    TMVAssert(x.ct() == NonConj);

    int n=A.size();
    int lda=A.isrm()?A.stepi():A.stepj();
    int xs=x.step();
    if (A.iscm() && A.isconj()) {
#ifdef CBLAS
      BLASNAME(ztrmv) (BLASRM BLASCH_LO, BLASCH_CT,
	  A.isunit()?BLASCH_U:BLASCH_NU,
	  BLASV(n),BLASP(A.cptr()),BLASV(lda),
	  BLASP(x.ptr()),BLASV(xs) BLAS1 BLAS1 BLAS1);
#else
      x.ConjugateSelf();
      BLASNAME(ztrmv) (BLASCM A.iscm()?BLASCH_UP:BLASCH_LO,
	  A.iscm()?BLASCH_NT:BLASCH_T, A.isunit()?BLASCH_U:BLASCH_NU,
	  BLASV(n),BLASP(A.cptr()),BLASV(lda),
	  BLASP(x.ptr()),BLASV(xs) BLAS1 BLAS1 BLAS1);
      x.ConjugateSelf();
#endif
    } else {
      BLASNAME(ztrmv) (BLASCM A.iscm()?BLASCH_UP:BLASCH_LO,
	  A.iscm()?BLASCH_NT:A.isconj()?BLASCH_CT:BLASCH_T, 
	  A.isunit()?BLASCH_U:BLASCH_NU,
	  BLASV(n),BLASP(A.cptr()),BLASV(lda),
	  BLASP(x.ptr()),BLASV(xs) BLAS1 BLAS1 BLAS1);
    }
  }
  template <> inline void BlasMultEqMV(
      const GenLowerTriMatrix<complex<double> >& A,
      const VectorView<complex<double> >& x)
  {
    TMVAssert(A.size() == x.size());
    TMVAssert(x.size() > 0);
    TMVAssert(x.ct() == NonConj);

    int n=A.size();
    int lda=A.isrm()?A.stepi():A.stepj();
    int xs=x.step();
    if (A.iscm() && A.isconj()) {
#ifdef CBLAS
      BLASNAME(ztrmv) (BLASRM BLASCH_UP, BLASCH_CT,
	  A.isunit()?BLASCH_U:BLASCH_NU,
	  BLASV(n),BLASP(A.cptr()),BLASV(lda),
	  BLASP(x.ptr()),BLASV(xs) BLAS1 BLAS1 BLAS1);
#else
      x.ConjugateSelf();
      BLASNAME(ztrmv) (BLASCM A.iscm()?BLASCH_LO:BLASCH_UP,
	  A.iscm()?BLASCH_NT:BLASCH_T, A.isunit()?BLASCH_U:BLASCH_NU,
	  BLASV(n),BLASP(A.cptr()),BLASV(lda),
	  BLASP(x.ptr()),BLASV(xs) BLAS1 BLAS1 BLAS1);
      x.ConjugateSelf();
#endif
    } else {
      BLASNAME(ztrmv) (BLASCM A.iscm()?BLASCH_LO:BLASCH_UP,
	  A.iscm()?BLASCH_NT:A.isconj()?BLASCH_CT:BLASCH_T, 
	  A.isunit()?BLASCH_U:BLASCH_NU,
	  BLASV(n),BLASP(A.cptr()),BLASV(lda),
	  BLASP(x.ptr()),BLASV(xs) BLAS1 BLAS1 BLAS1);
    }
  }
#ifndef NOFLOAT
  template <> inline void BlasMultEqMV( 
      const GenUpperTriMatrix<float>& A, const VectorView<float>& x)
  {
    TMVAssert(A.size() == x.size());
    TMVAssert(x.size() > 0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(x.ct() == NonConj);

    int n=A.size();
    int lda=A.isrm()?A.stepi():A.stepj();
    int xs=x.step();
    BLASNAME(strmv) (BLASCM A.iscm()?BLASCH_UP:BLASCH_LO,
	A.iscm()?BLASCH_NT:BLASCH_T, A.isunit()?BLASCH_U:BLASCH_NU,
	BLASV(n),BLASP(A.cptr()),BLASV(lda),
	BLASP(x.ptr()),BLASV(xs) BLAS1 BLAS1 BLAS1);
  }
  template <> inline void BlasMultEqMV( 
      const GenLowerTriMatrix<float>& A, const VectorView<float>& x)
  {
    TMVAssert(A.size() == x.size());
    TMVAssert(x.size() > 0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(x.ct() == NonConj);

    int n=A.size();
    int lda=A.isrm()?A.stepi():A.stepj();
    int xs=x.step();
    BLASNAME(strmv) (BLASCM A.iscm()?BLASCH_LO:BLASCH_UP,
	A.iscm()?BLASCH_NT:BLASCH_T, A.isunit()?BLASCH_U:BLASCH_NU,
	BLASV(n),BLASP(A.cptr()),BLASV(lda),
	BLASP(x.ptr()),BLASV(xs) BLAS1 BLAS1 BLAS1);
  }
  template <> inline void BlasMultEqMV(
      const GenUpperTriMatrix<complex<float> >& A,
      const VectorView<complex<float> >& x)
  {
    TMVAssert(A.size() == x.size());
    TMVAssert(x.size() > 0);
    TMVAssert(x.ct() == NonConj);

    int n=A.size();
    int lda=A.isrm()?A.stepi():A.stepj();
    int xs=x.step();
    if (A.iscm() && A.isconj()) {
#ifdef CBLAS
      BLASNAME(ctrmv) (BLASRM BLASCH_LO, BLASCH_CT,
	  A.isunit()?BLASCH_U:BLASCH_NU,
	  BLASV(n),BLASP(A.cptr()),BLASV(lda),
	  BLASP(x.ptr()),BLASV(xs) BLAS1 BLAS1 BLAS1);
#else
      x.ConjugateSelf();
      BLASNAME(ctrmv) (BLASCM A.iscm()?BLASCH_UP:BLASCH_LO,
	  A.iscm()?BLASCH_NT:BLASCH_T, A.isunit()?BLASCH_U:BLASCH_NU,
	  BLASV(n),BLASP(A.cptr()),BLASV(lda),
	  BLASP(x.ptr()),BLASV(xs) BLAS1 BLAS1 BLAS1);
      x.ConjugateSelf();
#endif
    } else {
      BLASNAME(ctrmv) (BLASCM A.iscm()?BLASCH_UP:BLASCH_LO,
	  A.iscm()?BLASCH_NT:A.isconj()?BLASCH_CT:BLASCH_T, 
	  A.isunit()?BLASCH_U:BLASCH_NU,
	  BLASV(n),BLASP(A.cptr()),BLASV(lda),
	  BLASP(x.ptr()),BLASV(xs) BLAS1 BLAS1 BLAS1);
    }
  }
  template <> inline void BlasMultEqMV(
      const GenLowerTriMatrix<complex<float> >& A,
      const VectorView<complex<float> >& x)
  {
    TMVAssert(A.size() == x.size());
    TMVAssert(x.size() > 0);
    TMVAssert(x.ct() == NonConj);

    int n=A.size();
    int lda=A.isrm()?A.stepi():A.stepj();
    int xs=x.step();
    if (A.iscm() && A.isconj()) {
#ifdef CBLAS
      BLASNAME(ctrmv) (BLASRM BLASCH_UP, BLASCH_CT,
	  A.isunit()?BLASCH_U:BLASCH_NU,
	  BLASV(n),BLASP(A.cptr()),BLASV(lda),
	  BLASP(x.ptr()),BLASV(xs) BLAS1 BLAS1 BLAS1);
#else
      x.ConjugateSelf();
      BLASNAME(ctrmv) (BLASCM A.iscm()?BLASCH_LO:BLASCH_UP,
	  A.iscm()?BLASCH_NT:BLASCH_T, A.isunit()?BLASCH_U:BLASCH_NU,
	  BLASV(n),BLASP(A.cptr()),BLASV(lda),
	  BLASP(x.ptr()),BLASV(xs) BLAS1 BLAS1 BLAS1);
      x.ConjugateSelf();
#endif
    } else {
      BLASNAME(ctrmv) (BLASCM A.iscm()?BLASCH_LO:BLASCH_UP,
	  A.iscm()?BLASCH_NT:A.isconj()?BLASCH_CT:BLASCH_T, 
	  A.isunit()?BLASCH_U:BLASCH_NU,
	  BLASV(n),BLASP(A.cptr()),BLASV(lda),
	  BLASP(x.ptr()),BLASV(xs) BLAS1 BLAS1 BLAS1);
    }
  }
#endif // FLOAT
#endif // BLAS

  template <class T, class Ta> inline void MultEqMV(
      const GenUpperTriMatrix<Ta>& A, const VectorView<T>& x)
  {
#ifdef XDEBUG
    Vector<T> x0 = x;
    Matrix<T> A0 = A;
    Vector<T> x2 = A0 * x0;
#endif
    TMVAssert(A.size() == x.size());
    TMVAssert(x.size() > 0);
    TMVAssert(x.step() == 1);

    if (x.isconj()) MultEqMV(A.Conjugate(),x.Conjugate());
    else {
#ifdef BLAS
      if (IsComplex(T()) && IsReal(Ta()))
	BlasMultEqMV(A,x);
      else if ((A.isrm() && A.stepi()>0) || (A.iscm() && A.stepj()>0)) 
	BlasMultEqMV(A,x);
      else {
	if (A.isunit()) {
	  UpperTriMatrix<Ta,UnitDiag,RowMajor> A2(A);
	  BlasMultEqMV(A2,x);
	} else {
	  UpperTriMatrix<Ta,NonUnitDiag,RowMajor> A2(A);
	  BlasMultEqMV(A2,x);
	}
      }
#else
      NonBlasMultEqMV(A,x);
#endif
    }
#ifdef XDEBUG
    if (Norm(x-x2) > 0.001*(Norm(A0)*Norm(x0))) {
      cerr<<"TriMultEqMV: \n";
      cerr<<"A = "<<A.cptr()<<"  "<<Type(A)<<"  "<<A0<<endl;
      cerr<<"x = "<<x.cptr()<<"  "<<Type(x)<<" step "<<x.step()<<"  "<<x0<<endl;
      cerr<<"-> x = "<<x<<endl;
      cerr<<"x2 = "<<x2<<endl;
      abort();
    }
#endif
  }

  template <class T, class Ta> inline void MultEqMV(
      const GenLowerTriMatrix<Ta>& A, const VectorView<T>& x)
  {
#ifdef XDEBUG
    Vector<T> x0 = x;
    Matrix<T> A0 = A;
    Vector<T> x2 = A0 * x0;
#endif
    TMVAssert(A.size() == x.size());
    TMVAssert(x.size() > 0);
    TMVAssert(x.step() == 1);

    if (x.isconj()) MultEqMV(A.Conjugate(),x.Conjugate());
    else {
#ifdef BLAS
      if (IsComplex(T()) && IsReal(Ta()))
	BlasMultEqMV(A,x);
      else if ( (A.isrm() && A.stepi()>0) || (A.iscm() && A.stepj()>0) )
	BlasMultEqMV(A,x);
      else {
	if (A.isunit()) {
	  LowerTriMatrix<Ta,UnitDiag,RowMajor> A2(A);
	  BlasMultEqMV(A2,x);
	} else {
	  LowerTriMatrix<Ta,NonUnitDiag,RowMajor> A2(A);
	  BlasMultEqMV(A2,x);
	}
      }
#else
      NonBlasMultEqMV(A,x);
#endif
    }

#ifdef XDEBUG
    if (Norm(x-x2) > 0.001*(Norm(A0)*Norm(x0))) {
      cerr<<"TriMultEqMV: \n";
      cerr<<"A = "<<A.cptr()<<"  "<<Type(A)<<"  "<<A0<<endl;
      cerr<<"x = "<<x.cptr()<<"  "<<Type(x)<<" step "<<x.step()<<"  "<<x0<<endl;
      cerr<<"-> x = "<<x<<endl;
      cerr<<"x2 = "<<x2<<endl;
      abort();
    }
#endif
  }

  template <class T, class Ta, class Tx> void MultMV(
      const T alpha, const GenUpperTriMatrix<Ta>& A,
      const GenVector<Tx>& x, const T beta, const VectorView<T>& y)
    // y = alpha * A * x + beta * y    
  { 
#ifdef XDEBUG
    Vector<T> x0 = x;
    Vector<T> y0 = y;
    Matrix<T> A0 = A;
    Vector<T> y2 = beta*y0+alpha*A0*x0;
#endif
    TMVAssert(A.size() == x.size());
    TMVAssert(A.size() == y.size());

    if (y.size() > 0) {
      if (alpha==T(0)) {
	y *= beta;
      } else if (beta == T(0) && y.step() == 1) {
	y = x;
	MultEqMV(A,y);
	y *= alpha;
      }
      else if (alpha != T(1) || x.step() != 1 || y.step() != 1 ||
	  ((x.SameStorageAs(y) || y.Real().ptr()==A.Real().cptr()) && 
	   beta != T(1))) {
	Vector<T> xx = alpha*x;
	MultEqMV(A,xx.View());
	AddVV(T(1),xx,beta,y);
      } else {
#ifdef BLAS
	Vector<T> xx = alpha*x;
	MultEqMV(A,xx.View());
	AddVV(T(1),xx,beta,y);
#else
	y *= beta;
	AddMultMV(A,x,y);
#endif
      }
    } 
#ifdef XDEBUG
    if (Norm(y-y2) > 0.001*(abs(alpha)*Norm(A0)*Norm(x0)+
	  (beta==T(0)?RealType(T)(0):abs(beta)*Norm(y0)))) {
      cerr<<"TriMultMV: alpha,beta = "<<alpha<<"  "<<beta<<endl;
      cerr<<"A = "<<Type(A)<<"  "<<A0<<endl;
      cerr<<"x = "<<Type(x)<<" step "<<x.step()<<"  "<<x0<<endl;
      cerr<<"y = "<<Type(y)<<" step "<<y.step()<<"  "<<y0<<endl;
      cerr<<"-> y = "<<y<<endl;
      cerr<<"y2 = "<<y2<<endl;
      abort();
    }
#endif
  }

  template <class T, class Ta, class Tx> void MultMV(
      const T alpha, const GenLowerTriMatrix<Ta>& A,
      const GenVector<Tx>& x, const T beta, const VectorView<T>& y)
    // y = alpha * A * x + beta * y    
  { 
#ifdef XDEBUG
    Vector<T> y0 = y;
    Vector<T> x0 = x;
    Matrix<T> A0 = A;
    Vector<T> y2 = beta*y0 + alpha*A0*x0;
#endif

    TMVAssert(A.size() == x.size());
    TMVAssert(A.size() == y.size());

    if (y.size() > 0) {
      if (alpha==T(0)) {
	y *= beta;
      } else if (beta == T(0) && y.step() == 1) {
	y = x;
	MultEqMV(A,y);
	if (alpha != T(1)) y *= alpha;
      } else if (alpha != T(1) || x.step() != 1 || y.step() != 1 ||
	  ((x.SameStorageAs(y) || y.ptr()==(const T*)(A.cptr())) && 
	   beta != T(1))) {
	Vector<T> xx = alpha*x;
	MultEqMV(A,xx.View());
	AddVV(T(1),xx,beta,y);
      } else {
#ifdef BLAS
	Vector<T> xx = alpha*x;
	MultEqMV(A,xx.View());
	AddVV(T(1),xx,beta,y);
#else
	y *= beta;
	AddMultMV(A,x,y);
#endif
      }
    }
#ifdef XDEBUG
    if (Norm(y-y2) > 0.001*(abs(alpha)*Norm(A0)*Norm(x0)+
	  (beta==T(0)?RealType(T)(0):abs(beta)*Norm(y0)))) {
      cerr<<"TriMultMV: alpha,beta = "<<alpha<<"  "<<beta<<endl;
      cerr<<"A = "<<A.cptr()<<"  "<<Type(A)<<"  "<<A0<<endl;
      cerr<<"x = "<<x.cptr()<<"  "<<Type(x)<<" step "<<x.step()<<"  "<<x0<<endl;
      cerr<<"y = "<<y.cptr()<<"  "<<Type(y)<<" step "<<y.step()<<"  "<<y0<<endl;
      cerr<<"-> y = "<<y<<endl;
      cerr<<"y2 = "<<y2<<endl;
      abort();
    }
#endif
  }

#define InstFile "TMV_TriMatrixArith_A.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


