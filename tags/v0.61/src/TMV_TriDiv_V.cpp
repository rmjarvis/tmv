///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// The Template Matrix/Vector Library for C++ was created by Mike Jarvis     //
// Copyright (C) 2008                                                        //
//                                                                           //
// The project is hosted at http://sourceforge.net/projects/tmv-cpp/         //
// where you can find the current version and current documention.           //
//                                                                           //
// For concerns or problems with the software, Mike may be contacted at      //
// mike_jarvis@users.sourceforge.net                                         //
//                                                                           //
// This program is free software; you can redistribute it and/or             //
// modify it under the terms of the GNU General Public License               //
// as published by the Free Software Foundation; either version 2            //
// of the License, or (at your option) any later version.                    //
//                                                                           //
// This program is distributed in the hope that it will be useful,           //
// but WITHOUT ANY WARRANTY; without even the implied warranty of            //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             //
// GNU General Public License for more details.                              //
//                                                                           //
// You should have received a copy of the GNU General Public License         //
// along with this program in the file LICENSE.                              //
//                                                                           //
// If not, write to:                                                         //
// The Free Software Foundation, Inc.                                        //
// 51 Franklin Street, Fifth Floor,                                          //
// Boston, MA  02110-1301, USA.                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////



#include "TMV_Blas.h"
#include "TMV_TriDiv.h"
#include "TMV_TriMatrix.h"
#include "TMV_Vector.h"

//#define XDEBUG

#ifdef XDEBUG
#include "TMV_MatrixArith.h"
#include "TMV_VectorArith.h"
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#endif

namespace tmv {

  //
  // TriLDivEq V
  //

  template <bool rm, bool ca, bool ua, class T, class Ta> 
    static void DoRowTri_LDivEq(
	const GenUpperTriMatrix<Ta>& A, const VectorView<T>& b)
    {
      // Solve A x = y  where A is an upper triangle matrix
      //cout<<"Row Upper\n";
      TMVAssert(b.step()==1);
      TMVAssert(A.size() == b.size());
      TMVAssert(b.size() > 0);
      TMVAssert(b.ct() == NonConj);
      TMVAssert(rm == A.isrm());
      TMVAssert(ca == A.isconj());
      TMVAssert(ua == A.isunit());

      const int N = A.size();

      const int sj = (rm?1:A.stepj());
      const int ds = A.stepi()+sj;
      const Ta* Aii = A.cptr() + (ua ? N-2 : N-1)*ds;
      T* bi = b.ptr() + (ua ? N-2 : N-1);

      if (!ua) {
	if (*Aii==Ta(0)) 
	  throw SingularUpperTriMatrix<Ta>(A);
#ifdef TMVFLDEBUG
	TMVAssert(bi >= b.first);
	TMVAssert(bi < b.last);
#endif
	*bi /= (ca ? CONJ(*Aii) : *Aii);
	Aii -= ds;
	--bi;
      }
      if (N==1) return;

      for(int i=N-1,len=1; i>0; --i,++len,Aii-=ds,--bi) {
	// Actual row being done is i-1, not i

	// *bi -= A.row(i,i+1,N) * b.SubVector(i+1,N);
	const T* bj = bi+1;
	const Ta* Aij = Aii + sj;
	for(int j=len;j>0;--j,++bj,(rm?++Aij:Aij+=sj)) {
#ifdef TMVFLDEBUG
	  TMVAssert(bi >= b.first);
	  TMVAssert(bi < b.last);
#endif
	  *bi -= (*bj) * (ca ? CONJ(*Aij) : *Aij);
	}

	if (!ua) {
	  if (*Aii==Ta(0)) 
	    throw SingularUpperTriMatrix<Ta>(A);
#ifdef TMVFLDEBUG
	  TMVAssert(bi >= b.first);
	  TMVAssert(bi < b.last);
#endif
	  *bi /= (ca ? CONJ(*Aii) : *Aii);
	}
      }
    }

  template <bool rm, class T, class Ta> static inline void RowTri_LDivEq(
      const GenUpperTriMatrix<Ta>& A, const VectorView<T>& b)
  {
    if (A.isconj())
      if (A.isunit())
	DoRowTri_LDivEq<rm,true,true>(A,b);
      else
	DoRowTri_LDivEq<rm,true,false>(A,b);
    else
      if (A.isunit())
	DoRowTri_LDivEq<rm,false,true>(A,b);
      else
	DoRowTri_LDivEq<rm,false,false>(A,b);
  }

  template <bool cm, bool ca, bool ua, class T, class Ta> 
    static void DoColTri_LDivEq(
	const GenUpperTriMatrix<Ta>& A, const VectorView<T>& b)
    {
      //cout<<"colmajor upper\n";
      // Solve A x = y  where A is an upper triangle matrix
      TMVAssert(b.step()==1);
      TMVAssert(A.size() == b.size());
      TMVAssert(b.size() > 0);
      TMVAssert(b.ct() == NonConj);
      TMVAssert(cm == A.iscm());
      TMVAssert(ca == A.isconj());
      TMVAssert(ua == A.isunit());

      const int N = A.size();

      const int si = (cm ? 1 : A.stepi());
      const int sj = A.stepj();
      const int ds = si+sj;
      const Ta* A0j = A.cptr()+(N-1)*sj;
      const Ta* Ajj = (ua ? 0 : A0j+(N-1)*si); // if unit, this isn't used.
      T*const b0 = b.ptr();
      T* bj = b0 + N-1;

      for(int j=N-1; j>0; --j,--bj,A0j-=sj) {
	if (*bj != T(0)) {
	  if (!ua) {
	    if (*Ajj==Ta(0)) 
	      throw SingularUpperTriMatrix<Ta>(A);
#ifdef TMVFLDEBUG
	    TMVAssert(bj >= b.first);
	    TMVAssert(bj < b.last);
#endif
	    *bj /= (ca ? CONJ(*Ajj) : *Ajj);
	    Ajj-=ds;
	  }

	  // b.SubVector(0,j) -= *bj * A.col(j,0,j);
	  T* bi = b0;
	  const Ta* Aij = A0j;
	  for(int i=j;i>0;--i,++bi,(cm?++Aij:Aij+=si)) {
#ifdef TMVFLDEBUG
	    TMVAssert(bi >= b.first);
	    TMVAssert(bi < b.last);
#endif
	    *bi -= *bj * (ca ? CONJ(*Aij) : *Aij);
	  }
	}
	else if (!ua) Ajj -= ds;
      }
      if (!ua && *bj != T(0)) {
	if (*Ajj==Ta(0)) 
	  throw SingularUpperTriMatrix<Ta>(A);
#ifdef TMVFLDEBUG
	TMVAssert(bj >= b.first);
	TMVAssert(bj < b.last);
#endif
	*bj /= (ca?CONJ(*Ajj):*Ajj);
      } 
    }

  template <bool cm, class T, class Ta> static inline void ColTri_LDivEq(
      const GenUpperTriMatrix<Ta>& A, const VectorView<T>& b)
  {
    if (A.isconj())
      if (A.isunit())
	DoColTri_LDivEq<cm,true,true>(A,b);
      else
	DoColTri_LDivEq<cm,true,false>(A,b);
    else
      if (A.isunit())
	DoColTri_LDivEq<cm,false,true>(A,b);
      else
	DoColTri_LDivEq<cm,false,false>(A,b);
  }

  template <bool rm, bool ca, bool ua, class T, class Ta> 
    static void DoRowTri_LDivEq(
	const GenLowerTriMatrix<Ta>& A, const VectorView<T>& b)
    {
      // Solve A x = y  where A is a lower triangle matrix
      TMVAssert(b.step()==1);
      TMVAssert(A.size() == b.size());
      TMVAssert(b.size() > 0);
      TMVAssert(b.ct() == NonConj);
      TMVAssert(rm == A.isrm());
      TMVAssert(ca == A.isconj());
      TMVAssert(ua == A.isunit());

      const int N = A.size();

      const int sj = (rm ? 1 : A.stepj());
      const int si = A.stepi();

      const Ta* Ai0 = A.cptr();
      T* b0 = b.ptr();

      if (!ua) {
	if (*Ai0==Ta(0)) 
	  throw SingularLowerTriMatrix<Ta>(A);
#ifdef TMVFLDEBUG
	TMVAssert(b0 >= b.first);
	TMVAssert(b0 < b.last);
#endif
	*b0 /= (ca ? CONJ(*Ai0) : *Ai0);
      }

      T* bi = b0+1;
      Ai0 += si;
      for(int i=1,len=1;i<N;++i,++len,++bi,Ai0+=si) {
	// *bi -= A.row(i,0,i) * b.SubVector(0,i);
	const Ta* Aij = Ai0;
	const T* bj = b0;
	for(int j=len;j>0;--j,++bj,(rm?++Aij:Aij+=sj)){
#ifdef TMVFLDEBUG
	  TMVAssert(bi >= b.first);
	  TMVAssert(bi < b.last);
#endif
	  *bi -= (*bj) * (ca ? CONJ(*Aij) : *Aij);
	}
	if (!ua) {
	  // Aij is Aii after the above for loop
	  if (*Aij==Ta(0)) 
	    throw SingularLowerTriMatrix<Ta>(A);
#ifdef TMVFLDEBUG
	  TMVAssert(bi >= b.first);
	  TMVAssert(bi < b.last);
#endif
	  *bi /= (ca ? CONJ(*Aij) : *Aij);
	}
      }
    }

  template <bool rm, class T, class Ta> static inline void RowTri_LDivEq(
      const GenLowerTriMatrix<Ta>& A, const VectorView<T>& b)
  {
    if (A.isconj())
      if (A.isunit())
	DoRowTri_LDivEq<rm,true,true>(A,b);
      else
	DoRowTri_LDivEq<rm,true,false>(A,b);
    else
      if (A.isunit())
	DoRowTri_LDivEq<rm,false,true>(A,b);
      else
	DoRowTri_LDivEq<rm,false,false>(A,b);
  }

  template <bool cm, bool ca, bool ua, class T, class Ta> 
    static void DoColTri_LDivEq(
	const GenLowerTriMatrix<Ta>& A, const VectorView<T>& b)
    {
      // Solve A x = y  where A is a lower triangle matrix
      TMVAssert(b.step()==1);
      TMVAssert(A.size() == b.size());
      TMVAssert(b.size() > 0);
      TMVAssert(b.ct() == NonConj);
      TMVAssert(cm == A.iscm());
      TMVAssert(ca == A.isconj());
      TMVAssert(ua == A.isunit());

      const int N = A.size();

      const int si = (cm ? 1 : A.stepi());
      const int ds = A.stepj()+si;
      const Ta* Ajj = A.cptr();
      T* bj = b.ptr();

      for(int j=0,len=N-1;len>0;++j,--len,++bj,Ajj+=ds) if (*bj != T(0)) {
	if (!ua) {
	  if (*Ajj==Ta(0)) 
	    throw SingularLowerTriMatrix<Ta>(A);
#ifdef TMVFLDEBUG
	  TMVAssert(bj >= b.first);
	  TMVAssert(bj < b.last);
#endif
	  *bj /= (ca ? CONJ(*Ajj) : *Ajj);
	}
	// b.SubVecotr(j+1,N) -= *bj * A.col(j,j+1,N);
	T* bi = bj+1;
	const Ta* Aij = Ajj+si;
	for(int i=len;i>0;--i,++bi,(cm?++Aij:Aij+=si)){
#ifdef TMVFLDEBUG
	  TMVAssert(bi >= b.first);
	  TMVAssert(bi < b.last);
#endif
	  *bi -= *bj * (ca ? CONJ(*Aij) : *Aij);
	}
      }
      if (!ua && *bj != T(0)) {
	if (*Ajj==Ta(0)) 
	  throw SingularLowerTriMatrix<Ta>(A);
#ifdef TMVFLDEBUG
	TMVAssert(bj >= b.first);
	TMVAssert(bj < b.last);
#endif
	*bj /= (ca ? CONJ(*Ajj) : *Ajj);
      } 
    }

  template <bool cm, class T, class Ta> static inline void ColTri_LDivEq(
      const GenLowerTriMatrix<Ta>& A, const VectorView<T>& b)
  {
    if (A.isconj())
      if (A.isunit())
	DoColTri_LDivEq<cm,true,true>(A,b);
      else
	DoColTri_LDivEq<cm,true,false>(A,b);
    else
      if (A.isunit())
	DoColTri_LDivEq<cm,false,true>(A,b);
      else
	DoColTri_LDivEq<cm,false,false>(A,b);
  }

  template <class T, class Ta> static inline void DoTri_LDivEq(
      const GenUpperTriMatrix<Ta>& A, const VectorView<T>& b)
  {
    if (A.isrm()) RowTri_LDivEq<true>(A,b);
    else if (A.iscm()) ColTri_LDivEq<true>(A,b);
    else RowTri_LDivEq<false>(A,b); 
  }

  template <class T, class Ta> static inline void DoTri_LDivEq(
      const GenLowerTriMatrix<Ta>& A, const VectorView<T>& b)
  {
    if (A.isrm()) RowTri_LDivEq<true>(A,b);
    else if (A.iscm()) ColTri_LDivEq<true>(A,b);
    else RowTri_LDivEq<false>(A,b); 
  }

  template <class T, class Ta> static void NonBlasTri_LDivEq(
      const GenUpperTriMatrix<Ta>& A, const VectorView<T>& b)
  {
    //cout<<"Upper LDivEq vect\n";
    TMVAssert(A.size() == b.size());
    TMVAssert(b.size()>0);
    TMVAssert(b.ct() == NonConj);

    if (b.step() == 1) {
      int i2 = b.size();
      for(const T* b2 = b.cptr()+i2-1; i2>0 && *b2==T(0); --i2,--b2);
      if (i2==0) return;
      else
	DoTri_LDivEq(A.SubTriMatrix(0,i2),b.SubVector(0,i2));
    } else {
      Vector<T> bb = b;
      NonBlasTri_LDivEq(A,bb.View());
      b = bb;
    }
  }

  template <class T, class Ta> static void NonBlasTri_LDivEq(
      const GenLowerTriMatrix<Ta>& A, const VectorView<T>& b)
  {
    //cout<<"Lower LDivEq vect\n";
    TMVAssert(A.size() == b.size());
    TMVAssert(b.size()>0);
    TMVAssert(b.ct() == NonConj);

    if (b.step() == 1) {
      const int N = b.size();
      int i1 = 0;
      for(const T* b1 = b.cptr(); i1<N && *b1==T(0); ++i1,++b1);
      if (i1==N) return;
      else
	DoTri_LDivEq(A.SubTriMatrix(i1,N),b.SubVector(i1,N));
    } else {
      Vector<T> bb = b;
      NonBlasTri_LDivEq(A,bb.View());
      b = bb;
    }
  }

#ifdef BLAS
  template <class T, class Ta> static inline void BlasTri_LDivEq(
      const GenUpperTriMatrix<Ta>& A, const VectorView<T>& b)
  { NonBlasTri_LDivEq(A,b); }
  template <class T, class Ta> static inline void BlasTri_LDivEq(
      const GenLowerTriMatrix<Ta>& A, const VectorView<T>& b)
  { NonBlasTri_LDivEq(A,b); }
#ifdef INST_DOUBLE
  template <> void BlasTri_LDivEq(
      const GenUpperTriMatrix<double>& A, const VectorView<double>& b)
  {
    TMVAssert(A.size() == b.size());
    TMVAssert(b.size()>0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(b.ct() == NonConj);
    int n=A.size();
    int lda = A.isrm()?A.stepi():A.stepj();
    int bs = b.step();
    BLASNAME(dtrsv) (BLASCM A.iscm()?BLASCH_UP:BLASCH_LO,
	A.iscm()?BLASCH_NT:BLASCH_T, A.isunit()?BLASCH_U:BLASCH_NU,
	BLASV(n),BLASP(A.cptr()),BLASV(lda),BLASP(b.ptr()),BLASV(bs)
	BLAS1 BLAS1 BLAS1);
  }
  template <> void BlasTri_LDivEq(
      const GenLowerTriMatrix<double>& A, const VectorView<double>& b)
  {
    TMVAssert(A.size() == b.size());
    TMVAssert(b.size()>0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(b.ct() == NonConj);
    int n=A.size();
    int lda = A.isrm()?A.stepi():A.stepj();
    int bs = b.step();
    BLASNAME(dtrsv) (BLASCM A.iscm()?BLASCH_LO:BLASCH_UP,
	A.iscm()?BLASCH_NT:BLASCH_T, A.isunit()?BLASCH_U:BLASCH_NU,
	BLASV(n),BLASP(A.cptr()),BLASV(lda),BLASP(b.ptr()),BLASV(bs)
	BLAS1 BLAS1 BLAS1);
  }
  template <> void BlasTri_LDivEq(
      const GenUpperTriMatrix<std::complex<double> >& A,
      const VectorView<std::complex<double> >& b)
  {
    TMVAssert(A.size() == b.size());
    TMVAssert(b.size()>0);
    TMVAssert(b.ct() == NonConj);

    int n=A.size();
    int lda = A.isrm()?A.stepi():A.stepj();
    int bs = b.step();
    if (A.iscm() && A.isconj()) {
#ifdef CBLAS
      BLASNAME(ztrsv) (BLASRM BLASCH_LO,BLASCH_CT,
	  A.isunit()?BLASCH_U:BLASCH_NU,
	  BLASV(n),BLASP(A.cptr()),BLASV(lda),BLASP(b.ptr()),BLASV(bs)
	  BLAS1 BLAS1 BLAS1);
#else
      b.ConjugateSelf();
      BLASNAME(ztrsv) (BLASCM BLASCH_UP,BLASCH_NT,
	  A.isunit()?BLASCH_U:BLASCH_NU,
	  BLASV(n),BLASP(A.cptr()),BLASV(lda),BLASP(b.ptr()),BLASV(bs)
	  BLAS1 BLAS1 BLAS1);
      b.ConjugateSelf();
#endif
    } else {
      BLASNAME(ztrsv) (BLASCM A.iscm()?BLASCH_UP:BLASCH_LO,
	  A.iscm()?BLASCH_NT:A.isconj()?BLASCH_CT:BLASCH_T,
	  A.isunit()?BLASCH_U:BLASCH_NU,
	  BLASV(n),BLASP(A.cptr()),BLASV(lda),BLASP(b.ptr()),BLASV(bs)
	  BLAS1 BLAS1 BLAS1);
    }
  }
  template <> void BlasTri_LDivEq(
      const GenLowerTriMatrix<std::complex<double> >& A,
      const VectorView<std::complex<double> >& b)
  {
    TMVAssert(A.size() == b.size());
    TMVAssert(b.size()>0);
    TMVAssert(b.ct() == NonConj);

    int n=A.size();
    int lda = A.isrm()?A.stepi():A.stepj();
    int bs = b.step();
    if (A.iscm() && A.isconj()) {
#ifdef CBLAS
      BLASNAME(ztrsv) (BLASRM BLASCH_UP,BLASCH_CT,
	  A.isunit()?BLASCH_U:BLASCH_NU,
	  BLASV(n),BLASP(A.cptr()),BLASV(lda),BLASP(b.ptr()),BLASV(bs)
	  BLAS1 BLAS1 BLAS1);
#else
      b.ConjugateSelf();
      BLASNAME(ztrsv) (BLASCM BLASCH_LO,BLASCH_NT,
	  A.isunit()?BLASCH_U:BLASCH_NU,
	  BLASV(n),BLASP(A.cptr()),BLASV(lda),BLASP(b.ptr()),BLASV(bs)
	  BLAS1 BLAS1 BLAS1);
      b.ConjugateSelf();
#endif
    } else {
      BLASNAME(ztrsv) (BLASCM A.iscm()?BLASCH_LO:BLASCH_UP,
	  A.iscm()?BLASCH_NT:A.isconj()?BLASCH_CT:BLASCH_T,
	  A.isunit()?BLASCH_U:BLASCH_NU,
	  BLASV(n),BLASP(A.cptr()),BLASV(lda),BLASP(b.ptr()),BLASV(bs)
	  BLAS1 BLAS1 BLAS1);
    }
  }
  template <> void BlasTri_LDivEq(
      const GenUpperTriMatrix<double>& A,
      const VectorView<std::complex<double> >& b)
  {
    TMVAssert(A.size() == b.size());
    TMVAssert(b.size()>0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(b.ct() == NonConj);
    int n=A.size();
    int lda = A.isrm()?A.stepi():A.stepj();
    int bs = 2*b.step();
    BLASNAME(dtrsv) (BLASCM A.iscm()?BLASCH_UP:BLASCH_LO,
	A.iscm()?BLASCH_NT:BLASCH_T, A.isunit()?BLASCH_U:BLASCH_NU,
	BLASV(n),BLASP(A.cptr()),BLASV(lda),BLASP((double*)b.ptr()),BLASV(bs)
	BLAS1 BLAS1 BLAS1);
    BLASNAME(dtrsv) (BLASCM A.iscm()?BLASCH_UP:BLASCH_LO,
	A.iscm()?BLASCH_NT:BLASCH_T, A.isunit()?BLASCH_U:BLASCH_NU,
	BLASV(n),BLASP(A.cptr()),BLASV(lda),BLASP((double*)b.ptr()+1),BLASV(bs)
	BLAS1 BLAS1 BLAS1);
  }
  template <> void BlasTri_LDivEq(
      const GenLowerTriMatrix<double>& A, 
      const VectorView<std::complex<double> >& b)
  {
    TMVAssert(A.size() == b.size());
    TMVAssert(b.size()>0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(b.ct() == NonConj);
    int n=A.size();
    int lda = A.isrm()?A.stepi():A.stepj();
    int bs = 2*b.step();
    BLASNAME(dtrsv) (BLASCM A.iscm()?BLASCH_LO:BLASCH_UP,
	A.iscm()?BLASCH_NT:BLASCH_T, A.isunit()?BLASCH_U:BLASCH_NU,
	BLASV(n),BLASP(A.cptr()),BLASV(lda),BLASP((double*)b.ptr()),BLASV(bs)
	BLAS1 BLAS1 BLAS1);
    BLASNAME(dtrsv) (BLASCM A.iscm()?BLASCH_LO:BLASCH_UP,
	A.iscm()?BLASCH_NT:BLASCH_T, A.isunit()?BLASCH_U:BLASCH_NU,
	BLASV(n),BLASP(A.cptr()),BLASV(lda),BLASP((double*)b.ptr()+1),BLASV(bs)
	BLAS1 BLAS1 BLAS1);
  }
#endif
#ifdef INST_FLOAT
  template <> void BlasTri_LDivEq(
      const GenUpperTriMatrix<float>& A, const VectorView<float>& b)
  {
    TMVAssert(A.size() == b.size());
    TMVAssert(b.size()>0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(b.ct() == NonConj);
    int n=A.size();
    int lda = A.isrm()?A.stepi():A.stepj();
    int bs = b.step();
    BLASNAME(strsv) (BLASCM A.iscm()?BLASCH_UP:BLASCH_LO,
	A.iscm()?BLASCH_NT:BLASCH_T, A.isunit()?BLASCH_U:BLASCH_NU,
	BLASV(n),BLASP(A.cptr()),BLASV(lda),BLASP(b.ptr()),BLASV(bs)
	BLAS1 BLAS1 BLAS1);
  }
  template <> void BlasTri_LDivEq(
      const GenLowerTriMatrix<float>& A, const VectorView<float>& b)
  {
    TMVAssert(A.size() == b.size());
    TMVAssert(b.size()>0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(b.ct() == NonConj);
    int n=A.size();
    int lda = A.isrm()?A.stepi():A.stepj();
    int bs = b.step();
    BLASNAME(strsv) (BLASCM A.iscm()?BLASCH_LO:BLASCH_UP,
	A.iscm()?BLASCH_NT:BLASCH_T, A.isunit()?BLASCH_U:BLASCH_NU,
	BLASV(n),BLASP(A.cptr()),BLASV(lda),BLASP(b.ptr()),BLASV(bs)
	BLAS1 BLAS1 BLAS1);
  }
  template <> void BlasTri_LDivEq(
      const GenUpperTriMatrix<std::complex<float> >& A,
      const VectorView<std::complex<float> >& b)
  {
    TMVAssert(A.size() == b.size());
    TMVAssert(b.size()>0);
    TMVAssert(b.ct() == NonConj);

    int n=A.size();
    int lda = A.isrm()?A.stepi():A.stepj();
    int bs = b.step();
    if (A.iscm() && A.isconj()) {
#ifdef CBLAS
      BLASNAME(ctrsv) (BLASRM BLASCH_LO,BLASCH_CT,
	  A.isunit()?BLASCH_U:BLASCH_NU,
	  BLASV(n),BLASP(A.cptr()),BLASV(lda),BLASP(b.ptr()),BLASV(bs)
	  BLAS1 BLAS1 BLAS1);
#else
      b.ConjugateSelf();
      BLASNAME(ctrsv) (BLASCM BLASCH_UP,BLASCH_NT,
	  A.isunit()?BLASCH_U:BLASCH_NU,
	  BLASV(n),BLASP(A.cptr()),BLASV(lda),BLASP(b.ptr()),BLASV(bs)
	  BLAS1 BLAS1 BLAS1);
      b.ConjugateSelf();
#endif
    } else {
      BLASNAME(ctrsv) (BLASCM A.iscm()?BLASCH_UP:BLASCH_LO,
	  A.iscm()?BLASCH_NT:A.isconj()?BLASCH_CT:BLASCH_T,
	  A.isunit()?BLASCH_U:BLASCH_NU,
	  BLASV(n),BLASP(A.cptr()),BLASV(lda),BLASP(b.ptr()),BLASV(bs)
	  BLAS1 BLAS1 BLAS1);
    }
  }
  template <> void BlasTri_LDivEq(
      const GenLowerTriMatrix<std::complex<float> >& A,
      const VectorView<std::complex<float> >& b)
  {
    TMVAssert(A.size() == b.size());
    TMVAssert(b.size()>0);
    TMVAssert(b.ct() == NonConj);

    int n=A.size();
    int lda = A.isrm()?A.stepi():A.stepj();
    int bs = b.step();
    if (A.iscm() && A.isconj()) {
#ifdef CBLAS
      BLASNAME(ctrsv) (BLASRM BLASCH_UP,BLASCH_CT,
	  A.isunit()?BLASCH_U:BLASCH_NU,
	  BLASV(n),BLASP(A.cptr()),BLASV(lda),BLASP(b.ptr()),BLASV(bs)
	  BLAS1 BLAS1 BLAS1);
#else
      b.ConjugateSelf();
      BLASNAME(ctrsv) (BLASCM BLASCH_LO,BLASCH_NT,
	  A.isunit()?BLASCH_U:BLASCH_NU,
	  BLASV(n),BLASP(A.cptr()),BLASV(lda),BLASP(b.ptr()),BLASV(bs)
	  BLAS1 BLAS1 BLAS1);
      b.ConjugateSelf();
#endif
    } else {
      BLASNAME(ctrsv) (BLASCM A.iscm()?BLASCH_LO:BLASCH_UP,
	  A.iscm()?BLASCH_NT:A.isconj()?BLASCH_CT:BLASCH_T,
	  A.isunit()?BLASCH_U:BLASCH_NU,
	  BLASV(n),BLASP(A.cptr()),BLASV(lda),BLASP(b.ptr()),BLASV(bs)
	  BLAS1 BLAS1 BLAS1);
    }
  }
  template <> void BlasTri_LDivEq(
      const GenUpperTriMatrix<float>& A,
      const VectorView<std::complex<float> >& b)
  {
    TMVAssert(A.size() == b.size());
    TMVAssert(b.size()>0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(b.ct() == NonConj);
    int n=A.size();
    int lda = A.isrm()?A.stepi():A.stepj();
    int bs = 2*b.step();
    BLASNAME(strsv) (BLASCM A.iscm()?BLASCH_UP:BLASCH_LO,
	A.iscm()?BLASCH_NT:BLASCH_T, A.isunit()?BLASCH_U:BLASCH_NU,
	BLASV(n),BLASP(A.cptr()),BLASV(lda),BLASP((float*)b.ptr()),BLASV(bs)
	BLAS1 BLAS1 BLAS1);
    BLASNAME(strsv) (BLASCM A.iscm()?BLASCH_UP:BLASCH_LO,
	A.iscm()?BLASCH_NT:BLASCH_T, A.isunit()?BLASCH_U:BLASCH_NU,
	BLASV(n),BLASP(A.cptr()),BLASV(lda),BLASP((float*)b.ptr()+1),BLASV(bs)
	BLAS1 BLAS1 BLAS1);
  }
  template <> void BlasTri_LDivEq(
      const GenLowerTriMatrix<float>& A, 
      const VectorView<std::complex<float> >& b)
  {
    TMVAssert(A.size() == b.size());
    TMVAssert(b.size()>0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(b.ct() == NonConj);
    int n=A.size();
    int lda = A.isrm()?A.stepi():A.stepj();
    int bs = 2*b.step();
    BLASNAME(strsv) (BLASCM A.iscm()?BLASCH_LO:BLASCH_UP,
	A.iscm()?BLASCH_NT:BLASCH_T, A.isunit()?BLASCH_U:BLASCH_NU,
	BLASV(n),BLASP(A.cptr()),BLASV(lda),BLASP((float*)b.ptr()),BLASV(bs)
	BLAS1 BLAS1 BLAS1);
    BLASNAME(strsv) (BLASCM A.iscm()?BLASCH_LO:BLASCH_UP,
	A.iscm()?BLASCH_NT:BLASCH_T, A.isunit()?BLASCH_U:BLASCH_NU,
	BLASV(n),BLASP(A.cptr()),BLASV(lda),BLASP((float*)b.ptr()+1),BLASV(bs)
	BLAS1 BLAS1 BLAS1);
  }
#endif // FLOAT
#endif // BLAS

  template <class T, class Ta> void Tri_LDivEq(
      const GenUpperTriMatrix<Ta>& A, const VectorView<T>& b)
  {
    TMVAssert(b.size() == A.size());
#ifdef XDEBUG
    Matrix<Ta> A0(A);
    Vector<T> b0(b);
#endif
    if (b.size() > 0) {
      if (b.isconj()) Tri_LDivEq(A.Conjugate(),b.Conjugate());
      else 
#ifdef BLAS
	if (IsComplex(T()) && IsReal(Ta()))
	  BlasTri_LDivEq(A,b);
	else if ( !((A.isrm() && A.stepi()>0) || (A.iscm() && A.stepj()>0)) ) {
	  if (A.isunit()) {
	    UpperTriMatrix<Ta,UnitDiag,ColMajor> AA = A;
	    BlasTri_LDivEq(AA,b);
	  } else {
	    UpperTriMatrix<Ta,NonUnitDiag,ColMajor> AA = A;
	    BlasTri_LDivEq(AA,b);
	  }
	}
	else BlasTri_LDivEq(A,b);
#else
	NonBlasTri_LDivEq(A,b);
#endif
    }
#ifdef XDEBUG
    Vector<T> b2 = A0*b;
    if (Norm(b2-b0) > 0.001*Norm(A0)*Norm(b0)) {
      cerr<<"Tri_LDivEq: v/Upper\n";
      cerr<<"A = "<<Type(A)<<"  "<<A0<<endl;
      cerr<<"b = "<<Type(b)<<"  "<<b0<<endl;
      cerr<<"Done: b = "<<b<<endl;
      cerr<<"A*b = "<<b2<<endl;
      abort();
    }
#endif
  }

  template <class T, class Ta> void Tri_LDivEq(
      const GenLowerTriMatrix<Ta>& A, const VectorView<T>& b)
  {
    TMVAssert(b.size() == A.size());
#ifdef XDEBUG
    Matrix<Ta> A0(A);
    Vector<T> b0(b);
#endif
    if (b.size() > 0) {
      if (b.isconj()) Tri_LDivEq(A.Conjugate(),b.Conjugate());
      else 
#ifdef BLAS
	if (IsComplex(T()) && IsReal(Ta()))
	  BlasTri_LDivEq(A,b);
	else if ( !((A.isrm() && A.stepi()>0) || (A.iscm() && A.stepj()>0)) ) {
	  if (A.isunit()) {
	    LowerTriMatrix<Ta,UnitDiag,ColMajor> AA = A;
	    BlasTri_LDivEq(AA,b);
	  } else {
	    LowerTriMatrix<Ta,NonUnitDiag,ColMajor> AA = A;
	    BlasTri_LDivEq(AA,b);
	  }
	}
	else BlasTri_LDivEq(A,b);
#else
	NonBlasTri_LDivEq(A,b);
#endif
    }
#ifdef XDEBUG
    Vector<T> b2 = A0*b;
    if (Norm(b2-b0) > 0.001*Norm(A0)*Norm(b0)) {
      cerr<<"Tri_LDivEq: v/Lower\n";
      cerr<<"A = "<<Type(A)<<"  "<<A<<endl;
      cerr<<"b = "<<Type(b)<<"  "<<b0<<endl;
      cerr<<"Done: b = "<<b<<endl;
      cerr<<"A*b = "<<b2<<endl;
      abort();
    }
#endif
  }

#define InstFile "TMV_TriDiv_V.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


