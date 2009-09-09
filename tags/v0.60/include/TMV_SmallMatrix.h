///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// The Template Matrix/Vector Library for C++ was created by Mike Jarvis     //
// Copyright (C) 2007                                                        //
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
// along with this program in the file gpl.txt.                              //
//                                                                           //
// If not, write to:                                                         //
// The Free Software Foundation, Inc.                                        //
// 51 Franklin Street, Fifth Floor,                                          //
// Boston, MA  02110-1301, USA.                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


//---------------------------------------------------------------------------
//
// This file defines the TMV SmallMatrix class.
//
// Constructors:
//
//    SmallMatrix<T,M,N,stor,I>()
//        Makes a SmallMatrix with column size = M and row size = N
//        with _uninitialized_ values
//
//    SmallMatrix<T,M,N,stor,I>(T x)
//        Makes a SmallMatrix of size n with all values = x
//
//    SmallMatrix<T,M,N,stor,I>(const vector<vector<T> >& m)
//        Makes a SmallMatrix with a_ij = m[i][j]
//
//    SmallMatrix<T,M,N,stor,I>(const T* m)
//    SmallMatrix<T,M,N,stor,I>(const vector<T>& m)
//        Make a SmallMatrix which copies the elements of m.
//
//    SmallMatrixViewOf<M,N,stor>(T* m)
//    SmallMatrixViewOf<M,N,stor>(const T* m)
//        Returns a MatrixView of the elements in m, using the actual
//        elements m for the storage.  
//
// A SmallMatrix _is a_ GenMatrix.  ie. it inherits from GenMatrix.
// So everything you can do with a Matrix, you can also do with a SmallMatrix.
// The advantage is that many access functions are doable at compile time,
// since M,N,stor are template parameters.
//
// Another advantage, which may instead be a disadvantage in some 
// circumstances, is taht every calculation is done inline.  This can
// speed up some calculations, especially for small matrices.  
// But it may take significantly longer to compile, depending on what
// calculations you are doing with them.


#ifndef TMV_SmallMatrix_H
#define TMV_SmallMatrix_H

#include "TMV_Matrix.h"
#include "TMV_SmallVector.h"

namespace tmv {

#define MIN(M,N) (M<N ? M : N)

  template <class T, size_t M, size_t N, int Si, int Sj, bool C> 
    class GenSmallMatrix;
  template <class T, size_t M, size_t N, int Si, int Sj, bool C=false, IndexStyle I=CStyle> 
    class ConstSmallMatrixView;
  template <class T, size_t M, size_t N, int Si, int Sj, bool C=false, IndexStyle I=CStyle> 
    class SmallMatrixView;
  template <class T, size_t M, size_t N, StorageType S=RowMajor, IndexStyle I=CStyle> 
    class SmallMatrix;

  template <class T1, class T2, size_t M, size_t N, int Si1, int Sj1, int Si2, int Sj2, bool C1> 
    inline void Copy(
	const GenSmallMatrix<T1,M,N,Si1,Sj1,C1>& m1, 
	const SmallMatrixView<T2,M,N,Si2,Sj2,false>& m2);
  template <class T1, class T2, size_t M, size_t N, int Si1, int Sj1, int Si2, int Sj2, bool C1> 
    inline void Copy(
	const GenSmallMatrix<T1,M,N,Si1,Sj1,C1>& m1, 
	const SmallMatrixView<T2,M,N,Si2,Sj2,true>& m2)
    { Copy(m1.Conjugate(),m2.Conjugate()); }
  template <class T1, class T2, size_t M, size_t N, int Si1, int Sj1, bool C1> 
    inline void Copy(
	const GenSmallMatrix<T1,M,N,Si1,Sj1,C1>& m1, 
	T2* m2p, const int Si2, const int Sj2);

  template <class T, size_t M, size_t N> struct AssignableToSmallMatrix :
    virtual public AssignableToMatrix<T>
  {
    virtual void AssignTom(
	const SmallMatrixView<RealType(T),M,N,N,1,false>&) const = 0;
    virtual void AssignTom(
	const SmallMatrixView<RealType(T),M,N,1,M,false>&) const = 0;
    inline void AssignTom(
	const SmallMatrixView<RealType(T),M,N,N,1,true>& m2) const 
    { AssignTom(m2.Conjugate()); }
    inline void AssignTom(
	const SmallMatrixView<RealType(T),M,N,1,M,true>& m2) const
    { AssignTom(m2.Conjugate()); }
    virtual void AssignTom(
	const SmallMatrixView<ComplexType(T),M,N,N,1,true>&) const = 0;
    virtual void AssignTom(
	const SmallMatrixView<ComplexType(T),M,N,1,M,true>&) const = 0;
    virtual void AssignTom(
	const SmallMatrixView<ComplexType(T),M,N,N,1,false>&) const = 0;
    virtual void AssignTom(
	const SmallMatrixView<ComplexType(T),M,N,1,M,false>&) const = 0;
    virtual void DoAssignTom(RealType(T)* m0p, 
	const int Si, const int Sj) const = 0;
    virtual void DoAssignTom(ComplexType(T)* m0p, 
	const int Si, const int Sj, const bool C) const = 0;
    template <int Si, int Sj, bool C> inline void AssignTom(
	const SmallMatrixView<RealType(T),M,N,Si,Sj,C>& m0) const
    { TMVAssert(IsReal(T())); DoAssignTom(m0.ptr(),Si,Sj); }
    template <int Si, int Sj, bool C> inline void AssignTom(
	const SmallMatrixView<ComplexType(T),M,N,Si,Sj,C>& m0) const
    { DoAssignTom(m0.ptr(),Si,Sj,C); }
  };

  // In general, I want the base matrix for GenSmallMatrix to be
  // GenMatrix, but it's convenient to be able to switch it to BasseMatrix
  // for testing purposes.   (Basically, to make sure I overload everything
  // correctly and not ever rely on the GenMatrix version of a function.)
#define GENBASE

  template <class T, size_t M, size_t N, int Si, int Sj, bool C> 
    class GenSmallMatrix : 
#ifdef GENBASE
      public GenMatrix<T>,
#else
      public BaseMatrix<T>,
      private DivHelper<T>,
#endif
      virtual public AssignableToSmallMatrix<T,M,N>
    {

      public:

	//
	// Constructors
	//

#ifdef GENBASE
	inline GenSmallMatrix() {}
	inline GenSmallMatrix(const GenSmallMatrix<T,M,N,Si,Sj,C>&) {}
#else
	inline GenSmallMatrix() : inst(0) {}
	inline GenSmallMatrix(const GenSmallMatrix<T,M,N,Si,Sj,C>&) : inst(0) {}
#endif
	virtual inline ~GenSmallMatrix() {}

	//
	// Access Functions
	//

#ifdef GENBASE
	using GenMatrix<T>::row;
	using GenMatrix<T>::col;
	using GenMatrix<T>::diag;
#endif
	inline ConstSmallVectorView<T,N,Sj,C> row(size_t i) const 
	{ 
	  TMVAssert(i<M);
	  return ConstSmallVectorView<T,N,Sj,C>(cptr()+int(i)*Si); 
	}

	inline ConstSmallVectorView<T,M,Si,C> col(size_t j) const
	{
	  TMVAssert(j<N);
	  return ConstSmallVectorView<T,M,Si,C>(cptr()+int(j)*Sj); 
	}

	inline ConstSmallVectorView<T,MIN(N,M),Si+Sj,C> diag() const
	{
	  return ConstSmallVectorView<T,MIN(N,M),Si+Sj,C>(cptr());
	}

#ifndef GENBASE
	inline ConstVectorView<T> diag(int i) const
	{ return RegView().diag(i); }

	inline ConstVectorView<T> row(int i, size_t j1, size_t j2) const
	{ return RegView().row(i,j1,j2); }

	inline ConstVectorView<T> col(int j, size_t i1, size_t i2) const
	{ return RegView().col(j,i1,i2); }

	inline ConstVectorView<T> diag(int i, size_t j1, size_t j2) const
	{ return RegView().diag(i,j1,j2); }
#endif

	inline T operator()(size_t i, size_t j) const
	{ 
	  TMVAssert(i<M);
	  TMVAssert(j<N);
	  return cref(i,j);
	}

	inline ConstSmallVectorView<T,N,Sj,C> operator[](size_t i) const
	{ 
	  TMVAssert(i<M);
	  return row(i); 
	}

	inline void AssignTom(
	    const SmallMatrixView<RealType(T),M,N,N,1,false>& m2) const
	{ TMVAssert(IsReal(T())); Copy(*this,m2); }
	inline void AssignTom(
	    const SmallMatrixView<RealType(T),M,N,1,M,false>& m2) const
	{ TMVAssert(IsReal(T())); Copy(*this,m2); }
	inline void AssignTom(
	    const SmallMatrixView<ComplexType(T),M,N,N,1,true>& m2) const 
	{ Copy(*this,m2); }
	inline void AssignTom(
	    const SmallMatrixView<ComplexType(T),M,N,1,M,true>& m2) const
	{ Copy(*this,m2); }
	inline void AssignTom(
	    const SmallMatrixView<ComplexType(T),M,N,N,1,false>& m2) const
	{ Copy(*this,m2); }
	inline void AssignTom(
	    const SmallMatrixView<ComplexType(T),M,N,1,M,false>& m2) const
	{ Copy(*this,m2); }
	inline void DoAssignTom(RealType(T)* m2p, 
	    const int Si2, const int Sj2) const
	{ TMVAssert(IsReal(T())); Copy(*this,m2p,Si2,Sj2); }
	inline void DoAssignTom(ComplexType(T)* m2p, 
	    const int Si2, const int Sj2, const bool C2) const
	{ 
	  if (C2) Copy(Conjugate(),m2p,Si2,Sj2); 
	  else Copy(*this,m2p,Si2,Sj2); 
	}
	inline void AssignToM(const MatrixView<RealType(T)>& m0) const
	{
	  TMVAssert(m0.colsize() == M);
	  TMVAssert(m0.rowsize() == N);
	  TMVAssert(IsReal(T()));
	  DoAssignTom(m0.ptr(),m0.stepi(),m0.stepj()); 
	}
	inline void AssignToM(const MatrixView<ComplexType(T)>& m0) const
	{
	  TMVAssert(m0.colsize() == M);
	  TMVAssert(m0.rowsize() == N);
	  if (m0.isconj()) Copy(Conjugate(),m0.ptr(),m0.stepi(),m0.stepj());
	  else Copy(*this,m0.ptr(),m0.stepi(),m0.stepj());
	}

	//
	// SubMatrix
	//

#ifdef GENBASE
	using GenMatrix<T>::SubMatrix;
	using GenMatrix<T>::SubVector;
	using GenMatrix<T>::ColPair;
	using GenMatrix<T>::RowPair;
	using GenMatrix<T>::Cols;
	using GenMatrix<T>::Rows;
#else
	inline ConstMatrixView<T> SubMatrix(
	    int i1, int i2, int j1, int j2) const
	{ return RegView().SubMatrix(i1,i2,j1,j2); }

	inline ConstMatrixView<T> SubMatrix(
	    int i1, int i2, int j1, int j2, int istep, int jstep) const
	{ return RegView().SubMatrix(i1,i2,j1,j2,istep,jstep); }

	inline ConstVectorView<T> SubVector(
	    int i, int j, int istep, int jstep, size_t s) const
	{ return RegView().SubVector(i,j,istep,jstep,s); }

	inline ConstMatrixView<T> ColPair(int j1, int j2) const
	{ return RegView().ColPair(j1,j2); }

	inline ConstMatrixView<T> RowPair(int i1, int i2) const
	{ return RegView().RowPair(i1,i2); }

	inline ConstMatrixView<T> Cols(int j1, int j2) const
	{ return RegView().Cols(j1,j2); }

	inline ConstMatrixView<T> Rows(int i1, int i2) const
	{ return RegView().Rows(i1,i2); }
#endif

	//
	// Views
	//

	inline ConstMatrixView<T> RegView() const
	{ return ConstMatrixView<T>(cptr(),M,N,Si,Sj,stor(),ct()); }

	inline ConstSmallMatrixView<T,M,N,Si,Sj,C> View() const
	{ return ConstSmallMatrixView<T,M,N,Si,Sj,C>(cptr()); }

	inline ConstSmallMatrixView<T,N,M,Sj,Si,C> Transpose() const
	{ return ConstSmallMatrixView<T,N,M,Sj,Si,C>(cptr()); }

	inline ConstSmallMatrixView<T,M,N,Si,Sj,!C> Conjugate() const
	{ return ConstSmallMatrixView<T,M,N,Si,Sj,!C>(cptr()); }

	inline ConstSmallMatrixView<T,N,M,Sj,Si,!C> Adjoint() const
	{ return ConstSmallMatrixView<T,N,M,Sj,Si,!C>(cptr()); }

	inline bool CanLinearize() const
	{ return ( (Si==1 && Sj==int(M)) || (Sj==1 && Si==int(N)) ); }

	inline ConstSmallVectorView<T,M*N,1,C> ConstLinearView() const
	{
	  TMVAssert(CanLinearize());
	  return ConstSmallVectorView<T,M*N,1,C>(cptr());
	}

	inline SmallMatrixView<T,M,N,Si,Sj,C> NonConst() const
	{
	  return SmallMatrixView<T,M,N,Si,Sj,C>(const_cast<T*>(cptr())
	      FIRSTLAST1(cptr(),M>0?row(M-1).end().GetP():cptr()));
	}

	//
	// Functions of Matrix
	//

	inline T Trace() const
	{ return diag().SumElements(); }

	inline RealType(T) Norm() const 
	{ return NormF(); }

	inline RealType(T) NormF() const
	{ return tmv::SQRT(NormSq()); }

	// NormF()^2
	inline RealType(T) NormSq() const
	{ 
	  if (CanLinearize()) return ConstLinearView().NormSq();
	  else {
	    RealType(T) sum(0);
	    if (Sj < Si) 
	      for(size_t i=0;i<M;++i) sum += row(i).NormSq();
	    else
	      for(size_t j=0;j<N;++j) sum += col(j).NormSq();
	    return sum;
	  }
	}

	// 1-Norm = max_j (sum_i |a_ij|)
	inline RealType(T) Norm1() const
	{
	  RealType(T) max(0);
	  for(size_t j=0;j<N;++j) {
	    RealType(T) temp = col(j).Norm1();
	    if (temp > max) max = temp;
	  }
	  return max;
	}

	// inf-Norm = max_i (sum_j |a_ij|)
	inline RealType(T) NormInf() const
	{ return Transpose().Norm1(); }

	// = max_i,j (|a_ij|)
	inline RealType(T) MaxAbsElement() const
	{
	  if (CanLinearize()) return ConstLinearView().MaxAbsElement();
	  else {
	    RealType(T) max(0);
	    if (Sj < Si) 
	      for(size_t i=0;i<M;++i) {
		RealType(T) temp = row(i).NormInf();
		if (temp > max) max = temp;
	      }
	    else
	      for(size_t j=0;j<N;++j) {
		RealType(T) temp = col(j).NormInf();
		if (temp > max) max = temp;
	      }
	    return max;
	  }
	}
#ifndef GENBASE
	inline T Det() const
	{ return DivHelper<T>::Det(); }
	inline RealType(T) Norm2() const
	{ return DivHelper<T>::Norm2(); }
	inline bool Singular() const
	{ return DivHelper<T>::Singular(); }
	inline RealType(T) Condition() const
	{ return DivHelper<T>::Condition(); }
#endif

	// 
	// Division Control
	//

#ifdef GENBASE
	using GenMatrix<T>::Inverse;
#else
	inline QuotXM<T,T> Inverse() const
	{ 
	  if (!inst) inst.reset(new Matrix<T>(*this));
	  return QuotXM<T>(T(1),*inst);
	}
	using DivHelper<T>::Inverse;
	using DivHelper<T>::InverseATA;
	inline void Inverse(const MatrixView<T>& minv) const
	{ DivHelper<T>::Inverse(minv); }
	inline void InverseATA(const MatrixView<T>& minv) const
	{ DivHelper<T>::InverseATA(minv); }
#endif

	template <class T2, int Si2, int Sj2> inline void Inverse(
	    const SmallMatrixView<T2,N,M,Si2,Sj2>& minv) const
	{ GenMatrix<T>::Inverse(minv.RegView()); }

	template <class T2, StorageType S2, IndexStyle I2> inline void Inverse(
	    SmallMatrix<T2,N,M,S2,I2>& minv) const
	{ GenMatrix<T>::Inverse(minv.RegView()); }

#ifndef GENBASE
	using DivHelper<T>::DivideInPlace;
	using DivHelper<T>::SaveDiv;
	using DivHelper<T>::SetDiv;
	using DivHelper<T>::UnSetDiv;
	using DivHelper<T>::ReSetDiv;
	using DivHelper<T>::GetDiv;
	using DivHelper<T>::CheckDecomp;
	inline void DivideUsing(DivType dt) const
	{
	  TMVAssert(dt == LU || dt == QR || dt == QRP || SV_Type(dt));
	  DivHelper<T>::DivideUsing(dt);
	}
	inline const LUDiv<T>& LUD() const
	{
	  TMVAssert(GetDiv());
	  TMVAssert(dynamic_cast<const LUDiv<T>*>(GetDiv()));
	  return *dynamic_cast<const LUDiv<T>*>(GetDiv());
	}
	inline const QRDiv<T>& QRD() const
	{
	  TMVAssert(GetDiv());
	  TMVAssert(dynamic_cast<const QRDiv<T>*>(GetDiv()));
	  return *dynamic_cast<const QRDiv<T>*>(GetDiv());
	}
	inline const QRPDiv<T>& QRPD() const
	{
	  TMVAssert(GetDiv());
	  TMVAssert(dynamic_cast<const QRPDiv<T>*>(GetDiv()));
	  return *dynamic_cast<const QRPDiv<T>*>(GetDiv());
	}
	inline const SVDiv<T>& SVD() const
	{
	  TMVAssert(GetDiv());
	  TMVAssert(dynamic_cast<const SVDiv<T>*>(GetDiv()));
	  return *dynamic_cast<const SVDiv<T>*>(GetDiv());
	}
	using DivHelper<T>::LDivEq;
	using DivHelper<T>::RDivEq;
	using DivHelper<T>::LDiv;
	using DivHelper<T>::RDiv;

	inline auto_ptr<BaseMatrix<T> > NewCopy() const
	{ return RegView().NewCopy(); }
	inline auto_ptr<BaseMatrix<T> > NewView() const
	{ return RegView().NewView(); }
	inline auto_ptr<BaseMatrix<T> > NewTranspose() const 
	{ return RegView().NewTranspose(); }
	inline auto_ptr<BaseMatrix<T> > NewConjugate() const
	{ return RegView().NewConjugate(); }
	inline auto_ptr<BaseMatrix<T> > NewAdjoint() const
	{ return RegView().NewAdjoint(); }
	inline auto_ptr<BaseMatrix<T> > NewInverse() const
	{ return RegView().NewInverse(); }
#endif

	//
	// I/O
	//

	inline void Write(std::ostream& os) const
	{
	  os << M << "  " << N << std::endl;
	  for(size_t i=0;i<M;++i) {
	    os << "( ";
	    typename SmallVectorView<T,N,Sj,C>::const_iterator it =
	      row(i).begin();
	    for(size_t k=N;k>0;--k,++it) os << ' '<<*it<<' ';
	    os << " )\n";
	  }
	}
	inline void Write(std::ostream& os, RealType(T) thresh) const
	{
	  os << M << "  " << N << std::endl;
	  for(size_t i=0;i<M;++i) {
	    os << "( ";
	    typename SmallVectorView<T,N,Sj,C>::const_iterator it =
	      row(i).begin();
	    for(size_t k=N;k>0;--k,++it) 
	      os << ' '<<(ABS(*it)<thresh ? T(0) : *it)<<' ';
	    os << " )\n";
	  }
	}

	virtual const T* cptr() const = 0;
	inline size_t colsize() const { return M; }
	inline size_t rowsize() const { return N; }
	inline int stepi() const { return Si; }
	inline int stepj() const { return Sj; }
	inline bool isrm() const { return Sj==1; }
	inline bool iscm() const { return Si==1; }
	inline bool isconj() const { return C && IsComplex(T()); }
	inline StorageType stor() const
	{ return isrm() ? RowMajor : iscm() ? ColMajor : NoMajor; }
	inline ConjType ct() const 
	{ return isconj() ? Conj : NonConj; }
	inline size_t ls() const
	{ return CanLinearize() ? M*N : 0; }

      protected :

	inline T cref(size_t i, size_t j) const
	{
	  TMVAssert(i<M && j<N);
	  const T* mi = cptr() + int(i)*Si + int(j)*Sj;
	  return isconj() ? CONJ(*mi) : *mi;
	}

#ifdef GENBASE
	using GenMatrix<T>::NewDivider;
#else
	inline void NewDivider() const { RegView().NewDivider(); }
	inline const BaseMatrix<T>& GetMatrix() const { return *this; }
#endif

      private :

	inline void operator=(const GenMatrix<T>&) { TMVAssert(FALSE); }

#ifndef GENBASE
	mutable auto_ptr<Matrix<T> > inst;
#endif

    }; // GenMatrix

  template <class T, size_t M, size_t N, int Si, int Sj, bool C, IndexStyle I> 
    class ConstSmallMatrixView : 
      public GenSmallMatrix<T,M,N,Si,Sj,C>
    {
      public :

	inline ConstSmallMatrixView(
	    const ConstSmallMatrixView<T,M,N,Si,Sj,C,I>& rhs) :
	  itsm(rhs.itsm) {}

	inline ConstSmallMatrixView(const GenSmallMatrix<T,M,N,Si,Sj,C>& rhs) :
	  itsm(rhs.cptr()) {}

	inline ConstSmallMatrixView(const GenMatrix<T>& rhs) :
	  itsm(rhs.cptr()) 
	{
	  TMVAssert(rhs.colsize() == M);
	  TMVAssert(rhs.rowsize() == N);
	  TMVAssert(rhs.stepi() == Si);
	  TMVAssert(rhs.stepj() == Sj);
	  TMVAssert(rhs.isconj() == C);
	}

	inline ConstSmallMatrixView(const T* _m) : itsm(_m) {}

	virtual inline ~ConstSmallMatrixView() {}

	virtual inline const T* cptr() const { return itsm; }

      private :

	const T*const itsm;

	inline void operator=(const ConstSmallMatrixView<T,M,N,Si,Sj,C,I>&) 
	{ TMVAssert(FALSE); }

    }; // ConstSmallMatrixView

  template <class T, size_t M, size_t N, int Si, int Sj, bool C> 
    class ConstSmallMatrixView<T,M,N,Si,Sj,C,FortranStyle> : 
    public ConstSmallMatrixView<T,M,N,Si,Sj,C,CStyle>
    {
      public :

	inline ConstSmallMatrixView(
	    const ConstSmallMatrixView<T,M,N,Si,Sj,C>& rhs) :
	  ConstSmallMatrixView<T,M,N,Si,Sj,C,CStyle>(rhs) {}

	inline ConstSmallMatrixView(const GenSmallMatrix<T,M,N,Si,Sj,C>& rhs) :
	  ConstSmallMatrixView<T,M,N,Si,Sj,C,CStyle>(rhs) {}

	inline ConstSmallMatrixView(const GenMatrix<T>& rhs) :
	  ConstSmallMatrixView<T,M,N,Si,Sj,C,CStyle>(rhs) {}

	inline ConstSmallMatrixView(const T* _m) :
	  ConstSmallMatrixView<T,M,N,Si,Sj,C,CStyle>(_m) {}

	virtual inline ~ConstSmallMatrixView() {}


	// 
	// Access
	//

	inline T operator()(size_t i, size_t j)
	{ 
	  TMVAssert(i>0 && i<=M);
	  TMVAssert(j>0 && j<=N);
	  return GenSmallMatrix<T,M,N,Si,Sj,C>::cref(i-1,j-1);
	}

	inline ConstSmallVectorView<T,N,Sj,C,FortranStyle> row(size_t i) const 
	{ 
	  TMVAssert(i>0 && i<=M);
	  return GenSmallMatrix<T,M,N,Si,Sj,C>::row(i-1);
	}

	inline ConstSmallVectorView<T,M,Si,C,FortranStyle> col(size_t j) const
	{
	  TMVAssert(j>0 && j<=N);
	  return GenSmallMatrix<T,M,N,Si,Sj,C>::col(j-1);
	}

	inline ConstSmallVectorView<T,MIN(N,M),Si+Sj,C,FortranStyle> 
	  diag() const
	{ return GenSmallMatrix<T,M,N,Si,Sj,C>::diag(); }

	inline ConstVectorView<T,FortranStyle> diag(int i) const
	{
	  TMVAssert(i<=int(N) && i>=-int(M));
	  return RegView().diag(i); 
	}

	inline ConstVectorView<T,FortranStyle> row(
	    size_t i, size_t j1, size_t j2) const
	{ 
	  TMVAssert(i>0 && i<=M);
	  TMVAssert(j1>0 && j1<=j2 && j2<=N);
	  return RegView().row(i,j1,j2); 
	}

	inline ConstVectorView<T,FortranStyle> col(
	    size_t j, size_t i1, size_t i2) const
	{
	  TMVAssert(j>0 && j<=N);
	  TMVAssert(i1>0 && i1<=i2 && i2<=M);
	  return RegView().col(j,i1,i2); 
	}

	inline ConstVectorView<T,FortranStyle> diag(
	    int i, size_t j1, size_t j2) const
	{
	  TMVAssert(i<=int(N) && i>=-int(M));
	  return RegView().diag(i,j1,j2); 
	}

	inline ConstSmallVectorView<T,N,Sj,C,FortranStyle> operator[](
	    size_t i) const
	{ return row(i); }

	//
	// SubMatrix
	//

	inline ConstMatrixView<T,FortranStyle> SubMatrix(
	    int i1, int i2, int j1, int j2) const
	{
	  TMVAssert(RegView().OKSubMatrix(i1,i2,j1,j2,1,1));
	  return RegView().SubMatrix(i1,i2,j1,j2);
	}

	inline ConstMatrixView<T,FortranStyle> SubMatrix(
	    int i1, int i2, int j1, int j2, int istep, int jstep) const
	{
	  TMVAssert(RegView().OKSubMatrix(i1,i2,j1,j2,istep,jstep));
	  return RegView().SubMatrix(i1,i2,j1,j2,istep,jstep);
	}

	inline ConstVectorView<T,FortranStyle> SubVector(
	    size_t i, size_t j, int istep, int jstep, size_t s) const
	{
	  TMVAssert(RegView().OKSubVector(i,j,istep,jstep,s));
	  return RegView().SubVector(i,j,istep,jstep,s);
	}

	inline ConstMatrixView<T,FortranStyle> ColPair(
	    size_t j1, size_t j2) const
	{
	  TMVAssert(j1 > 0 && j1 <= N);
	  TMVAssert(j2 > 0 && j2 <= N);
	  return RegView().ColPair(j1,j2);
	}

	inline ConstMatrixView<T,FortranStyle> RowPair(
	    size_t i1, size_t i2) const
	{
	  TMVAssert(i1 > 0 && i1 <= M);
	  TMVAssert(i2 > 0 && i2 <= M);
	  return RegView().RowPair(i1,i2);
	}

	inline ConstMatrixView<T,FortranStyle> Cols(
	    size_t j1, size_t j2) const
	{
	  TMVAssert(j1 > 0 && j1 <= j2);
	  TMVAssert(j2 <= N);
	  return RegView().Cols(j1,j2);
	}

	inline ConstMatrixView<T,FortranStyle> Rows(
	    size_t i1, size_t i2) const
	{
	  TMVAssert(i1 > 0 && i1 <= i2);
	  TMVAssert(i2 <= M);
	  return RegView().Rows(i1,i2);
	}

	inline ConstMatrixView<RealType(T),FortranStyle> Real() const
	{ return RegView().Real(); }

	inline ConstMatrixView<RealType(T),FortranStyle> Imag() const
	{ return RegView().Imag(); }


	//
	// Views
	//

	inline ConstSmallMatrixView<T,M,N,Si,Sj,C,FortranStyle> View() const
	{ return *this; }

	inline ConstMatrixView<T,FortranStyle> RegView() const
	{ 
	  return ConstMatrixView<T,FortranStyle>(
	      this->cptr(),M,N,Si,Sj,this->stor(),this->ct());
	}

	inline ConstSmallMatrixView<T,N,M,Sj,Si,C,FortranStyle> 
	  Transpose() const
	{ return GenSmallMatrix<T,M,N,Si,Sj,C>::Transpose(); }

	inline ConstSmallMatrixView<T,M,N,Si,Sj,!C,FortranStyle> 
	  Conjugate() const
	{ return GenSmallMatrix<T,M,N,Si,Sj,C>::Conjugate(); }

	inline ConstSmallMatrixView<T,N,M,Sj,Si,!C,FortranStyle> Adjoint() const
	{ return GenSmallMatrix<T,M,N,Si,Sj,C>::Adjoint(); }

	inline ConstSmallVectorView<T,M*N,1,C,FortranStyle> 
	  ConstLinearView() const
	{ return GenSmallMatrix<T,M,N,Si,Sj,C>::ConstLinearView(); }

	inline SmallMatrixView<T,M,N,Si,Sj,C,FortranStyle> NonConst() const
	{ return GenSmallMatrix<T,M,N,Si,Sj,C>::NonConst(); }

      private :

	inline void operator=(
	    const ConstSmallMatrixView<T,M,N,Si,Sj,C,FortranStyle>&) 
	{ TMVAssert(FALSE); }

    }; // FortranStyle ConstSmallMatrixView

  template <class T, size_t M, size_t N, int Si, int Sj, bool C, IndexStyle I> 
    class SmallMatrixView : 
      public GenSmallMatrix<T,M,N,Si,Sj,C>
    {

      public:

	//
	// Constructors
	//

	inline SmallMatrixView(const SmallMatrixView<T,M,N,Si,Sj,C,I>& rhs) : 
	  itsm(rhs.itsm) DEFFIRSTLAST(rhs.first,rhs.last)
	{ TMVAssert(I==CStyle); }

	inline SmallMatrixView(T* _m PARAMFIRSTLAST(T) ) :
	  itsm(_m) DEFFIRSTLAST(_first,_last)
	{ TMVAssert(I==CStyle); }

	inline SmallMatrixView(const MatrixView<T>& rhs) :
	  itsm(rhs.ptr()) DEFFIRSTLAST(rhs.first,rhs.last)
	{
	  TMVAssert(rhs.colsize() == M);
	  TMVAssert(rhs.rowsize() == N);
	  TMVAssert(rhs.stepi() == Si);
	  TMVAssert(rhs.stepj() == Sj);
	  TMVAssert(rhs.isconj() == C);
	  TMVAssert(I==CStyle);
	}

	virtual inline ~SmallMatrixView() {} 

	//
	// Op=
	//

	inline const SmallMatrixView<T,M,N,Si,Sj,C,I>& operator=(
	    const SmallMatrixView<T,M,N,Si,Sj,C,I>& m2) const
	{ Copy(m2,*this); return *this; }

	inline const SmallMatrixView<T,M,N,Si,Sj,C,I>& operator=(
	    const GenSmallMatrix<RealType(T),M,N,Si,Sj,C>& m2) const
	{ m2.AssignTom(*this); return *this; }

	inline const SmallMatrixView<T,M,N,Si,Sj,C,I>& operator=(
	    const GenSmallMatrix<ComplexType(T),M,N,Si,Sj,C>& m2) const
	{
	  TMVAssert(IsComplex(T())); 
	  m2.AssignTom(*this); 
	  return *this; 
	}

	template <class T2, int Si2, int Sj2, bool C2>
	  inline const SmallMatrixView<T,M,N,Si,Sj,C,I>& operator=(
	      const GenSmallMatrix<T2,M,N,Si2,Sj2,C2>& m2) const
	  { Copy(m2,*this); return *this; }

	inline const SmallMatrixView<T,M,N,Si,Sj,C,I>& operator=(
	    const GenMatrix<RealType(T)>& m2) const
	{ 
	  TMVAssert(m2.colsize() == M);
	  TMVAssert(m2.rowsize() == N);
	  m2.AssignToM(RegView()); 
	  return *this; 
	}

	inline const SmallMatrixView<T,M,N,Si,Sj,C,I>& operator=(
	    const GenMatrix<ComplexType(T)>& m2) const
	{
	  TMVAssert(m2.colsize() == M);
	  TMVAssert(m2.rowsize() == N);
	  TMVAssert(IsComplex(T()));
	  m2.AssignToM(RegView()); 
	  return *this; 
	}

	template <class T2>
	  inline const SmallMatrixView<T,M,N,Si,Sj,C,I>& operator=(
	      const GenMatrix<T2>& m2) const
	  {
	    TMVAssert(m2.colsize() == M);
	    TMVAssert(m2.rowsize() == N);
	    Copy(m2,RegView()); 
	    return *this; 
	  }

	inline const SmallMatrixView<T,M,N,Si,Sj,C,I>& operator=(T x) const 
	{ return SetToIdentity(x); }

	inline const SmallMatrixView<T,M,N,Si,Sj,C,I>& operator=(
	    const AssignableToSmallMatrix<RealType(T),M,N>& m2) const
	{ m2.AssignTom(*this); return *this; }

	inline const SmallMatrixView<T,M,N,Si,Sj,C,I>& operator=(
	    const AssignableToSmallMatrix<ComplexType(T),M,N>& m2) const
	{ m2.AssignTom(*this); return *this; }

	//
	// Access
	//

	inline RefType(T) operator()(size_t i,size_t j) const 
	{ 
	  TMVAssert(I==CStyle);
	  TMVAssert(i<M);
	  TMVAssert(j<N);
	  return ref(i,j); 
	}

	inline SmallVectorView<T,N,Sj,C,I> operator[](size_t i) const 
	{ 
	  TMVAssert(I==CStyle);
	  TMVAssert(i<M);
	  return row(i); 
	}

	typedef RefType(T) reference;

	inline SmallVectorView<T,N,Sj,C,I> row(size_t i) const
	{
	  TMVAssert(I==CStyle);
	  TMVAssert(i<M);
	  return SmallVectorView<T,N,Sj,C,I>(ptr()+int(i)*Si FIRSTLAST );
	}

	inline SmallVectorView<T,M,Si,C,I> col(size_t j) const
	{
	  TMVAssert(I==CStyle);
	  TMVAssert(j<N);
	  return SmallVectorView<T,M,Si,C,I>(ptr()+int(j)*Sj FIRSTLAST ); 
	}

	inline SmallVectorView<T,MIN(N,M),Si+Sj,C,I> diag() const
	{
	  TMVAssert(I==CStyle);
	  return SmallVectorView<T,MIN(N,M),Si+Sj,C,I>(ptr() FIRSTLAST);
	}

	inline VectorView<T,I> diag(int i) const
	{ return RegView().diag(i); }

	inline VectorView<T,I> row(size_t i, size_t j1, size_t j2) const 
	{ return RegView().row(i,j1,j2); }

	inline VectorView<T,I> col(size_t j, size_t i1, size_t i2) const
	{ return RegView().col(j,i1,i2); }

	inline VectorView<T,I> diag(int i, size_t j1, size_t j2) const
	{ return RegView().diag(i,j1,j2); }

	//
	// Modifying Functions
	//

	inline const SmallMatrixView<T,M,N,Si,Sj,C,I>& Zero() const 
	{ return SetAllTo(T(0)); }

	inline const SmallMatrixView<T,M,N,Si,Sj,C,I>& Clip(
	    RealType(T) thresh) const
	{
	  TMVAssert(I==CStyle);
	  if (this->CanLinearize()) LinearView().Clip(thresh);
	  else if (Sj < Si)
	    for(size_t i=0;i<M;++i) row(i).Clip(thresh);
	  else
	    for(size_t j=0;j<N;++j) col(j).Clip(thresh);
	  return *this;
	}

	inline const SmallMatrixView<T,M,N,Si,Sj,C,I>& SetAllTo(T x) const
	{
	  TMVAssert(I==CStyle);
	  if (this->CanLinearize()) LinearView().SetAllTo(x);
	  else if (Sj < Si)
	    for(size_t i=0;i<M;++i) row(i).SetAllTo(x);
	  else
	    for(size_t j=0;j<N;++j) col(j).SetAllTo(x);
	  return *this;
	}

	inline const SmallMatrixView<T,M,N,Si,Sj,C,I>& TransposeSelf() const
	{
	  TMVAssert(I==CStyle);
	  TMVAssert(M == N);
	  for(size_t i=1;i<N;++i)
	    for(size_t j=0;j<i;++j)
	      std::swap(ref(i,j),ref(j,i));
	  return *this;
	}

	inline const SmallMatrixView<T,M,N,Si,Sj,C,I>& ConjugateSelf() const
	{
	  TMVAssert(I==CStyle);
	  if (IsComplex(T())) 
	    if (this->CanLinearize()) LinearView().ConjugateSelf();
	    else if (Sj < Si)
	      for(size_t i=0;i<M;++i) row(i).ConjugateSelf();
	    else
	      for(size_t j=0;j<N;++j) col(j).ConjugateSelf();
	  return *this;
	}

	inline const SmallMatrixView<T,M,N,Si,Sj,C,I>& SetToIdentity(
	    T x=T(1)) const
	{
	  TMVAssert(I==CStyle);
	  TMVAssert(M == N);
	  Zero(); diag().SetAllTo(x);
	  return *this;
	}

	inline const SmallMatrixView<T,M,N,Si,Sj,C,I>& SwapRows(
	    size_t i1, size_t i2) const
	{
	  TMVAssert(I==CStyle);
	  TMVAssert(i1 < M && i2 < M);
	  if (i1!=i2) Swap(row(i1),row(i2));
	  return *this;
	}

	inline const SmallMatrixView<T,M,N,Si,Sj,C,I>& SwapCols(
	    size_t j1, size_t j2) const
	{
	  TMVAssert(I==CStyle);
	  TMVAssert(j1 < N && j2 < N);
	  if (j1!=j2) Swap(col(j1),col(j2));
	  return *this;
	}

	inline const SmallMatrixView<T,M,N,Si,Sj,C,I>& PermuteRows(
	    const size_t* p, size_t i1, size_t i2) const
	{
	  TMVAssert(I==CStyle);
	  TMVAssert(i2<=M);
	  TMVAssert(i1<=i2);
	  const size_t* pi = p+i1;
	  for(size_t i=i1;i<i2;++i,++pi) {
	    TMVAssert(*pi < M);
	    SwapRows(i,*pi);
	  }
	  return *this;
	}

	inline const SmallMatrixView<T,M,N,Si,Sj,C,I>& PermuteRows(
	    const size_t* p) const
	{ return PermuteRows(p,0,M); }

	inline const SmallMatrixView<T,M,N,Si,Sj,C,I>& PermuteCols(
	    const size_t* p, size_t j1, size_t j2) const
	{ Transpose().PermuteRows(p,j1,j2); return *this; }

	inline const SmallMatrixView<T,M,N,Si,Sj,C,I>& PermuteCols(
	    const size_t* p) const
	{ return PermuteCols(p,0,N); }

	const SmallMatrixView<T,M,N,Si,Sj,C,I>& ReversePermuteRows(
	    const size_t* p, size_t i1, size_t i2) const
	{
	  TMVAssert(I==CStyle);
	  TMVAssert(i2<=M);
	  TMVAssert(i1<=i2);
	  const size_t* pi = p+i2;
	  for(size_t i=i1;i<i2;) {
	    --i; --pi;
	    TMVAssert(*pi < M);
	    SwapRows(i,*pi);
	  }
	  return *this;
	}

	inline const SmallMatrixView<T,M,N,Si,Sj,C,I>& ReversePermuteRows(
	    const size_t* p) const
	{ return ReversePermuteRows(p,0,M); }

	inline const SmallMatrixView<T,M,N,Si,Sj,C,I>& ReversePermuteCols(
	    const size_t* p, size_t j1, size_t j2) const
	{ Transpose().ReversePermuteRows(p,j1,j2); return *this; }

	inline const SmallMatrixView<T,M,N,Si,Sj,C,I>& ReversePermuteCols(
	    const size_t* p) const
	{ return ReversePermuteCols(p,0,N); }

	//
	// SubMatrix
	//

	inline MatrixView<T,I> SubMatrix(int i1, int i2, int j1, int j2) const
	{
	  TMVAssert(RegView().OKSubMatrix(i1,i2,j1,j2,1,1));
	  return RegView().SubMatrix(i1,i2,j1,j2);
	}

	inline MatrixView<T,I> SubMatrix(
	    int i1, int i2, int j1, int j2, int istep, int jstep) const
	{
	  TMVAssert(RegView().OKSubMatrix(i1,i2,j1,j2,istep,jstep));
	  return RegView().SubMatrix(i1,i2,j1,j2,istep,jstep);
	}

	inline VectorView<T,I> SubVector(
	    size_t i, size_t j, int istep, int jstep, size_t size) const
	{
	  TMVAssert(RegView().OKSubVector(i,j,istep,jstep,size));
	  return RegView().SubVector(i,j,istep,jstep,size);
	}

	inline MatrixView<T,I> ColPair(size_t j1, size_t j2) const
	{
	  TMVAssert(j1<N && j2<N);
	  return RegView().ColPair(j1,j2);
	}

	inline MatrixView<T,I> RowPair(size_t i1, size_t i2) const
	{
	  TMVAssert(i1<M && i2<M);
	  return RegView().RowPair(i1,i2);
	}

	inline MatrixView<T,I> Cols(size_t j1, size_t j2) const
	{
	  TMVAssert(j1<=j2);
	  TMVAssert(j2<=N);
	  return RegView().Cols(j1,j2);
	}

	inline MatrixView<T,I> Rows(size_t i1, size_t i2) const
	{
	  TMVAssert(i1<=i2);
	  TMVAssert(i2<=M);
	  return RegView().Rows(i1,i2);
	}

	inline MatrixView<RealType(T)> Real() const
	{ return RegView().Real(); }

	inline MatrixView<RealType(T)> Imag() const
	{ return RegView().Imag(); }


	//
	// Views
	//

	inline SmallMatrixView<T,M,N,Si,Sj,C,I> View() const
	{ return *this; }

	inline MatrixView<T,I> RegView() const
	{
	  return MatrixView<T,I>(ptr(),M,N,Si,Sj,this->stor(),
	      this->ct(),(this->CanLinearize()?M*N:0) FIRSTLAST);
	}

	inline SmallMatrixView<T,N,M,Sj,Si,C,I> Transpose() const
	{ return SmallMatrixView<T,N,M,Sj,Si,C,I>(ptr() FIRSTLAST); }

	inline SmallMatrixView<T,M,N,Si,Sj,!C,I> Conjugate() const
	{ return SmallMatrixView<T,M,N,Si,Sj,!C,I>(ptr() FIRSTLAST); }

	inline SmallMatrixView<T,N,M,Sj,Si,!C,I> Adjoint() const
	{ return SmallMatrixView<T,N,M,Sj,Si,!C,I>(ptr() FIRSTLAST); }

	inline SmallVectorView<T,M*N,1,C,I> LinearView() const
	{
	  TMVAssert(this->CanLinearize());
	  return SmallVectorView<T,M*N,1,C,I>(ptr() FIRSTLAST );
	}


	//
	// I/O
	//

	inline void Read(std::istream& is) const
	{ RegView().Read(); }

	inline const T* cptr() const { return itsm; }
	inline T* ptr() const { return itsm; }

      protected:

	inline RefType(T) ref(size_t i, size_t j) const
	{
	  TMVAssert(I==CStyle);
	  TMVAssert(i<M && j<N);
	  T* mi = ptr() + int(i)*Si + int(j)*Sj;
	  return REF(mi,this->ct());
	}

      private:

	T*const itsm;

#ifdef TMVFLDEBUG
      public:
	const T*const first;
	const T*const last;
#endif

    }; // SmallMatrixView

  template <class T, size_t M, size_t N, int Si, int Sj, bool C> 
    class SmallMatrixView<T,M,N,Si,Sj,C,FortranStyle> : 
    public SmallMatrixView<T,M,N,Si,Sj,C,CStyle>
    {

      public:

	//
	// Constructors
	//

	inline SmallMatrixView(
	    const SmallMatrixView<T,M,N,Si,Sj,C,CStyle>& rhs) : 
	  SmallMatrixView<T,M,N,Si,Sj,C,CStyle>(rhs) {}

	inline SmallMatrixView(
	    const SmallMatrixView<T,M,N,Si,Sj,C,FortranStyle>& rhs) : 
	  SmallMatrixView<T,M,N,Si,Sj,C,CStyle>(rhs) {}

	inline SmallMatrixView(const MatrixView<T>& rhs) : 
	  SmallMatrixView<T,M,N,Si,Sj,C,CStyle>(rhs) {}

	inline SmallMatrixView(T* _m PARAMFIRSTLAST(T) ) :
	  SmallMatrixView<T,M,N,Si,Sj,C,CStyle>(_m FIRSTLAST1(_first,_last) ) {}

	virtual inline ~SmallMatrixView() {} 

	//
	// Op=
	//

	inline const SmallMatrixView<T,M,N,Si,Sj,C,FortranStyle>& operator=(
	    const SmallMatrixView<T,M,N,Si,Sj,C,FortranStyle>& m2) const
	{ SmallMatrixView<T,M,N,Si,Sj,C,CStyle>::operator=(m2); return *this; }

	inline const SmallMatrixView<T,M,N,Si,Sj,C,FortranStyle>& operator=(
	    const GenSmallMatrix<RealType(T),M,N,Si,Sj,C>& m2) const
	{ SmallMatrixView<T,M,N,Si,Sj,C,CStyle>::operator=(m2); return *this; }

	inline const SmallMatrixView<T,M,N,Si,Sj,C,FortranStyle>& operator=(
	    const GenSmallMatrix<ComplexType(T),M,N,Si,Sj,C>& m2) const
	{ SmallMatrixView<T,M,N,Si,Sj,C,CStyle>::operator=(m2); return *this; }

	template <class T2, int Si2, int Sj2, bool C2>
	  inline const SmallMatrixView<T,M,N,Si,Sj,C,FortranStyle>& operator=(
	      const GenSmallMatrix<T2,M,N,Si2,Sj2,C2>& m2) const
	{ SmallMatrixView<T,M,N,Si,Sj,C,CStyle>::operator=(m2); return *this; }

	inline const SmallMatrixView<T,M,N,Si,Sj,C,FortranStyle>& operator=(
	    const GenMatrix<RealType(T)>& m2) const
	{ SmallMatrixView<T,M,N,Si,Sj,C,CStyle>::operator=(m2); return *this; }

	inline const SmallMatrixView<T,M,N,Si,Sj,C,FortranStyle>& operator=(
	    const GenMatrix<ComplexType(T)>& m2) const
	{ SmallMatrixView<T,M,N,Si,Sj,C,CStyle>::operator=(m2); return *this; }

	template <class T2>
	  inline const SmallMatrixView<T,M,N,Si,Sj,C,FortranStyle>& operator=(
	      const GenMatrix<T2>& m2) const
	{ SmallMatrixView<T,M,N,Si,Sj,C,CStyle>::operator=(m2); return *this; }

	inline const SmallMatrixView<T,M,N,Si,Sj,C,FortranStyle>& operator=(
	    T x) const 
	{ SmallMatrixView<T,M,N,Si,Sj,C,CStyle>::operator=(x); return *this; }

	inline const SmallMatrixView<T,M,N,Si,Sj,C,FortranStyle>& operator=(
	    const AssignableToSmallMatrix<RealType(T),M,N>& m2) const
	{ SmallMatrixView<T,M,N,Si,Sj,C,CStyle>::operator=(m2); return *this; }

	inline const SmallMatrixView<T,M,N,Si,Sj,C,FortranStyle>& operator=(
	    const AssignableToSmallMatrix<ComplexType(T),M,N>& m2) const
	{ SmallMatrixView<T,M,N,Si,Sj,C,CStyle>::operator=(m2); return *this; }


	//
	// Access
	//

	inline RefType(T) operator()(size_t i,size_t j) const 
	{ 
	  TMVAssert(i > 0 && i <= M);
	  TMVAssert(j > 0 && j <= N);
	  return ref(i-1,j-1); 
	}

	inline SmallVectorView<T,N,Sj,C,FortranStyle> operator[](
	    size_t i) const 
	{ 
	  TMVAssert(i>0 && i<=M);
	  return SmallMatrixView<T,M,N,Si,Sj,C,CStyle>::row(i-1);
	}

	inline SmallVectorView<T,N,Sj,C,FortranStyle> row(size_t i) const
	{
	  TMVAssert(i>0 && i<=M);
	  return SmallMatrixView<T,M,N,Si,Sj,C,CStyle>::row(i-1);
	}

	inline SmallVectorView<T,M,Si,C,FortranStyle> col(size_t j) const
	{
	  TMVAssert(j>0 && j<=N);
	  return SmallMatrixView<T,M,N,Si,Sj,C,CStyle>::col(j-1);
	}

	inline SmallVectorView<T,MIN(M,N),Si+Sj,C,FortranStyle> diag() const
	{ return SmallMatrixView<T,M,N,Si,Sj,C,CStyle>::diag(); }

	inline VectorView<T,FortranStyle> diag(int i) const
	{ return RegView().diag(i); }

	inline VectorView<T,FortranStyle> row(
	    size_t i, size_t j1, size_t j2) const 
	{
	  TMVAssert(i>0 && i<=M);
	  TMVAssert(j1>0 && j1<=j2 && j2<=N);
	  return RegView().row(i,j1,j2); 
	}

	inline VectorView<T,FortranStyle> col(
	    size_t j, size_t i1, size_t i2) const
	{
	  TMVAssert(j>0 && j<=N);
	  TMVAssert(i1>0 && i1<=i2 && i2<=M);
	  return RegView().col(j,i1,i2); 
	}

	inline VectorView<T,FortranStyle> diag(
	    int i, size_t j1, size_t j2) const
	{
	  TMVAssert(i>-int(M) && i<=int(N));
	  return RegView().diag(i,j1,j2); 
	}

	//
	// Modifying Functions
	//

	inline const SmallMatrixView<T,M,N,Si,Sj,C,FortranStyle>& Zero() const 
	{ SmallMatrixView<T,M,N,Si,Sj,C,CStyle>::Zero(); return *this; }

	inline const SmallMatrixView<T,M,N,Si,Sj,C,FortranStyle>& Clip(
	    RealType(T) thresh) const
	{ SmallMatrixView<T,M,N,Si,Sj,C,CStyle>::Clip(thresh); return *this; }

	inline const SmallMatrixView<T,M,N,Si,Sj,C,FortranStyle>& SetAllTo(
	    T x) const
	{ SmallMatrixView<T,M,N,Si,Sj,C,CStyle>::SetAllTo(x); return *this; }

	inline const SmallMatrixView<T,M,N,Si,Sj,C,FortranStyle>& 
	  TransposeSelf() const
	  { 
	    SmallMatrixView<T,M,N,Si,Sj,C,CStyle>::TransposeSelf();
	    return *this; 
	  }

	inline const SmallMatrixView<T,M,N,Si,Sj,C,FortranStyle>& 
	  ConjugateSelf() const
	  { 
	    SmallMatrixView<T,M,N,Si,Sj,C,CStyle>::ConjugateSelf();
	    return *this;
	  }

	inline const SmallMatrixView<T,M,N,Si,Sj,C,FortranStyle>& SetToIdentity(
	    T x=T(1)) const
	{
	  SmallMatrixView<T,M,N,Si,Sj,C,CStyle>::SetToIdentity(x);
	  return *this;
	}

	inline const SmallMatrixView<T,M,N,Si,Sj,C,FortranStyle>& SwapRows(
	    size_t i1, size_t i2) const
	{ 
	  TMVAssert(i1 > 0 && i1 <= M);
	  TMVAssert(i2 > 0 && i2 <= M);
	  if (i1 != i2)
	    SmallMatrixView<T,M,N,Si,Sj,C,CStyle>::SwapRows(i1-1,i2-1); 
	  return *this; 
	}

	inline const SmallMatrixView<T,M,N,Si,Sj,C,FortranStyle>& SwapCols(
	    size_t j1, size_t j2) const
	{ 
	  TMVAssert(j1 > 0 && j1 <= N);
	  TMVAssert(j2 > 0 && j2 <= N);
	  if (j1 != j2)
	    SmallMatrixView<T,M,N,Si,Sj,C,CStyle>::SwapCols(j1-1,j2-1); 
	  return *this; 
	}

	inline const SmallMatrixView<T,M,N,Si,Sj,C,FortranStyle>& PermuteRows(
	    const size_t* p, size_t i1, size_t i2) const
	{
	  TMVAssert(i1>0);
	  SmallMatrixView<T,M,N,Si,Sj,C,CStyle>::PermuteRows(p,i1-1,i2);
	  return *this; 
	}

	inline const SmallMatrixView<T,M,N,Si,Sj,C,FortranStyle>& PermuteRows(
	    const size_t* p) const
	{ SmallMatrixView<T,M,N,Si,Sj,C,CStyle>::PermuteRows(p); return *this; }

	inline const SmallMatrixView<T,M,N,Si,Sj,C,FortranStyle>& PermuteCols(
	    const size_t* p, size_t j1, size_t j2) const
	{ Transpose().PermuteRows(p,j1,j2); return *this; }

	inline const SmallMatrixView<T,M,N,Si,Sj,C,FortranStyle>& PermuteCols(
	    const size_t* p) const
	{ Transpose().PermuteRows(p); return *this; }

	inline const SmallMatrixView<T,M,N,Si,Sj,C,FortranStyle>& 
	  ReversePermuteRows(const size_t* p, size_t i1, size_t i2) const
	  { 
	    TMVAssert(i1>0);
	    SmallMatrixView<T,M,N,Si,Sj,C,CStyle>::ReversePermuteRows(
		p,i1-1,i2);
	    return *this;
	  }

	inline const SmallMatrixView<T,M,N,Si,Sj,C,FortranStyle>& 
	  ReversePermuteRows(const size_t* p) const
	  {
	    SmallMatrixView<T,M,N,Si,Sj,C,CStyle>::ReversePermuteRows(p);
	    return *this;
	  }

	inline const SmallMatrixView<T,M,N,Si,Sj,C,FortranStyle>& 
	  ReversePermuteCols(
	      const size_t* p, size_t j1, size_t j2) const
	  { Transpose().ReversePermuteRows(p,j1,j2); return *this; }

	inline const SmallMatrixView<T,M,N,Si,Sj,C,FortranStyle>& 
	  ReversePermuteCols(
	      const size_t* p) const
	  { Transpose().ReversePermuteRows(p); return *this; }


	//
	// SubMatrix
	//

	inline MatrixView<T,FortranStyle> SubMatrix(
	    int i1, int i2, int j1, int j2) const
	{ return RegView().SubMatrix(i1,i2,j1,j2); }

	inline MatrixView<T,FortranStyle> SubMatrix(
	    int i1, int i2, int j1, int j2, int istep, int jstep) const
	{ return RegView().SubMatrix(i1,i2,j1,j2,istep,jstep); }

	inline VectorView<T,FortranStyle> SubVector(
	    size_t i, size_t j, int istep, int jstep, size_t s) const
	{ return RegView().SubVector(i,j,istep,jstep,s); }

	inline MatrixView<T,FortranStyle> ColPair(size_t j1, size_t j2) const
	{ return RegView().ColPair(j1,j2); }

	inline MatrixView<T,FortranStyle> RowPair(size_t i1, size_t i2) const
	{ return RegView().RowPair(i1,i2); }

	inline MatrixView<T,FortranStyle> Cols(size_t j1, size_t j2) const
	{ return RegView().Cols(j1,j2); }

	inline MatrixView<T,FortranStyle> Rows(size_t i1, size_t i2) const
	{ return RegView().Rows(i1,i2); }

	inline MatrixView<RealType(T),FortranStyle> Real() const
	{ return RegView().Real(); }

	inline MatrixView<RealType(T),FortranStyle> Imag() const
	{ return RegView().Imag(); }


	//
	// Views
	//

	inline SmallMatrixView<T,M,N,Si,Sj,C,FortranStyle> View() const
	{ return SmallMatrixView<T,M,N,Si,Sj,C,CStyle>::View(); }

	inline MatrixView<T,FortranStyle> RegView() const
	{ 
	  return MatrixView<T,FortranStyle>(this->ptr(),M,N,Si,Sj,this->stor(),
	      this->ct(),this->CanLinearize()?M*N:0 FIRSTLAST );
	}

	inline SmallMatrixView<T,N,M,Sj,Si,C,FortranStyle> Transpose() const
	{ return SmallMatrixView<T,M,N,Si,Sj,C,CStyle>::Transpose(); }

	inline SmallMatrixView<T,M,N,Si,Sj,!C,FortranStyle> Conjugate() const
	{ return SmallMatrixView<T,M,N,Si,Sj,C,CStyle>::Conjugate(); }

	inline SmallMatrixView<T,N,M,Sj,Si,!C,FortranStyle> Adjoint() const
	{ return SmallMatrixView<T,M,N,Si,Sj,C,CStyle>::Adjoint(); }

	inline SmallVectorView<T,M*N,1,C,FortranStyle> LinearView() const
	{ return SmallMatrixView<T,M,N,Si,Sj,C,CStyle>::LinearView(); }

      protected:

	using SmallMatrixView<T,M,N,Si,Sj,C,CStyle>::ref;

    }; // FortranStyle MatrixView

#define Si (S==RowMajor?int(N):int(1))
#define Sj (S==RowMajor?int(1):int(M))
  template <class T, size_t M, size_t N, StorageType S, IndexStyle I> 
    class SmallMatrix : 
      public GenSmallMatrix<T,M,N,Si,Sj,false>
    {

      public:

	//
	// Constructors
	//

#define NEW_SIZE \
	itsm(new T[M*N]) DEFFIRSTLAST(itsm.get(),itsm.get()+M*N)

	inline SmallMatrix() : NEW_SIZE
	{
	  TMVAssert(S==RowMajor || S==ColMajor);
#ifdef TMVDEBUG
	  SetAllTo(T(888));
#endif
	}

	explicit inline SmallMatrix(T x) : NEW_SIZE
	{ 
	  TMVAssert(S==RowMajor || S==ColMajor);
	  SetAllTo(x);
	}

	inline SmallMatrix(const T* vv) : NEW_SIZE
	{ 
	  TMVAssert(S==RowMajor || S==ColMajor);
	  memmove(itsm.get(),vv,M*N*sizeof(T));
	}

	inline SmallMatrix(const std::vector<T>& vv) : NEW_SIZE
	{ 
	  TMVAssert(S==RowMajor || S==ColMajor);
	  TMVAssert(vv.size() == M*N);
          T* vi = itsm.get();
	  typename std::vector<T>::const_iterator vvi = vv.begin();
	  for(size_t i=M*N;i>0;--i) *vi = *vvi;
	}

	explicit inline SmallMatrix(
	    const std::vector<std::vector<T> >& vv) : NEW_SIZE
	{ 
	  TMVAssert(S==RowMajor || S==ColMajor);
	  if (S == RowMajor) {
	    T* vi = itsm.get();
	    for(size_t i=0;i<M;++i,++vi) {
	      TMVAssert(vv[i].size() == N);
	      typename std::vector<T>::const_iterator vvi = vv[i].begin();
	      for(size_t j=0;j<N;++j,++vi,++vvi) {
		*vi = *vvi;
	      }
	    }
	  } else {
	    for(size_t i=0;i<M;++i) {
	      T* vi = itsm.get()+i;
	      typename std::vector<T>::const_iterator vvi = vv[i].begin();
	      for(size_t j=0;j<N;++j,vi+=Sj,++vvi) {
		*vi = *vvi;
	      }
	    }
	  }
	}

	inline SmallMatrix(const SmallMatrix<T,M,N,S,I>& rhs) : NEW_SIZE
	{ 
	  TMVAssert(S==RowMajor || S==ColMajor);
	  memmove(itsm.get(),rhs.itsm.get(),M*N*sizeof(T));
	}

	template <IndexStyle I2> inline SmallMatrix(
	    const SmallMatrix<T,M,N,S,I2>& rhs) : NEW_SIZE
	{ 
	  TMVAssert(S==RowMajor || S==ColMajor);
	  memmove(itsm.get(),rhs.itsm.get(),M*N*sizeof(T));
	}

	template <int Si2, int Sj2, bool C2> inline SmallMatrix(
	    const GenSmallMatrix<RealType(T),M,N,Si2,Sj2,C2>& m2) : NEW_SIZE
	{ 
	  TMVAssert(S==RowMajor || S==ColMajor);
	  m2.AssignTom(View());
	}

	template <int Si2, int Sj2, bool C2> inline SmallMatrix(
	    const GenSmallMatrix<ComplexType(T),M,N,Si2,Sj2,C2>& m2) : NEW_SIZE
	{ 
	  TMVAssert(S==RowMajor || S==ColMajor);
	  TMVAssert(IsComplex(T()));
	  m2.AssignTom(View());
	}

	template <IndexStyle I2> inline SmallMatrix(
	    const Matrix<T,S,I2>& rhs) : NEW_SIZE
	{ 
	  TMVAssert(S==RowMajor || S==ColMajor);
	  TMVAssert(rhs.colsize() == M);
	  TMVAssert(rhs.rowsize() == N);
	  memmove(itsm.get(),rhs.itsm.get(),M*N*sizeof(T));
	}

	template <class T2, int Si2, int Sj2, bool C2> 
	  inline SmallMatrix(const GenSmallMatrix<T2,M,N,Si2,Sj2,C2>& rhs) :
	  NEW_SIZE
	{ 
	  TMVAssert(IsComplex(T()) || IsReal(T2()));
	  TMVAssert(S==RowMajor || S==ColMajor);
	  Copy(rhs,View()); 
	}

	inline SmallMatrix(const GenMatrix<RealType(T)>& m2) : NEW_SIZE
	{ 
	  TMVAssert(S==RowMajor || S==ColMajor);
	  TMVAssert(m2.colsize() == M);
	  TMVAssert(m2.rowsize() == N);
	  m2.AssignToM(RegView());
	}

	inline SmallMatrix(const GenMatrix<ComplexType(T)>& m2) : NEW_SIZE
	{ 
	  TMVAssert(S==RowMajor || S==ColMajor);
	  TMVAssert(IsComplex(T()));
	  TMVAssert(m2.colsize() == M);
	  TMVAssert(m2.rowsize() == N);
	  m2.AssignToM(RegView());
	}

	template <class T2> inline SmallMatrix(const GenMatrix<T2>& m2) :
	  NEW_SIZE
	{ 
	  TMVAssert(IsComplex(T()) || IsReal(T2()));
	  TMVAssert(S==RowMajor || S==ColMajor);
	  TMVAssert(m2.colsize() == M);
	  TMVAssert(m2.rowsize() == N);
	  Copy(m2,RegView()); 
	}

	inline SmallMatrix(
	    const AssignableToSmallMatrix<RealType(T),M,N>& m2) : NEW_SIZE
	{ 
	  TMVAssert(S==RowMajor || S==ColMajor);
	  m2.AssignTom(View());
	}

	inline SmallMatrix(
	    const AssignableToSmallMatrix<ComplexType(T),M,N>& m2) : NEW_SIZE
	{ 
	  TMVAssert(S==RowMajor || S==ColMajor);
	  TMVAssert(IsComplex(T()));
	  m2.AssignTom(View());
	}

#undef NEW_SIZE

	virtual inline ~SmallMatrix() {}

	//
	// Op=
	//

	inline SmallMatrix<T,M,N,S,I>& operator=(
	    const SmallMatrix<T,M,N,S,I>& m2)
	{ 
	  if (&m2 != this) memmove(itsm.get(),m2.itsm.get(),M*N*sizeof(T));
	  return *this; 
	}

	template <class T2, int Si2, int Sj2, bool C2> 
	  inline SmallMatrix<T,M,N,S,I>& operator=(
	      const GenSmallMatrix<T2,M,N,Si2,Sj2,C2>& m2)
	  { 
	    TMVAssert(IsComplex(T()) || IsReal(T2()));
	    Copy(m2,View()); 
	    return *this; 
	  }

	inline SmallMatrix<T,M,N,S,I>& operator=(T x) 
	{ return SetToIdentity(x); }

	inline SmallMatrix<T,M,N,S,I>& operator=(
	    const GenMatrix<RealType(T)>& m2)
	{
	  TMVAssert(m2.colsize() == M);
	  TMVAssert(m2.rowsize() == N);
	  m2.AssignToM(RegView());
	  return *this;
	}

	inline SmallMatrix<T,M,N,S,I>& operator=(
	    const GenMatrix<ComplexType(T)>& m2)
	{
	  TMVAssert(m2.colsize() == M);
	  TMVAssert(m2.rowsize() == N);
	  TMVAssert(IsComplex(T()));
	  m2.AssignToM(RegView());
	  return *this;
	}

	template <class T2> inline SmallMatrix<T,M,N,S,I>& operator=(
	    const GenMatrix<T2>& m2)
	{
	  TMVAssert(m2.colsize() == M);
	  TMVAssert(m2.rowsize() == N);
	  Copy(m2,RegView());
	  return *this;
	}

	inline SmallMatrix<T,M,N,S,I>& operator=(
	    const AssignableToSmallMatrix<RealType(T),M,N>& m2)
	{ m2.AssignTom(View()); return *this; }

	inline SmallMatrix<T,M,N,S,I>& operator=(
	    const AssignableToSmallMatrix<ComplexType(T),M,N>& m2)
	{
	  TMVAssert(IsComplex(T()));
	  m2.AssignTom(View());
	  return *this;
	}


	//
	// Access
	//

	inline T operator()(size_t i,size_t j) const
	{ return View()(i,j); }

	inline T& operator()(size_t i,size_t j) 
	{ 
	  if (I==CStyle) {
	    TMVAssert(i<M && j<N);
	    return ref(i,j);
	  } else {
	    TMVAssert(i>0 && i<=M && j>0 && j<=N);
	    return ref(i-1,j-1);
	  }
	}

	inline ConstSmallVectorView<T,N,Sj,false,I> row(size_t i) const 
	{ return View().row(i); }

	inline ConstVectorView<T,I> row(size_t i, size_t j1, size_t j2) const 
	{ return View().row(i,j1,j2); }

	inline ConstSmallVectorView<T,N,Sj,false,I> operator[](size_t i) const
	{ return View().row(i); }

	inline ConstSmallVectorView<T,M,Si,false,I> col(size_t j) const
	{ return View().col(j); }

	inline ConstVectorView<T,I> col(size_t j, size_t i1, size_t i2) const
	{ return View().col(j,i1,i2); }

	inline ConstSmallVectorView<T,MIN(M,N),Si+Sj,false,I> diag() const
	{ return View().diag(); }

	inline ConstVectorView<T,I> diag(int i) const
	{ return View().diag(i); }

	inline ConstVectorView<T,I> diag(int i, size_t j1, size_t j2) const
	{ return View().diag(i,j1,j2); }

	inline SmallVectorView<T,N,Sj,false,I> row(size_t i)
	{ return View().row(i); }

	inline VectorView<T,I> row(size_t i, size_t j1, size_t j2)
	{ return View().row(i,j1,j2); }

	inline SmallVectorView<T,N,Sj,false,I> operator[](size_t i)
	{ return View().row(i); }

	inline SmallVectorView<T,M,Si,false,I> col(size_t j)
	{ return View().col(j); }

	inline VectorView<T,I> col(size_t j, size_t i1, size_t i2)
	{ return View().col(j,i1,i2); }

	inline SmallVectorView<T,MIN(M,N),Si+Sj,false,I> diag()
	{ return View().diag(); }

	inline VectorView<T,I> diag(int i)
	{ return View().diag(i); }

	inline VectorView<T,I> diag(int i, size_t j1, size_t j2) 
	{ return View().diag(i,j1,j2); }

	//
	// Modifying Functions
	//

	inline SmallMatrix<T,M,N,S,I>& Zero() 
	{ LinearView().Zero(); return *this; }

	inline SmallMatrix<T,M,N,S,I>& Clip(RealType(T) thresh)
	{ LinearView().Clip(thresh); return *this; }

	inline SmallMatrix<T,M,N,S,I>& SetAllTo(T x) 
	{ LinearView().SetAllTo(x); return *this; }

	inline SmallMatrix<T,M,N,S,I>& TransposeSelf() 
	{ LinearView().TransposeSelf(); return *this; }

	inline SmallMatrix<T,M,N,S,I>& ConjugateSelf() 
	{ LinearView().ConjugateSelf(); return *this; }

	inline SmallMatrix<T,M,N,S,I>& SetToIdentity(T x=T(1)) 
	{ View().SetToIdentity(x); return *this; }

	inline SmallMatrix<T,M,N,S,I>& SwapRows(size_t i1, size_t i2)
	{ View().SwapRows(i1,i2); return *this; }

	inline SmallMatrix<T,M,N,S,I>& SwapCols(size_t j1, size_t j2)
	{ View().SwapCols(j1,j2); return *this; }

	inline SmallMatrix<T,M,N,S,I>& PermuteRows(
	    const size_t* p, size_t i1, size_t i2)
	{ View().PermuteRows(p,i1,i2); return *this; }

	inline SmallMatrix<T,M,N,S,I>& PermuteRows(const size_t* p)
	{ View().PermuteRows(p); return *this; }

	inline SmallMatrix<T,M,N,S,I>& PermuteCols(
	    const size_t* p, size_t j1, size_t j2)
	{ View().PermuteCols(p,j1,j2); return *this; }

	inline SmallMatrix<T,M,N,S,I>& PermuteCols(const size_t* p)
	{ View().PermuteCols(p); return *this; }

	inline SmallMatrix<T,M,N,S,I>& ReversePermuteRows(
	    const size_t* p, size_t i1, size_t i2)
	{ View().ReversePermuteRows(p,i1,i2); return *this; }

	inline SmallMatrix<T,M,N,S,I>& ReversePermuteRows(const size_t* p)
	{ View().ReversePermuteRows(p); return *this; }

	inline SmallMatrix<T,M,N,S,I>& ReversePermuteCols(
	    const size_t* p, size_t j1, size_t j2)
	{ View().ReversePermuteCols(p,j1,j2); return *this; }

	inline SmallMatrix<T,M,N,S,I>& ReversePermuteCols(const size_t* p)
	{ View().ReversePermuteCols(p); return *this; }

	//
	// SubMatrix
	//

	inline ConstMatrixView<T,I> SubMatrix(
	    int i1, int i2, int j1, int j2) const
	{ return View().SubMatrix(i1,i2,j1,j2); }

	inline ConstMatrixView<T,I> SubMatrix(
	    int i1, int i2, int j1, int j2, int istep, int jstep) const
	{ return View().SubMatrix(i1,i2,j1,j2,istep,jstep); }

	inline ConstVectorView<T,I> SubVector(
	    size_t i, size_t j, int istep, int jstep, size_t s) const
	{ return View().SubVector(i,j,istep,jstep,s); }

	inline ConstMatrixView<T,I> ColPair(size_t j1, size_t j2) const
	{ return View().ColPair(j1,j2); }

	inline ConstMatrixView<T,I> RowPair(size_t i1, size_t i2) const
	{ return View().RowPair(i1,i2); }

	inline ConstMatrixView<T,I> Cols(size_t j1, size_t j2) const
	{ return View().Cols(j1,j2); }

	inline ConstMatrixView<T,I> Rows(size_t i1, size_t i2) const
	{ return View().Rows(i1,i2); }

	inline ConstMatrixView<RealType(T)> Real() const
	{ return View().Real(); }

	inline ConstMatrixView<RealType(T)> Imag() const
	{ return View().Imag(); }

	inline MatrixView<T,I> SubMatrix(int i1, int i2, int j1, int j2)
	{ return View().SubMatrix(i1,i2,j1,j2); }

	inline MatrixView<T,I> SubMatrix(
	    int i1, int i2, int j1, int j2, int istep, int jstep) 
	{ return View().SubMatrix(i1,i2,j1,j2,istep,jstep); }

	inline VectorView<T,I> SubVector(
	    size_t i, size_t j, int istep, int jstep, size_t s) 
	{ return View().SubVector(i,j,istep,jstep,s); }

	inline MatrixView<T,I> ColPair(size_t j1, size_t j2) 
	{ return View().ColPair(j1,j2); }

	inline MatrixView<T,I> RowPair(size_t i1, size_t i2) 
	{ return View().RowPair(i1,i2); }

	inline MatrixView<T,I> Cols(size_t j1, size_t j2) 
	{ return View().Cols(j1,j2); }

	inline MatrixView<T,I> Rows(size_t i1, size_t i2) 
	{ return View().Rows(i1,i2); }

	inline MatrixView<RealType(T)> Real() 
	{ return View().Real(); }

	inline MatrixView<RealType(T)> Imag() 
	{ return View().Imag(); }


	//
	// Views
	//

	inline ConstSmallMatrixView<T,M,N,Si,Sj,false,I> View() const
	{ return ConstSmallMatrixView<T,M,N,Si,Sj,false,I>(cptr()); }

	inline ConstMatrixView<T,I> RegView() const
	{ 
	  return ConstMatrixView<T,I>(cptr(),M,N,Si,Sj,this->stor(),
	    this->ct(),M*N); 
	}

	inline ConstSmallMatrixView<T,N,M,Sj,Si,false,I> Transpose() const
	{ return ConstSmallMatrixView<T,N,M,Sj,Si,false,I>(cptr()); }

	inline ConstSmallMatrixView<T,M,N,Si,Sj,true,I> Conjugate() const
	{ return ConstSmallMatrixView<T,M,N,Si,Sj,true,I>(cptr()); }

	inline ConstSmallMatrixView<T,N,M,Sj,Si,true,I> Adjoint() const
	{ return ConstSmallMatrixView<T,N,M,Sj,Si,true,I>(cptr()); }

	inline ConstSmallVectorView<T,M*N,1,false,I> ConstLinearView() const
	{ return ConstSmallVectorView<T,M*N,1,false,I>(cptr()); }

	inline SmallMatrixView<T,M,N,Si,Sj,false,I> View()
	{ return SmallMatrixView<T,M,N,Si,Sj,false,I>(ptr() FIRSTLAST); }

	inline MatrixView<T,I> RegView()
	{ 
	  return MatrixView<T,I>(ptr(),M,N,Si,Sj,this->stor(),
	    this->ct(),M*N FIRSTLAST); 
	}

	inline SmallMatrixView<T,N,M,Sj,Si,false,I> Transpose()
	{ return SmallMatrixView<T,N,M,Sj,Si,false,I>(ptr() FIRSTLAST); }

	inline SmallMatrixView<T,M,N,Si,Sj,true,I> Conjugate()
	{ return SmallMatrixView<T,M,N,Si,Sj,true,I>(ptr() FIRSTLAST); }

	inline SmallMatrixView<T,N,M,Sj,Si,true,I> Adjoint()
	{ return SmallMatrixView<T,N,M,Sj,Si,true,I>(ptr() FIRSTLAST); }

	inline SmallVectorView<T,M*N,1,false,I> LinearView()
	{ return SmallVectorView<T,M*N,1,false,I>(ptr() FIRSTLAST); }

	inline const T* cptr() const { return itsm.get(); }
	inline T* ptr() { return itsm.get(); }

      protected :

	auto_array<T> itsm;

	inline T& ref(size_t i, size_t j)
	{ 
	  TMVAssert(i<M);
	  TMVAssert(j<N);
	  T*const mi = ptr() + int(i)*Si + int(j)*Sj;
#ifdef TMVFLDEBUG
	  TMVAssert(mi >= first);
	  TMVAssert(mi < last);
#endif
	  return *mi;
	}

#ifdef TMVFLDEBUG
      public:
	const T*const first;
	const T*const last;
#endif

    }; // SmallMatrix

#undef Si
#undef Sj

//---------------------------------------------------------------------------

  //
  // Special Creators: 
  //   RowVectorViewOf(v) = 1xn Matrix with v in only row - Same Storage
  //   ColVectorViewOf(v) = nx1 Matrix with v in only col - Same Storage
  //   MatrixViewOf(m,colsize,rowsize,storage) = MatrixView of m 
  //

  template <class T, size_t N, int S, bool C> 
    inline ConstSmallMatrixView<T,1,N,int(N)*S,S,C> RowVectorViewOf(
	const GenSmallVector<T,N,S,C>& v)
    { return ConstSmallMatrixView<T,1,N,int(N)*S,S,C>(v.cptr()); }

  template <class T, size_t N, int S, bool C, IndexStyle I> 
    inline ConstSmallMatrixView<T,1,N,int(N)*S,S,C,I> RowVectorViewOf(
	const ConstSmallVectorView<T,N,S,C,I>& v)
    { return ConstSmallMatrixView<T,1,N,int(N)*S,S,C,I>(v.cptr()); }

  template <class T, size_t N, IndexStyle I> 
    inline ConstSmallMatrixView<T,1,N,N,1,false,I> RowVectorViewOf(
	const SmallVector<T,N,I>& v)
    { return ConstSmallMatrixView<T,1,N,N,1,false,I>(v.cptr()); }

  template <class T, size_t N, int S, bool C, IndexStyle I> 
    inline SmallMatrixView<T,1,N,int(N)*S,S,false,I> RowVectorViewOf(
      const SmallVectorView<T,N,S,C,I>& v)
    {
      return SmallMatrixView<T,1,N,int(N)*S,S,false,I>(v.ptr()
	  FIRSTLAST1(v.ptr(),v.ptr()+v.size()) );
    }

  template <class T, size_t N, IndexStyle I> 
    inline SmallMatrixView<T,1,N,N,1,false,I> RowVectorViewOf(
      SmallVector<T,N,I>& v)
    {
      return SmallMatrixView<T,1,N,N,1,false,I>(v.ptr()
	  FIRSTLAST1(v.ptr(),v.ptr()+v.size()) );
    }

  template <class T, size_t M, int S, bool C> 
    inline ConstSmallMatrixView<T,M,1,S,int(M)*S,C> ColVectorViewOf(
	const GenSmallVector<T,M,S,C>& v)
    { return ConstSmallMatrixView<T,M,1,S,int(M)*S,C>(v.cptr()); }

  template <class T, size_t M, int S, bool C, IndexStyle I> 
    inline ConstSmallMatrixView<T,M,1,S,int(M)*S,C,I> ColVectorViewOf(
	const ConstSmallVectorView<T,M,S,C,I>& v)
    { return ConstSmallMatrixView<T,M,1,S,int(M)*S,C,I>(v.cptr()); }

  template <class T, size_t M, IndexStyle I> 
    inline ConstSmallMatrixView<T,M,1,1,M,false,I> ColVectorViewOf(
	const SmallVector<T,M,I>& v)
    { return ConstSmallMatrixView<T,M,1,1,M,false,I>(v.cptr()); }

  template <class T, size_t M, int S, bool C, IndexStyle I> 
    inline SmallMatrixView<T,M,1,S,int(M)*S,false,I> ColVectorViewOf(
      const SmallVectorView<T,M,S,C,I>& v)
    {
      return SmallMatrixView<T,M,1,S,int(M)*S,false,I>(v.ptr()
	  FIRSTLAST1(v.ptr(),v.ptr()+v.size()) );
    }

  template <class T, size_t M, IndexStyle I> 
    inline SmallMatrixView<T,M,1,1,M,false,I> ColVectorViewOf(
      SmallVector<T,M,I>& v)
    {
      return SmallMatrixView<T,M,1,1,M,false,I>(v.ptr()
	  FIRSTLAST1(v.ptr(),v.ptr()+v.size()) );
    }

#define Si (S==RowMajor?int(N):int(1))
#define Sj (S==RowMajor?int(1):int(M))
  template <size_t M, size_t N, StorageType S, IndexStyle I, class T> 
    inline SmallMatrixView<T,M,N,Si,Sj,false,I> SmallMatrixViewOf(T* m)
    {
      TMVAssert(S == RowMajor || S == ColMajor);
      return SmallMatrixView<T,M,N,Si,Sj,false,I>(m FIRSTLAST1(m,m+M*N));
    }
  template <size_t M, size_t N, StorageType S, IndexStyle I, class T> 
    inline ConstSmallMatrixView<T,M,N,Si,Sj,false,I> SmallMatrixViewOf(
	const T* m)
    {
      TMVAssert(S == RowMajor || S == ColMajor);
      return ConstSmallMatrixView<T,M,N,Si,Sj,false,I>(m);
    }
  template <size_t M, size_t N, StorageType S, class T> 
    inline SmallMatrixView<T,M,N,Si,Sj,false> SmallMatrixViewOf(T* m)
    {
      TMVAssert(S == RowMajor || S == ColMajor);
      return SmallMatrixView<T,M,N,Si,Sj,false>(m FIRSTLAST1(m,m+M*N));
    }
  template <size_t M, size_t N, StorageType S, class T> 
    inline ConstSmallMatrixView<T,M,N,Si,Sj,false> SmallMatrixViewOf(
	const T* m)
    {
      TMVAssert(S == RowMajor || S == ColMajor);
      return ConstSmallMatrixView<T,M,N,Si,Sj,false>(m);
    }
#undef Si
#undef Sj


  //
  // Copy Matrices
  //

  template <class T1, class T2, size_t M, size_t N, int Si1, int Sj1, int Si2, int Sj2, bool C1> 
    inline void Copy(
	const GenSmallMatrix<T1,M,N,Si1,Sj1,C1>& m1, 
	const SmallMatrixView<T2,M,N,Si2,Sj2,false>& m2)
    {
      TMVAssert(IsComplex(T2()) || IsReal(T1()));

      if (Si1==Si2 && Sj1==Sj2 && m1.CanLinearize())
	Copy(m1.ConstLinearView(),m2.LinearView());
      else if (Si2==1 || Si1==1) 
	for(size_t j=0;j<N;++j) Copy(m1.col(j),m2.col(j));
      else
	for(size_t i=0;i<M;++i) Copy(m1.row(i),m2.row(i));
    }
  template <class T1, class T2, size_t M, size_t N, int Si1, int Sj1, bool C1> 
    inline void Copy(
	const GenSmallMatrix<T1,M,N,Si1,Sj1,C1>& m1, 
	T2* m2p, const int Si2, const int Sj2)
    { 
      if (Si1 < Sj1) 
	if (Si2 == 1)
	  for(size_t j=0;j<N;j++,m2p+=Sj2) 
	    Copy(m1.col(j),SmallVectorViewOf<M>(m2p));
	else
	  for(size_t j=0;j<N;j++,m2p+=Sj2) 
	    Copy(m1.col(j),m2p,Si2);
      else
	if (Sj2 == 1)
	  for(size_t i=0;i<M;i++,m2p+=Si2) 
	    Copy(m1.row(i),SmallVectorViewOf<N>(m2p));
	else
	  for(size_t i=0;i<M;i++,m2p+=Si2) 
	    Copy(m1.row(i),m2p,Sj2);
    }

  //
  // Swap Matrices
  //

  template <class T, size_t M, size_t N, int Si1, int Sj1, int Si2, int Sj2, bool C1, bool C2> 
    inline void Swap(
	const SmallMatrixView<T,M,N,Si1,Sj1,C1>& m1, 
	const SmallMatrixView<T,M,N,Si2,Sj2,C2>& m2) 
    {
      if (Si1==Si2 && Sj1==Sj2 && m1.CanLinearize())
	Swap(m1.LinearView(),m2.LinearView());
      else if (Si2==1 || Si1==1) 
	for(size_t j=0;j<N;++j) Swap(m1.col(j),m2.col(j));
      else
	for(size_t i=0;i<M;++i) Swap(m1.row(i),m2.row(i));
    }

  template <class T, size_t M, size_t N, int Si1, int Sj1, bool C1, StorageType S2, IndexStyle I2> 
    inline void Swap(
	const SmallMatrixView<T,M,N,Si1,Sj1,C1>& m1, 
	SmallMatrix<T,M,N,S2,I2>& m2)
    { Swap(m1,m2.View()); }

  template <class T, size_t M, size_t N, StorageType S1, IndexStyle I1, int Si2, int Sj2, bool C2> 
    inline void Swap(
	SmallMatrix<T,M,N,S1,I1>& m1,
	const SmallMatrixView<T,M,N,Si2,Sj2,C2>& m2)
    { Swap(m1.View(),m2); }

  template <class T, size_t M, size_t N, StorageType S1, IndexStyle I1, StorageType S2, IndexStyle I2> 
    inline void Swap(SmallMatrix<T,M,N,S1,I1>& m1, SmallMatrix<T,M,N,S2,I2>& m2)
    { Swap(m1.View(),m2.View()); }

  //
  // Functions of Matrices:
  //

  template <class T, size_t M, size_t N, int Si, int Sj, bool C> 
    inline T Det(const GenSmallMatrix<T,M,N,Si,Sj,C>& m)
    { return m.Det(); }

  template <class T, size_t M, size_t N, int Si, int Sj, bool C> 
    inline T Trace(const GenSmallMatrix<T,M,N,Si,Sj,C>& m)
    { return m.Trace(); }

  template <class T, size_t M, size_t N, int Si, int Sj, bool C> 
    inline RealType(T) Norm(const GenSmallMatrix<T,M,N,Si,Sj,C>& m)
    { return m.Norm(); }

  template <class T, size_t M, size_t N, int Si, int Sj, bool C> 
    inline RealType(T) NormSq(const GenSmallMatrix<T,M,N,Si,Sj,C>& m)
    { return m.NormSq(); }

  template <class T, size_t M, size_t N, int Si, int Sj, bool C> 
    inline RealType(T) NormF(const GenSmallMatrix<T,M,N,Si,Sj,C>& m)
    { return m.NormF(); }

  template <class T, size_t M, size_t N, int Si, int Sj, bool C> 
    inline RealType(T) Norm1(const GenSmallMatrix<T,M,N,Si,Sj,C>& m)
    { return m.Norm1(); }

  template <class T, size_t M, size_t N, int Si, int Sj, bool C> 
    inline RealType(T) Norm2(const GenSmallMatrix<T,M,N,Si,Sj,C>& m)
    { return m.Norm2(); }

  template <class T, size_t M, size_t N, int Si, int Sj, bool C> 
    inline RealType(T) NormInf(const GenSmallMatrix<T,M,N,Si,Sj,C>& m)
    { return m.NormInf(); }

  template <class T, size_t M, size_t N, int Si, int Sj, bool C> 
    inline RealType(T) MaxAbsElement(const GenSmallMatrix<T,M,N,Si,Sj,C>& m)
    { return m.MaxAbsElement(); }

  template <class T, size_t M, size_t N, int Si, int Sj, bool C> 
    inline ConstSmallMatrixView<T,N,M,Sj,Si,C> Transpose(
	const GenSmallMatrix<T,M,N,Si,Sj,C>& m)
    { return m.Transpose(); }

  template <class T, size_t M, size_t N, int Si, int Sj, bool C, IndexStyle I> 
    inline ConstSmallMatrixView<T,N,M,Sj,Si,C,I> Transpose(
	const ConstSmallMatrixView<T,M,N,Si,Sj,C,I>& m)
    { return m.Transpose(); }

  template <class T, size_t M, size_t N, int Si, int Sj, bool C, IndexStyle I> 
    inline SmallMatrixView<T,N,M,Sj,Si,C,I> Transpose(
	const SmallMatrixView<T,M,N,Si,Sj,C,I>& m)
    { return m.Transpose(); }

  template <class T, size_t M, size_t N, int Si, int Sj, bool C> 
    inline ConstSmallMatrixView<T,M,N,Si,Sj,!C> Conjugate(
	const GenSmallMatrix<T,M,N,Si,Sj,C>& m)
    { return m.Conjugate(); }

  template <class T, size_t M, size_t N, int Si, int Sj, bool C, IndexStyle I> 
    inline ConstSmallMatrixView<T,M,N,Si,Sj,!C,I> Conjugate(
	const ConstSmallMatrixView<T,M,N,Si,Sj,C,I>& m)
    { return m.Conjugate(); }

  template <class T, size_t M, size_t N, int Si, int Sj, bool C, IndexStyle I> 
    inline SmallMatrixView<T,M,N,Si,Sj,!C,I> Conjugate(
	const SmallMatrixView<T,M,N,Si,Sj,C,I>& m)
    { return m.Conjugate(); }

  template <class T, size_t M, size_t N, int Si, int Sj, bool C> 
    inline ConstSmallMatrixView<T,N,M,Sj,Si,!C> Adjoint(
	const GenSmallMatrix<T,M,N,Si,Sj,C>& m)
    { return m.Adjoint(); }

  template <class T, size_t M, size_t N, int Si, int Sj, bool C, IndexStyle I> 
    inline ConstSmallMatrixView<T,N,M,Sj,Si,!C,I> Adjoint(
	const ConstSmallMatrixView<T,M,N,Si,Sj,C,I>& m)
    { return m.Adjoint(); }

  template <class T, size_t M, size_t N, int Si, int Sj, bool C, IndexStyle I> 
    inline SmallMatrixView<T,N,M,Sj,Si,!C,I> Adjoint(
	const SmallMatrixView<T,M,N,Si,Sj,C,I>& m)
    { return m.Adjoint(); }

#define Si (S==RowMajor?int(N):int(1))
#define Sj (S==RowMajor?int(1):int(M))
  template <class T, size_t M, size_t N, StorageType S, IndexStyle I> 
    inline ConstSmallMatrixView<T,N,M,Sj,Si,false> Transpose(
	const SmallMatrix<T,M,N,S,I>& m)
    { return m.Transpose(); }

  template <class T, size_t M, size_t N, StorageType S, IndexStyle I> 
    inline SmallMatrixView<T,N,M,Sj,Si,false,I> Transpose(
	SmallMatrix<T,M,N,S,I>& m)
    { return m.Transpose(); }

  template <class T, size_t M, size_t N, StorageType S, IndexStyle I> 
    inline ConstSmallMatrixView<T,M,N,Si,Sj,true,I> Conjugate(
	const SmallMatrix<T,M,N,S,I>& m)
    { return m.Conjugate(); }

  template <class T, size_t M, size_t N, StorageType S, IndexStyle I> 
    inline SmallMatrixView<T,M,N,Si,Sj,true,I> Conjugate(
	SmallMatrix<T,M,N,S,I>& m)
    { return m.Conjugate(); }

  template <class T, size_t M, size_t N, StorageType S, IndexStyle I> 
    inline ConstSmallMatrixView<T,N,M,Sj,Si,true,I> Adjoint(
	const SmallMatrix<T,M,N,S,I>& m)
    { return m.Adjoint(); }

  template <class T, size_t M, size_t N, StorageType S, IndexStyle I> 
    inline SmallMatrixView<T,N,M,Sj,Si,true,I> Adjoint(
	SmallMatrix<T,M,N,S,I>& m)
    { return m.Adjoint(); }
#undef Si
#undef Sj

  //
  // Matrix ==, != Matrix
  //

  template <class T1, class T2, size_t M, size_t N, int Si1, int Sj1, int Si2, int Sj2, bool C1, bool C2> 
    inline bool operator==(
	const GenSmallMatrix<T1,M,N,Si1,Sj1,C1>& m1,
	const GenSmallMatrix<T2,M,N,Si2,Sj2,C2>& m2)
    { 
      if (Si1==Si2 && Sj1==Sj2 && m1.CanLinearize())
	return m1.ConstLinearView() == m2.ConstLinearView();
      else if (Si2==1 || Si1==1) {
	for(size_t j=0;j<N;++j) if (m1.col(j) != m2.col(j)) return false;
	return true;
      }
      else {
	for(size_t i=0;i<M;++i) if (m1.row(i) != m2.row(i)) return false;
	return true;
      }
    }

  template <class T1, class T2, size_t M, size_t N, int Si1, int Sj1, int Si2, int Sj2, bool C1, bool C2> 
    inline bool operator!=(
	const GenSmallMatrix<T1,M,N,Si1,Sj1,C1>& m1,
	const GenSmallMatrix<T2,M,N,Si2,Sj2,C2>& m2)
    { return !(m1 == m2); }


  //
  // I/O
  //
 
  template <class T, size_t M, size_t N, int Si, int Sj, bool C> 
    inline std::istream& operator>>(
	std::istream& is, const SmallMatrixView<T,M,N,Si,Sj,C>& m)
    { is >> m.RegView(); return is; }

  template <class T, size_t M, size_t N, StorageType S, IndexStyle I> 
    inline std::istream& operator>>(
	std::istream& is, SmallMatrix<T,M,N,S,I>& m)
    { return is>>m.RegView(); }

  template <class T, size_t M, size_t N, StorageType S, IndexStyle I> 
    inline std::string Type(
	const SmallMatrix<T,M,N,S,I>& )
    { 
      std::ostringstream s;
      s << std::string("SmallMatrix<")<<Type(T())<<','<<M<<','<<N;
      s <<','<<Text(S)<<','<<Text(I)<<">"; 
      return s.str();
    }
  template <class T, size_t M, size_t N, int Si, int Sj, bool C> 
    inline std::string Type(
	const GenSmallMatrix<T,M,N,Si,Sj,C>&)
    { 
      std::ostringstream s;
      s << std::string("GenSmallMatrix<")<<Type(T())<<','<<M<<','<<N;
      s <<','<<Si<<','<<Sj<<',';
      s << (C ? "true" : "false") << ">";
      return s.str();
    }
  template <class T, size_t M, size_t N, int Si, int Sj, bool C, IndexStyle I> 
    inline std::string Type(
	const ConstSmallMatrixView<T,M,N,Si,Sj,C,I>&)
    { 
      std::ostringstream s;
      s << std::string("ConstSmallMatrixView<")<<Type(T())<<','<<M<<','<<N;
      s <<','<<Si<<','<<Sj<<',';
      s << (C ? "true" : "false") <<','<<Text(I)<< ">";
      return s.str();
    }
  template <class T, size_t M, size_t N, int Si, int Sj, bool C, IndexStyle I> 
    inline std::string Type(
	const SmallMatrixView<T,M,N,Si,Sj,C,I>&)
    { 
      std::ostringstream s;
      s << std::string("SmallMatrixView<")<<Type(T())<<','<<M<<','<<N;
      s <<','<<Si<<','<<Sj<<',';
      s << (C ? "true" : "false") <<','<<Text(I)<< ">";
      return s.str();
    }

} // namespace tmv

#endif
