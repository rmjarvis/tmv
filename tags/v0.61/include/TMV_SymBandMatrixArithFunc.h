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


#ifndef TMV_SymBandMatrixArithFunc_H
#define TMV_SymBandMatrixArithFunc_H

#include "TMV_BaseSymBandMatrix.h"
#include "TMV_BandMatrixArithFunc.h"

#define CT std::complex<T>

namespace tmv {

  // y (+)= alpha * A * x
  template <bool add, class T, class Ta, class Tx> void MultMV(
      const T alpha, const GenSymBandMatrix<Ta>& A, 
      const GenVector<Tx>& x, const VectorView<T>& y);

  // A = alpha * A
  template <class T> inline void MultXM(
      const T alpha, const SymBandMatrixView<T>& A)
  { MultXM(alpha,A.UpperBand()); }

  // B += alpha * A
  template <class T, class Ta> inline void AddMM(
      const T alpha, const GenSymBandMatrix<Ta>& A,
      const SymBandMatrixView<T>& B)
  { AddMM(alpha,A.UpperBand(),B.UpperBand()); }
  template <class T, class Ta> void AddMM(
      const T alpha, const GenSymBandMatrix<Ta>& A, const BandMatrixView<T>& B);
  template <class T, class Ta> inline void AddMM(
      const T alpha, const GenSymBandMatrix<Ta>& A, const MatrixView<T>& B)
  { AddMM(alpha,A,BandMatrixViewOf(B,A.nlo(),A.nhi())); }

  // C = alpha * A + beta * B
  template <class T, class Ta, class Tb> inline void AddMM(
      const T alpha, const GenSymBandMatrix<Ta>& A, 
      const T beta, const GenSymBandMatrix<Tb>& B,
      const SymBandMatrixView<T>& C)
  { AddMM(alpha,A.UpperBand(),beta,B.UpperBand(),C.UpperBand()); }
  template <class T, class Ta, class Tb> void AddMM(
      const T alpha, const GenSymBandMatrix<Ta>& A, 
      const T beta, const GenSymBandMatrix<Tb>& B, const BandMatrixView<T>& C);
  template <class T, class Ta, class Tb> void AddMM(
      const T alpha, const GenSymBandMatrix<Ta>& A, 
      const T beta, const GenSymBandMatrix<Tb>& B, const MatrixView<T>& C);
  template <class T, class Ta, class Tb> void AddMM(
      const T alpha, const GenSymBandMatrix<Ta>& A, 
      const T beta, const GenBandMatrix<Tb>& B, const BandMatrixView<T>& C);
  template <class T, class Ta, class Tb> void AddMM(
      const T alpha, const GenSymBandMatrix<Ta>& A, 
      const T beta, const GenMatrix<Tb>& B, const MatrixView<T>& C);
  template <class T, class Ta, class Tb> inline void AddMM(
      const T alpha, const GenBandMatrix<Ta>& A, 
      const T beta, const GenSymBandMatrix<Tb>& B, const BandMatrixView<T>& C)
  { AddMM(beta,B,alpha,A,C); }
  template <class T, class Ta, class Tb> inline void AddMM(
      const T alpha, const GenMatrix<Ta>& A, 
      const T beta, const GenSymBandMatrix<Tb>& B, const MatrixView<T>& C)
  { AddMM(beta,B,alpha,A,C); }


  // C (+)= alpha * A * B
  template <bool add, class T, class Ta, class Tb> void MultMM(
      const T alpha, const GenSymBandMatrix<Ta>& A, 
      const GenMatrix<Tb>& B, const MatrixView<T>& C);
  template <bool add, class T, class Ta, class Tb> inline void MultMM(
      const T alpha, const GenMatrix<Ta>& A, 
      const GenSymBandMatrix<Tb>& B, const MatrixView<T>& C)
  { MultMM<add>(alpha,B.Transpose(),A.Transpose(),C.Transpose()); }
  template <bool add, class T, class Ta, class Tb> void MultMM(
      const T alpha, const GenSymBandMatrix<Ta>& A, 
      const GenSymBandMatrix<Tb>& B, const BandMatrixView<T>& C);
  template <bool add, class T, class Ta, class Tb> void MultMM(
      const T alpha, const GenSymBandMatrix<Ta>& A, 
      const GenBandMatrix<Tb>& B, const BandMatrixView<T>& C);
  template <bool add, class T, class Ta, class Tb> inline void MultMM(
      const T alpha, const GenBandMatrix<Ta>& A, 
      const GenSymBandMatrix<Tb>& B, const BandMatrixView<T>& C)
  { MultMM<add>(alpha,B.Transpose(),A.Transpose(),C.Transpose()); }

  template <class T, class Ta> inline void ElementProd(
      const T alpha, const GenSymBandMatrix<Ta>& A,
      const SymBandMatrixView<T>& B)
  { ElementProd(alpha,A.UpperBand(),B.UpperBand()); }
  template <class T, class Ta, class Tb> inline void AddElementProd(
      const T alpha, const GenSymBandMatrix<Ta>& A,
      const GenSymBandMatrix<Tb>& B,
      const SymBandMatrixView<T>& C)
  { AddElementProd(alpha,A.UpperBand(),B.UpperBand(),C.UpperBand()); }
 
  template <class T> class SymBandMatrixComposite :
    public GenSymBandMatrix<T>
  {
    public :
      inline SymBandMatrixComposite() : itsm1(0), itsm(0) {}
      inline SymBandMatrixComposite(const SymBandMatrixComposite<T>&) :
	itsm1(0), itsm(0) {}
      virtual inline ~SymBandMatrixComposite() {}

      // Definitions are in TMV_MultsBV.cpp
      const T* cptr() const;
      int stepi() const;
      int stepj() const;
      int diagstep() const;

      inline SymType sym() const { return Sym; }
      inline UpLoType uplo() const { return Lower; }
      inline ConjType ct() const { return NonConj; }

    private:
      mutable auto_array<T> itsm1;
      mutable T* itsm;

      SymBandMatrixComposite<T>& operator=(const SymBandMatrixComposite<T>&)
      { TMVAssert(FALSE); return *this; }
  };


  template <class T> class SymBandMatrixComposite<std::complex<T> > :
    public BandMatrixComposite<std::complex<T> >,
    virtual public AssignableToSymBandMatrix<std::complex<T> >
  {
    public :
      // Need to respecify that size() is pure virtual here, 
      // since GenBandMatrix defines size for compaibility
      // with AssignableToDiagMatrix, etc.
      virtual size_t size() const = 0;

      inline size_t colsize() const { return size(); }
      inline size_t rowsize() const { return size(); }
      inline int nhi() const { return this->nlo(); }
      inline void AssignToS(const SymMatrixView<RealType(T)>& m0) const
      {
	TMVAssert(m0.size() == size());
	TMVAssert(m0.sym() == this->sym());
	this->AssignTosB(SymBandMatrixViewOf(m0,this->nlo()));
	if (int(m0.size()) > this->nlo()+1)
	  BandMatrixViewOf(m0.UpperTri()).Diags(this->nlo()+1,m0.size()).Zero();
      }
      inline void AssignToS(const SymMatrixView<ComplexType(T)>& m0) const
      {
	TMVAssert(m0.size() == size());
	TMVAssert(m0.sym() == this->sym());
	this->AssignTosB(SymBandMatrixViewOf(m0,this->nlo()));
	if (int(m0.size()) > this->nlo()+1)
	  BandMatrixViewOf(m0.UpperTri()).Diags(this->nlo()+1,m0.size()).Zero();
      }

      inline void AssignToM(const MatrixView<RealType(T)>& m0) const
      { BandMatrixComposite<std::complex<T> >::AssignToM(m0); }
      inline void AssignToM(const MatrixView<ComplexType(T)>& m0) const
      { BandMatrixComposite<std::complex<T> >::AssignToM(m0); }
      inline void AssignToB(const BandMatrixView<RealType(T)>& m0) const
      { BandMatrixComposite<std::complex<T> >::AssignToB(m0); }
      inline void AssignToB(const BandMatrixView<ComplexType(T)>& m0) const
      { BandMatrixComposite<std::complex<T> >::AssignToB(m0); }
  };

  // Specialize allowed complex combinations:
  template <bool add, class T, class Ta, class Tx> inline void MultMV(
      const T alpha, const GenSymBandMatrix<Ta>& A, 
      const GenVector<Tx>& x, const VectorView<CT>& y)
  { MultMV<add>(CT(alpha),A,x,y); }

  template <class T> inline void MultXM(
      const T alpha, const SymBandMatrixView<CT>& A)
  { MultXM(CT(alpha),A); }

  template <class T, class Ta> inline void AddMM(
      const T alpha, const GenSymBandMatrix<Ta>& A,
      const SymBandMatrixView<CT>& B)
  { AddMM(CT(alpha),A,B); }
  template <class T, class Ta> inline void AddMM(
      const T alpha, const GenSymBandMatrix<Ta>& A, 
      const BandMatrixView<CT>& B)
  { AddMM(CT(alpha),A,B); }
  template <class T, class Ta> inline void AddMM(
      const T alpha, const GenSymBandMatrix<Ta>& A, const MatrixView<CT>& B)
  { AddMM(CT(alpha),A,B); }

  template <class T, class Ta, class Tb> inline void AddMM(
      const T alpha, const GenSymBandMatrix<Ta>& A, 
      const T beta, const GenSymBandMatrix<Tb>& B,
      const SymBandMatrixView<CT>& C)
  { AddMM(CT(alpha),A,CT(beta),B,C); }
  template <class T, class Ta, class Tb> inline void AddMM(
      const T alpha, const GenSymBandMatrix<Ta>& A, 
      const T beta, const GenSymBandMatrix<Tb>& B, const BandMatrixView<CT>& C)
  { AddMM(CT(alpha),A,CT(beta),B,C); }
  template <class T, class Ta, class Tb> inline void AddMM(
      const T alpha, const GenSymBandMatrix<Ta>& A, 
      const T beta, const GenSymBandMatrix<Tb>& B, const MatrixView<CT>& C)
  { AddMM(CT(alpha),A,CT(beta),B,C); }
  template <class T, class Ta, class Tb> inline void AddMM(
      const T alpha, const GenSymBandMatrix<Ta>& A, 
      const T beta, const GenBandMatrix<Tb>& B, const BandMatrixView<CT>& C)
  { AddMM(CT(alpha),A,CT(beta),B,C); }
  template <class T, class Ta, class Tb> inline void AddMM(
      const T alpha, const GenSymBandMatrix<Ta>& A, 
      const T beta, const GenMatrix<Tb>& B, const MatrixView<CT>& C)
  { AddMM(CT(alpha),A,CT(beta),B,C); }
  template <class T, class Ta, class Tb> inline void AddMM(
      const T alpha, const GenBandMatrix<Ta>& A, 
      const T beta, const GenSymBandMatrix<Tb>& B, const BandMatrixView<CT>& C)
  { AddMM(CT(alpha),A,CT(beta),B,C); }
  template <class T, class Ta, class Tb> inline void AddMM(
      const T alpha, const GenMatrix<Ta>& A, 
      const T beta, const GenSymBandMatrix<Tb>& B, const MatrixView<CT>& C)
  { AddMM(CT(alpha),A,CT(beta),B,C); }
  template <class T> inline void AddMM(
      const CT alpha, const GenSymBandMatrix<CT>& A, 
      const CT beta, const GenSymBandMatrix<T>& B, const MatrixView<CT>& C)
  { AddMM(beta,B,alpha,A,C); }
  template <class T> inline void AddMM(
      const CT alpha, const GenSymBandMatrix<CT>& A,
      const CT beta, const GenSymBandMatrix<T>& B, const BandMatrixView<CT>& C)
  { AddMM(beta,B,alpha,A,C); }

  template <bool add, class T, class Ta, class Tb> inline void MultMM(
      const T alpha, const GenSymBandMatrix<Ta>& A, 
      const GenMatrix<Tb>& B, const MatrixView<CT>& C)
  { MultMM<add>(CT(alpha),A,B,C); }
  template <bool add, class T, class Ta, class Tb> inline void MultMM(
      const T alpha, const GenMatrix<Ta>& A, 
      const GenSymBandMatrix<Tb>& B, const MatrixView<CT>& C)
  { MultMM<add>(CT(alpha),A,B,C); }
  template <bool add, class T, class Ta, class Tb> inline void MultMM(
      const T alpha, const GenSymBandMatrix<Ta>& A, 
      const GenSymBandMatrix<Tb>& B, const BandMatrixView<CT>& C)
  { MultMM<add>(CT(alpha),A,B,C); }
  template <bool add, class T, class Ta, class Tb> inline void MultMM(
      const T alpha, const GenSymBandMatrix<Ta>& A, 
      const GenBandMatrix<Tb>& B, const BandMatrixView<CT>& C)
  { MultMM<add>(CT(alpha),A,B,C); }
  template <bool add, class T, class Ta, class Tb> inline void MultMM(
      const T alpha, const GenBandMatrix<Ta>& A, 
      const GenSymBandMatrix<Tb>& B, const BandMatrixView<CT>& C)
  { MultMM<add>(CT(alpha),A,B,C); }

  template <class T, class Ta> inline void ElementProd(
      const T alpha, const GenSymBandMatrix<Ta>& A,
      const SymBandMatrixView<CT>& B)
  { ElementProd(CT(alpha),A,B); }
  template <class T, class Ta, class Tb> inline void AddElementProd(
      const T alpha, const GenSymBandMatrix<Ta>& A,
      const GenSymBandMatrix<Tb>& B, const SymBandMatrixView<CT>& C)
  { AddElementProd(CT(alpha),A,B,C); }

  // Specialize disallowed complex combinations:
  template <bool add, class T, class Ta, class Tb> inline void MultMV(
      const CT , const GenSymBandMatrix<Ta>& , 
      const GenVector<Tb>& , const VectorView<T>& )
  { TMVAssert(FALSE); }

  template <class T> inline void MultXM(const CT , 
      const SymBandMatrixView<T>& )
  { TMVAssert(FALSE); }

  template <class T, class Ta> inline void AddMM(
      const CT , const GenSymBandMatrix<Ta>& , const SymBandMatrixView<T>& )
  { TMVAssert(FALSE); }
  template <class T, class Ta> inline void AddMM(
      const CT , const GenSymBandMatrix<Ta>& , const BandMatrixView<T>& )
  { TMVAssert(FALSE); }
  template <class T, class Ta> inline void AddMM(
      const CT , const GenSymBandMatrix<Ta>& , const MatrixView<T>& )
  { TMVAssert(FALSE); }
  template <class T, class Ta, class Tb> inline void AddMM(
      const CT , const GenSymBandMatrix<Ta>& , 
      const CT , const GenSymBandMatrix<Tb>& , const SymBandMatrixView<T>& )
  { TMVAssert(FALSE); }
  template <class T, class Ta, class Tb> inline void AddMM(
      const CT , const GenSymBandMatrix<Ta>& , 
      const CT , const GenSymBandMatrix<Tb>& , const BandMatrixView<T>& )
  { TMVAssert(FALSE); }
  template <class T, class Ta, class Tb> inline void AddMM(
      const CT , const GenSymBandMatrix<Ta>& , 
      const CT , const GenSymBandMatrix<Tb>& , const MatrixView<T>& )
  { TMVAssert(FALSE); }
  template <class T, class Ta, class Tb> inline void AddMM(
      const CT , const GenSymBandMatrix<Ta>& , 
      const CT , const GenBandMatrix<Tb>& , const BandMatrixView<T>& )
  { TMVAssert(FALSE); }
  template <class T, class Ta, class Tb> inline void AddMM(
      const CT , const GenSymBandMatrix<Ta>& , 
      const CT , const GenMatrix<Tb>& , const MatrixView<T>& )
  { TMVAssert(FALSE); }
  template <class T, class Ta, class Tb> inline void AddMM(
      const CT , const GenBandMatrix<Ta>& , 
      const CT , const GenSymBandMatrix<Tb>& , const BandMatrixView<T>& )
  { TMVAssert(FALSE); }
  template <class T, class Ta, class Tb> inline void AddMM(
      const CT , const GenMatrix<Ta>& , 
      const CT , const GenSymBandMatrix<Tb>& , const MatrixView<T>& )
  { TMVAssert(FALSE); }

  template <bool add, class T, class Ta, class Tb> inline void MultMM(
      const CT , const GenSymBandMatrix<Ta>& , 
      const GenMatrix<Tb>& , const MatrixView<T>& )
  { TMVAssert(FALSE); }
  template <bool add, class T, class Ta, class Tb> inline void MultMM(
      const CT , const GenMatrix<Ta>& , 
      const GenSymBandMatrix<Tb>& , const MatrixView<T>& )
  { TMVAssert(FALSE); }
  template <bool add, class T, class Ta, class Tb> inline void MultMM(
      const CT , const GenSymBandMatrix<Ta>& , 
      const GenSymBandMatrix<Tb>& , const BandMatrixView<T>& )
  { TMVAssert(FALSE); }
  template <bool add, class T, class Ta, class Tb> inline void MultMM(
      const CT , const GenSymBandMatrix<Ta>& , 
      const GenBandMatrix<Tb>& , const BandMatrixView<T>& )
  { TMVAssert(FALSE); }
  template <bool add, class T, class Ta, class Tb> inline void MultMM(
      const CT , const GenBandMatrix<Ta>& , 
      const GenSymBandMatrix<Tb>& , const BandMatrixView<T>& )
  { TMVAssert(FALSE); }

} // namespace tmv

#undef CT

#endif
