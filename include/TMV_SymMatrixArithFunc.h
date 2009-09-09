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


#ifndef TMV_SymMatrixArithFunc_H
#define TMV_SymMatrixArithFunc_H

#include "TMV_MatrixArithFunc.h"

#define CT std::complex<T>

namespace tmv {

  // y (+)= alpha * A * x
  template <bool add, class T, class Ta, class Tx> void MultMV(
      const T alpha, const GenSymMatrix<Ta>& A, 
      const GenVector<Tx>& x, const VectorView<T>& y);

  // A = alpha * A
  template <class T> inline void MultXM(
      const T alpha, const SymMatrixView<T>& A)
  { MultXM(alpha,A.UpperTri()); }

  // B += alpha * A
  template <class T, class Ta> inline void AddMM(
      const T alpha, const GenSymMatrix<Ta>& A, const SymMatrixView<T>& B)
  { AddMM(alpha,A.UpperTri(),B.UpperTri()); }
  template <class T, class Ta> void AddMM(
      const T alpha, const GenSymMatrix<Ta>& A, const MatrixView<T>& B);
  // C = alpha * A + beta * B
  template <class T, class Ta, class Tb> inline void AddMM(
      const T alpha, const GenSymMatrix<Ta>& A, 
      const T beta, const GenSymMatrix<Tb>& B, const SymMatrixView<T>& C)
  { AddMM(alpha,A.UpperTri(),beta,B.UpperTri(),C.UpperTri()); }
  template <class T, class Ta, class Tb> void AddMM(
      const T alpha, const GenSymMatrix<Ta>& A, 
      const T beta, const GenSymMatrix<Tb>& B, const MatrixView<T>& C);
  template <class T, class Ta, class Tb> void AddMM(
      const T alpha, const GenSymMatrix<Ta>& A, 
      const T beta, const GenMatrix<Tb>& B, const MatrixView<T>& C);
  template <class T, class Ta, class Tb> inline void AddMM(
      const T alpha, const GenMatrix<Ta>& A, 
      const T beta, const GenSymMatrix<Tb>& B, const MatrixView<T>& C)
  { AddMM(beta,B,alpha,A,C); }

  // C (+)= alpha * A * B
  template <bool add, class T, class Ta, class Tb> void MultMM(
      const T alpha, const GenSymMatrix<Ta>& A, 
      const GenMatrix<Tb>& B, const MatrixView<T>& C);
  template <bool add, class T, class Ta, class Tb> inline void MultMM(
      const T alpha, const GenMatrix<Ta>& A, 
      const GenSymMatrix<Tb>& B, const MatrixView<T>& C)
  { MultMM<add>(alpha,B.Transpose(),A.Transpose(),C.Transpose()); }
  template <bool add, class T, class Ta, class Tb> void MultMM(
      const T alpha, const GenSymMatrix<Ta>& A, 
      const GenSymMatrix<Tb>& B, const MatrixView<T>& );

  // A (+)= alpha * (x^xT) (or x* if A is Herm)
  template <bool add, class T, class Tx> void Rank1Update(
      const T alpha, const GenVector<Tx>& x, const SymMatrixView<T>& A);
  // B (+)= alpha * (A * AT) (or At if B is Herm)
  template <bool add, class T, class Ta> void RankKUpdate(
      const T alpha, const GenMatrix<Ta>& A, const SymMatrixView<T>& B);
  template <bool add, class T, class Ta> void RankKUpdate(
      const T alpha, const GenLowerTriMatrix<Ta>& A, 
      const SymMatrixView<T>& B);
  template <bool add, class T, class Ta> void RankKUpdate(
      const T alpha, const GenUpperTriMatrix<Ta>& A, 
      const SymMatrixView<T>& B);

  // These two don't have += forms: they must called explicitly
  // A (+)= alpha * (x^y + y^x) (or x^y* + y^x* is A is Herm)
  template <bool add, class T, class Tx, class Ty> void Rank2Update(
      const T alpha, const GenVector<Tx>& x, const GenVector<Ty>& y, 
      const SymMatrixView<T>& A);
  // C (+)= alpha * (A * BT + B*AT) (or A*Bt + B*At if C is Herm)
  template <bool add, class T, class Ta, class Tb> void Rank2KUpdate(
      const T alpha, const GenMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const SymMatrixView<T>& C);

  // C (+)= alpha * A * B
  // This also needs to be called explicitly.
  // This is to prevent the programmer from doing this accidentally
  // when alpha * A * B is not necessarily symmetrix/hermitian.
  template <bool add, class T, class Ta, class Tb> void SymMultMM(
      const T alpha, const GenMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const SymMatrixView<T>& C);

  template <class T, class Ta> inline void ElementProd(
      const T alpha, const GenSymMatrix<Ta>& A, const SymMatrixView<T>& B)
  { ElementProd(alpha,A.UpperTri(),B.UpperTri()); }
  template <class T, class Ta, class Tb> inline void AddElementProd(
      const T alpha, const GenSymMatrix<Ta>& A, const GenSymMatrix<Tb>& B,
      const SymMatrixView<T>& C)
  { AddElementProd(alpha,A.UpperTri(),B.UpperTri(),C.UpperTri()); }
 
  template <class T> struct SymMatrixComposite :
    public MatrixComposite<T>,
    virtual public AssignableToSymMatrix<T>
  {
    inline size_t colsize() const { return this->size(); }
    inline size_t rowsize() const { return this->size(); }
  };

  // Specialize allowed complex combinations:
  template <bool add, class T, class Ta, class Tx> inline void MultMV(
      const T alpha, const GenSymMatrix<Ta>& A, 
      const GenVector<Tx>& x, const VectorView<CT>& y)
  { MultMV<add>(CT(alpha),A,x,y); }

  template <class T> inline void MultXM(
      const T alpha, const SymMatrixView<CT>& A)
  { MultXM(CT(alpha),A); }

  template <class T, class Ta> inline void AddMM(
      const T alpha, const GenSymMatrix<Ta>& A, const SymMatrixView<CT>& B)
  { AddMM(CT(alpha),A,B); }
  template <class T, class Ta> inline void AddMM(
      const T alpha, const GenSymMatrix<Ta>& A, const MatrixView<CT>& B)
  { AddMM(CT(alpha),A,B); }
  template <class T, class Ta, class Tb> inline void AddMM(
      const T alpha, const GenSymMatrix<Ta>& A, 
      const T beta, const GenSymMatrix<Tb>& B, const SymMatrixView<CT>& C)
  { AddMM(CT(alpha),A,CT(beta),B,C); }
  template <class T, class Ta, class Tb> inline void AddMM(
      const T alpha, const GenSymMatrix<Ta>& A, 
      const T beta, const GenSymMatrix<Tb>& B, const MatrixView<CT>& C)
  { AddMM(CT(alpha),A,CT(beta),B,C); }
  template <class T, class Ta, class Tb> inline void AddMM(
      const T alpha, const GenSymMatrix<Ta>& A, 
      const T beta, const GenMatrix<Tb>& B, const MatrixView<CT>& C)
  { AddMM(CT(alpha),A,CT(beta),B,C); }
  template <class T, class Ta, class Tb> inline void AddMM(
      const T alpha, const GenMatrix<Ta>& A, 
      const T beta, const GenSymMatrix<Tb>& B, const MatrixView<CT>& C)
  { AddMM(beta,B,alpha,A,C); }
  template <class T> inline void AddMM(
      const CT alpha, const GenSymMatrix<CT>& A, 
      const CT beta, const GenSymMatrix<T>& B, const MatrixView<CT>& C)
  { AddMM(beta,B,alpha,A,C); }

  template <bool add, class T, class Ta, class Tb> inline void MultMM(
      const T alpha, const GenSymMatrix<Ta>& A, 
      const GenMatrix<Tb>& B, const MatrixView<CT>& C)
  { MultMM<add>(CT(alpha),A,B,C); }
  template <bool add, class T, class Ta, class Tb> inline void MultMM(
      const T alpha, const GenMatrix<Ta>& A, 
      const GenSymMatrix<Tb>& B, const MatrixView<CT>& C)
  { MultMM<add>(CT(alpha),A,B,C); }
  template <bool add, class T, class Ta, class Tb> inline void MultMM(
      const T alpha, const GenSymMatrix<Ta>& A, 
      const GenSymMatrix<Tb>& B, const MatrixView<CT>& C)
  { MultMM<add>(CT(alpha),A,B,C); }

  template <bool add, class T, class Tx> inline void Rank1Update(
      const T alpha, const GenVector<Tx>& x, const SymMatrixView<CT>& A)
  { Rank1Update<add>(CT(alpha),x,A); }
  template <bool add, class T, class Ta> inline void RankKUpdate(
      const T alpha, const GenMatrix<Ta>& A, const SymMatrixView<CT>& B)
  { RankKUpdate<add>(CT(alpha),A,B); }
  template <bool add, class T, class Ta> inline void RankKUpdate(
      const T alpha, const GenLowerTriMatrix<Ta>& A, 
      const SymMatrixView<CT>& B)
  { RankKUpdate<add>(CT(alpha),A,B); }
  template <bool add, class T, class Ta> inline void RankKUpdate(
      const T alpha, const GenUpperTriMatrix<Ta>& A, 
      const SymMatrixView<CT>& B)
  { RankKUpdate<add>(CT(alpha),A,B); }

  template <bool add, class T, class Tx, class Ty> inline void Rank2Update(
      const T alpha, const GenVector<Tx>& x, const GenVector<Ty>& y, 
      const SymMatrixView<CT>& A)
  { Rank2Update<add>(CT(alpha),x,y,A); }
  template <bool add, class T, class Ta, class Tb> inline void Rank2KUpdate(
      const T alpha, const GenMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const SymMatrixView<CT>& C)
  { Rank2KUpdate<add>(CT(alpha),A,B,C); }

  template <bool add, class T, class Ta, class Tb> inline void SymMultMM(
      const T alpha, const GenMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const SymMatrixView<CT>& C)
  { SymMultMM<add>(CT(alpha),A,B,C); }

  template <class T, class Ta> inline void ElementProd(
      const T alpha, const GenSymMatrix<Ta>& A, const SymMatrixView<CT>& B)
  { ElementProd(CT(alpha),A,B); }
  template <class T, class Ta, class Tb> inline void AddElementProd(
      const T alpha, const GenSymMatrix<Ta>& A, const GenSymMatrix<Tb>& B,
      const SymMatrixView<CT>& C)
  { AddElementProd(CT(alpha),A,B,C); }
 

  // Specialize disallowed complex combinations:
  template <bool add, class T, class Ta, class Tb> inline void MultMV(
      const CT , const GenSymMatrix<Ta>& , 
      const GenVector<Tb>& , const VectorView<T>& )
  { TMVAssert(FALSE); }
  
  template <class T> inline void MultXM(const CT , 
      const SymMatrixView<T>& )
  { TMVAssert(FALSE); }

  template <class T, class Ta, class Tb> inline void AddMM(
      const CT , const GenSymMatrix<Ta>& , 
      const CT , const GenSymMatrix<Tb>& , const SymMatrixView<T>& )
  { TMVAssert(FALSE); }
  template <class T, class Ta, class Tb> inline void AddMM(
      const CT , const GenSymMatrix<Ta>& , 
      const CT , const GenSymMatrix<Tb>& , const MatrixView<T>& )
  { TMVAssert(FALSE); }
  template <class T, class Ta, class Tb> inline void AddMM(
      const CT , const GenSymMatrix<Ta>& , 
      const CT , const GenMatrix<Tb>& , const MatrixView<T>& )
  { TMVAssert(FALSE); }
  template <class T, class Ta, class Tb> inline void AddMM(
      const CT , const GenMatrix<Ta>& , 
      const CT , const GenSymMatrix<Tb>& , const MatrixView<T>& )
  { TMVAssert(FALSE); }

  template <bool add, class T, class Ta, class Tb> inline void MultMM(
      const CT , const GenSymMatrix<Ta>& , 
      const GenMatrix<Tb>& , const MatrixView<T>& )
  { TMVAssert(FALSE); }
  template <bool add, class T, class Ta, class Tb> inline void MultMM(
      const CT , const GenMatrix<Ta>& , 
      const GenSymMatrix<Tb>& , const MatrixView<T>& )
  { TMVAssert(FALSE); }
  template <bool add, class T, class Ta, class Tb> inline void MultMM(
      const CT , const GenSymMatrix<Ta>& , 
      const GenSymMatrix<Tb>& , const MatrixView<T>& )
  { TMVAssert(FALSE); }


} // namespace tmv

#undef CT

#endif
