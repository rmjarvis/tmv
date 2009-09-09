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


#ifndef TMV_SmallMatrixArithFunc_H
#define TMV_SmallMatrixArithFunc_H

#define CT std::complex<T>

namespace tmv {

  // y = alpha * A * x
  template <class Ta, class Tx, class Ty, int M, int N, int Sia, int Sja, bool Ca, int Sx, bool Cx, int Sy, bool Cy> 
    inline void MultMV_1(
	const GenSmallMatrix<Ta,M,N,Sia,Sja,Ca>& A, 
	const GenSmallVector<Tx,N,Sx,Cx>& x,
	const SmallVectorView<Ty,M,Sy,Cy>& y)
    {
      if (Sia == 1) {
	if (N==0) y.Zero();
	else {
	  y = x(0) * A.col(0);
	  for(int j=1;j<N;++j) y += x(j) * A.col(j);
	}
      } else {
	for(int i=0;i<M;++i) y(i) = A.row(i) * x; 
      }
    }
  template <class Ta, class Tx, class Ty, int M, int N, int Sia, int Sja, bool Ca, int Sx, bool Cx, int Sy, bool Cy> 
    inline void MultMV_m1(
	const GenSmallMatrix<Ta,M,N,Sia,Sja,Ca>& A, 
	const GenSmallVector<Tx,N,Sx,Cx>& x,
	const SmallVectorView<Ty,M,Sy,Cy>& y)
    {
      if (Sia == 1) {
	if (N==0) y.Zero();
	else {
	  y = (-x(0)) * A.col(0);
	  for(int j=1;j<N;++j) y += (-x(j)) * A.col(j);
	}
      } else {
	for(int i=0;i<M;++i) y(i) = -A.row(i) * x; 
      }
    }
  template <class T, class Ta, class Tx, class Ty, int M, int N, int Sia, int Sja, bool Ca, int Sx, bool Cx, int Sy, bool Cy> 
    inline void MultMV(
	const T alpha, const GenSmallMatrix<Ta,M,N,Sia,Sja,Ca>& A, 
	const GenSmallVector<Tx,N,Sx,Cx>& x,
	const SmallVectorView<Ty,M,Sy,Cy>& y)
    {
      if (Sia == 1) {
	MultMV_1(A,x,y);
	y *= alpha;
      } else {
	for(int i=0;i<M;++i) y(i) = alpha * (A.row(i) * x); 
      }
    }
  template <class T, class Tx, int M, int N, int Sia, int Sja, bool Ca, int Sx, bool Cx, int Sy, bool Cy> 
    inline void MultMV_1(
	const GenSmallMatrix<CT,M,N,Sia,Sja,Ca>& , 
	const GenSmallVector<Tx,N,Sx,Cx>& ,
	const SmallVectorView<T,M,Sy,Cy>& )
    { TMVAssert(FALSE); }
  template <class T, int M, int N, int Sia, int Sja, bool Ca, int Sx, bool Cx, int Sy, bool Cy> 
    inline void MultMV_1(
	const GenSmallMatrix<T,M,N,Sia,Sja,Ca>& , 
	const GenSmallVector<CT,N,Sx,Cx>& ,
	const SmallVectorView<T,M,Sy,Cy>& )
    { TMVAssert(FALSE); }
  template <class T, class Tx, int M, int N, int Sia, int Sja, bool Ca, int Sx, bool Cx, int Sy, bool Cy> 
    inline void MultMV_m1(
	const GenSmallMatrix<CT,M,N,Sia,Sja,Ca>& , 
	const GenSmallVector<Tx,N,Sx,Cx>& ,
	const SmallVectorView<T,M,Sy,Cy>& )
    { TMVAssert(FALSE); }
  template <class T, int M, int N, int Sia, int Sja, bool Ca, int Sx, bool Cx, int Sy, bool Cy> 
    inline void MultMV_m1(
	const GenSmallMatrix<T,M,N,Sia,Sja,Ca>& , 
	const GenSmallVector<CT,N,Sx,Cx>& ,
	const SmallVectorView<T,M,Sy,Cy>& )
    { TMVAssert(FALSE); }
  template <class T, class Ta, class Tx, int M, int N, int Sia, int Sja, bool Ca, int Sx, bool Cx, int Sy, bool Cy> 
    inline void MultMV(
	const CT , const GenSmallMatrix<Ta,M,N,Sia,Sja,Ca>& , 
	const GenSmallVector<Tx,N,Sx,Cx>& ,
	const SmallVectorView<T,M,Sy,Cy>& )
    { TMVAssert(FALSE); }
  template <class Ta, class Tx, class Ty, int M, int N, int Sia, int Sja, bool Ca, int Sx, bool Cx> 
    inline void MultMV_1(
	const GenSmallMatrix<Ta,M,N,Sia,Sja,Ca>& A, 
	const GenSmallVector<Tx,N,Sx,Cx>& x, Ty* yp, int Sy)
    { 
      if (Sia == 1) {
	if (N==0) 
	  for(int i=0,ii=0;i<M;++i,ii+=Sy) yp[ii] = Ty(0);
	else {
	  MultXV(x(0),A.col(0),yp,Sy);
	  for(int j=1;j<N;++j) 
	    AddVV(x(j),A.col(j),yp,Sy);
	}
      } else {
	for(int i=0,ii=0;i<M;++i,ii+=Sy) yp[ii] = A.row(i) * x; 
      }
    }
  template <class Ta, class Tx, class Ty, int M, int N, int Sia, int Sja, bool Ca, int Sx, bool Cx> 
    inline void MultMV_m1(
	const GenSmallMatrix<Ta,M,N,Sia,Sja,Ca>& A, 
	const GenSmallVector<Tx,N,Sx,Cx>& x, Ty* yp, int Sy)
    { 
      if (Sia == 1) {
	if (N==0) 
	  for(int i=0,ii=0;i<M;++i,ii+=Sy) yp[ii] = Ty(0);
	else {
	  MultXV(-x(0),A.col(0),yp,Sy);
	  for(int j=1;j<N;++j) 
	    AddVV(-x(j),A.col(j),yp,Sy);
	}
      } else {
	for(int i=0,ii=0;i<M;++i,ii+=Sy) yp[ii] = -(A.row(i) * x); 
      } 
    }
  template <class T, class Ta, class Tx, class Ty, int M, int N, int Sia, int Sja, bool Ca, int Sx, bool Cx> 
    inline void MultMV(
	const T alpha, const GenSmallMatrix<Ta,M,N,Sia,Sja,Ca>& A, 
	const GenSmallVector<Tx,N,Sx,Cx>& x, Ty* yp, int Sy)
    {
      if (Sia == 1) {
	MultMV_1(A,x,yp,Sy);
	MultXV<M>(alpha,yp,Sy);
      } else {
	for(int i=0,ii=0;i<M;++i,ii+=Sy) yp[ii] = alpha * (A.row(i) * x); 
      }
    }
  template <class T, class Tx, int M, int N, int Sia, int Sja, bool Ca, int Sx, bool Cx> 
    inline void MultMV_1(
	const GenSmallMatrix<CT,M,N,Sia,Sja,Ca>& , 
	const GenSmallVector<Tx,N,Sx,Cx>& , T* , int )
    { TMVAssert(FALSE); }
  template <class T, int M, int N, int Sia, int Sja, bool Ca, int Sx, bool Cx> 
    inline void MultMV_1(
	const GenSmallMatrix<T,M,N,Sia,Sja,Ca>& , 
	const GenSmallVector<CT,N,Sx,Cx>& , T* , int )
    { TMVAssert(FALSE); }
  template <class T, class Tx, int M, int N, int Sia, int Sja, bool Ca, int Sx, bool Cx> 
    inline void MultMV_m1(
	const GenSmallMatrix<CT,M,N,Sia,Sja,Ca>& , 
	const GenSmallVector<Tx,N,Sx,Cx>& , T* , int )
    { TMVAssert(FALSE); }
  template <class T, int M, int N, int Sia, int Sja, bool Ca, int Sx, bool Cx> 
    inline void MultMV_m1(
	const GenSmallMatrix<T,M,N,Sia,Sja,Ca>& , 
	const GenSmallVector<CT,N,Sx,Cx>& , T* , int )
    { TMVAssert(FALSE); }
  template <class T, class Ta, class Tx, int M, int N, int Sia, int Sja, bool Ca, int Sx, bool Cx> 
    inline void MultMV(
	const CT , const GenSmallMatrix<Ta,M,N,Sia,Sja,Ca>& , 
	const GenSmallVector<Tx,N,Sx,Cx>& , T* , int )
    { TMVAssert(FALSE); }

  // y += alpha * A * x
  template <class Ta, class Tx, class Ty, int M, int N, int Sia, int Sja, bool Ca, int Sx, bool Cx, int Sy, bool Cy> 
    inline void AddMultMV_1(
	const GenSmallMatrix<Ta,M,N,Sia,Sja,Ca>& A, 
	const GenSmallVector<Tx,N,Sx,Cx>& x,
	const SmallVectorView<Ty,M,Sy,Cy>& y)
    {
      if (Sia == 1)
	for(int j=0;j<N;++j) y += x(j) * A.col(j); 
      else
	for(int i=0;i<M;++i) y(i) += A.row(i) * x; 
    } 
  template <class Ta, class Tx, class Ty, int M, int N, int Sia, int Sja, bool Ca, int Sx, bool Cx, int Sy, bool Cy> 
    inline void AddMultMV_m1(
	const GenSmallMatrix<Ta,M,N,Sia,Sja,Ca>& A, 
	const GenSmallVector<Tx,N,Sx,Cx>& x,
	const SmallVectorView<Ty,M,Sy,Cy>& y)
    {
      if (Sia == 1)
	for(int j=0;j<N;++j) y -= x(j) * A.col(j); 
      else
	for(int i=0;i<M;++i) y(i) -= A.row(i) * x; 
    }
  template <class T, class Ta, class Tx, class Ty, int M, int N, int Sia, int Sja, bool Ca, int Sx, bool Cx, int Sy, bool Cy> 
    inline void AddMultMV(
	const T alpha, const GenSmallMatrix<Ta,M,N,Sia,Sja,Ca>& A, 
	const GenSmallVector<Tx,N,Sx,Cx>& x,
	const SmallVectorView<Ty,M,Sy,Cy>& y)
    {
      if (Sia == 1)
	for(int j=0;j<N;++j) y += (alpha*x(j)) * A.col(j); 
      else
	for(int i=0;i<M;++i) y(i) += alpha * (A.row(i) * x); 
    }
  template <class T, class Tx, int M, int N, int Sia, int Sja, bool Ca, int Sx, bool Cx, int Sy, bool Cy> 
    inline void AddMultMV_1(
	const GenSmallMatrix<CT,M,N,Sia,Sja,Ca>& , 
	const GenSmallVector<Tx,N,Sx,Cx>& ,
	const SmallVectorView<T,M,Sy,Cy>& )
    { TMVAssert(FALSE); }
  template <class T, int M, int N, int Sia, int Sja, bool Ca, int Sx, bool Cx, int Sy, bool Cy> 
    inline void AddMultMV_1(
	const GenSmallMatrix<T,M,N,Sia,Sja,Ca>& , 
	const GenSmallVector<CT,N,Sx,Cx>& ,
	const SmallVectorView<T,M,Sy,Cy>& )
    { TMVAssert(FALSE); }
  template <class T, class Tx, int M, int N, int Sia, int Sja, bool Ca, int Sx, bool Cx, int Sy, bool Cy> 
    inline void AddMultMV_m1(
	const GenSmallMatrix<CT,M,N,Sia,Sja,Ca>& , 
	const GenSmallVector<Tx,N,Sx,Cx>& ,
	const SmallVectorView<T,M,Sy,Cy>& )
    { TMVAssert(FALSE); }
  template <class T, int M, int N, int Sia, int Sja, bool Ca, int Sx, bool Cx, int Sy, bool Cy> 
    inline void AddMultMV_m1(
	const GenSmallMatrix<T,M,N,Sia,Sja,Ca>& , 
	const GenSmallVector<CT,N,Sx,Cx>& ,
	const SmallVectorView<T,M,Sy,Cy>& )
    { TMVAssert(FALSE); }
  template <class T, class Ta, class Tx, int M, int N, int Sia, int Sja, bool Ca, int Sx, bool Cx, int Sy, bool Cy> 
    inline void AddMultMV(
	const CT , const GenSmallMatrix<Ta,M,N,Sia,Sja,Ca>& , 
	const GenSmallVector<Tx,N,Sx,Cx>& ,
	const SmallVectorView<T,M,Sy,Cy>& )
    { TMVAssert(FALSE); }
  template <class Ta, class Tx, class Ty, int M, int N, int Sia, int Sja, bool Ca, int Sx, bool Cx> 
    inline void AddMultMV_1(
	const GenSmallMatrix<Ta,M,N,Sia,Sja,Ca>& A, 
	const GenSmallVector<Tx,N,Sx,Cx>& x, Ty* yp, int Sy)
    {
      if (Sia == 1)
	for(int j=0;j<N;++j) AddVV(x(j),A.col(j),yp,Sy); 
      else
	for(int i=0,ii=0;i<M;++i,ii+=Sy) yp[ii] += A.row(i) * x; 
    }
  template <class Ta, class Tx, class Ty, int M, int N, int Sia, int Sja, bool Ca, int Sx, bool Cx> 
    inline void AddMultMV_m1(
	const GenSmallMatrix<Ta,M,N,Sia,Sja,Ca>& A, 
	const GenSmallVector<Tx,N,Sx,Cx>& x, Ty* yp, int Sy)
    {
      if (Sia == 1)
	for(int j=0;j<N;++j) AddVV(-x(j),A.col(j),yp,Sy); 
      else
	for(int i=0,ii=0;i<M;++i,ii+=Sy) yp[ii] -= A.row(i) * x; 
    }
  template <class T, class Ta, class Tx, class Ty, int M, int N, int Sia, int Sja, bool Ca, int Sx, bool Cx> 
    inline void AddMultMV(
	const T alpha, const GenSmallMatrix<Ta,M,N,Sia,Sja,Ca>& A, 
	const GenSmallVector<Tx,N,Sx,Cx>& x, Ty* yp, int Sy)
    {
      if (Sia == 1)
	for(int j=0;j<N;++j) AddVV(alpha*x(j),A.col(j),yp,Sy); 
      else
	for(int i=0,ii=0;i<M;++i,ii+=Sy) yp[ii] += alpha * (A.row(i) * x);
    }
  template <class T, class Tx, int M, int N, int Sia, int Sja, bool Ca, int Sx, bool Cx> 
    inline void AddMultMV_1(
	const GenSmallMatrix<CT,M,N,Sia,Sja,Ca>& , 
	const GenSmallVector<Tx,N,Sx,Cx>& , T* , int )
    { TMVAssert(FALSE); }
  template <class T, int M, int N, int Sia, int Sja, bool Ca, int Sx, bool Cx> 
    inline void AddMultMV_1(
	const GenSmallMatrix<T,M,N,Sia,Sja,Ca>& , 
	const GenSmallVector<CT,N,Sx,Cx>& , T* , int )
    { TMVAssert(FALSE); }
  template <class T, class Tx, int M, int N, int Sia, int Sja, bool Ca, int Sx, bool Cx> 
    inline void AddMultMV_m1(
	const GenSmallMatrix<CT,M,N,Sia,Sja,Ca>& , 
	const GenSmallVector<Tx,N,Sx,Cx>& , T* , int )
    { TMVAssert(FALSE); }
  template <class T, int M, int N, int Sia, int Sja, bool Ca, int Sx, bool Cx> 
    inline void AddMultMV_m1(
	const GenSmallMatrix<T,M,N,Sia,Sja,Ca>& , 
	const GenSmallVector<CT,N,Sx,Cx>& , T* , int )
    { TMVAssert(FALSE); }
  template <class T, class Ta, class Tx, int M, int N, int Sia, int Sja, bool Ca, int Sx, bool Cx> 
    inline void AddMultMV(
	const CT , const GenSmallMatrix<Ta,M,N,Sia,Sja,Ca>& , 
	const GenSmallVector<Tx,N,Sx,Cx>& , T* , int )
    { TMVAssert(FALSE); }

  // A = alpha * A
  template <class T, class Ta, int M, int N, int Si, int Sj, bool C> 
    inline void MultXM(const T alpha, 
	const SmallMatrixView<Ta,M,N,Si,Sj,C>& A)
    {
      if (A.CanLinearize()) A.LinearView() *= alpha;
      else if (Si == 1) for(int j=0;j<N;++j) A.col(j) *= alpha;
      else for(int i=0;i<M;++i) A.row(i) *= alpha;
    }
  template <class T, int M, int N, int Si, int Sj, bool C> 
    inline void MultXM(const CT , const SmallMatrixView<T,M,N,Si,Sj,C>& )
    { TMVAssert(FALSE); }

  template <class T, class Ta, class Tb, int M, int N, int Si1, int Sj1, bool C1> 
    inline void MultXM(const T alpha, 
	const GenSmallMatrix<Ta,M,N,Si1,Sj1,C1>& A,
	Tb* Bp, int Si2, int Sj2)
    { for(int j=0;j<N;++j,Bp+=Sj2) MultXV(alpha,A.col(j),Bp,Si2); }
  template <class T, class Ta, int M, int N, int Si1, int Sj1, bool C1> 
    inline void MultXM(const CT , 
	const GenSmallMatrix<Ta,M,N,Si1,Sj1,C1>& , T* , int , int )
    { TMVAssert(FALSE); }

  // B += alpha * A
  template <class Ta, class Tb, int M, int N, int Sia, int Sja, bool Ca, int Sib, int Sjb, bool Cb> 
    inline void AddMM_1(const GenSmallMatrix<Ta,M,N,Sia,Sja,Ca>& A, 
	const SmallMatrixView<Tb,M,N,Sib,Sjb,Cb>& B)
    {
      if (Sia == Sib && Sja == Sjb && A.CanLinearize() && B.CanLinearize()) 
	B.LinearView() += A.ConstLinearView();
      else if (Sib == 1) 
	for(int j=0;j<N;++j) B.col(j) += A.col(j);
      else 
	for(int i=0;i<M;++i) B.row(i) += A.row(i);
    }
  template <class Ta, class Tb, int M, int N, int Sia, int Sja, bool Ca, int Sib, int Sjb, bool Cb> 
    inline void AddMM_m1(const GenSmallMatrix<Ta,M,N,Sia,Sja,Ca>& A, 
	const SmallMatrixView<Tb,M,N,Sib,Sjb,Cb>& B)
    {
      if (Sia == Sib && Sja == Sjb && A.CanLinearize() && B.CanLinearize()) 
	B.LinearView() -= A.ConstLinearView();
      else if (Sib == 1) 
	for(int j=0;j<N;++j) B.col(j) -= A.col(j);
      else
	for(int i=0;i<M;++i) B.row(i) -= A.row(i);
    }
  template <class T, class Ta, class Tb, int M, int N, int Sia, int Sja, bool Ca, int Sib, int Sjb, bool Cb> 
    inline void AddMM(const T alpha,
	const GenSmallMatrix<Ta,M,N,Sia,Sja,Ca>& A, 
	const SmallMatrixView<Tb,M,N,Sib,Sjb,Cb>& B)
    {
      if (Sia == Sib && Sja == Sjb && A.CanLinearize() && B.CanLinearize()) 
	B.LinearView() += alpha * A.ConstLinearView();
      else if (Sib == 1) 
	for(int j=0;j<N;++j) B.col(j) += alpha * A.col(j);
      else 
	for(int i=0;i<M;++i) B.row(i) += alpha * A.row(i);
    }
  template <class T, int M, int N, int Sia, int Sja, bool Ca, int Sib, int Sjb, bool Cb> 
    inline void AddMM_1(const GenSmallMatrix<CT,M,N,Sia,Sja,Ca>& , 
	const SmallMatrixView<T,M,N,Sib,Sjb,Cb>& )
    { TMVAssert(FALSE); }
  template <class T, int M, int N, int Sia, int Sja, bool Ca, int Sib, int Sjb, bool Cb> 
    inline void AddMM_m1(const GenSmallMatrix<CT,M,N,Sia,Sja,Ca>& , 
	const SmallMatrixView<T,M,N,Sib,Sjb,Cb>& )
    { TMVAssert(FALSE); }
  template <class T, class Ta, int M, int N, int Sia, int Sja, bool Ca, int Sib, int Sjb, bool Cb> 
    inline void AddMM(const CT , const GenSmallMatrix<Ta,M,N,Sia,Sja,Ca>& , 
	const SmallMatrixView<T,M,N,Sib,Sjb,Cb>& )
    { TMVAssert(FALSE); }

  template <class Ta, class Tb, int M, int N, int Sia, int Sja, bool Ca> 
    inline void AddMM_1(const GenSmallMatrix<Ta,M,N,Sia,Sja,Ca>& A, 
	Tb* Bp, int Sib, int Sjb)
    {
      if (Sib == 1) 
	for(int j=0;j<N;++j,Bp+=Sjb) 
	  AddVV_1(A.col(j),SmallVectorViewOf<M>(Bp));
      else if (Sjb == 1) 
	for(int i=0;i<M;++i,Bp+=Sib) 
	  AddVV_1(A.row(i),SmallVectorViewOf<N>(Bp));
      else 
	for(int i=0;i<M;++i,Bp+=Sib) AddVV_1(A.row(i),Bp,Sjb);
    }
  template <class Ta, class Tb, int M, int N, int Sia, int Sja, bool Ca> 
    inline void AddMM_m1(const GenSmallMatrix<Ta,M,N,Sia,Sja,Ca>& A, 
	Tb* Bp, int Sib, int Sjb)
    {
      if (Sib == 1) 
	for(int j=0;j<N;++j,Bp+=Sjb) 
	  AddVV_m1(A.col(j),SmallVectorViewOf<M>(Bp));
      else if (Sjb == 1) 
	for(int i=0;i<M;++i,Bp+=Sib) 
	  AddVV_m1(A.row(i),SmallVectorViewOf<N>(Bp));
      else 
	for(int i=0;i<M;++i,Bp+=Sib) AddVV_m1(A.row(i),Bp,Sjb);
    }
  template <class T, class Ta, class Tb, int M, int N, int Sia, int Sja, bool Ca> 
    inline void AddMM(const T alpha,
	const GenSmallMatrix<Ta,M,N,Sia,Sja,Ca>& A, 
	Tb* Bp, int Sib, int Sjb)
    {
      if (Sib == 1) 
	for(int j=0;j<N;++j,Bp+=Sjb) 
	  AddVV(alpha,A.col(j),SmallVectorViewOf<M>(Bp));
      else if (Sjb == 1) 
	for(int i=0;i<M;++i,Bp+=Sib) 
	  AddVV(alpha,A.row(i),SmallVectorViewOf<N>(Bp));
      else 
	for(int i=0;i<M;++i,Bp+=Sib) AddVV(alpha,A.row(i),Bp,Sjb);
    }
  template <class T, int M, int N, int Sia, int Sja, bool Ca> 
    inline void AddMM_1(const GenSmallMatrix<CT,M,N,Sia,Sja,Ca>& , 
	const T* , int , int )
    { TMVAssert(FALSE); }
  template <class T, int M, int N, int Sia, int Sja, bool Ca> 
    inline void AddMM_m1(const GenSmallMatrix<CT,M,N,Sia,Sja,Ca>& , 
	const T* , int , int )
    { TMVAssert(FALSE); }
  template <class T, class Ta, int M, int N, int Sia, int Sja, bool Ca> 
    inline void AddMM(const CT , const GenSmallMatrix<Ta,M,N,Sia,Sja,Ca>& , 
	const T* , int , int )
    { TMVAssert(FALSE); }

  // C = alpha * A + beta * B
  template <class Ta, class Tb, class Tc, int M, int N, int Sia, int Sja, bool Ca, int Sib, int Sjb, bool Cb, int Sic, int Sjc, bool Cc> 
    inline void AddMM_1_1(
	const GenSmallMatrix<Ta,M,N,Sia,Sja,Ca>& A, 
	const GenSmallMatrix<Tb,M,N,Sib,Sjb,Cb>& B,
	const SmallMatrixView<Tc,M,N,Sic,Sjc,Cc>& C)
    { Copy(A,C); AddMM_1(B,C); }
  template <class Ta, class Tb, class Tc, int M, int N, int Sia, int Sja, bool Ca, int Sib, int Sjb, bool Cb, int Sic, int Sjc, bool Cc> 
    inline void AddMM_1_m1(
	const GenSmallMatrix<Ta,M,N,Sia,Sja,Ca>& A, 
	const GenSmallMatrix<Tb,M,N,Sib,Sjb,Cb>& B,
	const SmallMatrixView<Tc,M,N,Sic,Sjc,Cc>& C)
    { Copy(A,C); AddMM_m1(B,C); }
  template <class T, class Ta, class Tb, class Tc, int M, int N, int Sia, int Sja, bool Ca, int Sib, int Sjb, bool Cb, int Sic, int Sjc, bool Cc> 
    inline void AddMM_1_x(
	const GenSmallMatrix<Ta,M,N,Sia,Sja,Ca>& A, 
	const T beta, const GenSmallMatrix<Tb,M,N,Sib,Sjb,Cb>& B,
	const SmallMatrixView<Tc,M,N,Sic,Sjc,Cc>& C)
    { MultXM(beta,B,C); AddMM_1(A,C); }
  template <class T, class Ta, class Tb, class Tc, int M, int N, int Sia, int Sja, bool Ca, int Sib, int Sjb, bool Cb, int Sic, int Sjc, bool Cc> 
    inline void AddMM_x_1(
	const T alpha, const GenSmallMatrix<Ta,M,N,Sia,Sja,Ca>& A, 
	const GenSmallMatrix<Tb,M,N,Sib,Sjb,Cb>& B,
	const SmallMatrixView<Tc,M,N,Sic,Sjc,Cc>& C)
    { MultXM(alpha,A,C); AddMM_1(B,C); }
  template <class T, class Ta, class Tb, class Tc, int M, int N, int Sia, int Sja, bool Ca, int Sib, int Sjb, bool Cb, int Sic, int Sjc, bool Cc> 
    inline void AddMM_x_m1(
	const T alpha, const GenSmallMatrix<Ta,M,N,Sia,Sja,Ca>& A, 
	const GenSmallMatrix<Tb,M,N,Sib,Sjb,Cb>& B,
	const SmallMatrixView<Tc,M,N,Sic,Sjc,Cc>& C)
    { MultXM(alpha,A,C); AddMM_m1(B,C); }
  template <class T, class Ta, class Tb, class Tc, int M, int N, int Sia, int Sja, bool Ca, int Sib, int Sjb, bool Cb, int Sic, int Sjc, bool Cc> 
    inline void AddMM(
	const T alpha, const GenSmallMatrix<Ta,M,N,Sia,Sja,Ca>& A, 
	const T beta, const GenSmallMatrix<Tb,M,N,Sib,Sjb,Cb>& B,
	const SmallMatrixView<Tc,M,N,Sic,Sjc,Cc>& C)
    { MultXM(alpha,A,C); AddMM(beta,B,C); }
  template <class T, class Tb, int M, int N, int Sia, int Sja, bool Ca, int Sib, int Sjb, bool Cb, int Sic, int Sjc, bool Cc> 
    inline void AddMM_1_1(
	const GenSmallMatrix<CT,M,N,Sia,Sja,Ca>& , 
	const GenSmallMatrix<Tb,M,N,Sib,Sjb,Cb>& ,
	const SmallMatrixView<T,M,N,Sic,Sjc,Cc>& )
    { TMVAssert(FALSE); }
  template <class T, int M, int N, int Sia, int Sja, bool Ca, int Sib, int Sjb, bool Cb, int Sic, int Sjc, bool Cc> 
    inline void AddMM_1_1(
	const GenSmallMatrix<T,M,N,Sia,Sja,Ca>& , 
	const GenSmallMatrix<CT,M,N,Sib,Sjb,Cb>& ,
	const SmallMatrixView<T,M,N,Sic,Sjc,Cc>& )
    { TMVAssert(FALSE); }
  template <class T, class Tb, int M, int N, int Sia, int Sja, bool Ca, int Sib, int Sjb, bool Cb, int Sic, int Sjc, bool Cc> 
    inline void AddMM_1_m1(
	const GenSmallMatrix<CT,M,N,Sia,Sja,Ca>& , 
	const GenSmallMatrix<Tb,M,N,Sib,Sjb,Cb>& ,
	const SmallMatrixView<T,M,N,Sic,Sjc,Cc>& )
    { TMVAssert(FALSE); }
  template <class T, int M, int N, int Sia, int Sja, bool Ca, int Sib, int Sjb, bool Cb, int Sic, int Sjc, bool Cc> 
    inline void AddMM_1_m1(
	const GenSmallMatrix<T,M,N,Sia,Sja,Ca>& , 
	const GenSmallMatrix<CT,M,N,Sib,Sjb,Cb>& ,
	const SmallMatrixView<T,M,N,Sic,Sjc,Cc>& )
    { TMVAssert(FALSE); }
  template <class T, class Ta, class Tb, int M, int N, int Sia, int Sja, bool Ca, int Sib, int Sjb, bool Cb, int Sic, int Sjc, bool Cc> 
    inline void AddMM_1_x(
	const GenSmallMatrix<Ta,M,N,Sia,Sja,Ca>& , 
	const CT , const GenSmallMatrix<Tb,M,N,Sib,Sjb,Cb>& ,
	const SmallMatrixView<T,M,N,Sic,Sjc,Cc>& )
    { TMVAssert(FALSE); }
  template <class T, class Ta, class Tb, int M, int N, int Sia, int Sja, bool Ca, int Sib, int Sjb, bool Cb, int Sic, int Sjc, bool Cc> 
    inline void AddMM_x_1(
	const CT , const GenSmallMatrix<Ta,M,N,Sia,Sja,Ca>& , 
	const GenSmallMatrix<Tb,M,N,Sib,Sjb,Cb>& ,
	const SmallMatrixView<T,M,N,Sic,Sjc,Cc>& )
    { TMVAssert(FALSE); }
  template <class T, class Ta, class Tb, int M, int N, int Sia, int Sja, bool Ca, int Sib, int Sjb, bool Cb, int Sic, int Sjc, bool Cc> 
    inline void AddMM_x_m1(
	const CT , const GenSmallMatrix<Ta,M,N,Sia,Sja,Ca>& , 
	const GenSmallMatrix<Tb,M,N,Sib,Sjb,Cb>& ,
	const SmallMatrixView<T,M,N,Sic,Sjc,Cc>& )
    { TMVAssert(FALSE); }
  template <class T, class Ta, class Tb, int M, int N, int Sia, int Sja, bool Ca, int Sib, int Sjb, bool Cb, int Sic, int Sjc, bool Cc> 
    inline void AddMM(
	const CT , const GenSmallMatrix<Ta,M,N,Sia,Sja,Ca>& , 
	const CT , const GenSmallMatrix<Tb,M,N,Sib,Sjb,Cb>& ,
	const SmallMatrixView<T,M,N,Sic,Sjc,Cc>& )
    { TMVAssert(FALSE); }
  template <class Ta, class Tb, class Tc, int M, int N, int Sia, int Sja, bool Ca, int Sib, int Sjb, bool Cb> 
    inline void AddMM_1_1(
	const GenSmallMatrix<Ta,M,N,Sia,Sja,Ca>& A, 
	const GenSmallMatrix<Tb,M,N,Sib,Sjb,Cb>& B,
	Tc* Cp, int Sic, int Sjc)
    { Copy(A,Cp,Sic,Sjc); AddMM_1(B,Cp,Sic,Sjc); }
  template <class Ta, class Tb, class Tc, int M, int N, int Sia, int Sja, bool Ca, int Sib, int Sjb, bool Cb> 
    inline void AddMM_1_m1(
	const GenSmallMatrix<Ta,M,N,Sia,Sja,Ca>& A, 
	const GenSmallMatrix<Tb,M,N,Sib,Sjb,Cb>& B,
	Tc* Cp, int Sic, int Sjc)
    { Copy(A,Cp,Sic,Sjc); AddMM_m1(B,Cp,Sic,Sjc); }
  template <class T, class Ta, class Tb, class Tc, int M, int N, int Sia, int Sja, bool Ca, int Sib, int Sjb, bool Cb> 
    inline void AddMM_1_x(
	const GenSmallMatrix<Ta,M,N,Sia,Sja,Ca>& A, 
	const T beta, const GenSmallMatrix<Tb,M,N,Sib,Sjb,Cb>& B,
	Tc* Cp, int Sic, int Sjc)
    { MultXM(beta,B,Cp,Sic,Sjc); AddMM_1(A,Cp,Sic,Sjc); }
  template <class T, class Ta, class Tb, class Tc, int M, int N, int Sia, int Sja, bool Ca, int Sib, int Sjb, bool Cb> 
    inline void AddMM_x_1(
	const T alpha, const GenSmallMatrix<Ta,M,N,Sia,Sja,Ca>& A, 
	const GenSmallMatrix<Tb,M,N,Sib,Sjb,Cb>& B,
	Tc* Cp, int Sic, int Sjc)
    { MultXM(alpha,A,Cp,Sic,Sjc); AddMM_1(B,Cp,Sic,Sjc); }
  template <class T, class Ta, class Tb, class Tc, int M, int N, int Sia, int Sja, bool Ca, int Sib, int Sjb, bool Cb> 
    inline void AddMM_x_m1(
	const T alpha, const GenSmallMatrix<Ta,M,N,Sia,Sja,Ca>& A, 
	const GenSmallMatrix<Tb,M,N,Sib,Sjb,Cb>& B,
	Tc* Cp, int Sic, int Sjc)
    { MultXM(alpha,A,Cp,Sic,Sjc); AddMM_m1(B,Cp,Sic,Sjc); }
  template <class T, class Ta, class Tb, class Tc, int M, int N, int Sia, int Sja, bool Ca, int Sib, int Sjb, bool Cb> 
    inline void AddMM(
	const T alpha, const GenSmallMatrix<Ta,M,N,Sia,Sja,Ca>& A, 
	const T beta, const GenSmallMatrix<Tb,M,N,Sib,Sjb,Cb>& B,
	Tc* Cp, int Sic, int Sjc)
    { MultXM(alpha,A,Cp,Sic,Sjc); AddMM(beta,B,Cp,Sic,Sjc); }
  template <class T, class Tb, int M, int N, int Sia, int Sja, bool Ca, int Sib, int Sjb, bool Cb> 
    inline void AddMM_1_1(
	const GenSmallMatrix<CT,M,N,Sia,Sja,Ca>& , 
	const GenSmallMatrix<Tb,M,N,Sib,Sjb,Cb>& ,
	T* , int , int )
    { TMVAssert(FALSE); }
  template <class T, int M, int N, int Sia, int Sja, bool Ca, int Sib, int Sjb, bool Cb> 
    inline void AddMM_1_1(
	const GenSmallMatrix<T,M,N,Sia,Sja,Ca>& , 
	const GenSmallMatrix<CT,M,N,Sib,Sjb,Cb>& ,
	T* , int , int )
    { TMVAssert(FALSE); }
  template <class T, class Tb, int M, int N, int Sia, int Sja, bool Ca, int Sib, int Sjb, bool Cb> 
    inline void AddMM_1_m1(
	const GenSmallMatrix<CT,M,N,Sia,Sja,Ca>& , 
	const GenSmallMatrix<Tb,M,N,Sib,Sjb,Cb>& ,
	T* , int , int )
    { TMVAssert(FALSE); }
  template <class T, int M, int N, int Sia, int Sja, bool Ca, int Sib, int Sjb, bool Cb> 
    inline void AddMM_1_m1(
	const GenSmallMatrix<T,M,N,Sia,Sja,Ca>& , 
	const GenSmallMatrix<CT,M,N,Sib,Sjb,Cb>& ,
	T* , int , int )
    { TMVAssert(FALSE); }
  template <class T, class Ta, class Tb, int M, int N, int Sia, int Sja, bool Ca, int Sib, int Sjb, bool Cb> 
    inline void AddMM_1_x(
	const GenSmallMatrix<Ta,M,N,Sia,Sja,Ca>& , 
	const CT , const GenSmallMatrix<Tb,M,N,Sib,Sjb,Cb>& ,
	T* , int , int )
    { TMVAssert(FALSE); }
  template <class T, class Ta, class Tb, int M, int N, int Sia, int Sja, bool Ca, int Sib, int Sjb, bool Cb> 
    inline void AddMM_x_1(
	const CT , const GenSmallMatrix<Ta,M,N,Sia,Sja,Ca>& , 
	const GenSmallMatrix<Tb,M,N,Sib,Sjb,Cb>& ,
	T* , int , int )
    { TMVAssert(FALSE); }
  template <class T, class Ta, class Tb, int M, int N, int Sia, int Sja, bool Ca, int Sib, int Sjb, bool Cb> 
    inline void AddMM_x_m1(
	const CT , const GenSmallMatrix<Ta,M,N,Sia,Sja,Ca>& , 
	const GenSmallMatrix<Tb,M,N,Sib,Sjb,Cb>& ,
	T* , int , int )
    { TMVAssert(FALSE); }
  template <class T, class Ta, class Tb, int M, int N, int Sia, int Sja, bool Ca, int Sib, int Sjb, bool Cb> 
    inline void AddMM(
	const CT , const GenSmallMatrix<Ta,M,N,Sia,Sja,Ca>& , 
	const CT , const GenSmallMatrix<Tb,M,N,Sib,Sjb,Cb>& ,
	T* , int , int )
    { TMVAssert(FALSE); }

  // C = alpha * A * B
  template <class Ta, class Tb, class Tc, int M, int N, int K, int Sia, int Sja, bool Ca, int Sib, int Sjb, bool Cb, int Sic, int Sjc, bool Cc> 
    inline void MultMM_1(
	const GenSmallMatrix<Ta,M,K,Sia,Sja,Ca>& A, 
	const GenSmallMatrix<Tb,K,N,Sib,Sjb,Cb>& B,
	const SmallMatrixView<Tc,M,N,Sic,Sjc,Cc>& C)
    {
      if (Sjc == 1) {
	for(int i=0;i<M;++i) 
	  MultMV_1(B.Transpose(),A.row(i),C.row(i));
      } else {
	for(int j=0;j<N;++j) 
	  MultMV_1(A,B.col(j),C.col(j));
      }
    }
  template <class Ta, class Tb, class Tc, int M, int N, int K, int Sia, int Sja, bool Ca, int Sib, int Sjb, bool Cb, int Sic, int Sjc, bool Cc> 
    inline void MultMM_m1(
	const GenSmallMatrix<Ta,M,K,Sia,Sja,Ca>& A, 
	const GenSmallMatrix<Tb,K,N,Sib,Sjb,Cb>& B,
	const SmallMatrixView<Tc,M,N,Sic,Sjc,Cc>& C)
    {
      if (Sjc == 1) {
	for(int i=0;i<M;++i) 
	  MultMV_m1(B.Transpose(),A.row(i),C.row(i));
      } else {
	for(int j=0;j<N;++j) 
	  MultMV_m1(A,B.col(j),C.col(j));
      }
    }
  template <class T, class Ta, class Tb, class Tc, int M, int N, int K, int Sia, int Sja, bool Ca, int Sib, int Sjb, bool Cb, int Sic, int Sjc, bool Cc> 
    inline void MultMM(const T alpha,
	const GenSmallMatrix<Ta,M,K,Sia,Sja,Ca>& A, 
	const GenSmallMatrix<Tb,K,N,Sib,Sjb,Cb>& B,
	const SmallMatrixView<Tc,M,N,Sic,Sjc,Cc>& C)
    {
      if (Sjc == 1) {
	for(int i=0;i<M;++i) 
	  MultMV(alpha,B.Transpose(),A.row(i),C.row(i));
      } else {
	for(int j=0;j<N;++j) 
	  MultMV(alpha,A,B.col(j),C.col(j));
      }
    }
  template <class T, class Tb, int M, int N, int K, int Sia, int Sja, bool Ca, int Sib, int Sjb, bool Cb, int Sic, int Sjc, bool Cc> 
    inline void MultMM_1(
	const GenSmallMatrix<CT,M,K,Sia,Sja,Ca>& , 
	const GenSmallMatrix<Tb,K,N,Sib,Sjb,Cb>& ,
	const SmallMatrixView<T,M,N,Sic,Sjc,Cc>& )
    { TMVAssert(FALSE); }
  template <class T, int M, int N, int K, int Sia, int Sja, bool Ca, int Sib, int Sjb, bool Cb, int Sic, int Sjc, bool Cc> 
    inline void MultMM_1(
	const GenSmallMatrix<T,M,K,Sia,Sja,Ca>& , 
	const GenSmallMatrix<CT,K,N,Sib,Sjb,Cb>& ,
	const SmallMatrixView<T,M,N,Sic,Sjc,Cc>& )
    { TMVAssert(FALSE); }
  template <class T, class Tb, int M, int N, int K, int Sia, int Sja, bool Ca, int Sib, int Sjb, bool Cb, int Sic, int Sjc, bool Cc> 
    inline void MultMM_m1(
	const GenSmallMatrix<CT,M,K,Sia,Sja,Ca>& , 
	const GenSmallMatrix<Tb,K,N,Sib,Sjb,Cb>& ,
	const SmallMatrixView<T,M,N,Sic,Sjc,Cc>& )
    { TMVAssert(FALSE); }
  template <class T, int M, int N, int K, int Sia, int Sja, bool Ca, int Sib, int Sjb, bool Cb, int Sic, int Sjc, bool Cc> 
    inline void MultMM_m1(
	const GenSmallMatrix<T,M,K,Sia,Sja,Ca>& , 
	const GenSmallMatrix<CT,K,N,Sib,Sjb,Cb>& ,
	const SmallMatrixView<T,M,N,Sic,Sjc,Cc>& )
    { TMVAssert(FALSE); }
  template <class T, class Ta, class Tb, int M, int N, int K, int Sia, int Sja, bool Ca, int Sib, int Sjb, bool Cb, int Sic, int Sjc, bool Cc> 
    inline void MultMM(const CT ,
	const GenSmallMatrix<Ta,M,K,Sia,Sja,Ca>& , 
	const GenSmallMatrix<Tb,K,N,Sib,Sjb,Cb>& ,
	const SmallMatrixView<T,M,N,Sic,Sjc,Cc>& )
    { TMVAssert(FALSE); }

  template <class Ta, class Tb, class Tc, int M, int N, int K, int Sia, int Sja, bool Ca, int Sib, int Sjb, bool Cb> 
    inline void MultMM_1(
	const GenSmallMatrix<Ta,M,K,Sia,Sja,Ca>& A, 
	const GenSmallMatrix<Tb,K,N,Sib,Sjb,Cb>& B,
	Tc* Cp, int Sic, int Sjc)
    {
      if (Sjc == 1) {
	for(int i=0;i<M;++i,Cp+=Sic) 
	  MultMV_1(B.Transpose(),A.row(i),SmallVectorViewOf<N>(Cp));
      } else if (Sic == 1) {
	for(int j=0;j<N;++j,Cp+=Sjc) 
	  MultMV_1(A,B.col(j),SmallVectorViewOf<M>(Cp));
      } else {
	for(int j=0;j<N;++j,Cp+=Sjc) 
	  MultMV_1(A,B.col(j),Cp,Sic);
      }
    }
  template <class Ta, class Tb, class Tc, int M, int N, int K, int Sia, int Sja, bool Ca, int Sib, int Sjb, bool Cb> 
    inline void MultMM_m1(
	const GenSmallMatrix<Ta,M,K,Sia,Sja,Ca>& A, 
	const GenSmallMatrix<Tb,K,N,Sib,Sjb,Cb>& B,
	Tc* Cp, int Sic, int Sjc)
    {
      if (Sjc == 1) {
	for(int i=0;i<M;++i,Cp+=Sic) 
	  MultMV_m1(B.Transpose(),A.row(i),SmallVectorViewOf<N>(Cp));
      } else if (Sic == 1) {
	for(int j=0;j<N;++j,Cp+=Sjc) 
	  MultMV_m1(A,B.col(j),SmallVectorViewOf<M>(Cp));
      } else {
	for(int j=0;j<N;++j,Cp+=Sjc) 
	  MultMV_m1(A,B.col(j),Cp,Sic);
      }
    }
  template <class T, class Ta, class Tb, class Tc, int M, int N, int K, int Sia, int Sja, bool Ca, int Sib, int Sjb, bool Cb> 
    inline void MultMM(const T alpha,
	const GenSmallMatrix<Ta,M,K,Sia,Sja,Ca>& A, 
	const GenSmallMatrix<Tb,K,N,Sib,Sjb,Cb>& B,
	Tc* Cp, int Sic, int Sjc)
    {
      if (Sjc == 1) {
	for(int i=0;i<M;++i,Cp+=Sic) 
	  MultMV(alpha,B.Transpose(),A.row(i),SmallVectorViewOf<N>(Cp));
      } else if (Sic == 1) {
	for(int j=0;j<N;++j,Cp+=Sjc) 
	  MultMV(alpha,A,B.col(j),SmallVectorViewOf<M>(Cp));
      } else {
	for(int j=0;j<N;++j,Cp+=Sjc) 
	  MultMV(alpha,A,B.col(j),Cp,Sic);
      }
    }
  template <class T, class Tb, int M, int N, int K, int Sia, int Sja, bool Ca, int Sib, int Sjb, bool Cb> 
    inline void MultMM_1(
	const GenSmallMatrix<CT,M,K,Sia,Sja,Ca>& , 
	const GenSmallMatrix<Tb,K,N,Sib,Sjb,Cb>& , T* , int , int )
    { TMVAssert(FALSE); }
  template <class T, int M, int N, int K, int Sia, int Sja, bool Ca, int Sib, int Sjb, bool Cb> 
    inline void MultMM_1(
	const GenSmallMatrix<T,M,K,Sia,Sja,Ca>& , 
	const GenSmallMatrix<CT,K,N,Sib,Sjb,Cb>& , T* , int , int )
    { TMVAssert(FALSE); }
  template <class T, class Tb, int M, int N, int K, int Sia, int Sja, bool Ca, int Sib, int Sjb, bool Cb> 
    inline void MultMM_m1(
	const GenSmallMatrix<CT,M,K,Sia,Sja,Ca>& , 
	const GenSmallMatrix<Tb,K,N,Sib,Sjb,Cb>& , T* , int , int )
    { TMVAssert(FALSE); }
  template <class T, int M, int N, int K, int Sia, int Sja, bool Ca, int Sib, int Sjb, bool Cb> 
    inline void MultMM_m1(
	const GenSmallMatrix<T,M,K,Sia,Sja,Ca>& , 
	const GenSmallMatrix<CT,K,N,Sib,Sjb,Cb>& , T* , int , int )
    { TMVAssert(FALSE); }
  template <class T, class Ta, class Tb, int M, int N, int K, int Sia, int Sja, bool Ca, int Sib, int Sjb, bool Cb> 
    inline void MultMM(const CT ,
	const GenSmallMatrix<Ta,M,K,Sia,Sja,Ca>& , 
	const GenSmallMatrix<Tb,K,N,Sib,Sjb,Cb>& , T* , int , int )
    { TMVAssert(FALSE); }

  // C += alpha * A * B
  template <class Ta, class Tb, class Tc, int M, int N, int K, int Sia, int Sja, bool Ca, int Sib, int Sjb, bool Cb, int Sic, int Sjc, bool Cc> 
    inline void AddMultMM_1(
	const GenSmallMatrix<Ta,M,K,Sia,Sja,Ca>& A, 
	const GenSmallMatrix<Tb,K,N,Sib,Sjb,Cb>& B,
	const SmallMatrixView<Tc,M,N,Sic,Sjc,Cc>& C)
    {
      if (Sjc == 1) {
	for(int i=0;i<M;++i) 
	  AddMultMV_1(B.Transpose(),A.row(i),C.row(i));
      } else {
	for(int j=0;j<N;++j) 
	  AddMultMV_1(A,B.col(j),C.col(j));
      }
    }
  template <class Ta, class Tb, class Tc, int M, int N, int K, int Sia, int Sja, bool Ca, int Sib, int Sjb, bool Cb, int Sic, int Sjc, bool Cc> 
    inline void AddMultMM_m1(
	const GenSmallMatrix<Ta,M,K,Sia,Sja,Ca>& A, 
	const GenSmallMatrix<Tb,K,N,Sib,Sjb,Cb>& B,
	const SmallMatrixView<Tc,M,N,Sic,Sjc,Cc>& C)
    {
      if (Sjc == 1) {
	for(int i=0;i<M;++i) 
	  AddMultMV_m1(B.Transpose(),A.row(i),C.row(i));
      } else {
	for(int j=0;j<N;++j) 
	  AddMultMV_m1(A,B.col(j),C.col(j));
      }
    }
  template <class T, class Ta, class Tb, class Tc, int M, int N, int K, int Sia, int Sja, bool Ca, int Sib, int Sjb, bool Cb, int Sic, int Sjc, bool Cc> 
    inline void AddMultMM(const T alpha,
	const GenSmallMatrix<Ta,M,K,Sia,Sja,Ca>& A, 
	const GenSmallMatrix<Tb,K,N,Sib,Sjb,Cb>& B,
	const SmallMatrixView<Tc,M,N,Sic,Sjc,Cc>& C)
    {
      if (Sjc == 1) {
	for(int i=0;i<M;++i) 
	  AddMultMV(alpha,B.Transpose(),A.row(i),C.row(i));
      } else {
	for(int j=0;j<N;++j) 
	  AddMultMV(alpha,A,B.col(j),C.col(j));
      }
    }
  template <class T, class Tb, int M, int N, int K, int Sia, int Sja, bool Ca, int Sib, int Sjb, bool Cb, int Sic, int Sjc, bool Cc> 
    inline void AddMultMM_1(
	const GenSmallMatrix<CT,M,K,Sia,Sja,Ca>& , 
	const GenSmallMatrix<Tb,K,N,Sib,Sjb,Cb>& ,
	const SmallMatrixView<T,M,N,Sic,Sjc,Cc>& )
    { TMVAssert(FALSE); }
  template <class T, int M, int N, int K, int Sia, int Sja, bool Ca, int Sib, int Sjb, bool Cb, int Sic, int Sjc, bool Cc> 
    inline void AddMultMM_1(
	const GenSmallMatrix<T,M,K,Sia,Sja,Ca>& , 
	const GenSmallMatrix<CT,K,N,Sib,Sjb,Cb>& ,
	const SmallMatrixView<T,M,N,Sic,Sjc,Cc>& )
    { TMVAssert(FALSE); }
  template <class T, class Tb, int M, int N, int K, int Sia, int Sja, bool Ca, int Sib, int Sjb, bool Cb, int Sic, int Sjc, bool Cc> 
    inline void AddMultMM_m1(
	const GenSmallMatrix<CT,M,K,Sia,Sja,Ca>& , 
	const GenSmallMatrix<Tb,K,N,Sib,Sjb,Cb>& ,
	const SmallMatrixView<T,M,N,Sic,Sjc,Cc>& )
    { TMVAssert(FALSE); }
  template <class T, int M, int N, int K, int Sia, int Sja, bool Ca, int Sib, int Sjb, bool Cb, int Sic, int Sjc, bool Cc> 
    inline void AddMultMM_m1(
	const GenSmallMatrix<T,M,K,Sia,Sja,Ca>& , 
	const GenSmallMatrix<CT,K,N,Sib,Sjb,Cb>& ,
	const SmallMatrixView<T,M,N,Sic,Sjc,Cc>& )
    { TMVAssert(FALSE); }
  template <class T, class Ta, class Tb, int M, int N, int K, int Sia, int Sja, bool Ca, int Sib, int Sjb, bool Cb, int Sic, int Sjc, bool Cc> 
    inline void AddMultMM(const CT ,
	const GenSmallMatrix<Ta,M,K,Sia,Sja,Ca>& , 
	const GenSmallMatrix<Tb,K,N,Sib,Sjb,Cb>& ,
	const SmallMatrixView<T,M,N,Sic,Sjc,Cc>& )
    { TMVAssert(FALSE); }

  // A = alpha * x * yT
  template <class Ta, class Tx, class Ty, int M, int N, int Sia, int Sja, bool Ca, int Sx, bool Cx, int Sy, bool Cy> 
    inline void Rank1Update_1(
	const GenSmallVector<Tx,M,Sx,Cx>& x, 
	const GenSmallVector<Ty,N,Sy,Cy>& y, 
	const SmallMatrixView<Ta,M,N,Sia,Sja,Ca>& A)
    {
      if (Sia == 1) {
	for(int j=0;j<N;++j) A.col(j) = y(j) * x;
      } else {
	for(int i=0;i<M;++i) A.row(i) = x(i) * y;
      }
    }
  template <class Ta, class Tx, class Ty, int M, int N, int Sia, int Sja, bool Ca, int Sx, bool Cx, int Sy, bool Cy> 
    inline void Rank1Update_m1(
	const GenSmallVector<Tx,M,Sx,Cx>& x, 
	const GenSmallVector<Ty,N,Sy,Cy>& y, 
	const SmallMatrixView<Ta,M,N,Sia,Sja,Ca>& A)
    {
      if (Sia == 1) {
	for(int j=0;j<N;++j) A.col(j) = (-y(j)) * x;
      } else {
	for(int i=0;i<M;++i) A.row(i) = (-x(i)) * y;
      }
    }
  template <class T, class Ta, class Tx, class Ty, int M, int N, int Sia, int Sja, bool Ca, int Sx, bool Cx, int Sy, bool Cy> 
    inline void Rank1Update(const T alpha,
	const GenSmallVector<Tx,M,Sx,Cx>& x, 
	const GenSmallVector<Ty,N,Sy,Cy>& y, 
	const SmallMatrixView<Ta,M,N,Sia,Sja,Ca>& A)
    {
      if (Sia == 1) {
	for(int j=0;j<N;++j) A.col(j) = (alpha*y(j)) * x;
      } else {
	for(int i=0;i<M;++i) A.row(i) = (alpha*x(i)) * y;
      }
    }
  template <class T, class Ty, int M, int N, int Sia, int Sja, bool Ca, int Sx, bool Cx, int Sy, bool Cy> 
    inline void Rank1Update_1(
	const GenSmallVector<CT,M,Sx,Cx>& , 
	const GenSmallVector<Ty,N,Sy,Cy>& , 
	const SmallMatrixView<T,M,N,Sia,Sja,Ca>& )
    { TMVAssert(FALSE); }
  template <class T, int M, int N, int Sia, int Sja, bool Ca, int Sx, bool Cx, int Sy, bool Cy> 
    inline void Rank1Update_1(
	const GenSmallVector<T,M,Sx,Cx>& , 
	const GenSmallVector<CT,N,Sy,Cy>& , 
	const SmallMatrixView<T,M,N,Sia,Sja,Ca>& )
    { TMVAssert(FALSE); }
  template <class T, class Ty, int M, int N, int Sia, int Sja, bool Ca, int Sx, bool Cx, int Sy, bool Cy> 
    inline void Rank1Update_m1(
	const GenSmallVector<CT,M,Sx,Cx>& , 
	const GenSmallVector<Ty,N,Sy,Cy>& , 
	const SmallMatrixView<T,M,N,Sia,Sja,Ca>& )
    { TMVAssert(FALSE); }
  template <class T, int M, int N, int Sia, int Sja, bool Ca, int Sx, bool Cx, int Sy, bool Cy> 
    inline void Rank1Update_m1(
	const GenSmallVector<T,M,Sx,Cx>& , 
	const GenSmallVector<CT,N,Sy,Cy>& , 
	const SmallMatrixView<T,M,N,Sia,Sja,Ca>& )
    { TMVAssert(FALSE); }
  template <class T, class Tx, class Ty, int M, int N, int Sia, int Sja, bool Ca, int Sx, bool Cx, int Sy, bool Cy> 
    inline void Rank1Update(const CT ,
	const GenSmallVector<Tx,M,Sx,Cx>& , 
	const GenSmallVector<Ty,N,Sy,Cy>& , 
	const SmallMatrixView<T,M,N,Sia,Sja,Ca>& )
    { TMVAssert(FALSE); }
  template <class Ta, class Tx, class Ty, int M, int N, int Sx, bool Cx, int Sy, bool Cy> 
    inline void Rank1Update_1(
	const GenSmallVector<Tx,M,Sx,Cx>& x, 
	const GenSmallVector<Ty,N,Sy,Cy>& y, 
	Ta* Ap, int Sia, int Sja)
    {
      if (Sia == 1) {
	for(int j=0;j<N;++j,Ap+=Sja) 
	  SmallVectorViewOf<M>(Ap) = y(j) * x;
      } else if (Sja == 1) {
	for(int i=0;i<M;++i,Ap+=Sia) 
	  SmallVectorViewOf<N>(Ap) = x(i) * y;
      } else {
	for(int i=0;i<M;++i,Ap+=Sia) MultXV(x(i),y,Ap,Sja);
      }
    }
  template <class Ta, class Tx, class Ty, int M, int N, int Sx, bool Cx, int Sy, bool Cy> 
    inline void Rank1Update_m1(
	const GenSmallVector<Tx,M,Sx,Cx>& x, 
	const GenSmallVector<Ty,N,Sy,Cy>& y, 
	Ta* Ap, int Sia, int Sja)
    {
      if (Sia == 1) {
	for(int j=0;j<N;++j,Ap+=Sja) 
	  SmallVectorViewOf<M>(Ap) = (-y(j)) * x;
      } else if (Sja == 1) {
	for(int i=0;i<M;++i,Ap+=Sia) 
	  SmallVectorViewOf<N>(Ap) = (-x(i)) * y;
      } else {
	for(int i=0;i<M;++i,Ap+=Sia) MultXV(-x(i),y,Ap,Sja);
      }
    }
  template <class T, class Ta, class Tx, class Ty, int M, int N, int Sx, bool Cx, int Sy, bool Cy> 
    inline void Rank1Update(const T alpha,
	const GenSmallVector<Tx,M,Sx,Cx>& x, 
	const GenSmallVector<Ty,N,Sy,Cy>& y, 
	Ta* Ap, int Sia, int Sja)
    {
      if (Sia == 1) {
	for(int j=0;j<N;++j,Ap+=Sja) 
	  SmallVectorViewOf<M>(Ap) = (alpha*y(j)) * x;
      } else if (Sja == 1) {
	for(int i=0;i<M;++i,Ap+=Sia) 
	  SmallVectorViewOf<N>(Ap) = (alpha*x(i)) * y;
      } else {
	for(int i=0;i<M;++i,Ap+=Sia) MultXV(alpha*x(i),y,Ap,Sja);
      }
    }
  template <class T, class Ty, int M, int N, int Sx, bool Cx, int Sy, bool Cy> 
    inline void Rank1Update_1(const GenSmallVector<CT,M,Sx,Cx>& , 
	const GenSmallVector<Ty,N,Sy,Cy>& , T* , int , int )
    { TMVAssert(FALSE); }
  template <class T, int M, int N, int Sx, bool Cx, int Sy, bool Cy> 
    inline void Rank1Update_1(const GenSmallVector<T,M,Sx,Cx>& , 
	const GenSmallVector<CT,N,Sy,Cy>& , T* , int , int )
    { TMVAssert(FALSE); }
  template <class T, class Ty, int M, int N, int Sx, bool Cx, int Sy, bool Cy> 
    inline void Rank1Update_m1(const GenSmallVector<CT,M,Sx,Cx>& , 
	const GenSmallVector<Ty,N,Sy,Cy>& , T* , int , int )
    { TMVAssert(FALSE); }
  template <class T, int M, int N, int Sx, bool Cx, int Sy, bool Cy> 
    inline void Rank1Update_m1(const GenSmallVector<T,M,Sx,Cx>& , 
	const GenSmallVector<CT,N,Sy,Cy>& , T* , int , int )
    { TMVAssert(FALSE); }
  template <class T, class Tx, class Ty, int M, int N, int Sx, bool Cx, int Sy, bool Cy> 
    inline void Rank1Update(const CT , const GenSmallVector<Tx,M,Sx,Cx>& , 
	const GenSmallVector<Ty,N,Sy,Cy>& , T* , int , int )
    { TMVAssert(FALSE); }

  // A += alpha * x * yT
  template <class Ta, class Tx, class Ty, int M, int N, int Sia, int Sja, bool Ca, int Sx, bool Cx, int Sy, bool Cy> 
    inline void AddRank1Update_1(
	const GenSmallVector<Tx,M,Sx,Cx>& x, 
	const GenSmallVector<Ty,N,Sy,Cy>& y, 
	const SmallMatrixView<Ta,M,N,Sia,Sja,Ca>& A)
    {
      if (Sia == 1) {
	for(int j=0;j<N;++j) A.col(j) += y(j) * x;
      } else {
	for(int i=0;i<M;++i) A.row(i) += x(i) * y;
      }
    }
  template <class Ta, class Tx, class Ty, int M, int N, int Sia, int Sja, bool Ca, int Sx, bool Cx, int Sy, bool Cy> 
    inline void AddRank1Update_m1(
	const GenSmallVector<Tx,M,Sx,Cx>& x, 
	const GenSmallVector<Ty,N,Sy,Cy>& y, 
	const SmallMatrixView<Ta,M,N,Sia,Sja,Ca>& A)
    {
      if (Sia == 1) {
	for(int j=0;j<N;++j) A.col(j) -= y(j) * x;
      } else {
	for(int i=0;i<M;++i) A.row(i) -= x(i) * y;
      }
    }
  template <class T, class Ta, class Tx, class Ty, int M, int N, int Sia, int Sja, bool Ca, int Sx, bool Cx, int Sy, bool Cy> 
    inline void AddRank1Update(const T alpha,
	const GenSmallVector<Tx,M,Sx,Cx>& x, 
	const GenSmallVector<Ty,N,Sy,Cy>& y, 
	const SmallMatrixView<Ta,M,N,Sia,Sja,Ca>& A)
    {
      if (Sia == 1) {
	for(int j=0;j<N;++j) A.col(j) += (alpha*y(j)) * x;
      } else {
	for(int i=0;i<M;++i) A.row(i) += (alpha*x(i)) * y;
      }
    }
  template <class T, class Ty, int M, int N, int Sia, int Sja, bool Ca, int Sx, bool Cx, int Sy, bool Cy> 
    inline void AddRank1Update_1(
	const GenSmallVector<CT,M,Sx,Cx>& , 
	const GenSmallVector<Ty,N,Sy,Cy>& , 
	const SmallMatrixView<T,M,N,Sia,Sja,Ca>& )
    { TMVAssert(FALSE); }
  template <class T, int M, int N, int Sia, int Sja, bool Ca, int Sx, bool Cx, int Sy, bool Cy> 
    inline void AddRank1Update_1(
	const GenSmallVector<T,M,Sx,Cx>& , 
	const GenSmallVector<CT,N,Sy,Cy>& , 
	const SmallMatrixView<T,M,N,Sia,Sja,Ca>& )
    { TMVAssert(FALSE); }
  template <class T, class Ty, int M, int N, int Sia, int Sja, bool Ca, int Sx, bool Cx, int Sy, bool Cy> 
    inline void AddRank1Update_m1(
	const GenSmallVector<CT,M,Sx,Cx>& , 
	const GenSmallVector<Ty,N,Sy,Cy>& , 
	const SmallMatrixView<T,M,N,Sia,Sja,Ca>& )
    { TMVAssert(FALSE); }
  template <class T, int M, int N, int Sia, int Sja, bool Ca, int Sx, bool Cx, int Sy, bool Cy> 
    inline void AddRank1Update_m1(
	const GenSmallVector<T,M,Sx,Cx>& , 
	const GenSmallVector<CT,N,Sy,Cy>& , 
	const SmallMatrixView<T,M,N,Sia,Sja,Ca>& )
    { TMVAssert(FALSE); }
  template <class T, class Tx, class Ty, int M, int N, int Sia, int Sja, bool Ca, int Sx, bool Cx, int Sy, bool Cy> 
    inline void AddRank1Update(const CT ,
	const GenSmallVector<Tx,M,Sx,Cx>& , 
	const GenSmallVector<Ty,N,Sy,Cy>& , 
	const SmallMatrixView<T,M,N,Sia,Sja,Ca>& )
    { TMVAssert(FALSE); }
  template <class Ta, class Tx, class Ty, int M, int N, int Sx, bool Cx, int Sy, bool Cy> 
    inline void AddRank1Update_1(
	const GenSmallVector<Tx,M,Sx,Cx>& x, 
	const GenSmallVector<Ty,N,Sy,Cy>& y, 
	Ta* Ap, int Sia, int Sja)
    {
      if (Sia == 1) {
	for(int j=0;j<N;++j,Ap+=Sja) 
	  SmallVectorViewOf<M>(Ap) += y(j) * x;
      } else if (Sja == 1) {
	for(int i=0;i<M;++i,Ap+=Sia) 
	  SmallVectorViewOf<N>(Ap) += x(i) * y;
      } else {
	for(int i=0;i<M;++i,Ap+=Sia) AddVV(x(i),y,Ap,Sja);
      }
    }
  template <class Ta, class Tx, class Ty, int M, int N, int Sx, bool Cx, int Sy, bool Cy> 
    inline void AddRank1Update_m1(
	const GenSmallVector<Tx,M,Sx,Cx>& x, 
	const GenSmallVector<Ty,N,Sy,Cy>& y, 
	Ta* Ap, int Sia, int Sja)
    {
      if (Sia == 1) {
	for(int j=0;j<N;++j,Ap+=Sja) 
	  SmallVectorViewOf<M>(Ap) -= y(j) * x;
      } else if (Sja == 1) {
	for(int i=0;i<M;++i,Ap+=Sia) 
	  SmallVectorViewOf<N>(Ap) -= x(i) * y;
      } else {
	for(int i=0;i<M;++i,Ap+=Sia) AddVV(-x(i),y,Ap,Sja);
      }
    }
  template <class T, class Ta, class Tx, class Ty, int M, int N, int Sx, bool Cx, int Sy, bool Cy> 
    inline void AddRank1Update(const T alpha,
	const GenSmallVector<Tx,M,Sx,Cx>& x, 
	const GenSmallVector<Ty,N,Sy,Cy>& y, 
	Ta* Ap, int Sia, int Sja)
    {
      if (N < M) {
	if (Sia == 1) {
	  for(int j=0;j<N;++j,Ap+=Sja) 
	    SmallVectorViewOf<M>(Ap) += (alpha*y(j)) * x;
	} else {
	  for(int j=0;j<N;++j,Ap+=Sja) AddVV(alpha*y(j),x,Ap,Sia);
	}
      } else if (M < N) {
	if (Sja == 1) {
	  for(int i=0;i<M;++i,Ap+=Sia) 
	    SmallVectorViewOf<N>(Ap) += (alpha*x(i)) * y;
	} else {
	  for(int i=0;i<M;++i,Ap+=Sia) AddVV(alpha*x(i),y,Ap,Sja);
	}
      } else {
	if (Sia == 1) {
	  for(int j=0;j<N;++j,Ap+=Sja) 
	    SmallVectorViewOf<M>(Ap) += (alpha*y(j)) * x;
	} else if (Sja == 1) {
	  for(int i=0;i<M;++i,Ap+=Sia) 
	    SmallVectorViewOf<N>(Ap) += (alpha*x(i)) * y;
	} else {
	  for(int i=0;i<M;++i,Ap+=Sia) AddVV(alpha*x(i),y,Ap,Sja);
	}
      }
    }
  template <class T, class Ty, int M, int N, int Sx, bool Cx, int Sy, bool Cy> 
    inline void AddRank1Update_1(const GenSmallVector<CT,M,Sx,Cx>& , 
	const GenSmallVector<Ty,N,Sy,Cy>& , T* , int , int )
    { TMVAssert(FALSE); }
  template <class T, int M, int N, int Sx, bool Cx, int Sy, bool Cy> 
    inline void AddRank1Update_1(const GenSmallVector<T,M,Sx,Cx>& , 
	const GenSmallVector<CT,N,Sy,Cy>& , T* , int , int )
    { TMVAssert(FALSE); }
  template <class T, class Ty, int M, int N, int Sx, bool Cx, int Sy, bool Cy> 
    inline void AddRank1Update_m1(const GenSmallVector<CT,M,Sx,Cx>& , 
	const GenSmallVector<Ty,N,Sy,Cy>& , T* , int , int )
    { TMVAssert(FALSE); }
  template <class T, int M, int N, int Sx, bool Cx, int Sy, bool Cy> 
    inline void AddRank1Update_m1(const GenSmallVector<T,M,Sx,Cx>& , 
	const GenSmallVector<CT,N,Sy,Cy>& , T* , int , int )
    { TMVAssert(FALSE); }
  template <class T, class Tx, class Ty, int M, int N, int Sx, bool Cx, int Sy, bool Cy> 
    inline void AddRank1Update(const CT , const GenSmallVector<Tx,M,Sx,Cx>& , 
	const GenSmallVector<Ty,N,Sy,Cy>& , T* , int , int )
    { TMVAssert(FALSE); }

  template <class T, class Ta, class Tb, int M, int N, int Sia, int Sja, bool Ca, int Sib, int Sjb, bool Cb> 
    inline void ElementProd(
	const T alpha, const GenSmallMatrix<Ta,M,N,Sia,Sja,Ca>& A, 
	const SmallMatrixView<Tb,M,N,Sib,Sjb,Cb>& B)
    {
      if (Sia == Sib && Sja == Sjb && A.CanLinearize() && B.CanLinearize()) {
	ElementProd(alpha,A.ConstLinearView(),B.LinearView());
      } else if (Sib == 1) {
	for(int j=0;j<N;++j) ElementProd(alpha,A.col(j),B.col(j));
      } else {
	for(int i=0;i<M;++i) ElementProd(alpha,A.row(i),B.row(i));
      }
    }
  template <class T, class Ta, int M, int N, int Sia, int Sja, bool Ca, int Sib, int Sjb, bool Cb> 
    inline void ElementProd(
	const CT , const GenSmallMatrix<Ta,M,N,Sia,Sja,Ca>& , 
	const SmallMatrixView<T,M,N,Sib,Sjb,Cb>& )
    { TMVAssert(FALSE); }

  template <class T, class Ta, class Tb, class Tc, int M, int N, int Sia, int Sja, bool Ca, int Sib, int Sjb, bool Cb, int Sic, int Sjc, bool Cc> 
    inline void AddElementProd(
	const T alpha, const GenSmallMatrix<Ta,M,N,Sia,Sja,Ca>& A,
	const GenSmallMatrix<Tb,M,N,Sib,Sjb,Cb>& B,
	const SmallMatrixView<Tc,M,N,Sic,Sjc,Cc>& C)
    {
      if (Sia == Sib && Sia == Sic && Sja == Sjb && Sja == Sjc && 
	  A.CanLinearize() && B.CanLinearize() && C.CanLinearize()) {
	AddElementProd(alpha,A.ConstLinearView(),B.ConstLinearView(),
	    C.LinearView());
      } else if (Sic == 1) {
	for(int j=0;j<N;++j) 
	  AddElementProd(alpha,A.col(j),B.col(j),C.col(j));
      } else {
	for(int i=0;i<M;++i) 
	  AddElementProd(alpha,A.row(i),B.row(i),C.row(i));
      }
    }
  template <class T, class Ta, class Tb, int M, int N, int Sia, int Sja, bool Ca, int Sib, int Sjb, bool Cb, int Sic, int Sjc, bool Cc> 
    inline void AddElementProd(
	const CT , const GenSmallMatrix<Ta,M,N,Sia,Sja,Ca>& ,
	const GenSmallMatrix<Tb,M,N,Sib,Sjb,Cb>& ,
	const SmallMatrixView<T,M,N,Sic,Sjc,Cc>& )
    { TMVAssert(FALSE); }

} // namespace tmv

#undef CT

#endif
