///////////////////////////////////////////////////////////////////////////////
// vim:et:ts=2:sw=2:ci:cino=f0,g0,t0,+0:
//                                                                           //
// The Template Matrix/Vector Library for C++ was created by Mike Jarvis     //
// Copyright (C) 1998 - 2009                                                 //
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

#ifndef TMV_Permute_H
#define TMV_Permute_H

#include "TMV_BaseMatrix_Rec.h"

namespace tmv {

  const int TMV_PERM_BLOCKSIZE = 32;

  //
  // PermuteRows
  //

  // TODO: This isn't optimized for SmallMatrices.  Also, might benefit
  // from multiple algorithms.  e.g. is the block thing better for both
  // rowmajor and colmajor?

  template <class M>
  inline void InlinePermuteRows(BaseMatrix_Rec_Mutable<M>& m, 
      const int*const p, int i1, int i2)
  {
    // This idea of doing the permutations a block at a time
    // is cribbed from the LAPack code.  It does speed things up 
    // quite a bit for large matrices.  On my machine where BLOCKSIZE=64
    // is optimal for most routines, blocks of 32 were optimal here.
    const int N = m.rowsize();
    const int Nx = N/TMV_PERM_BLOCKSIZE*TMV_PERM_BLOCKSIZE;
    if (Nx != 0) {
      for(int j=0;j<Nx;) {
        int j2 = j+TMV_PERM_BLOCKSIZE;
        const int* pi = p+i1;
        for(int i=i1;i<i2;++i,++pi) {
          TMVAssert(*pi < int(m.colsize()));
          m.CCols(j,j2).CSwapRows(i,*pi);
        }
        j = j2;
      }
    }
    if (Nx != N) {
      const int* pi = p+i1;
      for(int i=i1;i<i2;++i,++pi) {
        TMVAssert(*pi < int(m.colsize()));
        m.CCols(Nx,N).CSwapRows(i,*pi);
      }
    }
  }

  // Defined in TMV_Matrix.cpp
  template <class T>
  void InstPermuteRows(MatrixView<T> m, 
      const int*const p, int i1, int i2);

  template <class T>
  inline void InstPermuteRows(MatrixView<T,UNKNOWN,UNKNOWN,true> m, 
      const int*const p, int i1, int i2)
  { InstPermuteRows(m.Conjugate(),p,i1,i2); }

  template <bool inst, class M>
  struct CallPermuteRows // inst = true
  {
    static inline void call(M& m, const int*const p, int i1, int i2) 
    { InstPermuteRows(m.XView(),p,i1,i2); } 
  };
  template <class M>
  struct CallPermuteRows<false,M> // inst = false
  {
    static inline void call(M& m, const int*const p, int i1, int i2) 
    { InlinePermuteRows(m,p,i1,i2); } 
  };

  template <class M>
  inline void PermuteRows(BaseMatrix_Rec_Mutable<M>& m,
      const int*const p, int i1, int i2)
  {
    typedef typename M::value_type T;
    enum { inst = (
        Traits<T>::isinst &&
        (M::mrowmajor || M::mcolmajor) &&
        M::mcolsize == UNKNOWN &&
        M::mrowsize == UNKNOWN) };
    CallPermuteRows<inst,M>::call(m.mat(),p,i1,i2);
  }

  template <class M>
  inline void InlineReversePermuteRows(BaseMatrix_Rec_Mutable<M>& m, 
      const int*const p, int i1, int i2)
  {
    const int N = m.rowsize();
    const int Nx = N/TMV_PERM_BLOCKSIZE*TMV_PERM_BLOCKSIZE;
    if (Nx != 0) {
      for(int j=0;j<Nx;) {
        int j2 = j+TMV_PERM_BLOCKSIZE;
        const int* pi = p+i2;
        for(int i=i2;i>i1;) {
          --i; --pi;
          TMVAssert(*pi < int(m.colsize()));
          m.CCols(j,j2).CSwapRows(i,*pi);
        }
        j = j2;
      }
    }
    if (Nx != N) {
      const int* pi = p+i2;
      for(int i=i2;i>i1;) {
        --i; --pi;
        TMVAssert(*pi < int(m.colsize()));
        m.CCols(Nx,N).CSwapRows(i,*pi);
      }
    }
  }

  // Defined in TMV_Matrix.cpp
  template <class T>
  void InstReversePermuteRows(MatrixView<T> m, 
      const int*const p, int i1, int i2);

  template <class T>
  inline void InstReversePermuteRows(MatrixView<T,UNKNOWN,UNKNOWN,true> m, 
      const int*const p, int i1, int i2)
  { InstReversePermuteRows(m.Conjugate(),p,i1,i2); }

  template <bool inst, class M>
  struct CallReversePermuteRows // inst = true
  {
    static inline void call(M& m, const int*const p, int i1, int i2) 
    { InstReversePermuteRows(m.XView(),p,i1,i2); } 
  };
  template <class M>
  struct CallReversePermuteRows<false,M> // inst = false
  {
    static inline void call(M& m, const int*const p, int i1, int i2) 
    { InlineReversePermuteRows(m,p,i1,i2); } 
  };

  template <class M>
  inline void ReversePermuteRows(BaseMatrix_Rec_Mutable<M>& m,
      const int*const p, int i1, int i2)
  {
    typedef typename M::value_type T;
    enum { inst = (
        Traits<T>::isinst &&
        (M::mrowmajor || M::mcolmajor) &&
        M::mcolsize == UNKNOWN &&
        M::mrowsize == UNKNOWN) };
    CallReversePermuteRows<inst,M>::call(m.mat(),p,i1,i2);
  }

} // namespace tmv

#endif
