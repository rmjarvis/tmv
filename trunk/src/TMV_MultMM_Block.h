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

#ifndef TMV_MultMM_Block_H
#define TMV_MultMM_Block_H

#include "tmv/TMV_MultMM.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_SmallMatrix.h"
#include "tmv/TMV_MatrixIO.h"
#include "tmv/TMV_AddMM.h"

namespace tmv {

  const int XX = UNKNOWN;

#if defined(_OPENMP) && TMV_OPT >= 2

  // The block calculation is best done recursively to avoid cache misses 
  // as much as possible, but openmp doesn't do recursion very well,
  // so I just figure out the recursive order first and do a 
  // simple loop that openmp can divvy up appropriately.

  static void MakeTaskList(int M,int N, int i, int j,
      std::vector<int>& ilist, std::vector<int>& jlist, int& index)
  {
    if (M > N) {
      int M1 = M/2;
      MakeTaskList(M1,N,i,j,ilist,jlist,index);
      MakeTaskList(M-M1,N,i+M1,j,ilist,jlist,index);
    } else if (N > 1) {
      int N1 = N/2;
      MakeTaskList(M,N1,i,j,ilist,jlist,index);
      MakeTaskList(M,N-N1,i,j+N1,ilist,jlist,index);
    } else {
      ilist[index] = i;
      jlist[index] = j;
      index++;
    }
  }

  template <bool add, int ix, class T, class MAT1, class MAT2, class MAT3>
  static void DoMultMM(
      const Scaling<ix,T>& x, const MAT1& m1, const MAT2& m2, MAT3& m3)
  {
    TMVAssert(m1.colsize() == m3.colsize());
    TMVAssert(m1.rowsize() == m2.colsize());
    TMVAssert(m2.rowsize() == m3.rowsize());
    TMVAssert(m3.colsize() > 0);
    TMVAssert(m3.rowsize() > 0);
    TMVAssert(m1.rowsize() > 0);

    typedef typename MAT1::value_type T1;
    typedef typename MAT2::value_type T2;
    typedef typename MAT3::value_type T3;

    const int Si1 = MAT1::mstepi;
    const int Sj1 = MAT1::mstepj;
    const bool C1 = MAT1::mconj;

    const int Si2 = MAT2::mstepi;
    const int Sj2 = MAT2::mstepj;
    const bool C2 = MAT2::mconj;

    typedef RealType(T) RT;
    const Scaling<1,RT> one;

    // Notation:
    // N = full rowsize
    // NB = size of blocks in that direction.
    // Nb = number of full blocks
    // N1 = total size of all full blocks
    // N2 = remaining bit that doesn't fit into blocks
    // (likewize for M,K)

    const int M = m3.colsize();
    const int N = m3.rowsize();
    const int K = m1.rowsize();

    // Block sizes:
    //      MB x NB x KB
    // RRC:  4 x  4 x  4
    // RCC:  4 x  4 x 64
    // CRC: 16 x 16 x  4
    // CCC: 64 x  4 x  4
    const int MB = (
        MAT1::mrowmajor ? 
        ( MAT2::mrowmajor ? 4 : 4 ) : 
        ( MAT2::mrowmajor ? 16 : 64 ) );
    const int NB = (
        MAT1::mrowmajor ? 
        ( MAT2::mrowmajor ? 4 : 4 ) : 
        ( MAT2::mrowmajor ? 16 : 4 ) );
    const int KB = (
        MAT1::mrowmajor ? 
        ( MAT2::mrowmajor ? 4 : 64 ) : 
        ( MAT2::mrowmajor ? 4 : 4 ) );

    const int Mb = M/MB;
    const int Nb = N/NB;
    const int Kb = K/KB;

    const int M1 = Mb*MB;
    const int N1 = Nb*NB;
    const int K1 = Kb*KB;

    const int M2 = M%MB;
    const int N2 = N%NB;
    const int K2 = K%KB;

    if (Mb || Nb) {
      const int MbNb = Mb*Nb;
      std::vector<int> ilist(MbNb);
      std::vector<int> jlist(MbNb);
      if (MbNb) {
        int listindex=0;
        MakeTaskList(Mb,Nb,0,0,ilist,jlist,listindex);
      }
#ifdef _OPENMP
#pragma omp parallel
#endif
      {
        SmallMatrix<T3,MB,NB,ColMajor> mtemp;

        // Full MBxNB blocks in m3
        if (Mb && Nb) {
#ifdef _OPENMP
#pragma omp for nowait
#endif
          for(int ij=0;ij<int(MbNb);++ij) {
            const int i=ilist[ij];
            const int j=jlist[ij];
            const int ii = i*MB;
            const int jj = j*NB;

            if (Kb) {
              // Use SmallMatrixView explicitly, so we get the known size
              ConstSmallMatrixView<T1,MB,KB,Si1,Sj1,C1> m1a = (
                  m1.SubMatrix(ii,ii+MB,0,KB));
              ConstSmallMatrixView<T2,KB,NB,Si2,Sj2,C2> m2a = (
                  m2.SubMatrix(0,KB,jj,jj+NB));
              InlineMultMM<false>(one, m1a, m2a, mtemp);

              for (int k=1;k<Kb;++k) {
                const int kk = k*KB;
                ConstSmallMatrixView<T1,MB,KB,Si1,Sj1,C1> m1b = (
                    m1.SubMatrix(ii,ii+MB,kk,kk+KB));
                ConstSmallMatrixView<T2,KB,NB,Si2,Sj2,C2> m2b = (
                    m2.SubMatrix(kk,kk+KB,jj,jj+NB));
                InlineMultMM<true>(one, m1b, m2b, mtemp);
              }
              if (K2) {
                ConstSmallMatrixView<T1,MB,XX,Si1,Sj1,C1> m1c = (
                    m1.SubMatrix(ii,ii+MB,K1,K));
                ConstSmallMatrixView<T2,XX,NB,Si2,Sj2,C2> m2c = (
                    m2.SubMatrix(K1,K,jj,jj+NB));
                InlineMultMM<true>(one, m1c, m2c, mtemp);
              }
            } else {
              ConstSmallMatrixView<T1,MB,XX,Si1,Sj1,C1> m1a = (
                  m1.Rows(ii,ii+MB));
              ConstSmallMatrixView<T2,XX,NB,Si2,Sj2,C2> m2a = (
                  m2.Cols(jj,jj+NB));
              InlineMultMM<false>(one, m1a, m2a, mtemp);
            }
            if (add) m3.SubMatrix(ii,ii+MB,jj,jj+NB) += x * mtemp;
            else m3.SubMatrix(ii,ii+MB,jj,jj+NB) = x * mtemp;
          }
        } // Mb && Nb

        // MB x N2 partial blocks:
        if (Mb && N2) {
#ifdef _OPENMP
#pragma omp for nowait
#endif
          for(int i=0;i<int(Mb);++i) {
            const int ii = i*MB;

            SmallMatrixView<T3,MB,XX,1,MB> m3x = mtemp.Cols(0,N2);

            if (Kb) {
              ConstSmallMatrixView<T1,MB,KB,Si1,Sj1,C1> m1a = (
                  m1.SubMatrix(ii,ii+MB,0,KB));
              ConstSmallMatrixView<T2,KB,XX,Si2,Sj2,C2> m2a = (
                  m2.SubMatrix(0,KB,N1,N));
              InlineMultMM<false>(one, m1a, m2a, m3x);

              for (int k=1;k<Kb;++k) {
                const int kk = k*KB;
                ConstSmallMatrixView<T1,MB,KB,Si1,Sj1,C1> m1b = (
                    m1.SubMatrix(ii,ii+MB,kk,kk+KB));
                ConstSmallMatrixView<T2,KB,XX,Si2,Sj2,C2> m2b = (
                    m2.SubMatrix(kk,kk+KB,N1,N));
                InlineMultMM<true>(one, m1b, m2b, m3x);
              }
              if (K2) {
                ConstSmallMatrixView<T1,MB,XX,Si1,Sj1,C1> m1c = (
                    m1.SubMatrix(ii,ii+MB,K1,K));
                InlineMultMM<true>(one, m1c, m2.SubMatrix(K1,K,N1,N), m3x);
              }
            } else {
              ConstSmallMatrixView<T1,MB,XX,Si1,Sj1,C1> m1a = (
                  m1.Rows(ii,ii+MB));
              InlineMultMM<false>(one, m1a, m2.Cols(N1,N), m3x);
            }
            if (add) m3.SubMatrix(ii,ii+MB,N1,N) += x * m3x;
            else m3.SubMatrix(ii,ii+MB,N1,N) = x * m3x;
          }
        } // Mb && N2
        // M2 x NB partial blocks:
        if (M2 && Nb) {
#ifdef _OPENMP
#pragma omp for nowait
#endif
          for(int j=0;j<int(Nb);++j) {
            const int jj = j*NB;

            SmallMatrixView<T3,XX,NB,1,MB> m3x = mtemp.Rows(0,M2);

            if (Kb) {
              ConstSmallMatrixView<T1,XX,KB,Si1,Sj1,C1> m1a = (
                  m1.SubMatrix(M1,M,0,KB));
              ConstSmallMatrixView<T2,KB,NB,Si2,Sj2,C2> m2a = (
                  m2.SubMatrix(0,KB,jj,jj+NB));
              InlineMultMM<false>(one, m1a, m2a, m3x);

              for (int k=1;k<Kb;++k) {
                const int kk = k*KB;
                ConstSmallMatrixView<T1,XX,KB,Si1,Sj1,C1> m1b = (
                    m1.SubMatrix(M1,M,kk,kk+KB));
                ConstSmallMatrixView<T2,KB,NB,Si2,Sj2,C2> m2b = (
                    m2.SubMatrix(kk,kk+KB,jj,jj+NB));
                InlineMultMM<true>(one, m1b, m2b, m3x);
              }
              if (K2) {
                ConstSmallMatrixView<T2,XX,NB,Si2,Sj2,C2> m2c = (
                    m2.SubMatrix(K1,K,jj,jj+NB));
                InlineMultMM<true>(one, m1.SubMatrix(M1,M,K1,K), m2c, m3x);
              }
            } else {
              ConstSmallMatrixView<T2,XX,NB,Si2,Sj2,C2> m2a = (
                  m2.Cols(jj,jj+NB));
              InlineMultMM<false>(one, m1.Rows(M1,M), m2a, m3x);
            }
            if (add) m3.SubMatrix(M1,M,jj,jj+NB) += x * m3x;
            else m3.SubMatrix(M1,M,jj,jj+NB) = x * m3x;
          }
        } // M2 && Nb

        // Final M2 x N2 partial block
        if (M2 && N2) {
#ifdef _OPENMP
#pragma omp single
#endif
          {
            MatrixView<T3,1,MB> m3x = mtemp.SubMatrix(0,M2,0,N2);

            if (Kb) {
              ConstSmallMatrixView<T1,XX,KB,Si1,Sj1,C1> m1a = (
                  m1.SubMatrix(M1,M,0,KB));
              ConstSmallMatrixView<T2,KB,XX,Si2,Sj2,C2> m2a = (
                  m2.SubMatrix(0,KB,N1,N));
              InlineMultMM<false>(one, m1a, m2a, m3x);

              for (int k=1;k<Kb;++k) {
                const int kk = k*KB;
                ConstSmallMatrixView<T1,XX,KB,Si1,Sj1,C1> m1b = (
                    m1.SubMatrix(M1,M,kk,kk+KB));
                ConstSmallMatrixView<T2,KB,XX,Si2,Sj2,C2> m2b = (
                    m2.SubMatrix(kk,kk+KB,N1,N));
                InlineMultMM<true>(one, m1b, m2b, m3x);
              }
              if (K2) {
                InlineMultMM<true>(one, m1.SubMatrix(M1,M,K1,K),
                    m2.SubMatrix(K1,K,N1,N), m3x);
              }
            } else {
              InlineMultMM<false>(one, m1.Rows(M1,M), m2.Cols(N1,N), m3x);
            }
            if (add) m3.SubMatrix(M1,M,N1,N) += x * m3x;
            else m3.SubMatrix(M1,M,N1,N) = x * m3x;
          }
        } // M2 && N2
      } // end parallel
    } // Mb || Nb
    else if (Kb) {
      Matrix<T3,ColMajor> mtemp(m3.colsize(),m3.rowsize());

      ConstSmallMatrixView<T1,XX,KB,Si1,Sj1,C1> m1a = m1.Cols(0,KB);
      ConstSmallMatrixView<T2,KB,XX,Si2,Sj2,C2> m2a = m2.Rows(0,KB);
      InlineMultMM<false>(one, m1a, m2a, mtemp);

      for (int k=1;k<Kb;++k) {
        const int kk = k*KB;
        ConstSmallMatrixView<T1,XX,KB,Si1,Sj1,C1> m1b = m1.Cols(kk,kk+KB);
        ConstSmallMatrixView<T2,KB,XX,Si2,Sj2,C2> m2b = m2.Rows(kk,kk+KB);
        InlineMultMM<true>(one, m1b, m2b, mtemp);
      }
      if (K2) {
        InlineMultMM<true>(one, m1.Cols(K1,K), m2.Rows(K1,K), mtemp);
      }
      if (add) m3 += x * mtemp;
      else m3 = x * mtemp;
    }
    else {
      InlineMultMM<add>(x,m1,m2,m3);
    }
  }
#else
  template <bool add, int ix, class T, class MAT1, class MAT2, class MAT3>
  static void DoMultMM(
      const Scaling<ix,T>& x, const MAT1& m1, const MAT2& m2, MAT3& m3)
  { InlineMultMM<add>(x,m1,m2,m3); }
#endif

  template <bool add, class T, class M1, class M2, class M3>
  static void DoMultMM(const T x, const M1& m1, const M2& m2, M3& m3)
  {
    if (x == T(1))
      DoMultMM<add>(Scaling<1,T>(),m1,m2,m3);
    else if (x == T(-1))
      DoMultMM<add>(Scaling<-1,T>(),m1,m2,m3);
    else if (x == T(0))
    { if (!add) m3.Zero(); }
    else 
      DoMultMM<add>(Scaling<0,T>(x),m1,m2,m3);
  }

  template <bool add, class T, class M1, class M2, class M3>
  static void DoMultMM(
      const std::complex<T> x, const M1& m1, const M2& m2, M3& m3)
  {
    if (imag(x) == T(0)) {
      if (real(x) == T(1))
        DoMultMM<add>(Scaling<1,T>(),m1,m2,m3);
      else if (real(x) == T(-1))
        DoMultMM<add>(Scaling<-1,T>(),m1,m2,m3);
      else if (real(x) == T(0))
      { if (!add) m3.Zero(); }
      else
        DoMultMM<add>(Scaling<0,T>(real(x)),m1,m2,m3);
    } else 
      DoMultMM<add>(Scaling<0,std::complex<T> >(x),m1,m2,m3);
  }

} // namespace tmv

#endif

