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



#include "TMV_Blas.h"
#include "TMV_TriMatrix.h"
#include "TMV_TriDiv.h"
#include "TMV_TriMatrixArith.h"

//#define XDEBUG

#ifdef XDEBUG
#include "TMV_MatrixArith.h"
#include <iostream>
using std::cerr;
using std::endl;
#endif

namespace tmv {

#ifdef TMV_BLOCKSIZE
#define TRI_DIV_BLOCKSIZE TMV_BLOCKSIZE
#else
#define TRI_DIV_BLOCKSIZE 64
#endif
 
  template <bool unit, class T> inline void RecursiveInverse(
      const UpperTriMatrixView<T>& U)
  {
    TMVAssert(U.iscm() || U.isrm());
    TMVAssert(unit == U.isunit());

    const size_t N = U.size();
    const size_t nb = TRI_DIV_BLOCKSIZE;

    if (N == 1) {
      if (!unit) {
	T*const Uptr = U.ptr();
	if (*Uptr == T(0)) 
	  throw SingularUpperTriMatrix<T>(U);
	*Uptr = RealType(T)(1) / (*Uptr);
      }
    } else {
      size_t k = N/2;
      if (k > nb) k = k/nb*nb;

      UpperTriMatrixView<T> U00 = U.SubTriMatrix(0,k);
      MatrixView<T> U01 = U.SubMatrix(0,k,k,N);
      UpperTriMatrixView<T> U11 = U.SubTriMatrix(k,N);

      // U00 U01' + U01 U11' = 0
      // U00 U01' = -U01 U11'
      // U01' = -U00' U01 U11'

      RecursiveInverse<unit>(U00);
      RecursiveInverse<unit>(U11);
      U01 = -U00 * U01;
      U01 *= U11;
    }
  }

  template <class T> inline void NonLapInverse(const UpperTriMatrixView<T>& U)
  {
    try {
      if (U.isunit()) RecursiveInverse<true>(U);
      else RecursiveInverse<false>(U);
    }
    catch (Singular) {
      throw SingularUpperTriMatrix<T>(U);
    }
  }

#ifdef ALAP
  template <class T> inline void LapInverse(const UpperTriMatrixView<T>& m)
  { NonLapInverse(m); }
#ifdef INST_DOUBLE
  template <> inline void LapInverse(const UpperTriMatrixView<double>& m)
  {
    int n = m.size();
    int lda = m.iscm() ? m.stepj() : m.stepi();
    LAPNAME(dtrtri) (LAPCM m.iscm()?LAPCH_UP:LAPCH_LO,
	m.isunit()?LAPCH_U:LAPCH_NU,LAPV(n),LAPP(m.ptr()),LAPV(lda) LAPINFO
	LAP1 LAP1);
    LAP_Results("dtrtri");
  }
  template <> inline void LapInverse(
      const UpperTriMatrixView<std::complex<double> >& m)
  {
    int n = m.size();
    int lda = m.iscm() ? m.stepj() : m.stepi();
    LAPNAME(ztrtri) (LAPCM m.iscm()?LAPCH_UP:LAPCH_LO,
	m.isunit()?LAPCH_U:LAPCH_NU,LAPV(n),LAPP(m.ptr()),LAPV(lda) LAPINFO
	LAP1 LAP1);
    LAP_Results("ztrtri");
  }
#endif
#ifdef INST_FLOAT
  template <> inline void LapInverse(const UpperTriMatrixView<float>& m)
  {
    int n = m.size();
    int lda = m.iscm() ? m.stepj() : m.stepi();
    LAPNAME(strtri) (LAPCM m.iscm()?LAPCH_UP:LAPCH_LO,
	m.isunit()?LAPCH_U:LAPCH_NU,LAPV(n),LAPP(m.ptr()),LAPV(lda) LAPINFO
	LAP1 LAP1);
    LAP_Results("strtri");
  }
  template <> inline void LapInverse(
      const UpperTriMatrixView<std::complex<float> >& m)
  {
    int n = m.size();
    int lda = m.iscm() ? m.stepj() : m.stepi();
    LAPNAME(ctrtri) (LAPCM m.iscm()?LAPCH_UP:LAPCH_LO,
	m.isunit()?LAPCH_U:LAPCH_NU,LAPV(n),LAPP(m.ptr()),LAPV(lda) LAPINFO
	LAP1 LAP1);
    LAP_Results("ctrtri");
  }
#endif // FLOAT
#endif // ALAP

  template <class T> void Tri_Inverse(const UpperTriMatrixView<T>& U)
  {
#ifdef XDEBUG
    Matrix<T> U0(U);
#endif

    if (U.size() > 0) {
      if (!(U.iscm() || U.isrm())) {
	UpperTriMatrix<T> temp = U;
	Tri_Inverse(temp.View());
	U = temp;
      } else {
#ifdef ALAP
	LapInverse(U);
#else
	NonLapInverse(U);
#endif
      }
    }
#ifdef XDEBUG
    Matrix<T> eye = U*U0;
    if (Norm(eye-T(1)) > 0.0001*(Norm(U0)+Norm(U))) {
      cerr<<"UpperTriMatrix Inverse:\n";
      cerr<<"U = "<<tmv::Type(U)<<"  "<<U0<<endl;
      cerr<<"Uinv = "<<U<<endl;
      cerr<<"Uinv*U = "<<U*U0<<endl;
      cerr<<"U*Uinv = "<<U0*U<<endl;
      abort();
    }
#endif
  }

#define InstFile "TMV_TriDiv_D.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


