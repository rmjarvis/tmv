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
#ifndef TMV_SymBandCHDiv_H
#define TMV_SymBandCHDiv_H

#include "tmv/TMV_BaseSymBandMatrix.h"

namespace tmv {

  template <class T> 
  void LDL_Decompose(const SymBandMatrixView<T>& A);

  template <class T1, class T2> 
  void CH_LDivEq(const GenSymBandMatrix<T1>& L, const MatrixView<T2>& m);
  template <class T1, class T2> 
  void CH_RDivEq(const GenSymBandMatrix<T1>& L, const MatrixView<T2>& m);
  template <class T1, class T2> 
  void CH_Inverse(
      const GenSymBandMatrix<T1>& LLx, const SymMatrixView<T2>& sinv);

  template <class T1, class T2> 
  void LDL_LDivEq(const GenSymBandMatrix<T1>& L, const MatrixView<T2>& m);
  template <class T1, class T2> 
  void LDL_RDivEq(const GenSymBandMatrix<T1>& L, const MatrixView<T2>& m);
  template <class T1, class T2> 
  void LDL_Inverse(const GenSymBandMatrix<T1>& LLx,
      const SymMatrixView<T2>& sinv);

#ifndef NOTHROW
  template <class T> 
  class NonPosDefHermBandMatrix :
    public NonPosDef
  {
  public:
    mutable auto_ptr<HermBandMatrix<T> > A;

    inline NonPosDefHermBandMatrix(const GenSymBandMatrix<T>& _A) :
      NonPosDef("HermBandMatrix Cholesky decomposition"), 
      A(new HermBandMatrix<T>(_A)) {}
    inline NonPosDefHermBandMatrix(const NonPosDefHermBandMatrix<T>& rhs) :
      A(rhs.A) {}
    inline ~NonPosDefHermBandMatrix() throw() {}

    inline void Write(std::ostream& os) const throw()
    {
      NonPosDef::Write(os);
      os<<"The partially decomposed matrix is \n"<<*A<<std::endl;
    }
  };

  template <class T> 
  class NonPosDefHermBandMatrix2 :
    public NonPosDefHermBandMatrix<T>
  {
  public:
    mutable auto_ptr<HermBandMatrix<T> > A0;

    inline NonPosDefHermBandMatrix2(const GenSymBandMatrix<T>& _A,
        const GenSymBandMatrix<T>& _A0) :
      NonPosDefHermBandMatrix<T>(_A), A0(new HermBandMatrix<T>(_A0)) {}
    inline NonPosDefHermBandMatrix2(const NonPosDefHermBandMatrix2<T>& rhs) :
      NonPosDefHermBandMatrix<T>(rhs), A0(rhs.A0) {}
    inline ~NonPosDefHermBandMatrix2() throw() {}

    inline void Write(std::ostream& os) const throw()
    {
      NonPosDefHermBandMatrix<T>::Write(os);
      os<<"The original matrix was \n"<<*A0<<std::endl;
    }
  };

  template <class T> 
  class NonPosDefSymBandLDL :
    public NonPosDef
  {
  public:
    mutable auto_ptr<SymBandMatrix<T> > A;

    inline NonPosDefSymBandLDL(const GenSymBandMatrix<T>& _A) :
      NonPosDef("SymBandMatrix LDL decomposition."), 
      A(new SymBandMatrix<T>(_A)) {}
    inline NonPosDefSymBandLDL(
        const NonPosDefSymBandLDL<T>& rhs) : A(rhs.A) {}
    inline ~NonPosDefSymBandLDL() throw() {}

    inline void Write(std::ostream& os) const throw()
    {
      NonPosDef::Write(os);
      os<<"The partially decomposed matrix is \n"<<*A<<std::endl;
    }
  };

  template <class T> 
  class NonPosDefHermBandLDL :
    public NonPosDef
  {
  public:
    mutable auto_ptr<HermBandMatrix<T> > A;

    inline NonPosDefHermBandLDL(const GenSymBandMatrix<T>& _A) :
      NonPosDef("HermBandMatrix LDL decomposition."), 
      A(new HermBandMatrix<T>(_A)) {}
    inline NonPosDefHermBandLDL(
        const NonPosDefHermBandLDL<T>& rhs) : A(rhs.A) {}
    inline ~NonPosDefHermBandLDL() throw() {}

    inline void Write(std::ostream& os) const throw()
    {
      NonPosDef::Write(os);
      os<<"The partially decomposed matrix is \n"<<*A<<std::endl;
    }
  };
#endif

}

#endif
