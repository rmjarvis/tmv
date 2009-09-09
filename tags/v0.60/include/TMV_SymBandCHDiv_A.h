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


#ifndef TMV_SYMCHDiv_A_H
#define TMV_SYMCHDiv_A_H

namespace tmv {

  template <class T> class NonPosDefHermBandMatrix :
    public NonPosDef
  {
    public:
      mutable auto_ptr<HermBandMatrix<T> > A;

      inline NonPosDefHermBandMatrix(const GenSymBandMatrix<T>& _A) :
	NonPosDef("HermBandMatrix"), A(new HermBandMatrix<T>(_A)) {}
      inline NonPosDefHermBandMatrix(const NonPosDefHermBandMatrix<T>& rhs) :
	A(rhs.A) {}
      inline ~NonPosDefHermBandMatrix() throw() {}

      inline void Write(std::ostream& os) const throw()
      {
	NonPosDef::Write(os);
	os<<"The partially decomposed matrix is \n"<<*A<<std::endl;
      }
  };

  template <class T> class NonPosDefHermBandMatrix2 :
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


}

#endif
