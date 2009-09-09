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



namespace tmv {

  template <class T> class VectorReadError :
    public ReadError
  {
    public :
      int i;
      mutable auto_ptr<Vector<T> > v;
      char exp,got;
      size_t s;
      bool is, iseof, isbad;

      VectorReadError(int _i, const GenVector<T>& _v, std::istream& _is) :
	ReadError("Vector"),
	i(_i), v(new Vector<T>(_v)), exp(0), got(0), s(_v.size()),
	is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
      VectorReadError(std::istream& _is) :
	ReadError("Vector"),
	i(0), v(0), exp(0), got(0), s(0),
	is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
      VectorReadError(int _i, const GenVector<T>& _v, std::istream& _is,
	  char _e, char _g) :
	ReadError("Vector"),
	i(_i), v(new Vector<T>(_v)), exp(_e), got(_g), s(_v.size()),
	is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
      VectorReadError(const GenVector<T>& _v, std::istream& _is, size_t _s) :
	ReadError("Vector"),
	i(0), v(new Vector<T>(_v)), exp(0), got(0), s(_s),
	is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
      VectorReadError(const VectorReadError<T>& rhs) :
	ReadError("Vector"),
	i(rhs.i), v(rhs.v), exp(rhs.exp), got(rhs.got), s(rhs.s),
	is(rhs.is), iseof(rhs.iseof), isbad(rhs.isbad) {}

      virtual ~VectorReadError() throw() {}
      void Write(std::ostream& os) const throw();
  };

} // namespace mv


