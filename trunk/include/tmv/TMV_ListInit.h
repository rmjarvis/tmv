///////////////////////////////////////////////////////////////////////////////
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

#ifndef TMV_ListInit_H
#define TMV_ListInit_H

namespace tmv {

    class ListReadError : public ReadError
    {
    private :
        int n;

    public :
        inline ListReadError(int nleft) : n(nleft) {}
        inline void Write(std::ostream& os) const throw()
        {
            os<<"TMV Read Error: Reading from List initialization.\n";
            if (n == 0)
                os<<"List has more elements than expected.\n";
            else
                os<<"Reached end of list, but expecting "<<n<<
                    " more elements.\n";
        }
    };

    template <class T, class IT>
    class ListAssigner
    {
    public:
        inline ListAssigner(IT ptr, int nLeft ) :
            _ptr(ptr), _nLeft(nLeft), _isLast(true)
        {}

        inline ListAssigner(IT ptr, int nLeft, const T& x) :
            _ptr(ptr), _nLeft(nLeft), _isLast(false)
        {
            if (_nLeft == 0) throw ListReadError(0);
            TMVAssert((_nLeft > 0) && "Too many elements in ListInit");
            *_ptr++ = x;
            --_nLeft;
        }

        inline ListAssigner( const ListAssigner<T,IT>& rhs ) : 
            _ptr(rhs._ptr), _nLeft(rhs._nLeft), _isLast(true)
        { rhs._isLast = false; }

        inline ~ListAssigner()
        { if (_nLeft > 0 && _isLast) throw ListReadError(_nLeft); }

        inline ListAssigner<T,IT> operator,(const T& x)
        {
            if (_nLeft == 0) throw ListReadError(0);
            TMVAssert((_nLeft > 0) && "Too many elements in ListInit");
            *_ptr = x;
            _isLast = false;
            return ListAssigner<T,IT>(++_ptr,_nLeft-1);
        }

    protected:
        IT  _ptr;
        int _nLeft;
        mutable bool _isLast;
    };



} // namespace tmv

#endif
