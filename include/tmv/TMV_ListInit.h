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
        ptrdiff_t n;

    public :
        inline ListReadError(ptrdiff_t nLeft) : n(nLeft) {}
        inline void write(std::ostream& os) const throw()
        {
            os<<"TMV Error: Reading from list initialization.\n";
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
        ListAssigner(IT ptr, ptrdiff_t nLeft) :
            _ptr(ptr), _nLeft(nLeft), _isLast(true) 
        {}

        ListAssigner(IT ptr, ptrdiff_t nLeft, const T& x) : 
            _ptr(ptr), _nLeft(nLeft), _isLast(false)
        {
#ifdef TMV_DEBUG
            if (_nLeft == 0) throw ListReadError(0);
#endif
            TMVAssert(_nLeft > 0);
            *_ptr++ = x;
            --_nLeft;
        }

        ListAssigner(const ListAssigner<T,IT>& rhs) : 
            _ptr(rhs._ptr), _nLeft(rhs._nLeft), _isLast(true)
        { rhs._isLast = false; }

        ~ListAssigner()
        { 
#ifdef TMV_DEBUG
            if (_nLeft > 0 && _isLast) throw ListReadError(_nLeft); 
#endif
        }

        ListAssigner<T,IT> operator,(const T& x)
        {
#ifdef TMV_DEBUG
            if (_nLeft == 0) throw ListReadError(0);
#endif
            TMVAssert( _nLeft > 0 );
            *_ptr++ = x;
            _isLast = false;
            return ListAssigner<T,IT>(_ptr,_nLeft-1);
        }

        ListAssigner<T,IT> operator<<(const T& x)
        { return operator,(x); }

    protected:
        IT  _ptr;
        ptrdiff_t _nLeft;
        mutable bool _isLast;
    };



} // namespace tmv

#endif
