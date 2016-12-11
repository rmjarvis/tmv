///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// The Template Matrix/Vector Library for C++ was created by Mike Jarvis     //
// Copyright (C) 1998 - 2016                                                 //
// All rights reserved                                                       //
//                                                                           //
// The project is hosted at https://code.google.com/p/tmv-cpp/               //
// where you can find the current version and current documention.           //
//                                                                           //
// For concerns or problems with the software, Mike may be contacted at      //
// mike_jarvis17 [at] gmail.                                                 //
//                                                                           //
// This software is licensed under a FreeBSD license.  The file              //
// TMV_LICENSE should have bee included with this distribution.              //
// It not, you can get a copy from https://code.google.com/p/tmv-cpp/.       //
//                                                                           //
// Essentially, you can use this software however you want provided that     //
// you include the TMV_LICENSE file in any distribution that uses it.        //
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

    template <typename T, class IT>
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
