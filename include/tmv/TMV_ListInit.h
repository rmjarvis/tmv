
#ifndef TMV_ListInit_H
#define TMV_ListInit_H

namespace tmv {

    class ListReadError : 
        public ReadError
    {
    private :
        ptrdiff_t n;

    public :
        ListReadError(ptrdiff_t nleft) : n(nleft) {}
        void write(std::ostream& os) const throw()
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
        ListAssigner(IT ptr, ptrdiff_t nLeft ) :
            _ptr(ptr), _nLeft(nLeft), _isLast(true)
        {}

        ListAssigner(IT ptr, ptrdiff_t nLeft, const T& x) :
            _ptr(ptr), _nLeft(nLeft), _isLast(false)
        {
            if (_nLeft == 0) throw ListReadError(0);
            TMVAssert((_nLeft > 0) && "Too many elements in ListInit");
            *_ptr++ = x;
            --_nLeft;
        }

        ListAssigner( const ListAssigner<T,IT>& rhs ) : 
            _ptr(rhs._ptr), _nLeft(rhs._nLeft), _isLast(true)
        { rhs._isLast = false; }

        ~ListAssigner()
        { if (_nLeft > 0 && _isLast) throw ListReadError(_nLeft); }

        ListAssigner<T,IT> operator,(const T& x)
        {
            if (_nLeft == 0) throw ListReadError(0);
            TMVAssert((_nLeft > 0) && "Too many elements in ListInit");
            *_ptr = x;
            _isLast = false;
            return ListAssigner<T,IT>(++_ptr,_nLeft-1);
        }

    protected:
        IT  _ptr;
        ptrdiff_t _nLeft;
        mutable bool _isLast;
    };



} // namespace tmv

#endif
