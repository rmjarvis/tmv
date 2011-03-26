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


//-----------------------------------------------------------------------------
//
// This file defines the BasisVector class.
//
// A BasisVector is a vector whose element are all 0, except for one.
// Normally, this one element has a value of 1, but we allow for the 
// possibility of it being some other value.
// 
// Constructors:
//
//    BasisVector<T,A>(size_t n, int i)  
//        Makes a BasisVector of size n, with all values = 0, except for
//        v[i] = 1.
//
//    BasisVector<T,A>(size_t n, int i, T x)  
//        Makes a BasisVector of size n, with all values = 0, except for
//        v[i] = x.
//

#ifndef TMV_BasisVector_H
#define TMV_BasisVector_H

#include "TMV_BaseVector.h"

namespace tmv {

    template <class T, int A=0>
    class BasisVector;

    template <class T, int A>
    struct Traits<BasisVector<T,A> >
    {
        enum { okA = (
                A == CStyle || A == FortranStyle ) };

        typedef T value_type;

        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;
        enum { isreal = Traits<T>::isreal };
        enum { iscomplex = Traits<T>::iscomplex };

        typedef BasisVector<T,A> type;
        typedef Vector<T,A> calc_type;
        typedef const type& eval_type; 
        typedef calc_type copy_type;

        enum { _size = UNKNOWN }; 
        enum { _fort = Attrib<A>::fort };
        enum { _calc = true };
        enum { _checkalias = false };
    };

    template <class T, int A>
    class BasisVector : 
        public BaseVector<BasisVector<T,A> >
    {
    public:
        typedef BasisVector<T,A> type;

        enum { _size = Traits<type>::_size };
        enum { _fort = Traits<type>::_fort };
        enum { _calc = Traits<type>::_calc };
        enum { isreal = Traits<T>::isreal }; 
        enum { iscomplex = Traits<T>::iscomplex }; 

        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;

        //
        // Constructors
        //

        BasisVector(size_t n, int i, const T x=T(1)) : 
            itssize(n), itsindex(A==FortranStyle?i-1:i), itsval(x) 
            {
                TMVStaticAssert(Traits<type>::okA);
            }
        ~BasisVector() {}

        //
        // Auxilliary Functions
        //

        T cref(int i) const  { return i==itsindex ? itsval : T(0); }
        size_t size() const { return itssize; }
        int nElements() const { return 1; }
        template <class V2>
        void assignTo(BaseVector_Mutable<V2>& v2) const
        {
            TMVAssert(v2.size() == itssize);
            v2.setZero();
            v2.ref(itsindex) = itsval; 
        }

    protected :

        size_t itssize;
        int itsindex;
        T itsval;

    }; // BasisVector


    //
    // TMV_Text functions
    //

    template <class T, int A>
    static inline std::string TMV_Text(const BasisVector<T,A>& )
    {
        std::ostringstream s;
        s << "BasisVector<"<<TMV_Text(T());
        s <<","<<Attrib<A>::vtext()<<">";
        return s.str();
    }

} // namespace tmv

#endif
