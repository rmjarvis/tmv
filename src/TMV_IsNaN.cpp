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

// This needs to be math.h, not cmath, since sometimes
// cmath puts isnan in the std namespace and undefines the 
// regular isnan.  But sometimes it doesn't.
// So writing std::isnan is not portable.

#include "math.h"
#include "TMV_IsNaN.h"

namespace tmv {

    template <class T> 
    static bool DoIsNaN(T x)
    {
        return (x != x) || !(x * x >= 0);
    }

#ifdef TMV_INST_DOUBLE
    static bool DoIsNaN(double x)
    {
#ifdef isnan
        return isnan(x);
#else
        return (x != x) || !(x * x >= 0);
#endif
    }
#endif

#ifdef TMV_INST_FLOAT
    static bool DoIsNaN(float x)
    {
#ifdef isnan
        return isnan(x);
#else
        return (x != x) || !(x * x >= 0);
#endif
    }
#endif

#ifdef TMV_INST_LONGDOUBLE
    static bool DoIsNaN(long double x)
    {
#ifdef isnan
        return isnan(x);
#else
        return (x != x) || !(x * x >= 0);
#endif
    }
#endif

    template <class T> 
    bool IsNaN(T x)
    {
        return DoIsNaN(x);
    }

#define InstFile "TMV_IsNaN.inst"
#include "TMV_Inst.h"
#undef InstFile

}
