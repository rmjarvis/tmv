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

// This needs to be math.h, not cmath, since sometimes
// cmath puts isnan in the std namespace and undefines the 
// regular isnan.  But sometimes it doesn't.
// So writing std::isnan is not portable.

//#include <iostream>
#include "math.h"
#include "tmv/TMV_Base.h"
#include "TMV_IsNaN.h"

namespace tmv {

    template <class T> 
    static bool DoIsNaN(T x)
    {
        return (x != x) || !(x * x >= 0);
    }

#ifdef INST_DOUBLE
    static bool DoIsNaN(double x)
    {
#ifdef isnan
        return isnan(x);
#else
        return (x != x) || !(x * x >= 0);
#endif
    }
#endif

#ifdef INST_FLOAT
    static bool DoIsNaN(float x)
    {
#ifdef isnan
        return isnan(x);
#else
        return (x != x) || !(x * x >= 0);
#endif
    }
#endif

#ifdef INST_LONGDOUBLE
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
    bool isNaN(T x)
    { return DoIsNaN(x); }

#ifdef INST_INT
#undef INST_INT
#endif

#define InstFile "TMV_IsNaN.inst"
#include "TMV_Inst.h"
#undef InstFile

}
