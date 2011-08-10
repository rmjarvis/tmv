
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
