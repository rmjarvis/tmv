#ifndef TMV_IsNaN_H
#define TMV_IsNaN_H

// The isnan function/macro is not portable. 
// So I encapsulate all the non-portability issues with this in TMV_IsNaN.cpp.

namespace tmv {

    template <class T> 
    bool IsNaN(T x);

} // namespace mv

#endif
