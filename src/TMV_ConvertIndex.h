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


#ifndef TMV_ConvertIndex_H
#define TMV_ConvertIndex_H

#include <vector>

namespace tmv {

    template <class VI> 
    inline void ConvertIndexToPermute(ptrdiff_t n, const VI& newIndex, ptrdiff_t* P)
    {
        // newIndex[i]=j means value at original j location needs to go to i.
        std::vector<ptrdiff_t> currIndex(n);
        std::vector<ptrdiff_t> origIndex(n);
        for(ptrdiff_t i=0;i<n;++i) {
            currIndex[i] = i;
            origIndex[i] = i;
        } 
        // currIndex[i]=j means value at original i location is currently at j.
        // origIndex[j]=i means value at original i location is currently at j.
        for(ptrdiff_t i=0;i<n;++i) {
            ptrdiff_t ip = currIndex[newIndex[i]];
            P[i] = ip;
            if (i != ip) { 
                ptrdiff_t origi = origIndex[i];
                ptrdiff_t origip = origIndex[ip];
                currIndex[origi] = ip;
                currIndex[origip] = i;
                origIndex[i] = origip;
                origIndex[ip] = origi;
            } 
        } 
    }   

}   

#endif
