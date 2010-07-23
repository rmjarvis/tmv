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


#ifndef TMV_ConvertIndex_H
#define TMV_ConvertIndex_H

#include <vector>

namespace tmv {

    template <class VI> 
    inline void ConvertIndexToPermute(int n, const VI& newIndex, int* P)
    {
        // newIndex[i]=j means value at original j location needs to go to i.
        std::vector<int> currIndex(n);
        std::vector<int> origIndex(n);
        for(int i=0;i<n;++i) {
            currIndex[i] = i;
            origIndex[i] = i;
        } 
        // currIndex[i]=j means value at original i location is currently at j.
        // origIndex[j]=i means value at original i location is currently at j.
        for(int i=0;i<n;++i) {
            int ip = currIndex[newIndex[i]];
            P[i] = ip;
            if (i != ip) { 
                int origi = origIndex[i];
                int origip = origIndex[ip];
                currIndex[origi] = ip;
                currIndex[origip] = i;
                origIndex[i] = origip;
                origIndex[ip] = origi;
            } 
        } 
    }   

}   

#endif
