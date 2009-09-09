///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// The Template Matrix/Vector Library for C++ was created by Mike Jarvis     //
// Copyright (C) 2007                                                        //
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
// along with this program in the file gpl.txt.                              //
//                                                                           //
// If not, write to:                                                         //
// The Free Software Foundation, Inc.                                        //
// 51 Franklin Street, Fifth Floor,                                          //
// Boston, MA  02110-1301, USA.                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#ifndef TMV_VECTORAUX_H
#define TMV_VECTORAUX_H

#include <vector>

template <class VI> inline void ConvertIndexToPermute(
          size_t n, const VI& newindex, size_t* P)
{
  // newindex[i]=j means value at original j location needs to go to i.
  std::vector<size_t> currindex(n);
  std::vector<size_t> origindex(n);
  for(size_t i=0;i<n;++i) {
    currindex[i] = i;
    origindex[i] = i;
  } 
  // currindex[i]=j means value at original i location is currently at j.
  // origindex[j]=i means value at original i location is currently at j.
  for(size_t i=0;i<n;++i) {
    size_t ip = currindex[newindex[i]];
    P[i] = ip;
    if (i != ip) { 
      size_t origi = origindex[i];
      size_t origip = origindex[ip];
      currindex[origi] = ip;
      currindex[origip] = i;
      origindex[i] = origip;
      origindex[ip] = origi;
    } 
  } 
}   

#endif
