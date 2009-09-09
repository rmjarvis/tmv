///////////////////////////////////////////////////////////////////////////////
// vim:et:ts=2:sw=2:ci:cino=f0,g0,t0,+0:
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


#ifndef TMV_Prefetch_H
#define TMV_Prefetch_H

namespace tmv {

  // We have 4 functions that (might - depending on the compiler) 
  // prefetch data.  Two are for reading and two are for writing.
  // In each case the Multi version indicates that we are going to 
  // use the data multiple times, so try to keep it in cache.

  // The only supported prefetch so far is the GNU builtin prefetch.
  // TODO: Find out what the correlate on other compilers is.
  // I'm sure icpc must have something for this.

#ifdef __GNUC__
  // TODO: I'm not sure if the __GNUC__ guard is sufficient.  I might
  // need to check which version we are using.  I should find out
  // when gcc added this feature.
  static inline void Prefetch_Read(const void* p) 
  { __builtin_prefetch(p,0,0); }
  static inline void Prefetch_MultiRead(const void* p) 
  { __builtin_prefetch(p,0,3); }
  static inline void Prefetch_Write(void* p) 
  { __builtin_prefetch(p,1,0); }
  static inline void Prefetch_MultiWrite(void* p) 
  { __builtin_prefetch(p,1,3); }
#define TMV_PREFETCH
#else
  static inline void Prefetch_Read(const void* ) {}
  static inline void Prefetch_MultiRead(const void* ) {}
  static inline void Prefetch_Write(void* ) {}
  static inline void Prefetch_MultiWrite(void* ) {}
#endif

} // namespace tmv

#endif
