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


// This file is included at the end of each cpp file to instantiate
// the templates for each of the types desired.
// Add more types as desired using the format below:

#ifdef TMV_INST_DOUBLE
#define T double
#include InstFile
#undef T
#endif

#ifdef TMV_INST_FLOAT
#define T float
#include InstFile
#undef T
#endif

#ifdef TMV_INST_SKIP_BLAS
#undef TMV_INST_SKIP_BLAS
#endif

#ifdef TMV_INST_INT
// Define TISINT for any integer type.  e.g. long, short, etc.
#define TISINT 
#define T int
#include InstFile
#undef T
#endif

#ifdef TMV_INST_LONGDOUBLE
#define T long double
#include InstFile
#undef T
#endif

