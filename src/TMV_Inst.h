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


// This file is included at the end of each cpp file to instantiate
// the templates for each of the types desired.
// Add more types as desired using the format below:

#ifdef INST_DOUBLE
#define T double
#include InstFile
#undef T
#endif

#ifdef INST_FLOAT
#define T float
#include InstFile
#undef T
#endif

#ifdef INST_SKIP_BLAS
#undef INST_SKIP_BLAS
#endif

#ifdef INST_INT
// Define TISINT for any integer type.  e.g. long, short, etc.
#define TISINT 
#define T int
#include InstFile
#undef T
#undef TISINT
#endif

#ifdef INST_LONGDOUBLE
#define T long double
#include InstFile
#undef T
#endif

