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


#define RT TMV_RealType(T)
#define CT TMV_ComplexType(T)

#define DefDivEq(T) \
    inline void LDivEq(MatrixView<T> m) const \
{ \
    TMVAssert(colsize() == rowsize()); \
    TMVAssert(m.colsize() == colsize()); \
    doLDivEq(m);  \
} \
inline void RDivEq(MatrixView<T> m) const \
{ \
    TMVAssert(colsize() == rowsize()); \
    TMVAssert(m.rowsize() == rowsize()); \
    doRDivEq(m); \
} \
inline void makeInverse(MatrixView<T> minv) const \
{ \
    TMVAssert(minv.rowsize() == colsize()); \
    TMVAssert(minv.colsize() == rowsize()); \
    doMakeInverse(minv); \
} \


DefDivEq(RT);
DefDivEq(CT);
#undef DefDivEq

#define DefDiv(T1,T2) \
    inline void LDiv(const GenMatrix<T1>& m1, MatrixView<T2> m0) const \
{ \
    TMVAssert(m0.rowsize() == m1.rowsize()); \
    TMVAssert(m0.colsize() == rowsize()); \
    TMVAssert(m1.colsize() == colsize()); \
    doLDiv(m1,m0); \
} \
inline void RDiv(const GenMatrix<T1>& m1, MatrixView<T2> m0) const \
{ \
    TMVAssert(m0.colsize() == m1.colsize()); \
    TMVAssert(m1.rowsize() == rowsize()); \
    TMVAssert(m0.rowsize() == colsize()); \
    doRDiv(m1,m0); \
}
DefDiv(RT,RT);
DefDiv(RT,CT);
DefDiv(CT,CT);
#undef DefDiv

#undef RT
#undef CT

inline void makeInverseATA(MatrixView<T> minv) const
{
    if (colsize() < rowsize()) {
        TMVAssert(minv.rowsize() == colsize());
        TMVAssert(minv.colsize() == colsize());
    } else {
        TMVAssert(minv.rowsize() == rowsize());
        TMVAssert(minv.colsize() == rowsize());
    }
    doMakeInverseATA(minv);
}

