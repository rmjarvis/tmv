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


// Things that need to be #defined on entry:
// (The values for a normal Matrix+Matrix are given)
//
// SUMMM	SumMM
// GENMATRIX1	GenMatrix
// GENMATRIX2	GenMatrix
// PRODXM1	ProdXM
// PRODXM2	ProdXM

#ifndef X3
#define X3
#endif

#ifndef Y
#define Y
#endif

#ifndef CT
#define CT std::complex<T>
#endif
#ifndef CCT
#define CCT ConjRef<std::complex<T> >
#endif
#ifndef VCT
#define VCT VarConjRef<std::complex<T> >
#endif

#ifndef GETM1
#define GETM1 .getM1()
#endif

#ifndef GETM2
#define GETM2 .getM2()
#endif


// -(x*m+x*m)

template <typename T, typename T1, typename T2 Y>
inline SUMMM<T,T1,T2 X3> operator-(const SUMMM<T,T1,T2 X3>& smm)
{ return SUMMM<T,T1,T2 X3>(-smm.getX1(),smm GETM1,-smm.getX2, smm GETM2); }

// x*(x*m+x*m)

template <typename T, typename T1, typename T2 Y>
inline SUMMM<T,T1,T2 X3> operator*(const T x, const SUMMM<T,T1,T2 X3>& smm)
{ return SUMMM<T,T1,T2 X3>(smm.getX1()*x,smm GETM1,smm.getX2()*x, smm GETM2); }

template <typename T, typename T1, typename T2 Y>
inline SUMMM<CT,T1,T2 X3> operator*(const T x, const SUMMM<CT,T1,T2 X3>& smm)
{ return SUMMM<CT,T1,T2 X3>(smm.getX1()*x,smm GETM1,smm.getX2()*x, smm GETM2); }

template <typename T Y>
inline SUMMM<CT,T,T X3> operator*(const CT x, const SUMMM<T,T,T X3>& smm)
{ return SUMMM<CT,T,T X3>(smm.getX1()*x,smm GETM1,smm.getX2()*x, smm GETM2); }

template <typename T Y>
inline SUMMM<CT,T,T X3> operator*(const CCT x, const SUMMM<T,T,T X3>& smm)
{
    return SUMMM<CT,T,T X3>(
        smm.getX1()*CT(x),smm GETM1,smm.getX2()*CT(x),smm GETM2);
}

template <typename T Y>
inline SUMMM<CT,T,T X3> operator*(const VCT x, const SUMMM<T,T,T X3>& smm)
{
    return SUMMM<CT,T,T X3>(
        smm.getX1()*CT(x),smm GETM1,smm.getX2()*CT(x),smm GETM2);
}

template <typename T, typename T1, typename T2 Y>
inline SUMMM<CT,T1,T2 X3> operator*(const CCT x, const SUMMM<CT,T1,T2 X3>& smm)
{
    return SUMMM<CT,T1,T2 X3>(
        smm.getX1()*CT(x),smm GETM1,smm.getX2()*CT(x),smm GETM2);
}

template <typename T, typename T1, typename T2 Y>
inline SUMMM<CT,T1,T2 X3> operator*(const VCT x, const SUMMM<CT,T1,T2 X3>& smm)
{
    return SUMMM<CT,T1,T2 X3>(
        smm.getX1()*CT(x),smm GETM1,smm.getX2()*CT(x),smm GETM2);
}


// (x*m+x*m)*x

template <typename T, typename T1, typename T2 Y>
inline SUMMM<T,T1,T2 X3> operator*(const SUMMM<T,T1,T2 X3>& smm, const T x)
{ return SUMMM<T,T1,T2 X3>(smm.getX1()*x,smm GETM1,smm.getX2()*x, smm GETM2); }

template <typename T, typename T1, typename T2 Y>
inline SUMMM<CT,T1,T2 X3> operator*(const SUMMM<CT,T1,T2 X3>& smm, const T x)
{ return SUMMM<CT,T1,T2 X3>(smm.getX1()*x,smm GETM1,smm.getX2()*x, smm GETM2); }

template <typename T Y>
inline SUMMM<CT,T,T X3> operator*(const SUMMM<T,T,T X3>& smm, const CT x)
{ return SUMMM<CT,T,T X3>(smm.getX1()*x,smm GETM1,smm.getX2()*x, smm GETM2); }

template <typename T Y>
inline SUMMM<CT,T,T X3> operator*(const SUMMM<T,T,T X3>& smm, const CCT x)
{
    return SUMMM<CT,T,T X3>(
        smm.getX1()*CT(x),smm GETM1,smm.getX2()*CT(x),smm GETM2);
}

template <typename T Y>
inline SUMMM<CT,T,T X3> operator*(const SUMMM<T,T,T X3>& smm, const VCT x)
{
    return SUMMM<CT,T,T X3>(
        smm.getX1()*CT(x),smm GETM1,smm.getX2()*CT(x),smm GETM2);
}

template <typename T, typename T1, typename T2 Y>
inline SUMMM<CT,T1,T2 X3> operator*(const SUMMM<CT,T1,T2 X3>& smm, const CCT x)
{
    return SUMMM<CT,T1,T2 X3>(
        smm.getX1()*CT(x),smm GETM1,smm.getX2()*CT(x),smm GETM2);
}

template <typename T, typename T1, typename T2 Y>
inline SUMMM<CT,T1,T2 X3> operator*(const SUMMM<CT,T1,T2 X3>& smm, const VCT x)
{
    return SUMMM<CT,T1,T2 X3>(
        smm.getX1()*CT(x),smm GETM1,smm.getX2()*CT(x),smm GETM2);
}

// (x*m+x*m)/x

template <typename T, typename T1, typename T2 Y>
inline SUMMM<T,T1,T2 X3> operator/(const SUMMM<T,T1,T2 X3>& smm, const T x)
{
    return SUMMM<T,T1,T2 X3>(
        TMV_Divide(smm.getX1(),x),smm GETM1,TMV_Divide(smm.getX2(),x), smm GETM2);
}

template <typename T, typename T1, typename T2 Y>
inline SUMMM<CT,T1,T2 X3> operator/(const SUMMM<CT,T1,T2 X3>& smm, const T x)
{
    return SUMMM<CT,T1,T2 X3>(
        TMV_Divide(smm.getX1(),x),smm GETM1,TMV_Divide(smm.getX2(),x), smm GETM2);
}

template <typename T Y>
inline SUMMM<CT,T,T X3> operator/(const SUMMM<T,T,T X3>& smm, const CT x)
{
    return SUMMM<CT,T,T X3>(
        TMV_Divide(smm.getX1(),x),smm GETM1,TMV_Divide(smm.getX2(),x), smm GETM2);
}

template <typename T Y>
inline SUMMM<CT,T,T X3> operator/(const SUMMM<T,T,T X3>& smm, const CCT x)
{
    return SUMMM<CT,T,T X3>(
        TMV_Divide(smm.getX1(),CT(x)),smm GETM1,TMV_Divide(smm.getX2(),CT(x)),smm GETM2);
}

template <typename T Y>
inline SUMMM<CT,T,T X3> operator/(const SUMMM<T,T,T X3>& smm, const VCT x)
{
    return SUMMM<CT,T,T X3>(
        TMV_Divide(smm.getX1(),CT(x)),smm GETM1,TMV_Divide(smm.getX2(),CT(x)),smm GETM2);
}

template <typename T, typename T1, typename T2 Y>
inline SUMMM<CT,T1,T2 X3> operator/(const SUMMM<CT,T1,T2 X3>& smm, const CCT x)
{
    return SUMMM<CT,T1,T2 X3>(
        TMV_Divide(smm.getX1(),CT(x)),smm GETM1,TMV_Divide(smm.getX2(),CT(x)),smm GETM2);
}

template <typename T, typename T1, typename T2 Y>
inline SUMMM<CT,T1,T2 X3> operator/(const SUMMM<CT,T1,T2 X3>& smm, const VCT x)
{
    return SUMMM<CT,T1,T2 X3>(
        TMV_Divide(smm.getX1(),CT(x)),smm GETM1,TMV_Divide(smm.getX2(),CT(x)),smm GETM2);
}

#ifdef SUMMM_1_1

// -(x*m+x*m)

template <typename T, typename T1, typename T2 Y>
inline SUMMM_x_m1<T,T1,T2 X3> operator-(const SUMMM_1_1<T,T1,T2 X3>& smm)
{ return SUMMM<T,T1,T2 X3>(T(-1),smm GETM1,T(-1),smm GETM2); }

// x*(x*m+x*m)

template <typename T, typename T1, typename T2 Y>
inline SUMMM<T,T1,T2 X3> operator*(const T x, const SUMMM_1_1<T,T1,T2 X3>& smm)
{ return SUMMM<T,T1,T2 X3>(x,smm GETM1,x,smm GETM2); }

template <typename T, typename T1, typename T2 Y>
inline SUMMM<CT,T1,T2 X3> operator*(const T x, const SUMMM_1_1<CT,T1,T2 X3>& smm)
{ return SUMMM<CT,T1,T2 X3>(CT(x),smm GETM1,CT(x),smm GETM2); }

template <typename T Y>
inline SUMMM<CT,T,T X3> operator*(const CT x, const SUMMM_1_1<T,T,T X3>& smm)
{ return SUMMM<CT,T,T X3>(x,smm GETM1,x,smm GETM2); }

template <typename T Y>
inline SUMMM<CT,T,T X3> operator*(const CCT x, const SUMMM_1_1<T,T,T X3>& smm)
{ return SUMMM<CT,T,T X3>(CT(x),smm GETM1,CT(x),smm GETM2); }

template <typename T Y>
inline SUMMM<CT,T,T X3> operator*(const VCT x, const SUMMM_1_1<T,T,T X3>& smm)
{ return SUMMM<CT,T,T X3>(CT(x),smm GETM1,CT(x),smm GETM2); }

template <typename T, typename T1, typename T2 Y>
inline SUMMM<CT,T1,T2 X3> operator*(
    const CCT x, const SUMMM_1_1<CT,T1,T2 X3>& smm)
{ return SUMMM<CT,T1,T2 X3>(CT(x),smm GETM1,CT(x),smm GETM2); }

template <typename T, typename T1, typename T2 Y>
inline SUMMM<CT,T1,T2 X3> operator*(
    const VCT x, const SUMMM_1_1<CT,T1,T2 X3>& smm)
{ return SUMMM<CT,T1,T2 X3>(CT(x),smm GETM1,CT(x),smm GETM2); }


// (x*m+x*m)*x

template <typename T, typename T1, typename T2 Y>
inline SUMMM<T,T1,T2 X3> operator*(
    const SUMMM_1_1<T,T1,T2 X3>& smm, const T x)
{ return SUMMM<T,T1,T2 X3>(x,smm GETM1,x,smm GETM2); }

template <typename T, typename T1, typename T2 Y>
inline SUMMM<CT,T1,T2 X3> operator*(
    const SUMMM_1_1<CT,T1,T2 X3>& smm, const T x)
{ return SUMMM<CT,T1,T2 X3>(CT(x),smm GETM1,CT(x),smm GETM2); }

template <typename T Y>
inline SUMMM<CT,T,T X3> operator*(const SUMMM_1_1<T,T,T X3>& smm, const CT x)
{ return SUMMM<CT,T,T X3>(x,smm GETM1,x,smm GETM2); }

template <typename T Y>
inline SUMMM<CT,T,T X3> operator*(const SUMMM_1_1<T,T,T X3>& smm, const CCT x)
{ return SUMMM<CT,T,T X3>(CT(x),smm GETM1,CT(x),smm GETM2); }

template <typename T Y>
inline SUMMM<CT,T,T X3> operator*(const SUMMM_1_1<T,T,T X3>& smm, const VCT x)
{ return SUMMM<CT,T,T X3>(CT(x),smm GETM1,CT(x),smm GETM2); }

template <typename T, typename T1, typename T2 Y>
inline SUMMM<CT,T1,T2 X3> operator*(
    const SUMMM_1_1<CT,T1,T2 X3>& smm, const CCT x)
{ return SUMMM<CT,T1,T2 X3>(CT(x),smm GETM1,CT(x),smm GETM2); }

template <typename T, typename T1, typename T2 Y>
inline SUMMM<CT,T1,T2 X3> operator*(
    const SUMMM_1_1<CT,T1,T2 X3>& smm, const VCT x)
{ return SUMMM<CT,T1,T2 X3>(CT(x),smm GETM1,CT(x),smm GETM2); }

// (x*m+x*m)/x

template <typename T, typename T1, typename T2 Y>
inline SUMMM<T,T1,T2 X3> operator/(
    const SUMMM_1_1<T,T1,T2 X3>& smm, const T x)
{
    return SUMMM<T,T1,T2 X3>(
        TMV_RealType(T)(1)/x,smm GETM1,TMV_RealType(T)(1)/x,smm GETM2);
}

template <typename T, typename T1, typename T2 Y>
inline SUMMM<CT,T1,T2 X3> operator/(
    const SUMMM_1_1<CT,T1,T2 X3>& smm, const T x)
{
    return SUMMM<CT,T1,T2 X3>(
        TMV_InverseOf(x),smm GETM1,TMV_InverseOf(x),smm GETM2);
}

template <typename T Y>
inline SUMMM<CT,T,T X3> operator/(const SUMMM_1_1<T,T,T X3>& smm, const CT x)
{
    return SUMMM<CT,T,T X3>(
        TMV_InverseOf(x),smm GETM1,TMV_InverseOf(x),smm GETM2);
}

template <typename T Y>
inline SUMMM<CT,T,T X3> operator/(const SUMMM_1_1<T,T,T X3>& smm, const CCT x)
{
    return SUMMM<CT,T,T X3>(
        TMV_InverseOf(CT(x)),smm GETM1,TMV_InverseOf(CT(x)),smm GETM2);
}

template <typename T Y>
inline SUMMM<CT,T,T X3> operator/(const SUMMM_1_1<T,T,T X3>& smm, const VCT x)
{
    return SUMMM<CT,T,T X3>(
        TMV_InverseOf(CT(x)),smm GETM1,TMV_InverseOf(CT(x)),smm GETM2);
}

template <typename T, typename T1, typename T2 Y>
inline SUMMM<CT,T1,T2 X3> operator/(
    const SUMMM_1_1<CT,T1,T2 X3>& smm, const CCT x)
{
    return SUMMM<CT,T1,T2 X3>(
        TMV_InverseOf(CT(x)),smm GETM1,TMV_InverseOf(CT(x)),smm GETM2);
}

template <typename T, typename T1, typename T2 Y>
inline SUMMM<CT,T1,T2 X3> operator/(
    const SUMMM_1_1<CT,T1,T2 X3>& smm, const VCT x)
{
    return SUMMM<CT,T1,T2 X3>(
        TMV_InverseOf(CT(x)),smm GETM1,TMV_InverseOf(CT(x)),smm GETM2);
}


// SUMMM_1_m1

// -(x*m+x*m)

template <typename T, typename T1, typename T2 Y>
inline SUMMM_x_1<T,T1,T2 X3> operator-(const SUMMM_1_m1<T,T1,T2 X3>& smm)
{ return SUMMM<T,T1,T2 X3>(T(-1),smm GETM1,T(1),smm GETM2); }

// x*(x*m+x*m)

template <typename T, typename T1, typename T2 Y>
inline SUMMM<T,T1,T2 X3> operator*(
    const T x, const SUMMM_1_m1<T,T1,T2 X3>& smm)
{ return SUMMM<T,T1,T2 X3>(x,smm GETM1,-x,smm GETM2); }

template <typename T, typename T1, typename T2 Y>
inline SUMMM<CT,T1,T2 X3> operator*(
    const T x, const SUMMM_1_m1<CT,T1,T2 X3>& smm)
{ return SUMMM<CT,T1,T2 X3>(CT(x),smm GETM1,CT(-x),smm GETM2); }

template <typename T Y>
inline SUMMM<CT,T,T X3> operator*(const CT x, const SUMMM_1_m1<T,T,T X3>& smm)
{ return SUMMM<CT,T,T X3>(x,smm GETM1,-x,smm GETM2); }

template <typename T Y>
inline SUMMM<CT,T,T X3> operator*(const CCT x, const SUMMM_1_m1<T,T,T X3>& smm)
{ return SUMMM<CT,T,T X3>(CT(x),smm GETM1,-CT(x),smm GETM2); }

template <typename T Y>
inline SUMMM<CT,T,T X3> operator*(const VCT x, const SUMMM_1_m1<T,T,T X3>& smm)
{ return SUMMM<CT,T,T X3>(CT(x),smm GETM1,-CT(x),smm GETM2); }

template <typename T, typename T1, typename T2 Y>
inline SUMMM<CT,T1,T2 X3> operator*(
    const CCT x, const SUMMM_1_m1<CT,T1,T2 X3>& smm)
{ return SUMMM<CT,T1,T2 X3>(CT(x),smm GETM1,-CT(x),smm GETM2); }

template <typename T, typename T1, typename T2 Y>
inline SUMMM<CT,T1,T2 X3> operator*(
    const VCT x, const SUMMM_1_m1<CT,T1,T2 X3>& smm)
{ return SUMMM<CT,T1,T2 X3>(CT(x),smm GETM1,-CT(x),smm GETM2); }


// (x*m+x*m)*x

template <typename T, typename T1, typename T2 Y>
inline SUMMM<T,T1,T2 X3> operator*(
    const SUMMM_1_m1<T,T1,T2 X3>& smm, const T x)
{ return SUMMM<T,T1,T2 X3>(x,smm GETM1,-x,smm GETM2); }

template <typename T, typename T1, typename T2 Y>
inline SUMMM<CT,T1,T2 X3> operator*(
    const SUMMM_1_m1<CT,T1,T2 X3>& smm, const T x)
{ return SUMMM<CT,T1,T2 X3>(x,smm GETM1,CT(-x),smm GETM2); }

template <typename T Y>
inline SUMMM<CT,T,T X3> operator*(const SUMMM_1_m1<T,T,T X3>& smm, const CT x)
{ return SUMMM<CT,T,T X3>(x,smm GETM1,-x,smm GETM2); }

template <typename T Y>
inline SUMMM<CT,T,T X3> operator*(const SUMMM_1_m1<T,T,T X3>& smm, const CCT x)
{ return SUMMM<CT,T,T X3>(CT(x),smm GETM1,-CT(x),smm GETM2); }

template <typename T Y>
inline SUMMM<CT,T,T X3> operator*(const SUMMM_1_m1<T,T,T X3>& smm, const VCT x)
{ return SUMMM<CT,T,T X3>(CT(x),smm GETM1,-CT(x),smm GETM2); }

template <typename T, typename T1, typename T2 Y>
inline SUMMM<CT,T1,T2 X3> operator*(
    const SUMMM_1_m1<CT,T1,T2 X3>& smm, const CCT x)
{ return SUMMM<CT,T1,T2 X3>(CT(x),smm GETM1,-CT(x),smm GETM2); }

template <typename T, typename T1, typename T2 Y>
inline SUMMM<CT,T1,T2 X3> operator*(
    const SUMMM_1_m1<CT,T1,T2 X3>& smm, const VCT x)
{ return SUMMM<CT,T1,T2 X3>(CT(x),smm GETM1,-CT(x),smm GETM2); }

// (x*m+x*m)/x

template <typename T, typename T1, typename T2 Y>
inline SUMMM<T,T1,T2 X3> operator/(
    const SUMMM_1_m1<T,T1,T2 X3>& smm, const T x)
{
    return SUMMM<T,T1,T2 X3>(
        TMV_InverseOf(x),smm GETM1,TMV_InverseOf(-x),smm GETM2);
}

template <typename T, typename T1, typename T2 Y>
inline SUMMM<CT,T1,T2 X3> operator/(
    const SUMMM_1_m1<CT,T1,T2 X3>& smm, const T x)
{
    return SUMMM<CT,T1,T2 X3>(
        TMV_InverseOf(x),smm GETM1,TMV_InverseOf(-x),smm GETM2);
}

template <typename T Y>
inline SUMMM<CT,T,T X3> operator/(const SUMMM_1_m1<T,T,T X3>& smm, const CT x)
{
    return SUMMM<CT,T,T X3>(
        TMV_InverseOf(x),smm GETM1,TMV_InverseOf(-x),smm GETM2);
}

template <typename T Y>
inline SUMMM<CT,T,T X3> operator/(const SUMMM_1_m1<T,T,T X3>& smm, const CCT x)
{
    return SUMMM<CT,T,T X3>(
        TMV_InverseOf(CT(x)),smm GETM1,TMV_InverseOf(-CT(x)),smm GETM2);
}

template <typename T Y>
inline SUMMM<CT,T,T X3> operator/(const SUMMM_1_m1<T,T,T X3>& smm, const VCT x)
{
    return SUMMM<CT,T,T X3>(
        TMV_InverseOf(CT(x)),smm GETM1,TMV_InverseOf(-CT(x)),smm GETM2);
}

template <typename T, typename T1, typename T2 Y>
inline SUMMM<CT,T1,T2 X3> operator/(
    const SUMMM_1_m1<CT,T1,T2 X3>& smm, const CCT x)
{
    return SUMMM<CT,T1,T2 X3>(
        TMV_InverseOf(CT(x)),smm GETM1,TMV_InverseOf(-CT(x)),smm GETM2);
}

template <typename T, typename T1, typename T2 Y>
inline SUMMM<CT,T1,T2 X3> operator/(
    const SUMMM_1_m1<CT,T1,T2 X3>& smm, const VCT x)
{
    return SUMMM<CT,T1,T2 X3>(
        TMV_InverseOf(CT(x)),smm GETM1,TMV_InverseOf(-CT(x)),smm GETM2);
}

// SUMMM_1_m1

// -(x*m+x*m)

template <typename T, typename T1, typename T2 Y>
inline SUMMM<T,T1,T2 X3> operator-(const SUMMM_1_x<T,T1,T2 X3>& smm)
{ return SUMMM<T,T1,T2 X3>(T(-1),smm GETM1,-smm.getX2,smm GETM2); }

// x*(x*m+x*m)

template <typename T, typename T1, typename T2 Y>
inline SUMMM<T,T1,T2 X3> operator*(
    const T x, const SUMMM_1_x<T,T1,T2 X3>& smm)
{ return SUMMM<T,T1,T2 X3>(x,smm GETM1,smm.getX2()*x,smm GETM2); }

template <typename T, typename T1, typename T2 Y>
inline SUMMM<CT,T1,T2 X3> operator*(
    const T x, const SUMMM_1_x<CT,T1,T2 X3>& smm)
{ return SUMMM<CT,T1,T2 X3>(CT(x),smm GETM1,smm.getX2()*x,smm GETM2); }

template <typename T Y>
inline SUMMM<CT,T,T X3> operator*(const CT x, const SUMMM_1_x<T,T,T X3>& smm)
{ return SUMMM<CT,T,T X3>(x,smm GETM1,smm.getX2()*x,smm GETM2); }

template <typename T Y>
inline SUMMM<CT,T,T X3> operator*(const CCT x, const SUMMM_1_x<T,T,T X3>& smm)
{ return SUMMM<CT,T,T X3>(CT(x),smm GETM1,smm.getX2()*CT(x),smm GETM2); }

template <typename T Y>
inline SUMMM<CT,T,T X3> operator*(const VCT x, const SUMMM_1_x<T,T,T X3>& smm)
{ return SUMMM<CT,T,T X3>(CT(x),smm GETM1,smm.getX2()*CT(x),smm GETM2); }

template <typename T, typename T1, typename T2 Y>
inline SUMMM<CT,T1,T2 X3> operator*(
    const CCT x, const SUMMM_1_x<CT,T1,T2 X3>& smm)
{ return SUMMM<CT,T1,T2 X3>(CT(x),smm GETM1,smm.getX2()*CT(x),smm GETM2); }

template <typename T, typename T1, typename T2 Y>
inline SUMMM<CT,T1,T2 X3> operator*(
    const VCT x, const SUMMM_1_x<CT,T1,T2 X3>& smm)
{ return SUMMM<CT,T1,T2 X3>(CT(x),smm GETM1,smm.getX2()*CT(x),smm GETM2); }


// (x*m+x*m)*x

template <typename T, typename T1, typename T2 Y>
inline SUMMM<T,T1,T2 X3> operator*(
    const SUMMM_1_x<T,T1,T2 X3>& smm, const T x)
{ return SUMMM<T,T1,T2 X3>(x,smm GETM1,smm.getX2()*x,smm GETM2); }

template <typename T, typename T1, typename T2 Y>
inline SUMMM<CT,T1,T2 X3> operator*(
    const SUMMM_1_x<CT,T1,T2 X3>& smm, const T x)
{ return SUMMM<CT,T1,T2 X3>(CT(x),smm GETM1,smm.getX2()*x,smm GETM2); }

template <typename T Y>
inline SUMMM<CT,T,T X3> operator*(const SUMMM_1_x<T,T,T X3>& smm, const CT x)
{ return SUMMM<CT,T,T X3>(x,smm GETM1,smm.getX2()*x,smm GETM2); }

template <typename T Y>
inline SUMMM<CT,T,T X3> operator*(const SUMMM_1_x<T,T,T X3>& smm, const CCT x)
{ return SUMMM<CT,T,T X3>(CT(x),smm GETM1,smm.getX2()*CT(x),smm GETM2); }

template <typename T Y>
inline SUMMM<CT,T,T X3> operator*(const SUMMM_1_x<T,T,T X3>& smm, const VCT x)
{ return SUMMM<CT,T,T X3>(CT(x),smm GETM1,smm.getX2()*CT(x),smm GETM2); }

template <typename T, typename T1, typename T2 Y>
inline SUMMM<CT,T1,T2 X3> operator*(
    const SUMMM_1_x<CT,T1,T2 X3>& smm, const CCT x)
{ return SUMMM<CT,T1,T2 X3>(CT(x),smm GETM1,smm.getX2()*CT(x),smm GETM2); }

template <typename T, typename T1, typename T2 Y>
inline SUMMM<CT,T1,T2 X3> operator*(
    const SUMMM_1_x<CT,T1,T2 X3>& smm, const VCT x)
{ return SUMMM<CT,T1,T2 X3>(CT(x),smm GETM1,smm.getX2()*CT(x),smm GETM2); }

// (x*m+x*m)/x

template <typename T, typename T1, typename T2 Y>
inline SUMMM<T,T1,T2 X3> operator/(
    const SUMMM_1_x<T,T1,T2 X3>& smm, const T x)
{
    return SUMMM<T,T1,T2 X3>(
        TMV_InverseOf(x),smm GETM1,TMV_Divide(smm.getX2(),x),smm GETM2);
}

template <typename T, typename T1, typename T2 Y>
inline SUMMM<CT,T1,T2 X3> operator/(
    const SUMMM_1_x<CT,T1,T2 X3>& smm, const T x)
{
    return SUMMM<CT,T1,T2 X3>(
        TMV_InverseOf(x),smm GETM1,TMV_Divide(smm.getX2(),x),smm GETM2);
}

template <typename T Y>
inline SUMMM<CT,T,T X3> operator/(const SUMMM_1_x<T,T,T X3>& smm, const CT x)
{
    return SUMMM<CT,T,T X3>(
        TMV_InverseOf(x),smm GETM1,TMV_Divide(smm.getX2(),x),smm GETM2);
}

template <typename T Y>
inline SUMMM<CT,T,T X3> operator/(const SUMMM_1_x<T,T,T X3>& smm, const CCT x)
{
    return SUMMM<CT,T,T X3>(
        TMV_InverseOf(CT(x)),smm GETM1,TMV_Divide(smm.getX2(),CT(x)),smm GETM2);
}

template <typename T Y>
inline SUMMM<CT,T,T X3> operator/(const SUMMM_1_x<T,T,T X3>& smm, const VCT x)
{
    return SUMMM<CT,T,T X3>(
        TMV_InverseOf(CT(x)),smm GETM1,TMV_Divide(smm.getX2(),CT(x)),smm GETM2);
}

template <typename T, typename T1, typename T2 Y>
inline SUMMM<CT,T1,T2 X3> operator/(
    const SUMMM_1_x<CT,T1,T2 X3>& smm, const CCT x)
{
    return SUMMM<CT,T1,T2 X3>(
        TMV_InverseOf(CT(x)),smm GETM1,TMV_Divide(smm.getX2(),CT(x)),smm GETM2);
}

template <typename T, typename T1, typename T2 Y>
inline SUMMM<CT,T1,T2 X3> operator/(
    const SUMMM_1_x<CT,T1,T2 X3>& smm, const VCT x)
{
    return SUMMM<CT,T1,T2 X3>(
        TMV_InverseOf(CT(x)),smm GETM1,TMV_Divide(smm.getX2(),CT(x)),smm GETM2);
}

// SUMMM_1_m1

// -(x*m+x*m)

template <typename T, typename T1, typename T2 Y>
inline SUMMM_x_m1<T,T1,T2 X3> operator-(const SUMMM_x_1<T,T1,T2 X3>& smm)
{ return SUMMM_x_m1<T,T1,T2 X3>(-smm.getX1(),smm GETM1,T(-1),smm GETM2); }

// x*(x*m+x*m)

template <typename T, typename T1, typename T2 Y>
inline SUMMM<T,T1,T2 X3> operator*(const T x, const SUMMM_x_1<T,T1,T2 X3>& smm)
{ return SUMMM<T,T1,T2 X3>(smm.getX1()*x,smm GETM1,x,smm GETM2); }

template <typename T, typename T1, typename T2 Y>
inline SUMMM<CT,T1,T2 X3> operator*(
    const T x, const SUMMM_x_1<CT,T1,T2 X3>& smm)
{ return SUMMM<CT,T1,T2 X3>(smm.getX1()*x,smm GETM1,CT(x),smm GETM2); }

template <typename T Y>
inline SUMMM<CT,T,T X3> operator*(const CT x, const SUMMM_x_1<T,T,T X3>& smm)
{ return SUMMM<CT,T,T X3>(smm.getX1()*x,smm GETM1,x,smm GETM2); }

template <typename T Y>
inline SUMMM<CT,T,T X3> operator*(const CCT x, const SUMMM_x_1<T,T,T X3>& smm)
{ return SUMMM<CT,T,T X3>(smm.getX1()*CT(x),smm GETM1,CT(x),smm GETM2); }

template <typename T Y>
inline SUMMM<CT,T,T X3> operator*(const VCT x, const SUMMM_x_1<T,T,T X3>& smm)
{ return SUMMM<CT,T,T X3>(smm.getX1()*CT(x),smm GETM1,CT(x),smm GETM2); }

template <typename T, typename T1, typename T2 Y>
inline SUMMM<CT,T1,T2 X3> operator*(
    const CCT x, const SUMMM_x_1<CT,T1,T2 X3>& smm)
{ return SUMMM<CT,T1,T2 X3>(smm.getX1()*CT(x),smm GETM1,CT(x),smm GETM2); }

template <typename T, typename T1, typename T2 Y>
inline SUMMM<CT,T1,T2 X3> operator*(
    const VCT x, const SUMMM_x_1<CT,T1,T2 X3>& smm)
{ return SUMMM<CT,T1,T2 X3>(smm.getX1()*CT(x),smm GETM1,CT(x),smm GETM2); }


// (x*m+x*m)*x

template <typename T, typename T1, typename T2 Y>
inline SUMMM<T,T1,T2 X3> operator*(const SUMMM_x_1<T,T1,T2 X3>& smm, const T x)
{ return SUMMM<T,T1,T2 X3>(smm.getX1()*x,smm GETM1,x,smm GETM2); }

template <typename T, typename T1, typename T2 Y>
inline SUMMM<CT,T1,T2 X3> operator*(
    const SUMMM_x_1<CT,T1,T2 X3>& smm, const T x)
{ return SUMMM<CT,T1,T2 X3>(smm.getX1()*x,smm GETM1,CT(x),smm GETM2); }

template <typename T Y>
inline SUMMM<CT,T,T X3> operator*(const SUMMM_x_1<T,T,T X3>& smm, const CT x)
{ return SUMMM<CT,T,T X3>(smm.getX1()*x,smm GETM1,x,smm GETM2); }

template <typename T Y>
inline SUMMM<CT,T,T X3> operator*(const SUMMM_x_1<T,T,T X3>& smm, const CCT x)
{ return SUMMM<CT,T,T X3>(smm.getX1()*CT(x),smm GETM1,CT(x),smm GETM2); }

template <typename T Y>
inline SUMMM<CT,T,T X3> operator*(const SUMMM_x_1<T,T,T X3>& smm, const VCT x)
{ return SUMMM<CT,T,T X3>(smm.getX1()*CT(x),smm GETM1,CT(x),smm GETM2); }

template <typename T, typename T1, typename T2 Y>
inline SUMMM<CT,T1,T2 X3> operator*(
    const SUMMM_x_1<CT,T1,T2 X3>& smm, const CCT x)
{ return SUMMM<CT,T1,T2 X3>(smm.getX1()*CT(x),smm GETM1,CT(x),smm GETM2); }

template <typename T, typename T1, typename T2 Y>
inline SUMMM<CT,T1,T2 X3> operator*(
    const SUMMM_x_1<CT,T1,T2 X3>& smm, const VCT x)
{ return SUMMM<CT,T1,T2 X3>(smm.getX1()*CT(x),smm GETM1,CT(x),smm GETM2); }

// (x*m+x*m)/x

template <typename T, typename T1, typename T2 Y>
inline SUMMM<T,T1,T2 X3> operator/(const SUMMM_x_1<T,T1,T2 X3>& smm, const T x)
{
    return SUMMM<T,T1,T2 X3>(
        TMV_Divide(smm.getX1(),x),smm GETM1,TMV_InverseOf(x),smm GETM2);
}

template <typename T, typename T1, typename T2 Y>
inline SUMMM<CT,T1,T2 X3> operator/(const SUMMM_x_1<CT,T1,T2 X3>& smm, const T x)
{
    return SUMMM<CT,T1,T2 X3>(
        TMV_Divide(smm.getX1(),x),smm GETM1,TMV_InverseOf(x),smm GETM2);
}

template <typename T Y>
inline SUMMM<CT,T,T X3> operator/(const SUMMM_x_1<T,T,T X3>& smm, const CT x)
{
    return SUMMM<CT,T,T X3>(
        TMV_Divide(smm.getX1(),x),smm GETM1,TMV_InverseOf(x),smm GETM2);
}

template <typename T Y>
inline SUMMM<CT,T,T X3> operator/(const SUMMM_x_1<T,T,T X3>& smm, const CCT x)
{
    return SUMMM<CT,T,T X3>(
        TMV_Divide(smm.getX1(),CT(x)),smm GETM1,TMV_InverseOf(CT(x)),smm GETM2);
}

template <typename T Y>
inline SUMMM<CT,T,T X3> operator/(const SUMMM_x_1<T,T,T X3>& smm, const VCT x)
{
    return SUMMM<CT,T,T X3>(
        TMV_Divide(smm.getX1(),CT(x)),smm GETM1,TMV_InverseOf(CT(x)),smm GETM2);
}

template <typename T, typename T1, typename T2 Y>
inline SUMMM<CT,T1,T2 X3> operator/(
    const SUMMM_x_1<CT,T1,T2 X3>& smm, const CCT x)
{
    return SUMMM<CT,T1,T2 X3>(
        TMV_Divide(smm.getX1(),CT(x)),smm GETM1,TMV_InverseOf(CT(x)),smm GETM2);
}

template <typename T, typename T1, typename T2 Y>
inline SUMMM<CT,T1,T2 X3> operator/(
    const SUMMM_x_1<CT,T1,T2 X3>& smm, const VCT x)
{
    return SUMMM<CT,T1,T2 X3>(
        TMV_Divide(smm.getX1(),CT(x)),smm GETM1,TMV_InverseOf(CT(x)),smm GETM2);
}

// SUMMM_1_m1

// -(x*m+x*m)

template <typename T, typename T1, typename T2 Y>
inline SUMMM_x_1<T,T1,T2 X3> operator-(const SUMMM_x_m1<T,T1,T2 X3>& smm)
{ return SUMMM<T,T1,T2 X3>(-smm.getX1(),smm GETM1,T(1),smm GETM2); }

// x*(x*m+x*m)

template <typename T, typename T1, typename T2 Y>
inline SUMMM<T,T1,T2 X3> operator*(const T x, const SUMMM_x_m1<T,T1,T2 X3>& smm)
{ return SUMMM<T,T1,T2 X3>(smm.getX1()*x,smm GETM1,x,smm GETM2); }

template <typename T, typename T1, typename T2 Y>
inline SUMMM<CT,T1,T2 X3> operator*(
    const T x, const SUMMM_x_m1<CT,T1,T2 X3>& smm)
{ return SUMMM<CT,T1,T2 X3>(smm.getX1()*x,smm GETM1,CT(x),smm GETM2); }

template <typename T Y>
inline SUMMM<CT,T,T X3> operator*(const CT x, const SUMMM_x_m1<T,T,T X3>& smm)
{ return SUMMM<CT,T,T X3>(smm.getX1()*x,smm GETM1,x,smm GETM2); }

template <typename T Y>
inline SUMMM<CT,T,T X3> operator*(const CCT x, const SUMMM_x_m1<T,T,T X3>& smm)
{ return SUMMM<CT,T,T X3>(smm.getX1()*CT(x),smm GETM1,CT(x),smm GETM2); }

template <typename T Y>
inline SUMMM<CT,T,T X3> operator*(const VCT x, const SUMMM_x_m1<T,T,T X3>& smm)
{ return SUMMM<CT,T,T X3>(smm.getX1()*CT(x),smm GETM1,CT(x),smm GETM2); }

template <typename T, typename T1, typename T2 Y>
inline SUMMM<CT,T1,T2 X3> operator*(
    const CCT x, const SUMMM_x_m1<CT,T1,T2 X3>& smm)
{ return SUMMM<CT,T1,T2 X3>(smm.getX1()*CT(x),smm GETM1,CT(x),smm GETM2); }

template <typename T, typename T1, typename T2 Y>
inline SUMMM<CT,T1,T2 X3> operator*(
    const VCT x, const SUMMM_x_m1<CT,T1,T2 X3>& smm)
{ return SUMMM<CT,T1,T2 X3>(smm.getX1()*CT(x),smm GETM1,CT(x),smm GETM2); }


// (x*m+x*m)*x

template <typename T, typename T1, typename T2 Y>
inline SUMMM<T,T1,T2 X3> operator*(
    const SUMMM_x_m1<T,T1,T2 X3>& smm, const T x)
{ return SUMMM<T,T1,T2 X3>(smm.getX1()*x,smm GETM1,x,smm GETM2); }

template <typename T, typename T1, typename T2 Y>
inline SUMMM<CT,T1,T2 X3> operator*(
    const SUMMM_x_m1<CT,T1,T2 X3>& smm, const T x)
{ return SUMMM<CT,T1,T2 X3>(smm.getX1()*x,smm GETM1,CT(x),smm GETM2); }

template <typename T Y>
inline SUMMM<CT,T,T X3> operator*(const SUMMM_x_m1<T,T,T X3>& smm, const CT x)
{ return SUMMM<CT,T,T X3>(smm.getX1()*x,smm GETM1,x,smm GETM2); }

template <typename T Y>
inline SUMMM<CT,T,T X3> operator*(const SUMMM_x_m1<T,T,T X3>& smm, const CCT x)
{ return SUMMM<CT,T,T X3>(smm.getX1()*CT(x),smm GETM1,CT(x),smm GETM2); }

template <typename T Y>
inline SUMMM<CT,T,T X3> operator*(const SUMMM_x_m1<T,T,T X3>& smm, const VCT x)
{ return SUMMM<CT,T,T X3>(smm.getX1()*CT(x),smm GETM1,CT(x),smm GETM2); }

template <typename T, typename T1, typename T2 Y>
inline SUMMM<CT,T1,T2 X3> operator*(
    const SUMMM_x_m1<CT,T1,T2 X3>& smm, const CCT x)
{ return SUMMM<CT,T1,T2 X3>(smm.getX1()*CT(x),smm GETM1,CT(x),smm GETM2); }

template <typename T, typename T1, typename T2 Y>
inline SUMMM<CT,T1,T2 X3> operator*(
    const SUMMM_x_m1<CT,T1,T2 X3>& smm, const VCT x)
{ return SUMMM<CT,T1,T2 X3>(smm.getX1()*CT(x),smm GETM1,CT(x),smm GETM2); }

// (x*m+x*m)/x

template <typename T, typename T1, typename T2 Y>
inline SUMMM<T,T1,T2 X3> operator/(const SUMMM_x_m1<T,T1,T2 X3>& smm, const T x)
{
    return SUMMM<T,T1,T2 X3>(
        TMV_Divide(smm.getX1(),x),smm GETM1,TMV_InverseOf(-x),smm GETM2);
}

template <typename T, typename T1, typename T2 Y>
inline SUMMM<CT,T1,T2 X3> operator/(
    const SUMMM_x_m1<CT,T1,T2 X3>& smm, const T x)
{
    return SUMMM<CT,T1,T2 X3>(
        TMV_Divide(smm.getX1(),x),smm GETM1,TMV_InverseOf(-x),smm GETM2);
}

template <typename T Y>
inline SUMMM<CT,T,T X3> operator/(const SUMMM_x_m1<T,T,T X3>& smm, const CT x)
{
    return SUMMM<CT,T,T X3>(
        TMV_Divide(smm.getX1(),x),smm GETM1,TMV_InverseOf(-x),smm GETM2);
}

template <typename T Y>
inline SUMMM<CT,T,T X3> operator/(
    const SUMMM_x_m1<T,T,T X3>& smm, const CCT x)
{
    return SUMMM<CT,T,T X3>(
        TMV_Divide(smm.getX1(),CT(x)),smm GETM1,TMV_InverseOf(-CT(x)),smm GETM2);
}

template <typename T Y>
inline SUMMM<CT,T,T X3> operator/(const SUMMM_x_m1<T,T,T X3>& smm, const VCT x)
{
    return SUMMM<CT,T,T X3>(
        TMV_Divide(smm.getX1(),CT(x)),smm GETM1,TMV_InverseOf(-CT(x)),smm GETM2);
}

template <typename T, typename T1, typename T2 Y>
inline SUMMM<CT,T1,T2 X3> operator/(
    const SUMMM_x_m1<CT,T1,T2 X3>& smm, const CCT x)
{
    return SUMMM<CT,T1,T2 X3>(
        TMV_Divide(smm.getX1(),CT(x)),smm GETM1,TMV_InverseOf(-CT(x)),smm GETM2);
}

template <typename T, typename T1, typename T2 Y>
inline SUMMM<CT,T1,T2 X3> operator/(
    const SUMMM_x_m1<CT,T1,T2 X3>& smm, const VCT x)
{
    return SUMMM<CT,T1,T2 X3>(
        TMV_Divide(smm.getX1(),CT(x)),smm GETM1,TMV_InverseOf(-CT(x)),smm GETM2);
}
#undef SUMMM_1_1
#undef SUMMM_1_m1
#undef SUMMM_1_x
#undef SUMMM_x_m1
#undef SUMMM_x_1
#endif

#undef X3
#undef Y
#undef GETM1
#undef GETM2

