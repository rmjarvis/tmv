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


// Need to define the following with #define statements.
// (The given definition is for a DiagMatrix/Matrix.  Modify as
// appropriate for the various other matrices.)
// Note: Vector/Matrix uses the same format - no need to write that separately.
//
// #define GENMATRIX1 GenDiagMatrix
// #define GENMATRIX2 GenMatrix
// #define TQUOTMM TransientQuotMM
// #define TRQUOTMM TransientRQuotMM
// #define PRODXM1 ProdXM
// #define PRODXM2 ProdXM
// #define QUOTXM QuotXM

template <typename T, typename T1, typename T2>
class TQUOTMM;
template <typename T, typename T1, typename T2>
class TRQUOTMM;

// m/m

template <typename T>
inline TQUOTMM<T,T,T> operator/(
    const GENMATRIX1<T>& m1, const GENMATRIX2<T>& m2)
{
    return TQUOTMM<T,T,T>(T(1),
                          auto_ptr<Matrix<T,ColMajor> >(new Matrix<T,ColMajor>(m1)),
                          m2);
}

template <typename T>
inline TQUOTMM<CT,CT,T> operator/(
    const GENMATRIX1<CT>& m1, const GENMATRIX2<T>& m2)
{
    return TQUOTMM<CT,CT,T>(T(1),
                            auto_ptr<Matrix<CT,ColMajor> >(new Matrix<CT,ColMajor>(m1)),
                            m2);
}

template <typename T>
inline TQUOTMM<CT,T,CT> operator/(
    const GENMATRIX1<T>& m1, const GENMATRIX2<CT>& m2)
{
    return TQUOTMM<CT,T,CT>(CT(1),
                            auto_ptr<Matrix<T,ColMajor> >(new Matrix<T,ColMajor>(m1)),
                            m2);
}

// m%m

template <typename T>
inline TRQUOTMM<T,T,T> operator%(
    const GENMATRIX1<T>& m1, const GENMATRIX2<T>& m2)
{
    return TRQUOTMM<T,T,T>(T(1),
                           auto_ptr<Matrix<T,RowMajor> >(new Matrix<T,RowMajor>(m1)),
                           m2);
}

template <typename T>
inline TRQUOTMM<CT,CT,T> operator%(
    const GENMATRIX1<CT>& m1, const GENMATRIX2<T>& m2)
{
    return TRQUOTMM<CT,CT,T>(T(1),
                             auto_ptr<Matrix<CT,RowMajor> >(new Matrix<CT,RowMajor>(m1)),
                             m2);
}

template <typename T>
inline TRQUOTMM<CT,T,CT> operator%(
    const GENMATRIX1<T>& m1, const GENMATRIX2<CT>& m2)
{
    return TRQUOTMM<CT,T,CT>(CT(1),
                             auto_ptr<Matrix<T,RowMajor> >(new Matrix<T,RowMajor>(m1)),
                             m2);
}

// (x*m)/m

template <typename T, typename T1>
inline TQUOTMM<T,T1,T> operator/(
    const PRODXM1<T,T1>& pxm, const GENMATRIX2<T>& m)
{
    return TQUOTMM<T,T1,T>(pxm.getX(),
                           auto_ptr<Matrix<T1,ColMajor> >(new Matrix<T1,ColMajor>(pxm.getM())),
                           m);
}

template <typename T>
inline TQUOTMM<CT,T,CT> operator/(
    const PRODXM1<T,T>& pxm, const GENMATRIX2<CT>& m)
{
    return TQUOTMM<CT,T,CT>(CT(pxm.getX()),
                            auto_ptr<Matrix<T,ColMajor> >(new Matrix<T,ColMajor>(pxm.getM())),
                            m);
}

template <typename T, typename T1>
inline TQUOTMM<CT,T1,T> operator/(
    const PRODXM1<CT,T1>& pxm, const GENMATRIX2<T>& m)
{
    return TQUOTMM<CT,T1,T>(pxm.getX(),
                            auto_ptr<Matrix<T1,ColMajor> >(new Matrix<T1,ColMajor>(pxm.getM())),
                            m);
}

// (x*m)%m

template <typename T, typename T1>
inline TRQUOTMM<T,T1,T> operator%(
    const PRODXM1<T,T1>& pxm, const GENMATRIX2<T>& m)
{
    return TRQUOTMM<T,T1,T>(pxm.getX(),
                            auto_ptr<Matrix<T1,RowMajor> >(new Matrix<T1,RowMajor>(pxm.getM())),
                            m);
}

template <typename T>
inline TRQUOTMM<CT,T,CT> operator%(
    const PRODXM1<T,T>& pxm, const GENMATRIX2<CT>& m)
{
    return TRQUOTMM<CT,T,CT>(CT(pxm.getX()),
                             auto_ptr<Matrix<T,RowMajor> >(new Matrix<T,RowMajor>(pxm.getM())),
                             m);
}

template <typename T, typename T1>
inline TRQUOTMM<CT,T1,T> operator%(
    const PRODXM1<CT,T1>& pxm, const GENMATRIX2<T>& m)
{
    return TRQUOTMM<CT,T1,T>(pxm.getX(),
                             auto_ptr<Matrix<T1,RowMajor> >(new Matrix<T1,RowMajor>(pxm.getM())),
                             m);
}

// m/(x*m)

template <typename T, typename T2>
inline TQUOTMM<T,T,T2> operator/(
    const GENMATRIX1<T>& m, const PRODXM2<T,T2>& pxm)
{
    return TQUOTMM<T,T,T2>(TMV_InverseOf(pxm.getX()),
                           auto_ptr<Matrix<T,ColMajor> >(new Matrix<T,ColMajor>(m)),
                           pxm.getM());
}

template <typename T>
inline TQUOTMM<CT,CT,T> operator/(
    const GENMATRIX1<CT>& m, const PRODXM2<T,T>& pxm)
{
    return TQUOTMM<CT,CT,T>(TMV_InverseOf(pxm.getX()),
                            auto_ptr<Matrix<CT,ColMajor> >(new Matrix<CT,ColMajor>(m)),
                            pxm.getM());
}

template <typename T, typename T2>
inline TQUOTMM<CT,T,T2> operator/(
    const GENMATRIX1<T>& m, const PRODXM2<CT,T2>& pxm)
{
    return TQUOTMM<CT,T,T2>(TMV_InverseOf(pxm.getX()),
                            auto_ptr<Matrix<T,ColMajor> >(new Matrix<T,ColMajor>(m)),
                            pxm.getM());
}

// m%(x*m)

template <typename T, typename T2>
inline TRQUOTMM<T,T,T2> operator%(
    const GENMATRIX1<T>& m, const PRODXM2<T,T2>& pxm)
{
    return TRQUOTMM<T,T,T2>(TMV_InverseOf(pxm.getX()),
                            auto_ptr<Matrix<T,RowMajor> >(new Matrix<T,RowMajor>(m)),
                            pxm.getM());
}

template <typename T>
inline TRQUOTMM<CT,CT,T> operator%(
    const GENMATRIX1<CT>& m, const PRODXM2<T,T>& pxm)
{
    return TRQUOTMM<CT,CT,T>(TMV_InverseOf(pxm.getX()),
                             auto_ptr<Matrix<CT,RowMajor> >(new Matrix<CT,RowMajor>(m)),
                             pxm.getM());
}

template <typename T, typename T2>
inline TRQUOTMM<CT,T,T2> operator%(
    const GENMATRIX1<T>& m, const PRODXM2<CT,T2>& pxm)
{
    return TRQUOTMM<CT,T,T2>(TMV_InverseOf(pxm.getX()),
                             auto_ptr<Matrix<T,RowMajor> >(new Matrix<T,RowMajor>(m)),
                             pxm.getM());
}

// (x*m)/(x*m)

template <typename T, typename T1, typename T2>
inline TQUOTMM<T,T1,T2> operator/(
    const PRODXM1<T,T1>& pxm1, const PRODXM2<T,T2>& pxm2)
{
    return TQUOTMM<T,T1,T2>(TMV_Divide(pxm1.getX(),pxm2.getX()),
                            auto_ptr<Matrix<T1,ColMajor> >(new Matrix<T1,ColMajor>(pxm1.getM())),
                            pxm2.getM());
}

template <typename T, typename T1>
inline TQUOTMM<CT,T1,T> operator/(
    const PRODXM1<CT,T1>& pxm1, const PRODXM2<T,T>& pxm2)
{
    return TQUOTMM<CT,T1,T>(TMV_Divide(pxm1.getX(),pxm2.getX()),
                            auto_ptr<Matrix<T1,ColMajor> >(new Matrix<T1,ColMajor>(pxm1.getM())),
                            pxm2.getM());
}

template <typename T, typename T2>
inline TQUOTMM<CT,T,T2> operator/(
    const PRODXM1<T,T>& pxm1, const PRODXM2<CT,T2>& pxm2)
{
    return TQUOTMM<CT,T,T2>(TMV_Divide(pxm1.getX(),pxm2.getX()),
                            auto_ptr<Matrix<T,ColMajor> >(new Matrix<T,ColMajor>(pxm1.getM())),
                            pxm2.getM());
}

// (x*m)%(x*m)

template <typename T, typename T1, typename T2>
inline TRQUOTMM<T,T1,T2> operator%(
    const PRODXM1<T,T1>& pxm1, const PRODXM2<T,T2>& pxm2)
{
    return TRQUOTMM<T,T1,T2>(TMV_Divide(pxm1.getX(),pxm2.getX()),
                             auto_ptr<Matrix<T1,RowMajor> >(new Matrix<T1,RowMajor>(pxm1.getM())),
                             pxm2.getM());
}

template <typename T, typename T1>
inline TRQUOTMM<CT,T1,T> operator%(
    const PRODXM1<CT,T1>& pxm1, const PRODXM2<T,T>& pxm2)
{
    return TRQUOTMM<CT,T1,T>(TMV_Divide(pxm1.getX(),pxm2.getX()),
                             auto_ptr<Matrix<T1,RowMajor> >(new Matrix<T1,RowMajor>(pxm1.getM())),
                             pxm2.getM());
}

template <typename T, typename T2>
inline TRQUOTMM<CT,T,T2> operator%(
    const PRODXM1<T,T>& pxm1, const PRODXM2<CT,T2>& pxm2)
{
    return TRQUOTMM<CT,T,T2>(TMV_Divide(pxm1.getX(),pxm2.getX()),
                             auto_ptr<Matrix<T,RowMajor> >(new Matrix<T,RowMajor>(pxm1.getM())),
                             pxm2.getM());
}

// (x/m)*m

template <typename T, typename T2>
inline TQUOTMM<T,T,T2> operator*(
    const QUOTXM<T,T2>& qxm, const GENMATRIX1<T>& m)
{
    return TQUOTMM<T,T,T2>(qxm.getX(),
                           auto_ptr<Matrix<T,ColMajor> >(new Matrix<T,ColMajor>(m)),
                           qxm.getM());
}

template <typename T>
inline TQUOTMM<CT,CT,T> operator*(
    const QUOTXM<T,T>& qxm, const GENMATRIX1<CT>& m)
{
    return TQUOTMM<CT,CT,T>(qxm.getX(),
                            auto_ptr<Matrix<CT,ColMajor> >(new Matrix<CT,ColMajor>(m)),
                            qxm.getM());
}

template <typename T, typename T1>
inline TQUOTMM<CT,T,T1> operator*(
    const QUOTXM<CT,T1>& qxm, const GENMATRIX1<T>& m)
{
    return TQUOTMM<CT,T,T1>(qxm.getX(),
                            auto_ptr<Matrix<T,ColMajor> >(new Matrix<T,ColMajor>(m)),
                            qxm.getM());
}

// m*(x/m)

template <typename T, typename T2>
inline TRQUOTMM<T,T,T2> operator*(
    const GENMATRIX1<T>& m, const QUOTXM<T,T2>& qxm)
{
    return TRQUOTMM<T,T,T2>(qxm.getX(),
                            auto_ptr<Matrix<T,RowMajor> >(new Matrix<T,RowMajor>(m)),
                            qxm.getM());
}

template <typename T>
inline TRQUOTMM<CT,CT,T> operator*(
    const GENMATRIX1<CT>& m, const QUOTXM<T,T>& qxm)
{
    return TRQUOTMM<CT,CT,T>(qxm.getX(),
                             auto_ptr<Matrix<CT,RowMajor> >(new Matrix<CT,RowMajor>(m)),
                             qxm.getM());
}

template <typename T, typename T2>
inline TRQUOTMM<CT,T,T2> operator*(
    const GENMATRIX1<T>& m, const QUOTXM<CT,T2>& qxm)
{
    return TRQUOTMM<CT,T,T2>(qxm.getX(),
                             auto_ptr<Matrix<T,RowMajor> >(new Matrix<T,RowMajor>(m)),
                             qxm.getM());
}

// (x/m)*(x*m)

template <typename T, typename T1, typename T2>
inline TQUOTMM<T,T1,T2> operator*(
    const QUOTXM<T,T2>& qxm, const PRODXM1<T,T1>& pxm)
{
    return TQUOTMM<T,T1,T2>(pxm.getX()*qxm.getX(),
                            auto_ptr<Matrix<T1,ColMajor> >(new Matrix<T1,ColMajor>(pxm.getM())),
                            qxm.getM());
}

template <typename T, typename T1>
inline TQUOTMM<CT,T1,T> operator*(
    const QUOTXM<T,T>& qxm, const PRODXM1<CT,T1>& pxm)
{
    return TQUOTMM<CT,T1,T>(pxm.getX()*qxm.getX(),
                            auto_ptr<Matrix<T1,ColMajor> >(new Matrix<T1,ColMajor>(pxm.getM())),
                            qxm.getM());
}

template <typename T, typename T2>
inline TQUOTMM<CT,T,T2> operator*(
    const QUOTXM<CT,T2>& qxm, const PRODXM1<T,T>& pxm)
{
    return TQUOTMM<CT,T,T2>(pxm.getX()*qxm.getX(),
                            auto_ptr<Matrix<T,ColMajor> >(new Matrix<T,ColMajor>(pxm.getM())),
                            qxm.getM());
}

// (x*m)*(x/m)

template <typename T, typename T1, typename T2>
inline TRQUOTMM<T,T1,T2> operator*(
    const PRODXM1<T,T1>& pxm, const QUOTXM<T,T2>& qxm)
{
    return TRQUOTMM<T,T1,T2>(pxm.getX()*qxm.getX(),
                             auto_ptr<Matrix<T1,RowMajor> >(new Matrix<T1,RowMajor>(pxm.getM())),
                             qxm.getM());
}

template <typename T, typename T1>
inline TRQUOTMM<CT,T1,T> operator*(
    const PRODXM1<CT,T1>& pxm, const QUOTXM<T,T>& qxm)
{
    return TRQUOTMM<CT,T1,T>(pxm.getX()*qxm.getX(),
                             auto_ptr<Matrix<T1,RowMajor> >(new Matrix<T1,RowMajor>(pxm.getM())),
                             qxm.getM());
}

template <typename T, typename T2>
inline TRQUOTMM<CT,T,T2> operator*(
    const PRODXM1<T,T>& pxm, const QUOTXM<CT,T2>& qxm)
{
    return TRQUOTMM<CT,T,T2>(pxm.getX()*qxm.getX(),
                             auto_ptr<Matrix<T,RowMajor> >(new Matrix<T,RowMajor>(pxm.getM())),
                             qxm.getM());
}

