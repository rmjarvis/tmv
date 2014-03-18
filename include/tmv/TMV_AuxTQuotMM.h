///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// The Template Matrix/Vector Library for C++ was created by Mike Jarvis     //
// Copyright (C) 1998 - 2014                                                 //
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

template <class T, class T1, class T2>
class TQUOTMM;
template <class T, class T1, class T2> 
class TRQUOTMM;

// m/m

template <class T> 
inline TQUOTMM<T,T,T> operator/(
    const GENMATRIX1<T>& m1, const GENMATRIX2<T>& m2)
{
    std::auto_ptr<Matrix<T,ColMajor> > newm1(new Matrix<T,ColMajor>(m1));
    return TQUOTMM<T,T,T>(T(1),newm1,m2); 
}

template <class T> 
inline TQUOTMM<CT,CT,T> operator/(
    const GENMATRIX1<CT>& m1, const GENMATRIX2<T>& m2)
{
    std::auto_ptr<Matrix<CT,ColMajor> > newm1(new Matrix<CT,ColMajor>(m1));
    return TQUOTMM<CT,CT,T>(T(1),newm1,m2); 
}

template <class T> 
inline TQUOTMM<CT,T,CT> operator/(
    const GENMATRIX1<T>& m1, const GENMATRIX2<CT>& m2)
{
    std::auto_ptr<Matrix<T,ColMajor> > newm1(new Matrix<T,ColMajor>(m1));
    return TQUOTMM<CT,T,CT>(CT(1),newm1,m2); 
}

// m%m

template <class T> 
inline TRQUOTMM<T,T,T> operator%(
    const GENMATRIX1<T>& m1, const GENMATRIX2<T>& m2)
{
    std::auto_ptr<Matrix<T,RowMajor> > newm1(new Matrix<T,RowMajor>(m1));
    return TRQUOTMM<T,T,T>(T(1),newm1,m2); 
}

template <class T> 
inline TRQUOTMM<CT,CT,T> operator%(
    const GENMATRIX1<CT>& m1, const GENMATRIX2<T>& m2)
{
    std::auto_ptr<Matrix<CT,RowMajor> > newm1(new Matrix<CT,RowMajor>(m1));
    return TRQUOTMM<CT,CT,T>(T(1),newm1,m2); 
}

template <class T> 
inline TRQUOTMM<CT,T,CT> operator%(
    const GENMATRIX1<T>& m1, const GENMATRIX2<CT>& m2)
{
    std::auto_ptr<Matrix<T,RowMajor> > newm1(new Matrix<T,RowMajor>(m1));
    return TRQUOTMM<CT,T,CT>(CT(1),newm1,m2); 
}

// (x*m)/m

template <class T, class T1> 
inline TQUOTMM<T,T1,T> operator/(
    const PRODXM1<T,T1>& pxm, const GENMATRIX2<T>& m)
{ 
    std::auto_ptr<Matrix<T1,ColMajor> > newm1(new Matrix<T1,ColMajor>(pxm.getM()));
    return TQUOTMM<T,T1,T>(pxm.getX(),newm1,m); 
}

template <class T> 
inline TQUOTMM<CT,T,CT> operator/(
    const PRODXM1<T,T>& pxm, const GENMATRIX2<CT>& m)
{ 
    std::auto_ptr<Matrix<T,ColMajor> > newm1(new Matrix<T,ColMajor>(pxm.getM()));
    return TQUOTMM<CT,T,CT>(CT(pxm.getX()),newm1,m); 
}

template <class T, class T1> 
inline TQUOTMM<CT,T1,T> operator/(
    const PRODXM1<CT,T1>& pxm, const GENMATRIX2<T>& m)
{ 
    std::auto_ptr<Matrix<T1,ColMajor> > newm1(new Matrix<T1,ColMajor>(pxm.getM()));
    return TQUOTMM<CT,T1,T>(pxm.getX(),newm1,m); 
}

// (x*m)%m

template <class T, class T1> 
inline TRQUOTMM<T,T1,T> operator%(
    const PRODXM1<T,T1>& pxm, const GENMATRIX2<T>& m)
{ 
    std::auto_ptr<Matrix<T1,RowMajor> > newm1(new Matrix<T1,RowMajor>(pxm.getM()));
    return TRQUOTMM<T,T1,T>(pxm.getX(),newm1,m); 
}

template <class T> 
inline TRQUOTMM<CT,T,CT> operator%(
    const PRODXM1<T,T>& pxm, const GENMATRIX2<CT>& m)
{ 
    std::auto_ptr<Matrix<T,RowMajor> > newm1(new Matrix<T,RowMajor>(pxm.getM()));
    return TRQUOTMM<CT,T,CT>(CT(pxm.getX()),newm1,m); 
}

template <class T, class T1> 
inline TRQUOTMM<CT,T1,T> operator%(
    const PRODXM1<CT,T1>& pxm, const GENMATRIX2<T>& m)
{ 
    std::auto_ptr<Matrix<T1,RowMajor> > newm1(new Matrix<T1,RowMajor>(pxm.getM()));
    return TRQUOTMM<CT,T1,T>(pxm.getX(),newm1,m); 
}

// m/(x*m)

template <class T, class T2> 
inline TQUOTMM<T,T,T2> operator/(
    const GENMATRIX1<T>& m, const PRODXM2<T,T2>& pxm)
{ 
    std::auto_ptr<Matrix<T,ColMajor> > newm1(new Matrix<T,ColMajor>(m));
    return TQUOTMM<T,T,T2>(TMV_InverseOf(pxm.getX()),newm1,pxm.getM());
}

template <class T> 
inline TQUOTMM<CT,CT,T> operator/(
    const GENMATRIX1<CT>& m, const PRODXM2<T,T>& pxm)
{ 
    std::auto_ptr<Matrix<CT,ColMajor> > newm1(new Matrix<CT,ColMajor>(m));
    return TQUOTMM<CT,CT,T>(TMV_InverseOf(pxm.getX()),newm1,pxm.getM());
}

template <class T, class T2> 
inline TQUOTMM<CT,T,T2> operator/(
    const GENMATRIX1<T>& m, const PRODXM2<CT,T2>& pxm)
{ 
    std::auto_ptr<Matrix<T,ColMajor> > newm1(new Matrix<T,ColMajor>(m));
    return TQUOTMM<CT,T,T2>(TMV_InverseOf(pxm.getX()),newm1,pxm.getM());
}

// m%(x*m)

template <class T, class T2> 
inline TRQUOTMM<T,T,T2> operator%(
    const GENMATRIX1<T>& m, const PRODXM2<T,T2>& pxm)
{ 
    std::auto_ptr<Matrix<T,RowMajor> > newm1(new Matrix<T,RowMajor>(m));
    return TRQUOTMM<T,T,T2>(TMV_InverseOf(pxm.getX()),newm1,pxm.getM());
}

template <class T> 
inline TRQUOTMM<CT,CT,T> operator%(
    const GENMATRIX1<CT>& m, const PRODXM2<T,T>& pxm)
{ 
    std::auto_ptr<Matrix<CT,RowMajor> > newm1(new Matrix<CT,RowMajor>(m));
    return TRQUOTMM<CT,CT,T>(TMV_InverseOf(pxm.getX()),newm1,pxm.getM());
}

template <class T, class T2> 
inline TRQUOTMM<CT,T,T2> operator%(
    const GENMATRIX1<T>& m, const PRODXM2<CT,T2>& pxm)
{ 
    std::auto_ptr<Matrix<T,RowMajor> > newm1(new Matrix<T,RowMajor>(m));
    return TRQUOTMM<CT,T,T2>(TMV_InverseOf(pxm.getX()),newm1,pxm.getM());
}

// (x*m)/(x*m)

template <class T, class T1, class T2> 
inline TQUOTMM<T,T1,T2> operator/(
    const PRODXM1<T,T1>& pxm1, const PRODXM2<T,T2>& pxm2)
{ 
    std::auto_ptr<Matrix<T1,ColMajor> > newm1(new Matrix<T1,ColMajor>(pxm1.getM()));
    return TQUOTMM<T,T1,T2>(TMV_Divide(pxm1.getX(),pxm2.getX()),newm1,pxm2.getM());
}

template <class T, class T1> 
inline TQUOTMM<CT,T1,T> operator/(
    const PRODXM1<CT,T1>& pxm1, const PRODXM2<T,T>& pxm2)
{ 
    std::auto_ptr<Matrix<T1,ColMajor> > newm1(new Matrix<T1,ColMajor>(pxm1.getM()));
    return TQUOTMM<CT,T1,T>(TMV_Divide(pxm1.getX(),pxm2.getX()),newm1,pxm2.getM());
}

template <class T, class T2> 
inline TQUOTMM<CT,T,T2> operator/(
    const PRODXM1<T,T>& pxm1, const PRODXM2<CT,T2>& pxm2)
{ 
    std::auto_ptr<Matrix<T,ColMajor> > newm1(new Matrix<T,ColMajor>(pxm1.getM()));
    return TQUOTMM<CT,T,T2>(TMV_Divide(pxm1.getX(),pxm2.getX()),newm1,pxm2.getM());
}

// (x*m)%(x*m)

template <class T, class T1, class T2> 
inline TRQUOTMM<T,T1,T2> operator%(
    const PRODXM1<T,T1>& pxm1, const PRODXM2<T,T2>& pxm2)
{ 
    std::auto_ptr<Matrix<T1,RowMajor> > newm1(new Matrix<T1,RowMajor>(pxm1.getM()));
    return TRQUOTMM<T,T1,T2>(TMV_Divide(pxm1.getX(),pxm2.getX()),newm1,pxm2.getM());
}

template <class T, class T1> 
inline TRQUOTMM<CT,T1,T> operator%(
    const PRODXM1<CT,T1>& pxm1, const PRODXM2<T,T>& pxm2)
{ 
    std::auto_ptr<Matrix<T1,RowMajor> > newm1(new Matrix<T1,RowMajor>(pxm1.getM()));
    return TRQUOTMM<CT,T1,T>(TMV_Divide(pxm1.getX(),pxm2.getX()),newm1,pxm2.getM());
}

template <class T, class T2> 
inline TRQUOTMM<CT,T,T2> operator%(
    const PRODXM1<T,T>& pxm1, const PRODXM2<CT,T2>& pxm2)
{ 
    std::auto_ptr<Matrix<T,RowMajor> > newm1(new Matrix<T,RowMajor>(pxm1.getM()));
    return TRQUOTMM<CT,T,T2>(TMV_Divide(pxm1.getX(),pxm2.getX()),newm1,pxm2.getM());
}

// (x/m)*m

template <class T, class T2> 
inline TQUOTMM<T,T,T2> operator*(
    const QUOTXM<T,T2>& qxm, const GENMATRIX1<T>& m)
{ 
    std::auto_ptr<Matrix<T,ColMajor> > newm1(new Matrix<T,ColMajor>(m));
    return TQUOTMM<T,T,T2>(qxm.getX(),newm1,qxm.getM()); 
}

template <class T> 
inline TQUOTMM<CT,CT,T> operator*(
    const QUOTXM<T,T>& qxm, const GENMATRIX1<CT>& m)
{ 
    std::auto_ptr<Matrix<CT,ColMajor> > newm1(new Matrix<CT,ColMajor>(m));
    return TQUOTMM<CT,CT,T>(qxm.getX(),newm1,qxm.getM()); 
}

template <class T, class T1> 
inline TQUOTMM<CT,T,T1> operator*(
    const QUOTXM<CT,T1>& qxm, const GENMATRIX1<T>& m)
{ 
    std::auto_ptr<Matrix<T,ColMajor> > newm1(new Matrix<T,ColMajor>(m));
    return TQUOTMM<CT,T,T1>(qxm.getX(),newm1,qxm.getM()); 
}

// m*(x/m)

template <class T, class T2> 
inline TRQUOTMM<T,T,T2> operator*(
    const GENMATRIX1<T>& m, const QUOTXM<T,T2>& qxm)
{ 
    std::auto_ptr<Matrix<T,RowMajor> > newm1(new Matrix<T,RowMajor>(m));
    return TRQUOTMM<T,T,T2>(qxm.getX(),newm1,qxm.getM()); 
}

template <class T> 
inline TRQUOTMM<CT,CT,T> operator*(
    const GENMATRIX1<CT>& m, const QUOTXM<T,T>& qxm)
{ 
    std::auto_ptr<Matrix<CT,RowMajor> > newm1(new Matrix<CT,RowMajor>(m));
    return TRQUOTMM<CT,CT,T>(qxm.getX(),newm1,qxm.getM()); 
}

template <class T, class T2> 
inline TRQUOTMM<CT,T,T2> operator*(
    const GENMATRIX1<T>& m, const QUOTXM<CT,T2>& qxm)
{ 
    std::auto_ptr<Matrix<T,RowMajor> > newm1(new Matrix<T,RowMajor>(m));
    return TRQUOTMM<CT,T,T2>(qxm.getX(),newm1,qxm.getM()); 
}

// (x/m)*(x*m)

template <class T, class T1, class T2> 
inline TQUOTMM<T,T1,T2> operator*(
    const QUOTXM<T,T2>& qxm, const PRODXM1<T,T1>& pxm)
{ 
    std::auto_ptr<Matrix<T1,ColMajor> > newm1(new Matrix<T1,ColMajor>(pxm.getM()));
    return TQUOTMM<T,T1,T2>(pxm.getX()*qxm.getX(),newm1,qxm.getM());
}

template <class T, class T1> 
inline TQUOTMM<CT,T1,T> operator*(
    const QUOTXM<T,T>& qxm, const PRODXM1<CT,T1>& pxm)
{ 
    std::auto_ptr<Matrix<T1,ColMajor> > newm1(new Matrix<T1,ColMajor>(pxm.getM()));
    return TQUOTMM<CT,T1,T>(pxm.getX()*qxm.getX(),newm1,qxm.getM());
}

template <class T, class T2> 
inline TQUOTMM<CT,T,T2> operator*(
    const QUOTXM<CT,T2>& qxm, const PRODXM1<T,T>& pxm)
{ 
    std::auto_ptr<Matrix<T,ColMajor> > newm1(new Matrix<T,ColMajor>(pxm.getM()));
    return TQUOTMM<CT,T,T2>(pxm.getX()*qxm.getX(),newm1,qxm.getM());
}

// (x*m)*(x/m)

template <class T, class T1, class T2> 
inline TRQUOTMM<T,T1,T2> operator*(
    const PRODXM1<T,T1>& pxm, const QUOTXM<T,T2>& qxm)
{ 
    std::auto_ptr<Matrix<T1,RowMajor> > newm1(new Matrix<T1,RowMajor>(pxm.getM()));
    return TRQUOTMM<T,T1,T2>(pxm.getX()*qxm.getX(),newm1,qxm.getM());
}

template <class T, class T1> 
inline TRQUOTMM<CT,T1,T> operator*(
    const PRODXM1<CT,T1>& pxm, const QUOTXM<T,T>& qxm)
{ 
    std::auto_ptr<Matrix<T1,RowMajor> > newm1(new Matrix<T1,RowMajor>(pxm.getM()));
    return TRQUOTMM<CT,T1,T>(pxm.getX()*qxm.getX(),newm1,qxm.getM());
}

template <class T, class T2> 
inline TRQUOTMM<CT,T,T2> operator*(
    const PRODXM1<T,T>& pxm, const QUOTXM<CT,T2>& qxm)
{ 
    std::auto_ptr<Matrix<T,RowMajor> > newm1(new Matrix<T,RowMajor>(pxm.getM()));
    return TRQUOTMM<CT,T,T2>(pxm.getX()*qxm.getX(),newm1,qxm.getM());
}

