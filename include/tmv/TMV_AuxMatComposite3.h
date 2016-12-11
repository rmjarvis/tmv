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


// This file sets up the Composite QuotXM classes for Matrix types which
// returns a Matrix for the operations x/m or x%m.
// (Some do not: DiagMatrix for example, returns another DiagMatix.)
//
// Need to define the following with #define statements.
// (The given definition is for a Band Matrix.  Modify as
// appropriate for the various other matrices.)
//
// #define GENMATRIX2 GenBandMatrix
//
// #define QUOTXM QuotXB


//
// Scalar / Matrix
//

template <typename T, typename Tm>
class QUOTXM : public MatrixComposite<T>
{
public:
    inline QUOTXM(const T _x, const GENMATRIX2<Tm>& _m) : x(_x), m(_m) {}
    inline ptrdiff_t colsize() const { return m.rowsize(); }
    inline ptrdiff_t rowsize() const { return m.colsize(); }
    inline StorageType stor() const { return ColMajor; }
    inline T getX() const { return x; }
    inline const GENMATRIX2<Tm>& getM() const { return m; }
    inline void assignToM(MatrixView<TMV_RealType(T)> m0) const
    {
        TMVAssert(isReal(T()));
        TMVAssert(m0.colsize() == colsize() && m0.rowsize() == rowsize());
        m.makeInverse(m0);
        MultXM(x,m0);
    }
    inline void assignToM(MatrixView<TMV_ComplexType(T)> m0) const
    {
        TMVAssert(m0.colsize() == colsize() && m0.rowsize() == rowsize());
        m.makeInverse(m0);
        MultXM(x,m0);
    }
private:
    const T x;
    const GENMATRIX2<Tm>& m;
};

#include "tmv/TMV_AuxQuotXM.h"
#include "tmv/TMV_AuxQuotXMa.h"

