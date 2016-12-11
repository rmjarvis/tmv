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


//---------------------------------------------------------------------------
//
// This file defines the TMV BaseMatrix class.
//
// This base class defines some of the things that all
// matrices need to be able to do, as well as some of the
// arithmetic operations (those that return a Vector).
// This should be used as the base class for generic
// matrices as well as any special ones (eg. sparse,
// symmetric, etc.)
//
//

#ifndef TMV_BaseMatrix_H
#define TMV_BaseMatrix_H

#include "tmv/TMV_Base.h"
#include "tmv/TMV_BaseVector.h"
#include "tmv/TMV_IOStyle.h"

namespace tmv {

    template <typename T>
    class BaseMatrix;

    template <typename T>
    class GenMatrix;

    template <typename T, int A=0>
    class ConstMatrixView;

    template <typename T, int A=0>
    class MatrixView;

    template <typename T, int A=0>
    class Matrix;

    template <typename T, ptrdiff_t M, ptrdiff_t N, int A=0>
    class SmallMatrix;

    template <typename T, ptrdiff_t M, ptrdiff_t N>
    class SmallMatrixComposite;

    template <typename T>
    class Divider;

    template <typename T>
    struct AssignableToMatrix
    {
        typedef TMV_RealType(T) RT;
        typedef TMV_ComplexType(T) CT;

        virtual ptrdiff_t colsize() const = 0;
        virtual ptrdiff_t rowsize() const = 0;
        inline ptrdiff_t ncols() const
        { return rowsize(); }
        inline ptrdiff_t nrows() const
        { return colsize(); }
        inline bool isSquare() const
        { return colsize() == rowsize(); }

        virtual void assignToM(MatrixView<RT> m) const = 0;
        virtual void assignToM(MatrixView<CT> m) const = 0;

        virtual inline ~AssignableToMatrix() {}
    };

    template <typename T>
    class BaseMatrix : virtual public AssignableToMatrix<T>
    {
    public :
        typedef TMV_RealType(T) RT;

        //
        // Access Functions
        //

        using AssignableToMatrix<T>::colsize;
        using AssignableToMatrix<T>::rowsize;

        //
        // Functions of Matrix
        //

        virtual T det() const = 0;
        virtual RT logDet(T* sign=0) const = 0;
        virtual T trace() const = 0;
        virtual T sumElements() const = 0;
        virtual RT sumAbsElements() const = 0;
        virtual RT sumAbs2Elements() const = 0;

        virtual RT norm() const  = 0;
        virtual RT normSq(const RT scale = RT(1)) const = 0;
        virtual RT normF() const  = 0;
        virtual RT norm1() const = 0;
        virtual RT norm2() const  = 0;
        virtual RT doNorm2() const  = 0;
        virtual RT normInf() const = 0;
        virtual RT maxAbsElement() const = 0;
        virtual RT maxAbs2Element() const = 0;

        //
        // I/O: Write
        //

        virtual void write(const TMV_Writer& writer) const = 0;

        virtual inline ~BaseMatrix() {}

    }; // BaseMatrix

    template <typename T>
    class DivHelper : virtual public AssignableToMatrix<T>
    {
    public:

        typedef TMV_RealType(T) RT;

        //
        // Constructors
        //

        DivHelper();
        // Cannot do this inline, since need to delete pdiv,
        // and I only define DivImpl in BaseMatrix.cpp.
        virtual ~DivHelper();

        using AssignableToMatrix<T>::colsize;
        using AssignableToMatrix<T>::rowsize;

        T det() const
        {
            TMVAssert(rowsize() == colsize());
            return doDet();
        }

        RT logDet(T* sign) const
        {
            TMVAssert(rowsize() == colsize());
            return doLogDet(sign);
        }

        void makeInverse(MatrixView<T> minv) const
        {
            TMVAssert(minv.colsize() == rowsize());
            TMVAssert(minv.rowsize() == colsize());
            doMakeInverse(minv);
        }

        template <typename T1>
        inline void makeInverse(MatrixView<T1> minv) const
        {
            TMVAssert(minv.colsize() == rowsize());
            TMVAssert(minv.rowsize() == colsize());
            doMakeInverse(minv);
        }

        template <typename T1, int A>
        inline void makeInverse(Matrix<T1,A>& minv) const
        {
            TMVAssert(minv.colsize() == rowsize());
            TMVAssert(minv.rowsize() == colsize());
            doMakeInverse(minv.view());
        }

        inline void makeInverseATA(MatrixView<T> ata) const
        {
            TMVAssert(ata.colsize() ==
                      (rowsize() < colsize() ? rowsize() : colsize()));
            TMVAssert(ata.rowsize() ==
                      (rowsize() < colsize() ? rowsize() : colsize()));
            doMakeInverseATA(ata);
        }

        template <int A>
        inline void makeInverseATA(Matrix<T,A>& ata) const
        {
            TMVAssert(ata.colsize() ==
                      (rowsize() < colsize() ? rowsize() : colsize()));
            TMVAssert(ata.rowsize() ==
                      (rowsize() < colsize() ? rowsize() : colsize()));
            doMakeInverseATA(ata.view());
        }

        inline bool isSingular() const
        { return doIsSingular(); }

        inline RT norm2() const
        {
            TMVAssert(divIsSet() && getDivType() == SV);
            return doNorm2();
        }

        inline RT condition() const
        {
            TMVAssert(divIsSet() && getDivType() == SV);
            return doCondition();
        }

        // m^-1 * v -> v
        template <typename T1>
        inline void LDivEq(VectorView<T1> v) const
        {
            TMVAssert(colsize() == rowsize());
            TMVAssert(colsize() == v.size());
            doLDivEq(v);
        }

        template <typename T1>
        inline void LDivEq(MatrixView<T1> m) const
        {
            TMVAssert(colsize() == rowsize());
            TMVAssert(colsize() == m.colsize());
            doLDivEq(m);
        }

        // v * m^-1 -> v
        template <typename T1>
        inline void RDivEq(VectorView<T1> v) const
        {
            TMVAssert(colsize() == rowsize());
            TMVAssert(colsize() == v.size());
            doRDivEq(v);
        }

        template <typename T1>
        inline void RDivEq(MatrixView<T1> m) const
        {
            TMVAssert(colsize() == rowsize());
            TMVAssert(colsize() == m.rowsize());
            doRDivEq(m);
        }

        // m^-1 * v1 -> v0
        template <typename T1, typename T0>
        inline void LDiv(
            const GenVector<T1>& v1, VectorView<T0> v0) const
        {
            TMVAssert(rowsize() == v0.size());
            TMVAssert(colsize() == v1.size());
            doLDiv(v1,v0);
        }

        template <typename T1, typename T0>
        inline void LDiv(
            const GenMatrix<T1>& m1, MatrixView<T0> m0) const
        {
            TMVAssert(rowsize() == m0.colsize());
            TMVAssert(colsize() == m1.colsize());
            TMVAssert(m1.rowsize() == m0.rowsize());
            doLDiv(m1,m0);
        }

        // v1 * m^-1 -> v0
        template <typename T1, typename T0>
        inline void RDiv(
            const GenVector<T1>& v1, VectorView<T0> v0) const
        {
            TMVAssert(rowsize() == v1.size());
            TMVAssert(colsize() == v0.size());
            doRDiv(v1,v0);
        }

        template <typename T1, typename T0>
        inline void RDiv(
            const GenMatrix<T1>& m1, MatrixView<T0> m0) const
        {
            TMVAssert(rowsize() == m1.rowsize());
            TMVAssert(colsize() == m0.rowsize());
            TMVAssert(m1.colsize() == m0.colsize());
            doRDiv(m1,m0);
        }

        //
        // Division Control
        //

        void divideUsing(DivType dt) const;

        void divideInPlace() const;
        void dontDivideInPlace() const;
        void saveDiv() const;
        void dontSaveDiv() const;

        // setDiv is defined in the derived class.
        virtual void setDiv() const = 0;
        void unsetDiv() const;
        void resetDiv() const;

        DivType getDivType() const;
        bool divIsInPlace() const;
        bool divIsSaved() const;
        bool divIsSet() const;

        bool checkDecomp(std::ostream* fout=0) const;
        bool checkDecomp(const BaseMatrix<T>& m2, std::ostream* fout=0) const;

    protected :

        void doneDiv() const;
        const Divider<T>* getDiv() const;
        void resetDivType() const;

        // Two more that need to be defined in the derived class:
        virtual const BaseMatrix<T>& getMatrix() const = 0;

        mutable auto_ptr<Divider<T> > divider;
        mutable DivType divtype;

    private :

        DivHelper(const DivHelper<T>&);
        DivHelper<T>& operator=(const DivHelper<T>&);

        T doDet() const;
        RT doLogDet(T* sign) const;
        template <typename T1>
        void doMakeInverse(MatrixView<T1> minv) const;
        void doMakeInverseATA(MatrixView<T> minv) const;
        bool doIsSingular() const;
        RT doNorm2() const;
        RT doCondition() const;
        template <typename T1>
        void doLDivEq(VectorView<T1> v) const;
        template <typename T1>
        void doLDivEq(MatrixView<T1> m) const;
        template <typename T1>
        void doRDivEq(VectorView<T1> v) const;
        template <typename T1>
        void doRDivEq(MatrixView<T1> m) const;
        template <typename T1, typename T0>
        void doLDiv(
            const GenVector<T1>& v1, VectorView<T0> v0) const;
        template <typename T1, typename T0>
        void doLDiv(
            const GenMatrix<T1>& m1, MatrixView<T0> m0) const;
        template <typename T1, typename T0>
        void doRDiv(
            const GenVector<T1>& v1, VectorView<T0> v0) const;
        template <typename T1, typename T0>
        void doRDiv(
            const GenMatrix<T1>& m1, MatrixView<T0> m0) const;

    }; // DivHelper

    //
    // Functions of Matrices:
    //

    template <typename T>
    inline T Det(const BaseMatrix<T>& m)
    { return m.det(); }

    template <typename T>
    inline TMV_RealType(T) LogDet(const BaseMatrix<T>& m)
    { return m.logDet(); }

    template <typename T>
    inline T Trace(const BaseMatrix<T>& m)
    { return m.trace(); }

    template <typename T>
    inline T SumElements(const BaseMatrix<T>& m)
    { return m.sumElements(); }

    template <typename T>
    inline TMV_RealType(T) SumAbsElements(const BaseMatrix<T>& m)
    { return m.sumAbsElements(); }

    template <typename T>
    inline TMV_RealType(T) SumAbs2Elements(const BaseMatrix<T>& m)
    { return m.sumAbs2Elements(); }

    template <typename T>
    inline TMV_RealType(T) Norm(const BaseMatrix<T>& m)
    { return m.norm(); }

    template <typename T>
    inline TMV_RealType(T) NormSq(const BaseMatrix<T>& m)
    { return m.normSq(); }

    template <typename T>
    inline TMV_RealType(T) NormF(const BaseMatrix<T>& m)
    { return m.normF(); }

    template <typename T>
    inline TMV_RealType(T) Norm1(const BaseMatrix<T>& m)
    { return m.norm1(); }

    template <typename T>
    inline TMV_RealType(T) Norm2(const BaseMatrix<T>& m)
    { return m.norm2(); }

    template <typename T>
    inline TMV_RealType(T) NormInf(const BaseMatrix<T>& m)
    { return m.normInf(); }

    template <typename T>
    inline TMV_RealType(T) MaxAbsElement(const BaseMatrix<T>& m)
    { return m.maxAbsElement(); }

    template <typename T>
    inline TMV_RealType(T) MaxAbs2Element(const BaseMatrix<T>& m)
    { return m.maxAbs2Element(); }


    //
    // I/O
    //

    template <typename T>
    inline std::ostream& operator<<(
        const TMV_Writer& writer, const BaseMatrix<T>& m)
    { m.write(writer); return writer.getos(); }

    template <typename T>
    inline std::ostream& operator<<(
        std::ostream& os, const BaseMatrix<T>& m)
    { return os << IOStyle() << m; }


    template <typename T, int A>
    inline std::string TMV_Text(const Matrix<T,A>& )
    {
        return std::string("Matrix<") +
            TMV_Text(T()) + "," + Attrib<A>::text() + ">";
    }
    template <typename T>
    inline std::string TMV_Text(const GenMatrix<T>& )
    {
        return std::string("GenMatrix<") + TMV_Text(T()) + ">";
    }
    template <typename T, int A>
    inline std::string TMV_Text(const ConstMatrixView<T,A>& )
    {
        return std::string("ConstMatrixView<") +
            TMV_Text(T()) + "," + Attrib<A>::text() + ">";
    }
    template <typename T, int A>
    inline std::string TMV_Text(const MatrixView<T,A>& )
    {
        return std::string("MatrixView<") +
            TMV_Text(T()) + "," + Attrib<A>::text() + ">";
    }

} // namespace tmv

#endif
