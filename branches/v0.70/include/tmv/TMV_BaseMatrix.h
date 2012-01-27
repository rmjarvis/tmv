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

    template <class T> 
    class BaseMatrix;

    template <class T> 
    class GenMatrix;

    template <class T, int A=0>
    class ConstMatrixView;

    template <class T, int A=0>
    class MatrixView;

    template <class T, int A=0>
    class Matrix;

    template <class T, int M, int N, int A=0>
    class SmallMatrix;

    template <class T, int M, int N> 
    class SmallMatrixComposite;

    template <class T> 
    class Divider;

    template <class T> 
    struct AssignableToMatrix
    {
        typedef TMV_RealType(T) RT;
        typedef TMV_ComplexType(T) CT;

        virtual int colsize() const = 0;
        virtual int rowsize() const = 0;
        inline int ncols() const 
        { return rowsize(); }
        inline int nrows() const 
        { return colsize(); }
        inline bool isSquare() const 
        { return colsize() == rowsize(); }

        virtual void assignToM(MatrixView<RT> m) const = 0; 
        virtual void assignToM(MatrixView<CT> m) const = 0; 

        virtual inline ~AssignableToMatrix() {}
    };

    template <class T> 
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

    template <class T> 
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

        template <class T1> 
        inline void makeInverse(MatrixView<T1> minv) const
        {
            TMVAssert(minv.colsize() == rowsize());
            TMVAssert(minv.rowsize() == colsize());
            doMakeInverse(minv);
        }

        template <class T1, int A>
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
        template <class T1> 
        inline void LDivEq(VectorView<T1> v) const 
        {
            TMVAssert(colsize() == rowsize());
            TMVAssert(colsize() == v.size());
            doLDivEq(v);
        }

        template <class T1> 
        inline void LDivEq(MatrixView<T1> m) const 
        { 
            TMVAssert(colsize() == rowsize());
            TMVAssert(colsize() == m.colsize());
            doLDivEq(m);
        }

        // v * m^-1 -> v
        template <class T1> 
        inline void RDivEq(VectorView<T1> v) const 
        { 
            TMVAssert(colsize() == rowsize());
            TMVAssert(colsize() == v.size());
            doRDivEq(v);
        }

        template <class T1> 
        inline void RDivEq(MatrixView<T1> m) const 
        { 
            TMVAssert(colsize() == rowsize());
            TMVAssert(colsize() == m.rowsize());
            doRDivEq(m);
        }

        // m^-1 * v1 -> v0
        template <class T1, class T0> 
        inline void LDiv(
            const GenVector<T1>& v1, VectorView<T0> v0) const
        { 
            TMVAssert(rowsize() == v0.size());
            TMVAssert(colsize() == v1.size());
            doLDiv(v1,v0);
        }

        template <class T1, class T0> 
        inline void LDiv(
            const GenMatrix<T1>& m1, MatrixView<T0> m0) const
        { 
            TMVAssert(rowsize() == m0.colsize());
            TMVAssert(colsize() == m1.colsize());
            TMVAssert(m1.rowsize() == m0.rowsize());
            doLDiv(m1,m0);
        }

        // v1 * m^-1 -> v0
        template <class T1, class T0> 
        inline void RDiv(
            const GenVector<T1>& v1, VectorView<T0> v0) const
        { 
            TMVAssert(rowsize() == v1.size());
            TMVAssert(colsize() == v0.size());
            doRDiv(v1,v0);
        }

        template <class T1, class T0> 
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

        mutable std::auto_ptr<Divider<T> > divider;
        mutable DivType divtype;

    private :

        DivHelper(const DivHelper<T>&);
        DivHelper<T>& operator=(const DivHelper<T>&);

        T doDet() const;
        RT doLogDet(T* sign) const;
        template <class T1> 
        void doMakeInverse(MatrixView<T1> minv) const;
        void doMakeInverseATA(MatrixView<T> minv) const;
        bool doIsSingular() const;
        RT doNorm2() const;
        RT doCondition() const;
        template <class T1> 
        void doLDivEq(VectorView<T1> v) const;
        template <class T1> 
        void doLDivEq(MatrixView<T1> m) const;
        template <class T1> 
        void doRDivEq(VectorView<T1> v) const;
        template <class T1> 
        void doRDivEq(MatrixView<T1> m) const;
        template <class T1, class T0> 
        void doLDiv(
            const GenVector<T1>& v1, VectorView<T0> v0) const;
        template <class T1, class T0> 
        void doLDiv(
            const GenMatrix<T1>& m1, MatrixView<T0> m0) const;
        template <class T1, class T0> 
        void doRDiv(
            const GenVector<T1>& v1, VectorView<T0> v0) const;
        template <class T1, class T0> 
        void doRDiv(
            const GenMatrix<T1>& m1, MatrixView<T0> m0) const;

    }; // DivHelper

    //
    // Functions of Matrices:
    //

    template <class T> 
    inline T Det(const BaseMatrix<T>& m)
    { return m.det(); }

    template <class T> 
    inline TMV_RealType(T) LogDet(const BaseMatrix<T>& m)
    { return m.logDet(); }

    template <class T> 
    inline T Trace(const BaseMatrix<T>& m)
    { return m.trace(); }

    template <class T> 
    inline T SumElements(const BaseMatrix<T>& m)
    { return m.sumElements(); }

    template <class T> 
    inline TMV_RealType(T) SumAbsElements(const BaseMatrix<T>& m)
    { return m.sumAbsElements(); }

    template <class T> 
    inline TMV_RealType(T) SumAbs2Elements(const BaseMatrix<T>& m)
    { return m.sumAbs2Elements(); }

    template <class T> 
    inline TMV_RealType(T) Norm(const BaseMatrix<T>& m)
    { return m.norm(); }

    template <class T> 
    inline TMV_RealType(T) NormSq(const BaseMatrix<T>& m)
    { return m.normSq(); }

    template <class T> 
    inline TMV_RealType(T) NormF(const BaseMatrix<T>& m)
    { return m.normF(); }

    template <class T> 
    inline TMV_RealType(T) Norm1(const BaseMatrix<T>& m)
    { return m.norm1(); }

    template <class T> 
    inline TMV_RealType(T) Norm2(const BaseMatrix<T>& m)
    { return m.norm2(); }

    template <class T> 
    inline TMV_RealType(T) NormInf(const BaseMatrix<T>& m)
    { return m.normInf(); }

    template <class T> 
    inline TMV_RealType(T) MaxAbsElement(const BaseMatrix<T>& m)
    { return m.maxAbsElement(); }

    template <class T> 
    inline TMV_RealType(T) MaxAbs2Element(const BaseMatrix<T>& m)
    { return m.maxAbs2Element(); }


    //
    // I/O
    //

    template <class T>
    inline std::ostream& operator<<(
        const TMV_Writer& writer, const BaseMatrix<T>& m)
    { m.write(writer); return writer.getos(); }

    template <class T> 
    inline std::ostream& operator<<(
        std::ostream& os, const BaseMatrix<T>& m)
    { return os << IOStyle() << m; }


    template <class T, int A>
    inline std::string TMV_Text(const Matrix<T,A>& )
    {
        return std::string("Matrix<") +
            TMV_Text(T()) + "," + Attrib<A>::text() + ">";
    }
    template <class T> 
    inline std::string TMV_Text(const GenMatrix<T>& )
    {
        return std::string("GenMatrix<") + TMV_Text(T()) + ">";
    }
    template <class T, int A>
    inline std::string TMV_Text(const ConstMatrixView<T,A>& )
    {
        return std::string("ConstMatrixView<") +
            TMV_Text(T()) + "," + Attrib<A>::text() + ">";
    }
    template <class T, int A>
    inline std::string TMV_Text(const MatrixView<T,A>& )
    {
        return std::string("MatrixView<") +
            TMV_Text(T()) + "," + Attrib<A>::text() + ">";
    }

} // namespace tmv

#endif
