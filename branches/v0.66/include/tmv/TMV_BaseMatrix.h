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

namespace tmv {

    template <class T> 
    class BaseMatrix;

    template <class T> 
    class GenMatrix;

    template <class T, IndexStyle I=CStyle> 
    class ConstMatrixView;

    template <class T, IndexStyle I=CStyle> 
    class MatrixView;

    template <class T, StorageType S=ColMajor, IndexStyle I=CStyle> 
    class Matrix;

    template <class T, int M, int N, StorageType S=ColMajor, IndexStyle I=CStyle> 
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

        virtual size_t colsize() const = 0;
        virtual size_t rowsize() const = 0;
        inline size_t ncols() const 
        { return rowsize(); }
        inline size_t nrows() const 
        { return colsize(); }
        inline bool isSquare() const 
        { return colsize() == rowsize(); }
        TMV_DEPRECATED(bool IsSquare() const)
        { return isSquare(); }

        virtual void assignToM(const MatrixView<RT>& m) const = 0; 
        virtual void assignToM(const MatrixView<CT>& m) const = 0; 

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

        virtual void makeInverse(const MatrixView<T>& minv) const = 0;
        virtual void makeInverseATA(const MatrixView<T>& ata) const = 0;
        virtual bool isSingular() const = 0;
        virtual RT condition() const = 0;
        virtual RT doCondition() const = 0;

        TMV_DEPRECATED(virtual auto_ptr<BaseMatrix<T> > newCopy() const) = 0;
        TMV_DEPRECATED(virtual auto_ptr<BaseMatrix<T> > newView() const) = 0;
        TMV_DEPRECATED(virtual auto_ptr<BaseMatrix<T> > newTranspose() const) = 0;
        TMV_DEPRECATED(virtual auto_ptr<BaseMatrix<T> > newConjugate() const) = 0;
        TMV_DEPRECATED(virtual auto_ptr<BaseMatrix<T> > newAdjoint() const) = 0;
        TMV_DEPRECATED(virtual auto_ptr<BaseMatrix<T> > newInverse() const) = 0;

        TMV_DEPRECATED(T Det() const)
        { return det(); }
        TMV_DEPRECATED(RT LogDet(T* sign=0) const)
        { return logDet(sign); }
        TMV_DEPRECATED(T Trace() const)
        { return trace(); }
        TMV_DEPRECATED(RT Norm() const )
        { return norm(); }
        TMV_DEPRECATED(RT NormF() const)
        { return normF(); }
        TMV_DEPRECATED(RT NormSq(const RT scale = RT(1)) const)
        { return normSq(); }
        TMV_DEPRECATED(RT Norm1() const)
        { return norm1(); }
        TMV_DEPRECATED(RT NormInf() const)
        { return normInf(); }
        TMV_DEPRECATED(RT MaxAbsElement() const)
        { return maxAbsElement(); }
        TMV_DEPRECATED(bool Singular() const)
        { return isSingular(); }
        TMV_DEPRECATED(RT Norm2() const)
        { return norm2(); }
        TMV_DEPRECATED(RT Condition() const)
        { return condition(); }
        TMV_DEPRECATED(RT DoNorm2() const)
        { return doNorm2(); }
        TMV_DEPRECATED(RT DoCondition() const)
        { return doCondition(); }
        TMV_DEPRECATED(void Inverse(const MatrixView<T>& minv) const)
        { makeInverse(minv); }
        TMV_DEPRECATED(void InverseATA(const MatrixView<T>& ata) const)
        { makeInverseATA(ata); }
        TMV_DEPRECATED(auto_ptr<BaseMatrix<T> > NewCopy() const)
        { return newCopy(); }
        TMV_DEPRECATED(auto_ptr<BaseMatrix<T> > NewView() const)
        { return newView(); }
        TMV_DEPRECATED(auto_ptr<BaseMatrix<T> > NewTranspose() const)
        { return newTranspose(); }
        TMV_DEPRECATED(auto_ptr<BaseMatrix<T> > NewConjugate() const)
        { return newConjugate(); }
        TMV_DEPRECATED(auto_ptr<BaseMatrix<T> > NewAdjoint() const)
        { return newAdjoint(); }
        TMV_DEPRECATED(auto_ptr<BaseMatrix<T> > NewInverse() const)
        { return newInverse(); }


        // 
        // I/O: Write
        //

        virtual void write(std::ostream& os) const = 0;
        virtual void write(std::ostream& os, RT thresh) const = 0;

        TMV_DEPRECATED(void Write(std::ostream& os) const)
        { write(os); }
        TMV_DEPRECATED(void Write(std::ostream& os, RT thresh) const)
        { write(os,thresh); }

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

        inline DivHelper() : pdiv(0) {}
        inline DivHelper(const BaseMatrix<T>&) : pdiv(0) {}
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

        void makeInverse(const MatrixView<T>& minv) const
        {
            TMVAssert(minv.colsize() == rowsize());
            TMVAssert(minv.rowsize() == colsize());
            doMakeInverse(minv);
        }

        template <class T1> 
        inline void makeInverse(const MatrixView<T1>& minv) const
        {
            TMVAssert(minv.colsize() == rowsize());
            TMVAssert(minv.rowsize() == colsize());
            doMakeInverse(minv);
        }

        template <class T1, StorageType S, IndexStyle I> 
        inline void makeInverse(Matrix<T1,S,I>& minv) const
        {
            TMVAssert(minv.colsize() == rowsize());
            TMVAssert(minv.rowsize() == colsize());
            doMakeInverse(minv.view());
        }

        inline void makeInverseATA(const MatrixView<T>& ata) const
        { 
            TMVAssert(ata.colsize() == 
                      (rowsize() < colsize() ? rowsize() : colsize()));
            TMVAssert(ata.rowsize() == 
                      (rowsize() < colsize() ? rowsize() : colsize()));
            doMakeInverseATA(ata);
        }

        template <StorageType S, IndexStyle I> 
        inline void makeInverseATA(Matrix<T,S,I>& ata) const
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
        inline void LDivEq(const VectorView<T1>& v) const 
        {
            TMVAssert(colsize() == rowsize());
            TMVAssert(colsize() == v.size());
            doLDivEq(v);
        }

        template <class T1> 
        inline void LDivEq(const MatrixView<T1>& m) const 
        { 
            TMVAssert(colsize() == rowsize());
            TMVAssert(colsize() == m.colsize());
            doLDivEq(m);
        }

        // v * m^-1 -> v
        template <class T1> 
        inline void RDivEq(const VectorView<T1>& v) const 
        { 
            TMVAssert(colsize() == rowsize());
            TMVAssert(colsize() == v.size());
            doRDivEq(v);
        }

        template <class T1> 
        inline void RDivEq(const MatrixView<T1>& m) const 
        { 
            TMVAssert(colsize() == rowsize());
            TMVAssert(colsize() == m.rowsize());
            doRDivEq(m);
        }

        // m^-1 * v1 -> v0
        template <class T1, class T0> 
        inline void LDiv(
            const GenVector<T1>& v1, const VectorView<T0>& v0) const
        { 
            TMVAssert(rowsize() == v0.size());
            TMVAssert(colsize() == v1.size());
            doLDiv(v1,v0);
        }

        template <class T1, class T0> 
        inline void LDiv(
            const GenMatrix<T1>& m1, const MatrixView<T0>& m0) const
        { 
            TMVAssert(rowsize() == m0.colsize());
            TMVAssert(colsize() == m1.colsize());
            TMVAssert(m1.rowsize() == m0.rowsize());
            doLDiv(m1,m0);
        }

        // v1 * m^-1 -> v0
        template <class T1, class T0> 
        inline void RDiv(
            const GenVector<T1>& v1, const VectorView<T0>& v0) const
        { 
            TMVAssert(rowsize() == v1.size());
            TMVAssert(colsize() == v0.size());
            doRDiv(v1,v0);
        }

        template <class T1, class T0> 
        inline void RDiv(
            const GenMatrix<T1>& m1, const MatrixView<T0>& m0) const
        { 
            TMVAssert(rowsize() == m1.rowsize());
            TMVAssert(colsize() == m0.rowsize());
            TMVAssert(m1.colsize() == m0.colsize());
            doRDiv(m1,m0);
        }

        //
        // Division Control
        //

        void divideInPlace() const;
        void saveDiv() const;
        void divideUsing(DivType dt) const;
        void setDiv() const;
        bool divIsSet() const;
        void unsetDiv() const;
        void resetDiv() const;

        bool checkDecomp(std::ostream* fout=0) const;
        bool checkDecomp(const BaseMatrix<T>& m2, std::ostream* fout=0) const;

        TMV_DEPRECATED(void DivideInPlace() const)
        { divideInPlace(); }
        TMV_DEPRECATED(void SaveDiv() const)
        { saveDiv(); }
        TMV_DEPRECATED(void DivideUsing(DivType dt) const)
        { divideUsing(dt); }
        TMV_DEPRECATED(void SetDiv() const)
        { setDiv(); }
        TMV_DEPRECATED(bool DivIsSet() const)
        { return divIsSet(); }
        TMV_DEPRECATED(void UnSetDiv() const)
        { return unsetDiv(); }
        TMV_DEPRECATED(void ReSetDiv() const)
        { return resetDiv(); }
        TMV_DEPRECATED(bool CheckDecomp(std::ostream* fout=0) const)
        { return checkDecomp(fout); }
        TMV_DEPRECATED(bool CheckDecomp(
                const BaseMatrix<T>& m2, std::ostream* fout=0) const)
        { return checkDecomp(m2,fout); }


    protected :

        struct DivImpl;

        mutable DivImpl* pdiv;

        const Divider<T>* getDiv() const;
        void setDiv(Divider<T>*) const;
        DivType getDivType() const;
        bool isDivInPlace() const;
        void doneDiv() const;

        // This is why the divider stuff is implemented using private
        // inheritance.  NewDivider needs to be defined in the 
        // derived class.
        virtual void newDivider() const = 0;
        virtual const BaseMatrix<T>& getMatrix() const = 0;

    private :

        DivHelper(const DivHelper<T>&);
        DivHelper<T>& operator=(const DivHelper<T>&);

        void setupDiv() const;

        T doDet() const;
        RT doLogDet(T* sign) const;
        template <class T1> 
        void doMakeInverse(const MatrixView<T1>& minv) const;
        void doMakeInverseATA(const MatrixView<T>& minv) const;
        bool doIsSingular() const;
        RT doNorm2() const;
        RT doCondition() const;
        void doWrite(std::ostream& os) const;
        void doWrite(std::ostream& os, RT thresh) const;
        template <class T1> 
        void doLDivEq(const VectorView<T1>& v) const;
        template <class T1> 
        void doLDivEq(const MatrixView<T1>& m) const;
        template <class T1> 
        void doRDivEq(const VectorView<T1>& v) const;
        template <class T1> 
        void doRDivEq(const MatrixView<T1>& m) const;
        template <class T1, class T0> 
        void doLDiv(
            const GenVector<T1>& v1, const VectorView<T0>& v0) const;
        template <class T1, class T0> 
        void doLDiv(
            const GenMatrix<T1>& m1, const MatrixView<T0>& m0) const;
        template <class T1, class T0> 
        void doRDiv(
            const GenVector<T1>& v1, const VectorView<T0>& v0) const;
        template <class T1, class T0> 
        void doRDiv(
            const GenMatrix<T1>& m1, const MatrixView<T0>& m0) const;

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
    inline std::ostream& operator<<(std::ostream& os, const BaseMatrix<T>& m)
    { m.write(os); return os; }

    inline std::string TMV_Text(StorageType s)
    {
        return 
            s == RowMajor ? "RowMajor" :
            s == ColMajor ? "ColMajor" :
            s == DiagMajor ? "DiagMajor" :
            s == NoMajor ? "NoMajor" : 
            "Unknown";
    }

    template <class T, StorageType S, IndexStyle I> 
    inline std::string TMV_Text(const Matrix<T,S,I>& )
    {
        return std::string("Matrix<") +
            TMV_Text(T()) + "," + 
            TMV_Text(S) + "," +
            TMV_Text(I) + ">";
    }
    template <class T> 
    inline std::string TMV_Text(const GenMatrix<T>& m)
    {
        return std::string("GenMatrix<") +
            TMV_Text(T()) + "," +
            TMV_Text(m.stor()) + "," +
            TMV_Text(m.ct()) + ">";
    }
    template <class T, IndexStyle I> 
    inline std::string TMV_Text(const ConstMatrixView<T,I>& m)
    {
        return std::string("ConstMatrixView<") +
            TMV_Text(T()) + "," +
            TMV_Text(m.stor()) + "," +
            TMV_Text(I) + "," +
            TMV_Text(m.ct()) + ">";
    }
    template <class T, IndexStyle I> 
    inline std::string TMV_Text(const MatrixView<T,I>& m)
    {
        return std::string("MatrixView<") +
            TMV_Text(T()) + "," +
            TMV_Text(m.stor()) + "," +
            TMV_Text(I) + "," +
            TMV_Text(m.ct()) + ">";
    }

} // namespace tmv

#endif
