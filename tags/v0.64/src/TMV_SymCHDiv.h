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
#ifndef TMV_SymCHDiv_H
#define TMV_SymCHDiv_H

#include "tmv/TMV_BaseSymMatrix.h"

namespace tmv {

    template <class T1, class T2> 
    void CH_LDivEq(const GenSymMatrix<T1>& L, const MatrixView<T2>& m);

    template <class T1, class T2> 
    void CH_RDivEq(const GenSymMatrix<T1>& L, const MatrixView<T2>& m);

    template <class T1, class T2> 
    void CH_Inverse(const GenSymMatrix<T1>& LLx, const SymMatrixView<T2>& sinv);

#ifndef NOTHROW
    template <class T> 
    class NonPosDefHermMatrix : public NonPosDef
    {
    public:

        mutable auto_ptr<HermMatrix<T> > A;

        NonPosDefHermMatrix(const GenSymMatrix<T>& _A) :
            NonPosDef("HermMatrix Cholesky decmposition."), 
            A(new HermMatrix<T>(_A)) {}
        NonPosDefHermMatrix(const NonPosDefHermMatrix<T>& rhs) :
            A(rhs.A) {}
        ~NonPosDefHermMatrix() throw() {}

        void write(std::ostream& os) const throw()
        {
            NonPosDef::write(os);
            os<<"The partially decomposed matrix is \n"<<*A<<std::endl;
        }
    };

    template <class T> 
    class NonPosDefHermMatrix2 : public NonPosDefHermMatrix<T>
    {
    public:
        mutable auto_ptr<HermMatrix<T> > A0;

        NonPosDefHermMatrix2(
            const GenSymMatrix<T>& _A, const GenSymMatrix<T>& _A0) :
            NonPosDefHermMatrix<T>(_A), A0(new HermMatrix<T>(_A0)) {}
        NonPosDefHermMatrix2(const NonPosDefHermMatrix2<T>& rhs) :
            NonPosDefHermMatrix<T>(rhs), A0(rhs.A0) {}
        ~NonPosDefHermMatrix2() throw() {}

        void write(std::ostream& os) const throw()
        {
            NonPosDefHermMatrix<T>::write(os);
            os<<"The original matrix was \n"<<*A0<<std::endl;
        }
    };
#endif

}

#endif
