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

#include "TMV_Blas.h"
#include "tmv/TMV_MultMM.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_ProdXM.h"
#include "tmv/TMV_SumMM.h"
#include "tmv/TMV_Vector.h"

#ifdef BLAS
#include "TMV_MultMM_Blas.h"
#endif

namespace tmv {

    template <bool add, class T1, int C1, class T2, int C2, class T3>
    static void DoMultMM(
        const T3 x,
        const ConstMatrixView<T1,C1>& m1,
        const ConstMatrixView<T2,C2>& m2, MatrixView<T3,ColMajor> m3)
    {
        typedef typename Traits<T3>::real_type RT;
        if (x == RT(1))
            InlineMultMM<add>(Scaling<1,RT>(),m1,m2,m3);
        else if (x == RT(-1))
            InlineMultMM<add>(Scaling<-1,RT>(),m1,m2,m3);
        else if (x == RT(0))
            Maybe<!add>::zero(m3); 
        else if (TMV_IMAG(x) == RT(0))
            InlineMultMM<add>(Scaling<0,RT>(TMV_REAL(x)),m1,m2,m3);
        else
            InlineMultMM<add>(Scaling<0,T3>(x),m1,m2,m3);
    }

#ifdef BLAS
#ifdef TMV_INST_DOUBLE
    template <bool add>
    static void DoMultMM(
        const double x,
        const ConstMatrixView<double,ColMajor>& m1,
        const ConstMatrixView<double,ColMajor>& m2,
        MatrixView<double,ColMajor> m3)
    { BlasMultMM(x,m1,m2,add?1:0,m3); }
    template <bool add, class T1, int C1, class T2, int C2>
    static void DoMultMM(
        const std::complex<double> x,
        const ConstMatrixView<T1,C1>& m1,
        const ConstMatrixView<T2,C2>& m2,
        MatrixView<std::complex<double>,ColMajor> m3)
    { BlasMultMM(x,m1,m2,add?1:0,m3); }
#endif // TMV_INST_DOUBLE
#ifdef TMV_INST_FLOAT
    template <bool add>
    static void DoMultMM(
        const float x,
        const ConstMatrixView<float,ColMajor>& m1,
        const ConstMatrixView<float,ColMajor>& m2,
        MatrixView<float,ColMajor> m3)
    { BlasMultMM(x,m1,m2,add?1:0,m3); }
    template <bool add, class T1, int C1, class T2, int C2>
    static void DoMultMM(
        const std::complex<float> x,
        const ConstMatrixView<T1,C1>& m1,
        const ConstMatrixView<T2,C2>& m2,
        MatrixView<std::complex<float>,ColMajor> m3)
    { BlasMultMM(x,m1,m2,add?1:0,m3); }
#endif // TMV_INST_FLOAT
#endif // BLAS

    template <bool add, class T1, int C1, class T2, int C2, class T3>
    void DoInstMultMM(
        const T3 x,
        const ConstMatrixView<T1,C1>& m1,
        const ConstMatrixView<T2,C2>& m2, MatrixView<T3,ColMajor> m3)
    {
        //std::cout<<"Start MultMM_CCC:\n";
        //std::cout<<"x = "<<x<<std::endl;
        //std::cout<<"m1 = "<<TMV_Text(m1)<<"  "<<m1<<std::endl;
        //std::cout<<"m2 = "<<TMV_Text(m2)<<"  "<<m2<<std::endl;
        //std::cout<<"m3 = "<<TMV_Text(m3)<<"  "<<m3<<std::endl;
        //Matrix<T3> m3x = m3;
        //if (!add) m3x.setZero();
        //for(size_t j=0;j<m3.rowsize();++j) 
            //InstAddMultMV(x,m1.xView(),m2.col(j).xView(),m3x.col(j).xView());
        //std::cout<<"m3x = "<<m3x<<std::endl;

        if (m3.colsize() > 0 && m3.rowsize() > 0) {
            if (m1.rowsize() == 0) Maybe<!add>::zero(m3);
            else DoMultMM<add>(x,m1,m2,m3); 
        }
        //std::cout<<"m3 => "<<m3<<std::endl;
        //typename Traits<T3>::real_type normdiff = Norm(m3-m3x);
        //std::cout<<"Norm(diff) = "<<normdiff<<std::endl;
        //if (normdiff > 1.e-3) abort();
    }
 

#define InstFile "TMV_MultMM_CCC.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


