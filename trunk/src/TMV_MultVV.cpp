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
#include "tmv/TMV_MultVV.h"
#include "tmv/TMV_Vector.h"

namespace tmv {

    //
    // MultVV
    //

    template <class V1, class V2>
    static typename V2::value_type DoMultVV(const V1& v1, const V2& v2)
    {
        if (v1.step() == 1 && v2.step() == 1)
            return InlineMultVV(v1.unitView(),v2.unitView());
        else
            return InlineMultVV(v1,v2);
    }

#ifdef BLAS
#ifndef BLASNORETURN
#ifdef TMV_INST_DOUBLE
    double DoMultVV(
        const ConstVectorView<double>& v1, const ConstVectorView<double>& v2) 
    { 
        TMVAssert(v1.size()==v2.size());
        int n=v1.size();
        if (n == 0) return 0.;
        int s1=v1.step();
        int s2=v2.step();
        const double* v1p = v1.cptr();
        if (s1 < 0) v1p += (n-1)*s1;
        const double* v2p = v2.cptr();
        if (s2 < 0) v2p += (n-1)*s2;
        return BLASNAME(ddot) (
            BLASV(n),BLASP(v1p),BLASV(s1),BLASP(v2p),BLASV(s2));
    }
    std::complex<double> DoMultVV(
        const ConstVectorView<std::complex<double> >& v1, 
        const ConstVectorView<std::complex<double> >& v2) 
    {
        TMVAssert(v1.size()==v2.size());
        int n=v1.size();
        if (n == 0) return 0.;
        int s1=v1.step();
        int s2=v2.step();
        const std::complex<double>* v1p = v1.cptr();
        if (s1 < 0) v1p += (n-1)*s1;
        const std::complex<double>* v2p = v2.cptr();
        if (s2 < 0) v2p += (n-1)*s2;
        std::complex<double> res;
        BLASZDOTSET( res, BLASZDOTNAME(zdotu) (
                BLASZDOT1(BLASP(&res))
                BLASV(n),BLASP(v2p),BLASV(s2),
                BLASP(v1p),BLASV(s1)
                BLASZDOT2(BLASP(&res)) ));
        return res;
    }
    std::complex<double> DoMultVV(
        const ConstVectorView<std::complex<double>,UNKNOWN,true>& v1,
        const ConstVectorView<std::complex<double> >& v2) 
    {
        TMVAssert(v1.size()==v2.size());
        int n=v1.size();
        if (n == 0) return 0.;
        int s1=v1.step();
        int s2=v2.step();
        const std::complex<double>* v1p = v1.cptr();
        if (s1 < 0) v1p += (n-1)*s1;
        const std::complex<double>* v2p = v2.cptr();
        if (s2 < 0) v2p += (n-1)*s2;
        std::complex<double> res;
        BLASZDOTSET( res, BLASZDOTNAME(zdotc) (
                BLASZDOT1(BLASP(&res))
                BLASV(n),BLASP(v1p),BLASV(s1),
                BLASP(v2p),BLASV(s2)
                BLASZDOT2(BLASP(&res)) ));
        return res;
    }
#ifdef TMV_INST_MIX
    std::complex<double> DoMultVV(
        const ConstVectorView<double>& v1, 
        const ConstVectorView<std::complex<double> >& v2) 
    {
        TMVAssert(v1.size()==v2.size());
        int n=v1.size();
        if (n == 0) return 0.F;
        int s1=v1.step();
        int s2=2*v2.step();
        const double* v1p = v1.cptr();
        if (s1 < 0) v1p += (n-1)*s1;
        const std::complex<double>* v2p = v2.cptr();
        if (s2 < 0) v2p += (n-1)*v2.step();
        double resr = BLASNAME(ddot) (
            BLASV(n),BLASP(v1p),BLASV(s1),BLASP((double*)v2p),BLASV(s2));
        double resi = BLASNAME(ddot) (
            BLASV(n),BLASP(v1p),BLASV(s1),BLASP((double*)v2p+1),BLASV(s2));
        return std::complex<double>(resr,resi);
    }
#endif
#endif
#ifdef TMV_INST_FLOAT
    float DoMultVV(
        const ConstVectorView<float>& v1, const ConstVectorView<float>& v2) 
    {
        TMVAssert(v1.size()==v2.size());
        int n=v1.size();
        if (n == 0) return 0.F;
        int s1=v1.step();
        int s2=v2.step();
        const float* v1p = v1.cptr();
        if (s1 < 0) v1p += (n-1)*s1;
        const float* v2p = v2.cptr();
        if (s2 < 0) v2p += (n-1)*s2;
        return BLASNAME(sdot) (
            BLASV(n),BLASP(v1p),BLASV(s1),BLASP(v2p),BLASV(s2));
    }
    std::complex<float> DoMultVV(
        const ConstVectorView<std::complex<float> >& v1, 
        const ConstVectorView<std::complex<float> >& v2) 
    {
        TMVAssert(v1.size()==v2.size());
        int n=v1.size();
        if (n == 0) return 0.F;
        int s1=v1.step();
        int s2=v2.step();
        const std::complex<float>* v1p = v1.cptr();
        if (s1 < 0) v1p += (n-1)*s1;
        const std::complex<float>* v2p = v2.cptr();
        if (s2 < 0) v2p += (n-1)*s2;
        std::complex<float> res;
        BLASZDOTSET( res, BLASZDOTNAME(cdotu) (
                BLASZDOT1(BLASP(&res))
                BLASV(n),BLASP(v2p),BLASV(s2),
                BLASP(v1p),BLASV(s1)
                BLASZDOT2(BLASP(&res)) ));
        return res;
    }
    std::complex<float> DoMultVV(
        const ConstVectorView<std::complex<float>,UNKNOWN,true>& v2,
        const ConstVectorView<std::complex<float> >& v1) 
    {
        TMVAssert(v1.size()==v2.size());
        int n=v1.size();
        if (n == 0) return 0.F;
        int s1=v1.step();
        int s2=v2.step();
        const std::complex<float>* v1p = v1.cptr();
        if (s1 < 0) v1p += (n-1)*s1;
        const std::complex<float>* v2p = v2.cptr();
        if (s2 < 0) v2p += (n-1)*s2;
        std::complex<float> res;
        BLASZDOTSET( res, BLASZDOTNAME(cdotc) (
                BLASZDOT1(BLASP(&res))
                BLASV(n),BLASP(v1p),BLASV(s1),
                BLASP(v2p),BLASV(s2)
                BLASZDOT2(BLASP(&res)) ));
        return res;
    }
#ifdef TMV_INST_MIX
    std::complex<float> DoMultVV(
        const ConstVectorView<float>& v1, 
        const ConstVectorView<std::complex<float> >& v2) 
    {
        TMVAssert(v1.size()==v2.size());
        int n=v1.size();
        if (n == 0) return 0.F;
        int s1=v1.step();
        int s2=2*v2.step();
        const float* v1p = v1.cptr();
        if (s1 < 0) v1p += (n-1)*s1;
        const std::complex<float>* v2p = v2.cptr();
        if (s2 < 0) v2p += (n-1)*v2.step();
        float resr = BLASNAME(sdot) (
            BLASV(n),BLASP(v1p),BLASV(s1),BLASP((float*)v2p),BLASV(s2));
        float resi = BLASNAME(sdot) (
            BLASV(n),BLASP(v1p),BLASV(s1),BLASP((float*)v2p+1),BLASV(s2));
        return std::complex<float>(resr,resi);
    }
#endif
#endif
#endif // BLASNORETURN
#endif // BLAS

    template <class T1, bool C1, class T2>
    T2 InstMultVV(
        const ConstVectorView<T1,UNKNOWN,C1>& v1, const ConstVectorView<T2>& v2)
    { return DoMultVV(v1,v2); }


#define InstFile "TMV_MultVV.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


