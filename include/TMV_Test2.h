///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// The Template Matrix/Vector Library for C++ was created by Mike Jarvis     //
// Copyright (C) 2007                                                        //
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
// along with this program in the file gpl.txt.                              //
//                                                                           //
// If not, write to:                                                         //
// The Free Software Foundation, Inc.                                        //
// 51 Franklin Street, Fifth Floor,                                          //
// Boston, MA  02110-1301, USA.                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#define ALIASOK

template <class T> void TestBandMatrix();
template <class T> void TestAllBandDiv();
template <class T> void TestSymMatrix();
template <class T> void TestAllSymDiv();
template <class T> void TestSymBandMatrix();
template <class T> void TestAllSymBandDiv();

template <class T> void TestBandMatrixArith_A();
template <class T> void TestBandMatrixArith_B1();
template <class T> void TestBandMatrixArith_B2();
template <class T> void TestBandMatrixArith_C1();
template <class T> void TestBandMatrixArith_C2();
template <class T> void TestBandMatrixArith_D1();
template <class T> void TestBandMatrixArith_D2();
template <class T> void TestBandDiv_A(tmv::DivType dt);
template <class T> void TestBandDiv_B(tmv::DivType dt);
template <class T> void TestBandDiv_C(tmv::DivType dt);
template <class T> void TestBandDiv_D(tmv::DivType dt);

enum PosDefCode { PosDef, InDef, Sing };
inline std::string PDLabel(PosDefCode pdc)
{
  if (pdc == PosDef) return "Positive Definite";
  else if (pdc == InDef) return "Indefinite";
  else return "Singular";
}

template <class T> void TestSymMatrixArith_A();
template <class T> void TestSymMatrixArith_B1();
template <class T> void TestSymMatrixArith_B2();
template <class T> void TestSymMatrixArith_C1();
template <class T> void TestSymMatrixArith_C2();
template <class T> void TestSymMatrixArith_D1();
template <class T> void TestSymMatrixArith_D2();
template <class T> void TestSymMatrixArith_E1();
template <class T> void TestSymMatrixArith_E2();
template <class T> void TestSymDiv_A(tmv::DivType dt, PosDefCode pc);
template <class T> void TestSymDiv_B(tmv::DivType dt, PosDefCode pc);
template <class T> void TestSymDiv_C(tmv::DivType dt, PosDefCode pc);
template <class T> void TestSymDiv_D(tmv::DivType dt, PosDefCode pc);
template <class T> void TestSymDiv_E(tmv::DivType dt, PosDefCode pc);

template <class T> void TestSymBandMatrixArith_A();
template <class T> void TestSymBandMatrixArith_B1();
template <class T> void TestSymBandMatrixArith_B2();
template <class T> void TestSymBandMatrixArith_C1();
template <class T> void TestSymBandMatrixArith_C2();
template <class T> void TestSymBandMatrixArith_D1();
template <class T> void TestSymBandMatrixArith_D2();
template <class T> void TestSymBandMatrixArith_E1();
template <class T> void TestSymBandMatrixArith_E2();
template <class T> void TestSymBandMatrixArith_F1();
template <class T> void TestSymBandMatrixArith_F2();
template <class T> void TestSymBandDiv_A(tmv::DivType dt, PosDefCode pc);
template <class T> void TestSymBandDiv_B(tmv::DivType dt, PosDefCode pc);
template <class T> void TestSymBandDiv_C(tmv::DivType dt, PosDefCode pc);
template <class T> void TestSymBandDiv_D(tmv::DivType dt, PosDefCode pc);
template <class T> void TestSymBandDiv_E(tmv::DivType dt, PosDefCode pc);
template <class T> void TestSymBandDiv_F(tmv::DivType dt, PosDefCode pc);

