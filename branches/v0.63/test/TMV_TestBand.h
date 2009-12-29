// vim:et:ts=2:sw=2:ci:cino=f0,g0,t0,+0:
#ifndef TMV_TESTBAND_H
#define TMV_TESTBAND_H

#include "TMV_Test.h"
#include "TMV.h"
#include "TMV_Band.h"
#include "TMV_TestMatrixArith.h"

template <class T, tmv::StorageType S> void TestBasicBandMatrix();
template <class T> void TestBandMatrixArith_A();
template <class T> void TestBandMatrixArith_B();

#endif
