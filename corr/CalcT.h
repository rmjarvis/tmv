#ifndef CORRIO_H
#define CORRIO_H

#include "Corr.h"
#include "Form.h"
#include "MemDebug.h"

void CalcT(
    double s, std::complex<double> t, 
    std::complex<double>* T0, std::complex<double>* T1,
    std::complex<double>* T2, std::complex<double>* T3,
    double R1, double R2=0., double R3=0.);

double Tplus(double x);
double Tminus(double x);
double Tcross(double x);
double Splus(double r);
double Sminus(double r);

#endif
