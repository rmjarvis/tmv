// This program calculates the <M^3> statistic from the 
// correlation function.
//

#include "Corr.h"
#include "CorrIO.h"
#include "dbg.h"

#include <fstream>

#include <algorithm>

const double minsep = 1.;       // (arcsec) Minimum separation to consider
const double maxsep = 200.*60.;  // (arcsec) Maximum separation to consider
const double binsize = 0.05;     // Size of bins in logr, u, v

// Derived Constants:

// The number of bins needed given the binsize
const int nrbins = int(ceil(log(maxsep/minsep)/binsize));
const int nubins = int(ceil(1./binsize));
const int nvbins = 2*nubins;

bool XDEBUG = false;
std::ostream* dbgout=0;

int main()
{
    dbgout = &std::cout;
    //dbgout = new std::ofstream("e3_m3.debug");

    std::vector<std::vector<std::vector<EEEData> > > threeptdata(
        nrbins, std::vector<std::vector<EEEData> >(
            nubins,std::vector<EEEData>(nvbins)));
    std::ifstream e3in("e3.out");
    Read(e3in,threeptdata);

    std::ofstream m3out("m3.out");
    WriteM3(m3out,threeptdata);

#if 0
    std::ofstream m3out112("e3.m3112.out");
    WriteM3(m3out112,threeptdata,1.,1.,2.);

    std::ofstream m3out121("e3.m3121.out");
    WriteM3(m3out121,threeptdata,1.,2.,1.);

    std::ofstream m3out211("e3.m3211.out");
    WriteM3(m3out211,threeptdata,2.,1.,1.);

    std::ofstream m3out122("e3.m3122.out");
    WriteM3(m3out122,threeptdata,1.,2.,2.);

    std::ofstream m3out124("e3.m3124.out");
    WriteM3(m3out124,threeptdata,1.,2.,4.);

    std::ofstream m3out114("e3.m3114.out");
    WriteM3(m3out114,threeptdata,1.,1.,4.);

    std::ofstream m3out144("e3.m3144.out");
    WriteM3(m3out144,threeptdata,1.,4.,4.);
#endif

    if (dbgout && dbgout != &std::cout) 
    { delete dbgout; dbgout=0; }
    return 0;
}

