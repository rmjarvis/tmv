#include "Corr.h"
#include "CorrIO.h"
#include <vector>
#include <fstream>
#include <string>
#include <complex>
#include "dbg.h"

//#define ALLDEBUG
const double outputsize = 3.e3;
// Only output progress lines when size > outputsize.
// if ALLDEBUG is defined output all progress lines.

//#define XAssert(s) Assert(s)  
// XAssert adds extra (time-consuming) checks
#define XAssert(s)

#ifdef MEMDEBUG
AllocList* allocList;
#endif

int recursen=-1;
#include "Process2.h"

// Constants to set:
const double minsep = 1.;        // (arcsec) Minimum separation to consider
const double maxsep = 200.*60.;  // (arcsec) Maximum separation to consider
const double binsize = 0.05;     // Ratio of distances for consecutive bins
const double binslop = 1.0;      // Allowed slop on getting right distance bin
const double smoothscale = 2.;   // Smooth from r/smoothscale to r*smoothscale

// Derived Constants:
const double logminsep = log(minsep);
const double halfminsep = 0.5*minsep;
const int nbins = int(ceil(log(maxsep/minsep)/binsize));
const double b = binslop*binsize;
const double minsepsq = minsep*minsep;
const double maxsepsq = maxsep*maxsep;
const double bsq = b*b;

std::ostream* dbgout = 0;
bool XDEBUG = false;

// Switch commented line of these two to get more
// error checking during the execution

std::vector<int> bincount(nbins,0);

void DirectProcess11(
    std::vector<EEData>& data, const Cell& c1, const Cell& c2,
    const double dsq, const Position2D& r)
{
    std::complex<double> cr(r.GetX(),r.GetY());
    XAssert(c1.Size()+c2.Size() < sqrt(dsq)*b + 0.0001);
    XAssert(std::abs(c2.MeanPos() - c1.MeanPos() - cr) < 0.0001);
    XAssert(std::abs(dsq - std::norm(r)) < 0.0001);

    Assert(dsq > minsepsq);
    Assert(dsq < maxsepsq);

    std::complex<double> expm4iarg = conj(cr*cr)/dsq; // now expm2iarg
    expm4iarg *= expm4iarg; // now expm4iarg

    const double logr = log(dsq)/2.;
    Assert(logr >= logminsep);

    const int k = int(floor((logr - logminsep)/binsize));
    Assert(k >= 0); Assert(k<int(data.size()));

    const double ww = c1.Weight()*c2.Weight();
    const std::complex<double> e1 = c1.WE();
    const std::complex<double> e2 = c2.WE();

    const std::complex<double> ee = e1*e2*expm4iarg;
    const std::complex<double> eec = e1*conj(e2);

    const double npairs = c1.GetN()*c2.GetN();

    EEData& bindata = data[k];
    bindata.xiplus += eec;
    bindata.ximinus += ee;
    bindata.weight += ww;
    bindata.meanlogr += ww*logr;
    bindata.npair += npairs;

    bincount[k]++;
}

void FinalizeProcess(std::vector<EEData>& data, double vare)
{
    for(int i=0;i<nbins;i++) {
        EEData& bindata = data[i];
        double wt = data[i].weight;
        if (wt == 0.) 
            bindata.xiplus = bindata.ximinus =
                bindata.meanlogr = bindata.varxi = 0.;
        else {
            bindata.xiplus /= wt;
            bindata.ximinus /= wt;
            bindata.meanlogr /= wt;
            bindata.varxi = vare*vare/bindata.npair;
        }
        if (bindata.npair < 100.) bindata.meanlogr = logminsep+(i+0.5)*binsize; 
    }
}

int main(int argc, char* argv[])
{
#ifdef MEMDEBUG
    atexit(&DumpUnfreed);
#endif

    if (argc < 2 || argc > 3) myerror("Usage: corree filename [filename2]");

    dbgout = &std::cout;

    std::ifstream fin(argv[1]);

    double vare=0.;
    std::vector<CellData> celldata;
    Read(fin,celldata,vare);
    Cell wholefield(celldata);

    dbg<<"ngal = "<<celldata.size()<<std::endl;
    dbg<<"vare = "<<vare<<": sig_sn (per component) = "<<sqrt(vare)<<std::endl;
    dbg<<"nbins = "<<nbins<<": min,maxsep = "<<minsep<<','<<maxsep<<std::endl;

    std::vector<EEData> data(nbins);

    if (argc == 3) {
        std::ifstream fin2(argv[2]);
        double vare2=0.;
        std::vector<CellData> celldata2;
        Read(fin2,celldata2,vare2);
        Cell wholefield2(celldata2);
        Process11(data,minsep,maxsep,minsepsq,maxsepsq,wholefield,wholefield2);
        vare = sqrt(vare*vare2);
    } else {
        Process2(data,minsep,maxsep,minsepsq,maxsepsq,wholefield);
    }

    FinalizeProcess(data,vare);
    dbg<<"done processing\n";

    std::ofstream e2out("e2.out");
    std::ofstream m2out("m2.out");

    WriteEE(e2out,data);
    WriteM2(m2out,data);

    if (dbgout && dbgout != &std::cout) 
    { delete dbgout; dbgout=0; }
    return 0;
}

