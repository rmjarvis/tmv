#include "Cell.h"
#include "Form.h"
#include "Corr.h"
#include "CorrIO.h"
#include "dbg.h"

#include <fstream>
#include <string>

double outputsize = 1.e3;
//#define ALLDEBUG
//#define XAssert(s) Assert(s)
#define XAssert(s)

#ifdef MEMDEBUG
AllocList* allocList;
#endif

int recursen=-1;
#include "Process2.h"

//#define DOCROSSRAND

std::ostream* dbgout = 0;
bool XDEBUG = false;

// Constants to set:
const double minsep = 1.;        // (arcsec) Minimum separation to consider
const double maxsep = 200.*60.;  // (arcsec) Maximum separation to consider
const double binsize = 0.05;    // Ratio of distances for consecutive bins
const double binslop = 1.0;      // Allowed slop on getting right distance bin
const double smoothscale = 2.;   // Smooth from r/smoothscale to r*smoothscale

const int nthetabins = 32; // nbins for theta binning test
const double rmintheta = 40.*60.;  // min sep (arcsec) for theta binning test

// Derived Constants:
const double logminsep = log(minsep);
const double halfminsep = 0.5*minsep;
const int nbins = int(ceil(log(maxsep/minsep)/binsize));
const double b = 2.*binslop*binsize;
const int ktheta = int(floor((log(rmintheta)-logminsep)/binsize));
const double minsepsq = minsep*minsep;
const double maxsepsq = maxsep*maxsep;
const double bsq = b*b;

template <class T> inline T SQR(const T& x) { return x*x; }

void DirectProcess11(
    std::vector<NEData>& data, const NCell& c1, const Cell& c2,
    const double dsq, const Position2D& r)
{
#ifdef ALLDEBUG
    dbg<<std::string(recursen+1,'-')<<"Direct: d = "<<d<<std::endl;
#endif
    XAssert(c1.Size()+c2.Size() < sqrt(dsq)*b + 0.0001);
    XAssert(std::abs(c2.MeanPos() - c1.MeanPos() - r) < 0.0001);
    XAssert(std::abs(dsq - std::norm(r)) < 0.0001);

    Assert(dsq >= minsepsq);
    Assert(dsq < maxsepsq);

    std::complex<double> cr(r.GetX(),r.GetY());
    const std::complex<double> expm2iarg = SQR(conj(cr))/dsq;

    const double logr = log(dsq)/2.;
    Assert(logr >= logminsep);

    const int k = int(floor((logr - logminsep)/binsize));
    Assert(k >= 0); Assert(k<int(data.size()));

    const double nw = c1.GetN()*c2.Weight();
    const std::complex<double> net = -double(c1.GetN())*c2.WE()*expm2iarg;
    const double npairs = c1.GetN()*c2.GetN();

    NEData& crossbin = data[k];
    crossbin.meangammat += net;
    crossbin.weight += nw;
    crossbin.meanlogr += nw*logr;
    crossbin.npair += npairs;
}

void FinalizeProcess(std::vector<NEData>& data, double vare)
{
    for(int i=0;i<nbins;i++) {
        NEData& crossbin = data[i];
        double wt = crossbin.weight;
        if (wt == 0.) 
            crossbin.meangammat = crossbin.meanlogr = crossbin.vargammat = 0.;
        else {
            crossbin.meangammat /= wt;
            crossbin.meanlogr /= wt;
            crossbin.vargammat = vare/crossbin.npair;
        }
        if (crossbin.npair<100.) crossbin.meanlogr = logminsep+(i+0.5)*binsize;
    }
}

int main(int argc, char* argv[])
{
#ifdef MEMDEBUG
    atexit(&DumpUnfreed);
#endif

    if (argc < 3) myerror("Usage: corrne brightfile faintfile");

    dbgout = new std::ofstream("ne.debug");

    std::ifstream brightfin(argv[1]);
    std::ifstream faintfin(argv[2]);

    double vare;
    std::vector<NCellData> brightdata;
    std::vector<CellData> faintdata;

    dbg << "Read bright gals\n";
    Read(brightfin,brightdata);
    if (!brightfin) myerror("reading brightfile ",argv[1]);
    NCell brightfield(brightdata);

    dbg << "Read faint gals\n";
    Read(faintfin,faintdata,vare);
    if (!faintfin) myerror("reading faintfile ",argv[2]);

    Cell faintfield(faintdata);

    dbg << "bright size = "<<brightdata.size();
    dbg <<", faintsize = "<<faintdata.size()<<std::endl;

    dbg<<"vare = "<<vare<<": sig_sn (per component) = "<<sqrt(vare)<<std::endl;
    dbg<<"nbins = "<<nbins<<": min,maxsep = "<<minsep<<','<<maxsep<<std::endl;

    std::vector<NEData> crossdata(nbins);
    Process11(crossdata,minsep,maxsep,minsepsq,maxsepsq,brightfield,faintfield);
    FinalizeProcess(crossdata,vare);

    if (dbgout) {
        dbg<<"done process crossdata:\n";
        //WriteNE(*dbgout,crossdata);
    }

    std::ofstream x2out("x2.out");
    WriteNE(x2out,crossdata);

    if (dbgout && dbgout != &std::cout) {delete dbgout; dbgout=0;}
    return 0;
}

