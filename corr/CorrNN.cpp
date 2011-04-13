#include "Cell.h"
#include "Form.h"
#include "Corr.h"
#include "CorrIO.h"
#include "dbg.h"

#include <fstream>
#include <string>
#include <algorithm>

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
const double binsize = 0.05;     // Ratio of distances for consecutive bins
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

//std::vector<int> bincount(nbins,0);

std::vector<std::complex<double> > testthetanet(nthetabins,0.);
std::vector<double> testthetanw(nthetabins,0.);
std::vector<int> testthetan(nthetabins,0.);

void DirectProcess11(
    std::vector<double>& data, const NCell& c1, const NCell& c2,
    const double dsq, const Position2D& r)
{
#ifdef ALLDEBUG
    dbg<<std::string(recursen+1,'-')<<"Direct: d = "<<d<<std::endl;
#endif
    XAssert(c1.Size()+c2.Size() < sqrt(dsq)*b + 0.0001);
    XAssert(Dist(c2.MeanPos() - c1.MeanPos(),r) < 0.0001);
    XAssert(std::abs(dsq - std::norm(r)) < 0.0001);

    Assert(dsq >= minsepsq);
    Assert(dsq < maxsepsq);

    const double logr = log(dsq)/2.;
    Assert(logr >= logminsep);

    const int k = int(floor((logr - logminsep)/binsize));
    Assert(k >= 0); 
    Assert(k<int(data.size()));

    const double npairs = c1.GetN()*c2.GetN();

    data[k] += npairs;
}

int main(int argc, char* argv[])
{
#ifdef MEMDEBUG
    atexit(&DumpUnfreed);
#endif

    if (argc < 3) myerror("Usage: corrnn datafile randfiles");

    dbgout = new std::ofstream("nn.debug");

    std::ifstream fin(argv[1]);
    std::ifstream randlistfin(argv[2]);

    std::vector<NCellData> data;

    dbg << "Read gals\n";
    Read(fin,data);
    if (!fin) myerror("reading file ",argv[1]);
    NCell wholefield(data);

    dbg << "ngals = "<<data.size();

    dbg << "Read rand fields\n";
    int nrandfields;
    randlistfin >> nrandfields;
    if (!randlistfin) myerror("reading randlistfile ",argv[2]);
    if (!(nrandfields > 0)) myerror("no random fields");
    std::vector<std::vector<NCellData> > randdata(nrandfields);
    std::vector<NCell*> randfield(nrandfields);
    for(int i=0;i<nrandfields;i++) {
        std::string randfieldname;
        randlistfin >> randfieldname;
        if (!randlistfin) myerror("reading randlistfile ",argv[2]);
        std::ifstream randfin(randfieldname.c_str());
        Read(randfin,randdata[i]);
        if (!randfin) myerror("reading randfile ",randfieldname.c_str());
        randfield[i] = new NCell(randdata[i]);
    }

    dbg<<"nbins = "<<nbins<<": min,maxsep = "<<minsep<<','<<maxsep<<std::endl;

    std::vector<double> DD(nbins,0.), RR(nbins,0.), DR(nbins,0.);

    Process2(DD,minsep,maxsep,minsepsq,maxsepsq,wholefield);
    for(int k=0;k<nbins;k++) DD[k] *= 2;

    for(int i=0;i<nrandfields;i++) {
        dbg<<"rand: i = "<<i<<std::endl;
        Process2(RR,minsep,maxsep,minsepsq,maxsepsq,*randfield[i]);
        Process11(DR,minsep,maxsep,minsepsq,maxsepsq,wholefield,*randfield[i]);
    }
    for(int k=0;k<nbins;k++) RR[k] *= 2;
    int ndr = nrandfields;
    int nrr = nrandfields;

#ifdef DOCROSSRAND
    for(int i=0;i<nrandfields;i++) {
        for(int j=i+1;j<nrandfields;j++) {
            dbg<<"j = "<<j<<std::endl;
            Process11(RR,minsep,maxsep,minsepsq,maxsepsq,
                      *randfield[i],*randfield[j]);
            nrr++;
        }
    }
#endif
    dbg<<"done process DD,DR,RR data\n";

    for(int i=0;i<nbins;i++) {
        dbg<<"RR["<<i<<"] = "<<RR[i]<<", DD = "<<DD[i]<<", DR = "<<DR[i]<<std::endl;
        RR[i] /= nrr;
        DR[i] /= ndr;
    }
    dbg<<"done divide by number of realizations\n";

    std::ofstream fout("nn.out");
    WriteNN(fout,DD,DR,RR,nrr);
    dbg<<"done write"<<std::endl;

    if (dbgout && dbgout != &std::cout) 
    { delete dbgout; dbgout=0; }
    return 0;
}

