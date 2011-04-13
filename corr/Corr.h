#ifndef CORR_H
#define CORR_H

#include "Cell.h"
#include "MemDebug.h"

struct EEData {
    EEData() :
        xiplus(0.), ximinus(0.), varxi(0.), meanlogr(0.), weight(0.), npair(0.) 
    {}

    std::complex<double> xiplus,ximinus;
    double varxi, meanlogr, weight, npair;
};

struct EEEData {
    EEEData() :
        gam0(0.), gam1(0.), gam2(0.), gam3(0.), vargam(0.),
        meanr(0.), meanu(0.), meanv(0.), weight(0.), ntri(0.) 
    {}

    std::complex<double> gam0, gam1, gam2, gam3;
    double vargam, meanr, meanu, meanv, weight, ntri;
};

struct NEData {
    NEData() : 
        meangammat(0.), vargammat(0.), meanlogr(0.), weight(0.), npair(0.) 
    {}

    std::complex<double> meangammat;
    double vargammat, meanlogr, weight, npair;
};

extern const double minsep;
extern const double maxsep;
extern const double minsepsq;
extern const double maxsepsq;
extern const double binsize;
extern const double b;
extern const double bsq;

template <class CellType1, class CellType2>
inline void CalcSplit(
    bool& split1, bool& split2, const CellType1& c1, 
    const CellType2& c2, const double d3)
{
    // This function determines whether either input cell needs to be
    // split.  It is written as a template so that the second cell
    // can be either a Cell or a PairCell.  (The same rules apply.)
    // If you already know that c1 needs to be split, then split1 can
    // be input as true, and we only check c2.  (and vice versa)
    // In normal operation, both are input as false, and we check
    // whether they need to be split.  
    //
    // If (s1+s2)/d > b, then we need to split one or both.
    //
    // If s1 > b*d, then it doesn't matter what s2 is -- we should
    // definitely split c1.  
    // Likewise if s2 > b*d
    //
    // If neither of these conditions is true, then we test
    // to see if s1+s2 > b*d
    // If we're close to the threshold, it will generally be quicker
    // to only split the larger one.  But if both s1 and s2 are 
    // reasonably large (compared to b/d) then we will probably end 
    // up needing to split both, so go ahead and do so now.
    // This small vs. large test is quantified by the parameter
    // splitfactor.  I varied split factor with the 2-point
    // correlation code until it ran fastest.  The result is 
    // given above.  I don't know if this value is also best for
    // 3 point uses, but it's probably reasonably close.


    const double splitfactor = 0.585;
    // The split factor helps determine whether to split
    // both cells or just one when the factor (s1+s2)/d 
    // is too large.  
    // If s1+s2 > f*d*b then split both.
    // Otherwise just split the large Cell.
    // The value of f was determined empirically by seeing 
    // when the code ran fastest.  This may be specific
    // to the data I was testing it on, but I would guess
    // that this value is close enough to optimal for most
    // datasets.

    if (split1) {
        if (split2) {
            // do nothing
        } else {
            const double s2 = c2.Size();
            const double maxs = d3*b;
            split2 = s2 > maxs;
        }
    } else {
        if (split2) {
            const double s1 = c1.Size();
            const double maxs = d3*b;
            split1 = s1 > maxs;
        } else {
            const double s1 = c1.Size();
            const double s2 = c2.Size();
            const double maxs = d3*b;
            split1 = s1 > maxs;
            split2 = s2 > maxs;
            if (!split1 && !split2) {
                const double sum = s1+s2;
                if (sum > maxs) {
                    double modmax = splitfactor*maxs;
                    if (s1 > s2) {
                        split1 = true;
                        split2 = (s2 > modmax);
                    }
                    else {
                        split2 = true;
                        split1 = (s1 > modmax);
                    }
                }
            }
        }
    }
}

template <class CellType1, class CellType2>
inline void CalcSplitSq(
    bool& split1, bool& split2, const CellType1& c1, 
    const CellType2& c2, const double d3sq)
{
    const double splitfactorsq = 0.3422;

    // The same as above, but when we know the distance squared rather
    // than just the distance.  We get some speed up by saving the 
    // square roots in some parts of the code.
    if (split1) {
        if (split2) {
            // do nothing
        } else {
            const double s2sq = c2.SizeSq();
            const double maxssq = bsq*d3sq;
            split2 = s2sq > maxssq;
        }
    } else {
        if (split2) {
            const double s1sq = c1.SizeSq();
            const double maxssq = bsq*d3sq;
            split1 = s1sq > maxssq;
        } else {
            const double s1sq = c1.SizeSq();
            const double s2sq = c2.SizeSq();
            const double maxssq = bsq*d3sq;
            split1 = s1sq > maxssq;
            split2 = s2sq > maxssq;
            if (!split1 && !split2) {
                double sumsq = c1.Size()+c2.Size();
                sumsq *= sumsq;
                if (sumsq > maxssq) {
                    double modmax = splitfactorsq*maxssq;
                    if (s1sq > s2sq) {
                        split1 = true;
                        split2 = (s2sq > modmax);
                    }
                    else {
                        split2 = true;
                        split1 = (s1sq > modmax);
                    }
                }
            }
        }
    }
}

inline bool NoSplit(const Cell& c2, const Cell& c3, const double d1)
{
    static const double altb = b/(1.-b);
    // A debugging routine.  Usually of the form:
    // XAssert(NoSplit(c2,c3,d1))
    // This just asserts that the cells obey the non-splitting eqn:
    // (s1 + s2)/d < b
    // Technically we use altb = b/(1-b) which = b for small b.
    if (c2.Size() + c3.Size() < d1*altb+0.0001)
        return true;
    else {
        std::cerr<<c2.Size()<<" + "<<c3.Size()<<" > "<<
            d1<<" * "<<b/(1-b)<<std::endl;
        return false;
    }
}

inline bool Check(
    const Cell& c1, const Cell& c2, const Cell& c3,
    const double d1, const double d2, const double d3)
{
    // Checks that d1,d2,d3 are correct for the three Cells given.
    // Used as a debugging check.
    bool ok=true;
    if (Dist(c3.MeanPos(),c2.MeanPos())-d1 > 0.0001) 
    { std::cerr<<"d1\n"; ok = false; }
    if (Dist(c1.MeanPos(),c3.MeanPos())-d2 > 0.0001) 
    { std::cerr<<"d2\n"; ok = false; }
    if (Dist(c2.MeanPos(),c1.MeanPos())-d3 > 0.0001) 
    { std::cerr<<"d3\n"; ok = false; }
    if (d1 > d2+d3+0.0001) {std::cerr<<"sum d1\n"; ok = false; }
    if (d2 > d1+d3+0.0001) {std::cerr<<"sum d2\n"; ok = false; }
    if (d3 > d1+d2+0.0001) {std::cerr<<"sum d3\n"; ok = false; }
    return ok;
}

#endif
