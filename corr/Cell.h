#ifndef CELL_H
#define CELL_H

enum SplitMethod { MIDDLE, MEDIAN, MEAN };

#define INF 1.e200

#include <iostream>
#include <algorithm>
#include <complex>
#include <vector>
#include "dbg.h"
#include "Bounds.h"
#include "MemDebug.h"

const double PI = 3.141592653589793;
const double TWOPI = 2.*PI;
const double IOTA = 1.e-10;

struct NCellData 
{
    Position2D pos;

    NCellData() {}
    NCellData(const Position2D& _pos) : pos(_pos) {}
};

struct CellData : 
    public NCellData 
{
    std::complex<double> e;
    double w;

    CellData() {}
    CellData(const Position2D& _pos, const std::complex<double>& _e, 
             double _w) : 
        NCellData(_pos), e(_e), w(_w) {}
};

template <class CellType>
inline int DoCountLeaves(const CellType* c)
{
    if (c->Left()) {
        Assert(c->Right());
        return c->Left()->CountLeaves() + c->Right()->CountLeaves();
    } else return 1;
}

template <class CellType>
inline std::vector<const CellType*> DoGetAllLeaves(const CellType* c)
{
    std::vector<const CellType*> ret;
    if (c->Left()) {
        std::vector<const CellType*> temp = c->Left()->GetAllLeaves();
        ret.insert(ret.end(),temp.begin(),temp.end());
        Assert(c->Right());
        temp = c->Right()->GetAllLeaves();
        ret.insert(ret.end(),temp.begin(),temp.end()); 
    } else {
        Assert(!c->Right());
        ret.push_back(static_cast<const CellType*>(c));
    }
    return ret;
}

template <class DataType>
struct DataCompare 
{

    bool splitonx;

    DataCompare(bool s) : splitonx(s) {}
    bool operator()(const DataType& cd1, const DataType& cd2) const 
    {
        return (splitonx ?
                cd1.pos.GetX() < cd2.pos.GetX() :
                cd1.pos.GetY() < cd2.pos.GetY());
    }
};

template <class DataType>
struct DataCompareToValue 
{

    bool splitonx;
    double splitvalue;

    DataCompareToValue(bool s, double v) : splitonx(s), splitvalue(v) {}
    bool operator()(const DataType& cd) const 
    {
        return (splitonx ? cd.pos.GetX() : cd.pos.GetY()) < splitvalue;
    }
};

template <class DataType>
inline size_t SplitCell(
    std::vector<DataType>& vdata, SplitMethod sm, 
    size_t start, size_t end, const Position2D& meanpos)
{
    size_t mid;

    Bounds2D b;
    for(size_t i=start;i<end;i++) b += vdata[i].pos;

    bool splitonx = ((b.GetXMax()-b.GetXMin()) > (b.GetYMax() - b.GetYMin()));

    switch (sm) { // three different split methods
      case MIDDLE :
           { // Middle is the average of the min and max value of x or y
               double splitvalue = 
                   ( splitonx ?
                     (b.GetXMax()+b.GetXMin())/2. :
                     (b.GetYMax()+b.GetYMin())/2.);
               DataCompareToValue<DataType> comp(splitonx,splitvalue);
               typename std::vector<DataType>::iterator middle = 
                   std::partition(vdata.begin()+start,vdata.begin()+end,comp);
               mid = middle - vdata.begin();
           } break;
      case MEDIAN :
           { // Median is the point which divides the group into equal numbers
               DataCompare<DataType> comp(splitonx);
               mid = (start+end)/2;
               typename std::vector<DataType>::iterator middle =
                   vdata.begin()+mid;
               std::nth_element(
                   vdata.begin()+start,middle,vdata.begin()+end,comp);
           } break;
      case MEAN :
           { // Mean is the weighted average value of x or y
               double splitvalue = 
                   (splitonx ? meanpos.GetX() : meanpos.GetY());
               DataCompareToValue<DataType> comp(splitonx,splitvalue);
               typename std::vector<DataType>::iterator middle = 
                   std::partition(vdata.begin()+start,vdata.begin()+end,comp);
               mid = middle - vdata.begin();
           } break;
      default :
           myerror("Invalid SplitMethod");
    }

    Assert(mid > start);
    Assert(mid < end);
    return mid;
}

class NCell 
{

    // A Cell that only counts the number of galaxies.

public:

    NCell() : meanpos(),size(0.),ngals(0),left(0),right(0) {}
    NCell(std::vector<NCellData>& data, SplitMethod sm=MEAN,
          size_t start=0, size_t end=0);
    ~NCell() { if (left) delete left; if(right) delete right;}

    const Position2D& MeanPos() const { return meanpos; }
    double Size() const { return size; }
    double SizeSq() const { return sizesq; }
    double AllSize() const { return size; } 
    // For PairCell's AllSize is different
    int GetN() const { return ngals; }

    const NCell* Left() const { return left; }
    const NCell* Right() const { return right; }

    int CountLeaves() const { return DoCountLeaves(this); }
    std::vector<const NCell*> GetAllLeaves() const 
    { return DoGetAllLeaves(this); }

    typedef Position2D PosType;

protected:
    Position2D meanpos;
    double size;
    double sizesq;
    int ngals;

    NCell* left;
    NCell* right;
};

inline NCell::NCell(
    std::vector<NCellData>& vdata, SplitMethod sm,
    size_t start, size_t end) :
    meanpos(),size(0.),sizesq(0.),ngals(0),left(0),right(0)
{
    if (end == 0) end = vdata.size();
    Assert(vdata.size()>0);
    Assert(end <= vdata.size());
    Assert(end > start);

    ngals = end-start;

    if (end - start == 1) {
        meanpos = vdata[start].pos;
    } else {
        for(size_t i=start;i<end;i++) {
            meanpos += vdata[i].pos;
        }
        meanpos /= double(ngals);

        for(size_t i=start;i<end;i++) {
            double devsq = DistSq(vdata[i].pos,meanpos);
            if (devsq > sizesq) sizesq = devsq;
        }
        size = sqrt(sizesq);

        if (size > 0.) {
            size_t mid = SplitCell(vdata,sm,start,end,meanpos);

            left = new NCell(vdata,sm,start,mid);
            right = new NCell(vdata,sm,mid,end);
            if (!left || !right) 
                myerror("out of memory - cannot create new NCell");
        }
    }
}

class Cell : public NCell 
{

    // A Cell contains the accumulated data for a bunch of galaxies.
    // It is characterized primarily by a centroid and a size.
    // The centroid is simply the weighted centroid of all the galaxy positions.
    // The size is the maximum deviation of any one of these galaxies 
    // from the centroid.  That is, all galaxies fall within a radius
    // size from the centroid.
    // The structure also keeps track of some averages and sums about
    // the galaxies which are used in the correlation function calculations.

public:

    Cell(std::vector<CellData>& data, SplitMethod sm=MEAN,
         size_t start=0, size_t end=0);
    ~Cell() {}

    const std::complex<double>& WE() const { return we; }
    double ENorm() const { return enorm; }
    double Weight() const { return w; }

    const Cell* Left() const 
    { return static_cast<const Cell*>(left); }
    const Cell* Right() const 
    { return static_cast<const Cell*>(right); }

    std::vector<const Cell*> GetAllLeaves() const {return DoGetAllLeaves(this);}

protected:
    std::complex<double> we;
    double enorm;
    double w;
};

inline Cell::Cell(
    std::vector<CellData>& vdata, SplitMethod sm,
    size_t start, size_t end) : NCell(),we(0.),enorm(0.),w(0.)
{
    if (end == 0) end = vdata.size();
    Assert(vdata.size()>0);
    Assert(end <= vdata.size());
    Assert(end > start);

    ngals = end-start;

    if (end - start == 1) {
        const CellData& data = vdata[start];
        meanpos = data.pos;
        we = data.e*data.w;
        enorm = std::norm(data.e);
        w = data.w;
        Assert(w>0.);
    } else {
        for(size_t i=start;i<end;i++) {
            const CellData& data = vdata[i];
            meanpos += data.pos * data.w;
            we += data.e * data.w;
            w += data.w;
        }
        Assert(w>0.);
        enorm = std::norm(we/w);
        meanpos /= w;

        sizesq = 0.;
        for(size_t i=start;i<end;i++) {
            double devsq = DistSq(vdata[i].pos,meanpos);
            if (devsq > sizesq) sizesq = devsq;
        }
        size = sqrt(sizesq);

        if (size > 0.) {
            size_t mid = SplitCell(vdata,sm,start,end,meanpos);

            left = new Cell(vdata,sm,start,mid);
            right = new Cell(vdata,sm,mid,end);
            if (!left || !right) 
                myerror("out of memory - cannot create new Cell");
        }
    }
}

#endif
