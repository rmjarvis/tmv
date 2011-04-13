
template <class DataType, class CellType1, class CellType2>
void Process11(
    std::vector<DataType>& data, double minr, double maxr,
    double minrsq, double maxrsq,
    const CellType1& c1, const CellType2& c2)
{
    const double dsq = DistSq(c1.MeanPos(),c2.MeanPos());
    const double s1ps2 = c1.AllSize()+c2.AllSize();

    if (dsq < minrsq) if (sqrt(dsq)+s1ps2 < minr) return;
    if (dsq > maxrsq) if (sqrt(dsq)-s1ps2 > maxr) return;

#ifdef ALLDEBUG
    recursen++;
    dbg<<std::string(recursen,'.')<<"Start P11: "<<c1.MeanPos()<<" -- "<<c2.MeanPos()<<"   N = "<<c1.GetN()<<","<<c2.GetN()<<" d = "<<sqrt(dsq)<<std::endl;
#endif

    // See if need to split:
    bool split1=false, split2=false;
    CalcSplitSq(split1,split2,c1,c2,dsq);

    if (split1) {
        if (split2) {
            if (!c1.Left()) {
                dbg<<"minr = "<<minr<<", maxr = "<<maxr<<std::endl;
                dbg<<"minrsq = "<<minrsq<<", maxrsq = "<<maxrsq<<std::endl;
                dbg<<"c1.Size = "<<c1.Size()<<", c2.Size = "<<c2.Size()<<std::endl;
                dbg<<"c1.SizeSq = "<<c1.SizeSq()<<", c2.SizeSq = "<<c2.SizeSq()<<std::endl;
                dbg<<"c1.GetN = "<<c1.GetN()<<", c2.GetN = "<<c2.GetN()<<std::endl;
                dbg<<"c1.Pos = "<<c1.MeanPos()<<", c2.Pos = "<<c2.MeanPos()<<std::endl;
                dbg<<"dsq = "<<dsq<<", s1ps2 = "<<s1ps2<<std::endl;
            }
            Assert(c1.Left());
            Assert(c1.Right());
            Assert(c2.Left());
            Assert(c2.Right());
            Process11(data,minr,maxr,minrsq,maxrsq,*c1.Left(),*c2.Left());
            Process11(data,minr,maxr,minrsq,maxrsq,*c1.Left(),*c2.Right());
            Process11(data,minr,maxr,minrsq,maxrsq,*c1.Right(),*c2.Left());
            Process11(data,minr,maxr,minrsq,maxrsq,*c1.Right(),*c2.Right());
        } else {
            Assert(c1.Left());
            Assert(c1.Right());
            Process11(data,minr,maxr,minrsq,maxrsq,*c1.Left(),c2);
            Process11(data,minr,maxr,minrsq,maxrsq,*c1.Right(),c2);
        }
    } else {
        if (split2) {
            Assert(c2.Left());
            Assert(c2.Right());
            Process11(data,minr,maxr,minrsq,maxrsq,c1,*c2.Left());
            Process11(data,minr,maxr,minrsq,maxrsq,c1,*c2.Right());
        } else if (dsq > minrsq && dsq < maxrsq) {
            const typename CellType1::PosType& r = c2.MeanPos()-c1.MeanPos();
            XAssert(NoSplit(c1,c2,sqrt(dsq)));
            XAssert(Dist(c2.MeanPos() - c1.MeanPos(),r) < 0.0001);
            DirectProcess11(data,c1,c2,dsq,r);
        }
    }

#ifdef ALLDEBUG
    dbg<<std::string(recursen,'.')<<"Done P11\n";
    recursen--;
#endif
}

template <class DataType, class CellType>
void Process2(
    std::vector<DataType>& data, double minr, double maxr,
    double minrsq, double maxrsq, const CellType& c12)
{
    if (c12.Size() < minr/2.) return;

    recursen++;

    std::ostream* tempdbg = dbgout;
    if (dbgout) {
#ifndef ALLDEBUG
        if (c12.Size() < outputsize) dbgout = 0;
#endif
        dbg<<std::string(recursen,'.')<<"Start P2: size = "<<c12.Size()<<", center = "<<c12.MeanPos()<<"   N = "<<c12.GetN()<<std::endl;
    }

#ifdef ALLDEBUG
    dbg<<std::string(recursen,'.')<<"P2 ("<<c12.Size()<<") call Left P2\n";
#endif
    Assert(c12.Left());
    Assert(c12.Right());
    Process2(data,minr,maxr,minrsq,maxrsq,*c12.Left());

#ifdef ALLDEBUG
    dbg<<std::string(recursen,'.')<<"P2 ("<<c12.Size()<<") call Right P2\n";
#endif
    Process2(data,minr,maxr,minrsq,maxrsq,*c12.Right());

#ifdef ALLDEBUG
    dbg<<std::string(recursen,'.')<<"P2 ("<<c12.Size()<<") call P11\n";
#endif
    Process11(data,minr,maxr,minrsq,maxrsq,*c12.Left(),*c12.Right());

    if (tempdbg) {
        dbg<<std::string(recursen,'.')<<"Done P2: size = "<<c12.Size()<<", center = "<<c12.MeanPos()<<std::endl;
        dbgout = tempdbg;
    }
    recursen--;
}
