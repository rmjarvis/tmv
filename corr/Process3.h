
template <class DataType, class CellType1, class CellType2, class CellType3>
void ProcessV(
    std::vector<DataType>& vdata,
    const CellType1& c1, const CellType2& c2, const CellType3& c3,
    const double d3, const bool swap12, const double u);

template <class DataType, class CellType1, class CellType2, class CellType3>
void ProcessV(
    std::vector<DataType>& vdata,
    const CellType1& c1, const CellType2& c2, const CellType3& c3,
    const double d1, const double d2, const double d3,
    const bool swap12, const double u)
// Preconditions: s1+s2/d3 < b
//                s3/d2 < b/u
//                d3 < d2 < d1
{
    static const double sqrttwob = sqrt(2.*b);
#define altb (b/(1.-b))

#ifdef ALLDEBUG
    recursen++;
    dbg<<std::string(recursen,'.')<<"Start ProcessV "<<d1<<" >= "<<d2<<" >= "<<d3<<std::endl;
    dbg<<std::string(recursen,'.')<<"sizes = "<<c1.Size()<<" -- "<<c2.Size()<<" -- "<<c3.Size()<<std::endl;
#endif

    XAssert(NoSplit(c1,c2,d3));
    XAssert(Check(c1,c2,c3,d1,d2,d3));

    Assert(d2 >= d3);
    Assert(d1 >= d2);
    XAssert(std::abs(u-d3/d2) < altb);
    XAssert(c3.Size()/d2 < altb*d2/MAX(d3-altb*d2,1.e-10));

    // v is a bit more complicated than u:
    // Consider the angle bisector of d1,d2.  Call this b
    // Then let phi = the angle between d1 (or d2) and b
    // And let theta = the (smaller) angle between d3 and b
    // Then projecting d1,d2,d3 onto b, one finds:
    // d1 cos phi - d2 cos phi = d3 cos theta
    // v = (d1-d2)/d3 = cos theta / cos phi
    // Note that phi < 30 degrees, so cos phi won't make much
    // of a difference here.
    // So the biggest change in v from moving c3 is in theta:
    // dv = abs(dv/dtheta) dtheta = abs(sin theta)/cos phi (s3/b)
    // b cos phi > d2, so 
    // dv < b ==> s < m b cos phi / sqrt(1-v^2) < m d2 / sqrt(1-v^2)
    // But - just like with u where the denominator was u+m, not just u,
    // we need to modify this.  
    // If v = 1, dv = 1-cos(dtheta) = 1-cos(s3/b) ~= 1/2(s3/b)^2
    // So in this case, s3/b < sqrt(2m)
    // In general we require both to be true.

    double soverd = c3.Size()/d2;
    double v = (d1-d2)/d3;
    const bool split3 = c3.Size() > 0. && 
        ( soverd > sqrttwob || 
          soverd > b / sqrt(1-SQR(v)) );

    if (split3) {
        Assert(c3.Size()>0.);
        Assert(c3.Left());
        Assert(c3.Right());
        ProcessV(vdata,c1,c2,*c3.Left(),d3,swap12,u);
        ProcessV(vdata,c1,c2,*c3.Right(),d3,swap12,u);
    } else {
        if (c3.Size() == 0.) { // then v didn't get set above
            v = (d1-d2)/d3;
        }
        if (v > 0.99999) v = 0.99999;

        const typename CellType1::PosType r3 = c2.MeanPos() - c1.MeanPos();
        const typename CellType1::PosType r1 = c3.MeanPos() - c2.MeanPos();
        // This next thing only makes sense for 2D.
        // Need to think about this if I want to extend 3 point stuff into 3D.
        std::complex<double> cr3(r3.GetX(),r3.GetY());
        std::complex<double> cr1(r1.GetX(),r1.GetY());
        if (imag(cr1*conj(cr3)) < 0.) v = -v;

        int kv = int(floor((v+1.)/dv));
        //if (kv < 0) {Assert(kv==-1); kv=0;}
        //if (kv >= nvbins) {Assert(kv==nvbins); kv = nvbins-1;}
        Assert(kv >= 0); Assert(kv<int(vdata.size()));
        DirectProcessV(vdata[kv],c1,c2,c3,d1,d2,d3,swap12,u,v);

    }
#ifdef ALLDEBUG
    dbg<<std::string(recursen,'.')<<"Done PV\n";
    recursen--;
#endif
}

template <class DataType, class CellType1, class CellType2, class CellType3>
void ProcessV(
    std::vector<DataType>& vdata,
    const CellType1& c1, const CellType2& c2, const CellType3& c3,
    const double d3, const bool swap12, const double u)
// Preconditions: s1+s2/d3 < b
//                s3/d2 < b/u
{
#define altb (b/(1.-b))

    XAssert(NoSplit(c1,c2,d3));

    const double d1 = Dist(c3.MeanPos(),c2.MeanPos());
    const double d2 = Dist(c1.MeanPos(),c3.MeanPos());

    if (d3 > d2) return;
    if (d3 > d1) return;

#ifdef ALLDEBUG
    recursen++;
    dbg<<std::string(recursen,'.')<<"Start ProcessV1 "<<d1<<" >= "<<d2<<" >= "<<d3<<std::endl;
    dbg<<std::string(recursen,'.')<<"sizes = "<<c1.Size()<<" -- "<<c2.Size()<<" -- "<<c3.Size()<<std::endl;
#endif

    XAssert(Check(c1,c2,c3,d1,d2,d3));
    XAssert(std::abs(u-d3/d2) < altb);

    if (d1 >= d2) {
        ProcessV(vdata,c1,c2,c3,d1,d2,d3,swap12,u);
    } else {
        ProcessV(vdata,c2,c1,c3,d2,d1,d3,!swap12,u);
    } 
#ifdef ALLDEBUG
    dbg<<std::string(recursen,'.')<<"Done PV3\n";
    recursen--;
#endif
}

template <class DataType, class CellType1, class CellType2, class CellType3>
void ProcessU(
    std::vector<std::vector<DataType> >& uvdata,
    const CellType1& c1, const CellType2& c2, const CellType3& c3,
    const double d3, const bool swap12);

template <class DataType, class CellType1, class CellType2, class CellType3>
void ProcessU(
    std::vector<std::vector<DataType> >& uvdata,
    const CellType1& c1, const CellType2& c2, const CellType3& c3,
    const double d1, const double d2, const double d3, const bool swap12)
// Preconditions: s1+s2/d3 < b
//                d1 >= d2 >= d3
{
#ifdef ALLDEBUG
    recursen++;
    dbg<<std::string(recursen,'.')<<"Start ProcessU3 "<<d1<<" >= "<<d2<<" >= "<<d3<<std::endl;
    dbg<<std::string(recursen,'.')<<"sizes = "<<c1.Size()<<" -- "<<c2.Size()<<" -- "<<c3.Size()<<std::endl;
#endif

    XAssert(NoSplit(c1,c2,d3));
    XAssert(Check(c1,c2,c3,d1,d2,d3));

    Assert(d1 >= d2);
    Assert(d2 >= d3);

    // u = d3/d2
    // u_max = d3/(d2-s) = d3/d2 * 1/(1-s/d2)
    // du = u_max - u = u * s/d2 / (1-s/d2)
    // du < b ==> s/d2 < m / (u+m)
    // Note: u-u_min < max ==> s/d2 < m / (u-m)
    //       which is less stringent

    const double u = d3/d2;
    const bool split3 = c3.Size()*(u+b) > d2*b;

    if (split3) {
        Assert(c3.Size()>0);
        Assert(c3.Left());
        Assert(c3.Right());
        ProcessU(uvdata,c1,c2,*c3.Left(),d3,swap12);
        ProcessU(uvdata,c1,c2,*c3.Right(),d3,swap12);
    } else {
        int ku = int(floor(u/du));
        if (ku >= nubins) {Assert(ku==nubins); ku--;}
        std::vector<DataType>& vdata = uvdata[ku];

        ProcessV(vdata,c1,c2,c3,d1,d2,d3,swap12,u);
    }
#ifdef ALLDEBUG
    dbg<<std::string(recursen,'.')<<"Done PU3\n";
    recursen--;
#endif
}

template <class DataType, class CellType1, class CellType2, class CellType3>
void ProcessU(
    std::vector<std::vector<DataType> >& uvdata,
    const CellType1& c1, const CellType2& c2, const CellType3& c3,
    const double d3, const bool swap12)
// Preconditions: s1+s2/d3 < b
{
#ifdef ALLDEBUG
    recursen++;
    dbg<<std::string(recursen,'.')<<"Start ProcessU1 "<<Dist(c2.MeanPos(),c3.MeanPos())<<" , "<<Dist(c1.MeanPos(),c3.MeanPos())<<" , "<<d3<<std::endl;
    dbg<<std::string(recursen,'.')<<"sizes = "<<c1.Size()<<" -- "<<c2.Size()<<" -- "<<c3.Size()<<std::endl;
#endif
    const double d2 = Dist(c1.MeanPos(),c3.MeanPos());
    if (d2 < d3) {
        if (d2 + c3.Size() < d3) {
#ifdef ALLDEBUG
            recursen--;
#endif
            return;
        } else {
            Assert(c3.Size()>0);
            Assert(c3.Left());
            Assert(c3.Right());
            ProcessU(uvdata,c1,c2,*c3.Left(),d3,swap12);
            ProcessU(uvdata,c1,c2,*c3.Right(),d3,swap12);
        }
    } else {
        const double d1 = Dist(c3.MeanPos(),c2.MeanPos());
        if (d1 < d3) {
            if (d1 + c3.Size() < d3) {
#ifdef ALLDEBUG
                recursen--;
#endif
                return;
            } else {
                Assert(c3.Size()>0);
                Assert(c3.Left());
                Assert(c3.Right());
                ProcessU(uvdata,c1,c2,*c3.Left(),d3,swap12);
                ProcessU(uvdata,c1,c2,*c3.Right(),d3,swap12);
            }
        } else {

            XAssert(NoSplit(c1,c2,d3));
            XAssert(Check(c1,c2,c3,d1,d2,d3));

            if (d1 >= d2) {
                ProcessU(uvdata,c1,c2,c3,d1,d2,d3,swap12);
            } else {
                ProcessU(uvdata,c2,c1,c3,d2,d1,d3,!swap12);
            }
        }
    }
#ifdef ALLDEBUG
    dbg<<std::string(recursen,'.')<<"Done PU1\n";
    recursen--;
#endif
}

template <class DataType, class CellType1, class CellType2, class CellType3>
void Process111(
    std::vector<std::vector<std::vector<DataType> > >& data,
    double minr, double maxr,
    const CellType1& c1, const CellType2& c2, const CellType3& c3)
// Does all triangles with 1 point each in c1, c2, c3
// for which d3 is the smallest side
{
    static const double logminsep = log(minsep);

    const typename CellType1::PosType r3 = c2.MeanPos() - c1.MeanPos();
    const double d3 = Dist(c2.MeanPos(),c1.MeanPos());
    const double s1ps2 = c1.AllSize()+c2.AllSize();
    if (d3+s1ps2 < minr) return;
    if (d3-s1ps2 > maxr) return;

    recursen++;
    std::ostream* tempdbg = 0;
    if (dbgout) {
        tempdbg = dbgout;
#ifndef ALLDEBUG
        if (std::max(c1.Size(),c2.Size()) < outputsize) dbgout = 0;
#endif
        dbg<<std::string(recursen,'.')<<"Start Process111 d3 = "<<d3;
        dbg<<"  sizes = "<<c1.Size()<<" -- "<<c2.Size()<<std::endl;
    }

    bool split1 = false, split2 = false;
    CalcSplit(split1,split2,c1,c2,d3);

    if (split1) {
        if (split2) {
            // split 1,2
#ifdef ALLDEBUG
            dbg<<std::string(recursen,'.')<<"split 1, 2\n";
#endif
            Assert(c1.Left());
            Assert(c1.Right());
            Assert(c2.Left());
            Assert(c2.Right());
            Process111(data,minr,maxr,*c1.Left(),*c2.Left(),c3);
            Process111(data,minr,maxr,*c1.Left(),*c2.Right(),c3);
            Process111(data,minr,maxr,*c1.Right(),*c2.Left(),c3);
            Process111(data,minr,maxr,*c1.Right(),*c2.Right(),c3);
        } else {
            // split 1 only
#ifdef ALLDEBUG
            dbg<<std::string(recursen,'.')<<"split 1\n";
#endif
            Assert(c1.Left());
            Assert(c1.Right());
            Process111(data,minr,maxr,*c1.Left(),c2,c3);
            Process111(data,minr,maxr,*c1.Right(),c2,c3);
        }
    } else {
        if (split2) {
            // split 2 only
#ifdef ALLDEBUG
            dbg<<std::string(recursen,'.')<<"split 2\n";
#endif
            Assert(c2.Left());
            Assert(c2.Right());
            Process111(data,minr,maxr,c1,*c2.Left(),c3);
            Process111(data,minr,maxr,c1,*c2.Right(),c3);
        } else if (d3 <= maxr) {
            // don't split any
#ifdef ALLDEBUG
            dbg<<std::string(recursen,'.')<<"do uv\n";
#endif

            double logr = log(d3);
            const int kr = int(floor((logr-logminsep)/binsize));
            Assert(kr >= 0);
            Assert(kr < nrbins);

            std::vector<std::vector<DataType> >& datakr = data[kr];
            std::vector<std::vector<DataType> > tempuvdata(
                datakr.size(),
                std::vector<DataType>(datakr[0].size()));

            ProcessU(tempuvdata,c1,c2,c3,d3,false);

            Add(datakr,tempuvdata,c1,c2,r3,d3);

        } else {
#ifdef ALLDEBUG
            dbg<<std::string(recursen,'.')<<"d3 > maxsep: "<<d3<<" > "<<maxsep<<std::endl;
#endif
        }
    }

    if (tempdbg) {
        dbg<<std::string(recursen,'.')<<"Done P111\n";
        dbgout = tempdbg;
    }
    recursen--;
}

template <class DataType, class CellType12, class CellType3>
void Process21(
    std::vector<std::vector<std::vector<DataType> > >& data,
    double minr, double maxr,
    const CellType12& c12, const CellType3& c3)
// Does all triangles with 2 points in c12 and 3rd point in c3
// for which d3 is the smallest side
{
    const double minsize = minr/2.;
    if (c12.Size() < minsize) return;

    recursen++;
    std::ostream* tempdbg = 0;
    if (dbgout) {
        tempdbg = dbgout;
#ifndef ALLDEBUG
        if (c12.Size() < outputsize) dbgout = 0;
#endif
        dbg<<std::string(recursen,'.')<<"Start Process21: size = "<<c12.Size();
        dbg<<"  MeanPos = "<<c12.MeanPos()<<std::endl;
    }

    Assert(c12.Left());
    Assert(c12.Right());

    Process21(data,minr,maxr,*c12.Left(),c3);
    Process21(data,minr,maxr,*c12.Right(),c3);
    Process111(data,minr,maxr,*c12.Left(),*c12.Right(),c3);

    if (tempdbg) {
        dbg<<std::string(recursen,'.')<<"Done P21: size = "<<c12.Size()<<std::endl;
        dbgout = tempdbg;
    }
    recursen--;
}

template <class DataType, class CellType>
void Process3(
    std::vector<std::vector<std::vector<DataType> > >& data,
    double minr, double maxr, const CellType& c123)
// Does all triangles with 3 points in c123 for which d3 is the smallest side
{ 
    static const double sqrt3 = sqrt(3.);
    const double minsize = minr/sqrt3;
    if (c123.Size() < minsize) return;

    recursen++;
    std::ostream* tempdbg = 0;
    if (dbgout) {
        tempdbg = dbgout;
#ifndef ALLDEBUG
        if (c123.Size() < outputsize) dbgout = 0;
#endif
        dbg<<std::string(recursen,'.')<<"Start Process3: size = "<<c123.Size();
        dbg<<"  MeanPos = "<<c123.MeanPos()<<std::endl;
    }

    Assert(c123.Left());
    Assert(c123.Right());
    Process21(data,minr,maxr,*c123.Left(),*c123.Right()); 
    Process21(data,minr,maxr,*c123.Right(),*c123.Left()); 
    Process3(data,minr,maxr,*c123.Left()); 
    Process3(data,minr,maxr,*c123.Right()); 

    if (tempdbg) {
        dbg<<std::string(recursen,'.')<<"Done P3: size = "<<c123.Size()<<std::endl;
        dbgout = tempdbg;
    }
    recursen--;
}
