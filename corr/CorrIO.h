#ifndef CORRIO_H
#define CORRIO_H

#include "Corr.h"
#include "Form.h"
#include "MemDebug.h"
#include "Cell.h"
#include "Cell3D.h"

void Read(std::istream& fin, std::vector<CellData>& celldata, double& vare);
void Read(std::istream& fin, std::vector<NCellData>& celldata);
void Read(std::istream& fin, std::vector<TCell3DData>& celldata);
void Read(std::istream& fin, std::vector<NCell3DData>& celldata);

void WriteEEE(
    std::ostream& fout,
    const std::vector<std::vector<std::vector<EEEData> > >& data);
void Read(
    std::istream& fin,
    std::vector<std::vector<std::vector<EEEData> > >& data);
void WriteM3(
    std::ostream& fout,
    const std::vector<std::vector<std::vector<EEEData> > >& data,
    double k1=1., double k2=1., double k3=1.);

void WriteEE(std::ostream& fout, std::vector<EEData>& data);
void WriteM2(std::ostream& fout, std::vector<EEData>& data);

void WriteNorm(
    std::ostream& fout,
    const std::vector<NEData>& crossdata,
    const std::vector<EEData>& twoptdata,
    const std::vector<double>& dd, const std::vector<double>& dr,
    const std::vector<double>& rr);
void WriteNE(std::ostream& fout, const std::vector<NEData>& crossdata);
void WriteNN(
    std::ostream& fout,
    const std::vector<double>& dd, const std::vector<double>& dr, 
    const std::vector<double>& rr, double nrr);

void WriteNNN(
    std::ostream& fout,
    const std::vector<std::vector<std::vector<double> > >& ddd,
    const std::vector<std::vector<std::vector<double> > >& ddr,
    const std::vector<std::vector<std::vector<double> > >& drd,
    const std::vector<std::vector<std::vector<double> > >& rdd,
    const std::vector<std::vector<std::vector<double> > >& drr,
    const std::vector<std::vector<std::vector<double> > >& rdr,
    const std::vector<std::vector<std::vector<double> > >& rrd,
    const std::vector<std::vector<std::vector<double> > >& rrr);
void WriteNNN(
    std::ostream& fout,
    const std::vector<std::vector<std::vector<double> > >& ddd,
    const std::vector<std::vector<std::vector<double> > >& ddr,
    const std::vector<std::vector<std::vector<double> > >& drr,
    const std::vector<std::vector<std::vector<double> > >& rrr);
void WriteNNN(
    std::ostream& fout,
    const std::vector<std::vector<std::vector<double> > >& ddd,
    const std::vector<std::vector<std::vector<double> > >& rrr);

#endif
