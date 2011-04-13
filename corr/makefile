CC= g++-4
CFLAGS= -g3 -O2 -fno-strict-aliasing
BINDIR= .
ODIR= .

all : $(BINDIR)/corrnn3d $(BINDIR)/corrnn $(BINDIR)/corree $(BINDIR)/correeea $(BINDIR)/correeeb $(BINDIR)/corrne $(BINDIR)/corrnn $(BINDIR)/corr_norm $(BINDIR)/e3_m3 $(BINDIR)/corrnnna $(BINDIR)/corree_multi

$(BINDIR)/e3_m3: $(ODIR)/E3_M3.o $(ODIR)/CorrIO.o $(ODIR)/CalcT.o
	$(CC) $(CFLAGS) $(ODIR)/E3_M3.o $(ODIR)/CorrIO.o $(ODIR)/CalcT.o -lm -o $(BINDIR)/e3_m3

$(BINDIR)/corrnnna: $(ODIR)/CorrNNNa.o $(ODIR)/CorrIO.o $(ODIR)/CalcT.o
	$(CC) $(CFLAGS) $(ODIR)/CorrNNNa.o $(ODIR)/CorrIO.o $(ODIR)/CalcT.o -lm -o $(BINDIR)/corrnnna

$(BINDIR)/correeea: $(ODIR)/CorrEEEa.o $(ODIR)/CorrIO.o $(ODIR)/CalcT.o
	$(CC) $(CFLAGS) $(ODIR)/CorrEEEa.o $(ODIR)/CorrIO.o $(ODIR)/CalcT.o -lm -o $(BINDIR)/correeea

$(BINDIR)/correeeb: $(ODIR)/CorrEEEb.o $(ODIR)/CorrIO.o $(ODIR)/CalcT.o
	$(CC) $(CFLAGS) $(ODIR)/CorrEEEb.o $(ODIR)/CorrIO.o $(ODIR)/CalcT.o -lm -o $(BINDIR)/correeeb

$(BINDIR)/corrneea: $(ODIR)/CorrNEEa.o $(ODIR)/CorrIO.o $(ODIR)/CalcT.o
	$(CC) $(CFLAGS) $(ODIR)/CorrNEEa.o $(ODIR)/CorrIO.o $(ODIR)/CalcT.o -lm -o $(BINDIR)/corrneea

$(BINDIR)/corree: $(ODIR)/CorrEE.o $(ODIR)/CorrIO.o $(ODIR)/CalcT.o 
	$(CC) $(CFLAGS) $(ODIR)/CorrEE.o $(ODIR)/CorrIO.o $(ODIR)/CalcT.o -lm -o $(BINDIR)/corree

$(BINDIR)/corrne: $(ODIR)/CorrNE.o $(ODIR)/CorrIO.o $(ODIR)/CalcT.o 
	$(CC) $(CFLAGS) $(ODIR)/CorrNE.o $(ODIR)/CorrIO.o $(ODIR)/CalcT.o -lm -o $(BINDIR)/corrne

$(BINDIR)/corrnn: $(ODIR)/CorrNN.o $(ODIR)/CorrIO.o $(ODIR)/CalcT.o 
	$(CC) $(CFLAGS) $(ODIR)/CorrNN.o $(ODIR)/CorrIO.o $(ODIR)/CalcT.o -lm -o $(BINDIR)/corrnn

$(BINDIR)/corrnn3d: $(ODIR)/CorrNN3D.o $(ODIR)/CorrIO.o $(ODIR)/CalcT.o 
	$(CC) $(CFLAGS) $(ODIR)/CorrNN3D.o $(ODIR)/CorrIO.o $(ODIR)/CalcT.o -lm -o $(BINDIR)/corrnn3d

$(BINDIR)/corr_norm: $(ODIR)/Corr_Norm.o $(ODIR)/CorrIO.o $(ODIR)/CalcT.o 
	$(CC) $(CFLAGS) $(ODIR)/Corr_Norm.o $(ODIR)/CorrIO.o $(ODIR)/CalcT.o -lm -o $(BINDIR)/corr_norm

$(BINDIR)/corree_multi: $(ODIR)/CorrEE_multi.o $(ODIR)/CorrIO.o $(ODIR)/CalcT.o 
	$(CC) $(CFLAGS) $(ODIR)/CorrEE_multi.o $(ODIR)/CorrIO.o $(ODIR)/CalcT.o -lm -o $(BINDIR)/corree_multi

$(ODIR)/E3_M3.o: E3_M3.cpp CorrIO.h Corr.h Form.h
	$(CC) $(CFLAGS) -c E3_M3.cpp -o $(ODIR)/E3_M3.o

$(ODIR)/CorrNNNa.o: CorrNNNa.cpp Cell.h CorrIO.h Corr.h Process3.h
	$(CC) $(CFLAGS) -c CorrNNNa.cpp -o $(ODIR)/CorrNNNa.o

$(ODIR)/CorrEEEa.o: CorrEEEa.cpp Cell.h CorrIO.h Corr.h Process3.h
	$(CC) $(CFLAGS) -c CorrEEEa.cpp -o $(ODIR)/CorrEEEa.o

$(ODIR)/CorrEEEb.o: CorrEEEb.cpp Cell.h CorrIO.h Corr.h PairCell.h Process2.h
	$(CC) $(CFLAGS) -c CorrEEEb.cpp -o $(ODIR)/CorrEEEb.o

$(ODIR)/CorrEE.o: CorrEE.cpp Cell.h CorrIO.h Corr.h Process2.h
	$(CC) $(CFLAGS) -c CorrEE.cpp -o $(ODIR)/CorrEE.o

$(ODIR)/CorrNE.o: CorrNE.cpp Cell.h CorrIO.h Corr.h Process2.h
	$(CC) $(CFLAGS) -c CorrNE.cpp -o $(ODIR)/CorrNE.o

$(ODIR)/CorrNN.o: CorrNN.cpp Cell.h CorrIO.h Corr.h Process2.h
	$(CC) $(CFLAGS) -c CorrNN.cpp -o $(ODIR)/CorrNN.o

$(ODIR)/CorrNN3D.o: CorrNN3D.cpp Cell3D.h CorrIO.h Corr.h Process2.h
	$(CC) $(CFLAGS) -c CorrNN3D.cpp -o $(ODIR)/CorrNN3D.o

$(ODIR)/Corr_Norm.o: Corr_Norm.cpp Cell.h CorrIO.h Corr.h Process2.h
	$(CC) $(CFLAGS) -c Corr_Norm.cpp -o $(ODIR)/Corr_Norm.o

$(ODIR)/CorrEE_multi.o: CorrEE_multi.cpp Cell.h CorrIO.h Corr.h Process2.h
	$(CC) $(CFLAGS) -c CorrEE_multi.cpp -o $(ODIR)/CorrEE_multi.o

$(ODIR)/CorrIO.o: CorrIO.cpp Cell.h CorrIO.h Corr.h
	$(CC) $(CFLAGS) -c CorrIO.cpp -o $(ODIR)/CorrIO.o

$(ODIR)/CalcT.o: CalcT.cpp CalcT.h
	$(CC) $(CFLAGS) -c CalcT.cpp -o $(ODIR)/CalcT.o
