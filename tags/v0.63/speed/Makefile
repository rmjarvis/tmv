CC= g++-4
INCLUDE= -I../include
CFLAGS= $(INCLUDE) -O2 -DNDEBUG -fopenmp
LIBS= -L../lib -ltmv -lcblas
LIBFILE= ../lib/libtmv.a

tmvspeed : TMV_Speed_MultMM.cpp $(LIBFILE)
	$(CC) $(CFLAGS) TMV_Speed_MultMM.cpp -o tmvspeed $(LIBS)
