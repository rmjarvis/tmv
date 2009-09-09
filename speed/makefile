CC= icpc
INCLUDE= -I../include -I/home/mjarvis/eigen-2.0.0
#CFLAGS= $(INCLUDE) -DNDEBUG -g -O3 -DTMV_INLINE -DTMV_OPT=3 -fopenmp
CFLAGS= $(INCLUDE) -g -O -DTMV_INLINE -DTMV_OPT=3 -fopenmp
#CFLAGS= $(INCLUDE) -g -O -DTMV_INLINE -DTMV_OPT=3 -fopenmp -DTMV_MEM_TEST
#CFLAGS= $(INCLUDE) -DNDEBUG -g -O3
#CFLAGS= $(INCLUDE) -DNDEBUG -g -O3 -fopenmp
#CFLAGS= $(INCLUDE) -g -O
#LIB= -L../lib -ltmv_3 -llapack -lblas -lgfortran
#LIB= -L../lib -ltmv_mem -llapack -lblas -lgfortran -lpthread
LIB= -L../lib -llapack -lblas -lgfortran -lpthread

#all : tmvspeed3_always_make 
all : tmvspeed1 tmvspeed2 tmvspeed3

tmvspeed1_always_make :
	$(CC) $(CFLAGS) TMV_Speed_1.cpp -o tmvspeed1 $(LIB)

tmvspeed2_always_make :
	$(CC) $(CFLAGS) TMV_Speed_2.cpp -o tmvspeed2 $(LIB)

tmvspeed3_always_make :
	$(CC) $(CFLAGS) TMV_Speed_3.cpp -o tmvspeed3 $(LIB)

tmvspeed1 : TMV_Speed_1.cpp
	$(CC) $(CFLAGS) TMV_Speed_1.cpp -o tmvspeed1 $(LIB)

tmvspeed2 : TMV_Speed_2.cpp
	$(CC) $(CFLAGS) TMV_Speed_2.cpp -o tmvspeed2 $(LIB)

tmvspeed3 : TMV_Speed_3.cpp
	$(CC) $(CFLAGS) TMV_Speed_3.cpp -o tmvspeed3 $(LIB)

