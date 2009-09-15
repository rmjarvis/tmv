#INCLUDE= -I../include -I/home/mjarvis/eigen-2.0.0

#CC= icpc
#CFLAGS= $(INCLUDE) -g -O -w -DTMV_INLINE -DTMV_OPT=3 -openmp -wd279
#CFLAGS= $(INCLUDE) -g -fast -w -DTMV_INLINE -DTMV_OPT=3 -openmp -wd279
#CFLAGS= $(INCLUDE) -g -w -O3 -DTMV_INLINE -DTMV_OPT=3 -wd279
#CFLAGS= $(INCLUDE) -openmp -DNOBLAS -wd279
#LIB= -L../lib -ltmv -L/usr/local/intel/mkl/10.0.1.014/lib/32 -lmkl_ia32 -lguide -lpthread

INCLUDE= -I../include -I/Users/Mike/Astro/eigen-2.0.0
CC= g++-4
CFLAGS= $(INCLUDE) -DNDEBUG -g -O3 -DTMV_INLINE -DTMV_OPT=3 -fopenmp
#CFLAGS= $(INCLUDE) -g -O -DTMV_INLINE -DTMV_OPT=3 -fopenmp 
#CFLAGS= $(INCLUDE) -g -O -DTMV_INLINE -DTMV_OPT=3 -fopenmp -DTMV_MEM_TEST
#CFLAGS= $(INCLUDE) -DNDEBUG -g -O3
#CFLAGS= $(INCLUDE) -DNDEBUG -g -O3 -fopenmp
#CFLAGS= $(INCLUDE) -g -O

LIB= -L../lib -ltmv -llapack -lblas -lgfortran
#LIB= -L../lib -ltmv_3 -llapack -lblas -lgfortran
#LIB= -L../lib -ltmv_mem -llapack -lblas -lgfortran -lpthread
#LIB= -L../lib -llapack -lblas -lgfortran -lpthread

all : tmvspeed3_always_make 
#all : tmvspeed1 tmvspeed2 tmvspeed3

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

