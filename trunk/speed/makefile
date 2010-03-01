#EIGEN= -I/home/mjarvis/eigen-2.0.0
#EIGEN= -I/Users/Mike/Astro/eigen-2.0.0
#EIGEN= -I/Users/Mike/Astro/eigen-2.0.10
EIGEN= -I../eigen

#INCLUDE= -I../include $(EIGEN)

#CC= icpc
#CFLAGS= $(INCLUDE) -g -O -w -DTMV_INLINE -DTMV_OPT=3 -openmp -wd279
#CFLAGS= $(INCLUDE) -g -fast -w -DTMV_INLINE -DTMV_OPT=3 -openmp -wd279
#CFLAGS= $(INCLUDE) -g -w -O3 -DTMV_INLINE -DTMV_OPT=3 -wd279
#CFLAGS= $(INCLUDE) -openmp -DNOBLAS -wd279
#LIB= -L../lib -ltmv -L/usr/local/intel/mkl/10.0.1.014/lib/32 -lmkl_ia32 -lguide -lpthread

#INCLUDE= -I../../tags/v0.63/include $(EIGEN) -DTMV_V063
INCLUDE= -I../include $(EIGEN)
CC= g++-4
CFLAGS= $(INCLUDE) -DNDEBUG -g -O2 -DTMV_INLINE -DTMV_OPT=3
#CFLAGS= $(INCLUDE) -DNDEBUG -g -Wall -O2 -DTMV_INLINE -DTMV_OPT=3 -fopenmp
#CFLAGS= $(INCLUDE) -DNDEBUG -g -O2 -DTMV_INLINE -DTMV_OPT=3 -fopenmp
#CFLAGS= $(INCLUDE) -g -O -DTMV_INLINE -DTMV_OPT=3 -fopenmp -DTMV_MEM_TEST
#CFLAGS= $(INCLUDE) -DNDEBUG -g -O2
#CFLAGS= $(INCLUDE) -DNDEBUG -g -O2 -fopenmp
#CFLAGS= $(INCLUDE) -g -O

LIB= 
#LIB= -L../../tags/v0.63/lib -ltmv_symband -ltmv -llapack -lblas -lgfortran -lpthread
#LIB= -L../lib -ltmv -llapack -lblas -lgfortran -lpthread
#LIB= -L../lib -ltmv_3 -llapack -lblas -lgfortran
#LIB= -L../lib -ltmv_mem -llapack -lblas -lgfortran -lpthread
#LIB= -L../lib -llapack -lblas -lgfortran -lpthread

all : tmvspeed2x_always_make 
#all : tmvspeed1 tmvspeed2 tmvspeed3

becker : becker.cpp
	$(CC) $(CFLAGS) becker.cpp -o becker $(LIB)

tmvspeed1_always_make :
	$(CC) $(CFLAGS) TMV_Speed_1.cpp -o tmvspeed1 $(LIB)

tmvspeed1x_always_make :
	$(CC) $(CFLAGS) TMV_Speed_1x.cpp -o tmvspeed1x $(LIB)

tmvspeed2_always_make :
	$(CC) $(CFLAGS) TMV_Speed_2.cpp -o tmvspeed2 $(LIB)

tmvspeed2x_always_make :
	$(CC) $(CFLAGS) TMV_Speed_2x.cpp -o tmvspeed2x $(LIB)

tmvspeed3_always_make :
	$(CC) $(CFLAGS) TMV_Speed_3.cpp -o tmvspeed3 $(LIB)

tmvspeed1 : TMV_Speed_1.cpp
	$(CC) $(CFLAGS) TMV_Speed_1.cpp -o tmvspeed1 $(LIB)

tmvspeed1x : TMV_Speed_1x.cpp
	$(CC) $(CFLAGS) TMV_Speed_1x.cpp -o tmvspeed1x $(LIB)

tmvspeed2 : TMV_Speed_2.cpp
	$(CC) $(CFLAGS) TMV_Speed_2.cpp -o tmvspeed2 $(LIB)

tmvspeed2x : TMV_Speed_2x.cpp
	$(CC) $(CFLAGS) TMV_Speed_2x.cpp -o tmvspeed2x $(LIB)

tmvspeed3 : TMV_Speed_3.cpp
	$(CC) $(CFLAGS) TMV_Speed_3.cpp -o tmvspeed3 $(LIB)

