#EIGEN= -I/home/mjarvis/eigen-2.0.0
#EIGEN= -I/Users/Mike/Astro/eigen-2.0.12
EIGEN= -I../eigen
#EIGEN= -I/Users/Mike/Astro/eigen-default

#INCLUDE= -I../include $(EIGEN)

#CC= icpc
#CFLAGS= $(INCLUDE) -g3 -O -w -DTMV_INLINE -DTMV_OPT=3 -openmp -wd279
#CFLAGS= $(INCLUDE) -g3 -fast -w -DTMV_INLINE -DTMV_OPT=3 -openmp -wd279
#CFLAGS= $(INCLUDE) -g3 -w -O3 -DTMV_INLINE -DTMV_OPT=3 -wd279
#CFLAGS= $(INCLUDE) -openmp -DNOBLAS -wd279
#LIB= -L../lib -ltmv -L/usr/local/intel/mkl/10.0.1.014/lib/32 -lmkl_ia32 -lguide -lpthread

#INCLUDE= -I../../tags/v0.63/include $(EIGEN) -DTMV_V063
INCLUDE= -I../include $(EIGEN)
CC= g++-4
#CFLAGS= $(INCLUDE) -DNDEBUG -g3 -O2 -DTMV_INLINE -DTMV_OPT=3
#CFLAGS= $(INCLUDE) -DNDEBUG -g3 -O2 -DTMV_OPT=3 -Wall -fopenmp
#CFLAGS= $(INCLUDE) -DNDEBUG -g3 -O2 -DTMV_INLINE -DTMV_OPT=3 -fopenmp
#CFLAGS= $(INCLUDE) -DNDEBUG -g3 -O2 -DTMV_OPT=3 -fopenmp
#CFLAGS= $(INCLUDE) -g3 -O -DTMV_INLINE -DTMV_OPT=3 -fopenmp -DTMV_MEM_TEST
CFLAGS= $(INCLUDE) -DNDEBUG -g3 -O2 -fopenmp
#CFLAGS= $(INCLUDE) -g3 -O

#LIB= 
#LIB= -L../../tags/v0.63/lib -ltmv_symband -ltmv -llapack -lblas -lgfortran -lpthread
LIB= -L../lib -ltmv -llapack -lblas -lgfortran -lpthread
#LIB= -L../lib -ltmv_3 -llapack -lblas -lgfortran
#LIB= -L../lib -ltmv_mem -llapack -lblas -lgfortran -lpthread
#LIB= -L../lib -llapack -lblas -lgfortran -lpthread

all : tmvspeed4_always_make 
#all : tmvspeed1 tmvspeed2 tmvspeed3 tmvspeed4

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

tmvspeed4_always_make :
	$(CC) $(CFLAGS) -c TMV_Speed_4.cpp -o TMV_Speed_4.o
	$(CC) $(CFLAGS) TMV_Speed_4.o -o tmvspeed4 $(LIB)

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

tmvspeed4 : TMV_Speed_4.cpp
	$(CC) $(CFLAGS) TMV_Speed_4.cpp -o tmvspeed4 $(LIB)

