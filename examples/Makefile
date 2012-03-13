CC= g++
INCLUDE= -I../include
CFLAGS= $(INCLUDE) -O2 -DNDEBUG
TMVLINK := $(shell cat ../share/tmv-link)
SBTMVLINK = -ltmv_symband $(TMVLINK)
LIBS = $(TMVLINK) 
SBLIBS = $(SBTMVLINK) 

#all: vector matrix division bandmatrix symmatrix
all: vector matrix division 
	check_comments

Vector.o : Vector.cpp $(LIB)
	$(CC) $(CFLAGS) -c Vector.cpp -o Vector.o

vector : Vector.o 
	$(CC) Vector.o -o vector $(LIBS)

Matrix.o : Matrix.cpp $(LIB)
	$(CC) $(CFLAGS) -c Matrix.cpp -o Matrix.o

matrix : Matrix.o 
	$(CC) Matrix.o -o matrix $(LIBS)

Division.o : Division.cpp $(LIB)
	$(CC) $(CFLAGS) -c Division.cpp -o Division.o

division : Division.o 
	$(CC) Division.o -o division $(LIBS)

BandMatrix.o : BandMatrix.cpp $(LIB)
	$(CC) $(CFLAGS) -c BandMatrix.cpp -o BandMatrix.o

bandmatrix : BandMatrix.o 
	$(CC) BandMatrix.o -o bandmatrix $(LIBS)

SymMatrix.o : SymMatrix.cpp $(LIB)
	$(CC) $(CFLAGS) -c SymMatrix.cpp -o SymMatrix.o

symmatrix : SymMatrix.o 
	$(CC) SymMatrix.o -o symmatrix $(LIBS)

