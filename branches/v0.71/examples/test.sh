#!/bin/bash

grep "\/\/\!" Vector.cpp | sed s%.*\/\/\!\ %% > test1
vector > test2
diff test1 test2

grep "\/\/\!" Matrix.cpp | sed s%.*\/\/\!\ %% > test1
matrix > test2
diff test1 test2

grep "\/\/\!" Division.cpp | sed s%.*\/\/\!\ %% > test1
division > test2
diff test1 test2

grep "\/\/\!" BandMatrix.cpp | sed s%.*\/\/\!\ %% > test1
bandmatrix > test2
diff test1 test2

grep "\/\/\!" SymMatrix.cpp | sed s%.*\/\/\!\ %% > test1
symmatrix > test2
diff test1 test2
