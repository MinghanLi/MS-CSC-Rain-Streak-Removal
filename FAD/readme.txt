For any questions regarding to this package, please contact Jie Wang, jiewangustc@gmail.com.
This package is only for non-commercial use. All copyrights are reserved to the authors. Before you run the code, please read the following instructions carefully.-------------------------------------------------------------------------------------------------------
The code implements the algorithm proposed in:
"A Highly Parallel Algorithm for Isotropic Total Variation Models",  
Jie Wang, Qingyang Li, Sen Yang, Wei Fan, Peter Wonka and Jieping Ye,
ICML2014
-------------------------------------------------------------------------------------------------------
Please start with the "start.m" file.

The FAD.cpp implements the sequential version of the proposed algorithm, and pFAD_OMP implements
the parallel version with OpenMP. Both of the files have been compiled. 


In case you make revisions of the files and need to recompile them, you can follow the simple instructions below.

To compile FAD.cpp, please type:

mex -O FAD.cpp

To compile the pFAD.cpp, please type:

mex -O pFAD_OMP.cpp COMPFLAGS="/openmp $COMPFLAGS"

The speedup gained by pFAD.cpp is roughly linear to the number of cores.