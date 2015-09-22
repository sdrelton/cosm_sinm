# cosm_sinm
This repository contains MATLAB functions to compute the matrix cosine, sine, and both simultaneously.
The algorithms are taken from a recent paper on the subject by Al-Mohy, Higham, and Relton (see bottom).

## Details
```matlab
% Example
A = randn(10);
schur_fact = 0; % No Schur decomposition
C = cosm(A, schur_fact);
S = sinm(A, schur_fact);
[C, S] = cosmsinm(A, schur_fact);
```

Each function has two input arguments:

1. The matrix at which to compute the function.
2. Flag to determine the use of the Schur decomposition (optional).

The input matrix must be finite and square, whilst the second argument can be one of the following.
- 0 -- No Schur decomposition.
- 1 -- Use a real Schur decomposition where possible.
- 2 -- Use a complex Schur decomposition.

Using a Schur decomposition can be faster when A is large and nonnormal, and can also be more accurate.

You can check that the functions are working by using the test script. In MATLAB run the following.

```matlab
test_cosm_sinm
```

For more detail on the algorithm details and performance please see the following paper.

A. H. Al-Mohy, N. J. Higham, and Samuel D. Relton,
**New Algorithms for the Matrix Sine and Cosine Separately or Simultaneously.** 
_SIAM J. Sci. Comput._, 37(1), A456-A487, 2015.

[Open Access PDF](http://epubs.siam.org/doi/abs/10.1137/140973979)
