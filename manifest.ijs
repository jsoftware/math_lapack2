NB. lapack2 manifest

CAPTION=: 'LAPACK2'

DESCRIPTION=: 0 : 0
LAPACK (Linear Algebra Package) is a set of routines for solving systems of simultaneous linear equations, least-squares solutions of linear systems of equations, eigenvalue problems, and singular value problems. The associated matrix factorizations (LU, Cholesky, QR, SVD, Schur, generalized Schur) are also provided, as are related computations such as reordering of the Schur factorizations and estimating condition numbers.

This addon is a leaner version of another math/lapack addon which is no longer maintained.

Binary for Mac/iOS is provided by the veclib framework.

Binary for Linux, install liblapack3 (or similar) from your distro repository. If available, install libopenblas-base or libatlas3-base which provides an optimized version of BLAS.

For Windows, run getbin_jlapack2_'' to install the shared library.

Both Windows and Android binary provided here use reference BLAS.

Reference BLAS implementation may be orders of magnitude slower than optimized implementations. Build your own optimized BLAS if speed performance is critical.

See wiki page: code.jsoftware.com/wiki/Vocabulary/LAPACK
)

VERSION=: '1.0.11'

RELEASE=: ''

FOLDER=: 'math/lapack2'

FILES=: 0 : 0
lapack2.ijs
lib/readme.txt
example/
test/
)

PLATFORMS=: ''
