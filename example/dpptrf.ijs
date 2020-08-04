NB. subroutine dpptrf ( character                         UPLO,
NB.                     integer                           N,
NB.                     double precision, dimension( * )  AP,
NB.                     integer                           INFO
NB.                   )
NB.
NB. DPPTRF
NB.
NB. Purpose:
NB.
NB.      DPPTRF computes the Cholesky factorization of a real symmetric
NB.      positive definite matrix A stored in packed format.
NB.
NB.      The factorization has the form
NB.         A = U**T * U,  if UPLO = 'U', or
NB.         A = L  * L**T,  if UPLO = 'L',
NB.      where U is an upper triangular matrix and L is lower triangular.
NB.
NB. Parameters
NB.
NB.                             UPLO is CHARACTER*1
NB.     [in]     UPLO           = 'U':  Upper triangle of A is stored;
NB.                             = 'L':  Lower triangle of A is stored.
NB.
NB.                             N is INTEGER
NB.     [in]     N              The order of the matrix A.  N >= 0.
NB.
NB.                             AP is DOUBLE PRECISION array, dimension (N*(N+1)/2)
NB.                             On entry, the upper or lower triangle of the symmetric matrix
NB.                             A, packed columnwise in a linear array.  The j-th column of A
NB.                             is stored in the array AP as follows:
NB.                             if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
NB.     [in,out] AP             if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.
NB.                             See below for further details.
NB.
NB.                             On exit, if INFO = 0, the triangular factor U or L from the
NB.                             Cholesky factorization A = U**T*U or A = L*L**T, in the same
NB.                             storage format as A.
NB.
NB.                             INFO is INTEGER
NB.                             = 0:  successful exit
NB.     [out]    INFO           < 0:  if INFO = -i, the i-th argument had an illegal value
NB.                             > 0:  if INFO = i, the leading minor of order i is not
NB.                                   positive definite, and the factorization could not be
NB.                                   completed.

NB. subroutine dpptrs ( character                              UPLO,
NB.                     integer                                N,
NB.                     integer                                NRHS,
NB.                     double precision, dimension( * )       AP,
NB.                     double precision, dimension( ldb, * )  B,
NB.                     integer                                LDB,
NB.                     integer                                INFO
NB.                   )
NB.
NB. DPPTRS
NB.
NB. Purpose:
NB.
NB.      DPPTRS solves a system of linear equations A*X = B with a symmetric
NB.      positive definite matrix A in packed storage using the Cholesky
NB.      factorization A = U**T*U or A = L*L**T computed by DPPTRF.
NB.
NB. Parameters
NB.
NB.                             UPLO is CHARACTER*1
NB.     [in]     UPLO           = 'U':  Upper triangle of A is stored;
NB.                             = 'L':  Lower triangle of A is stored.
NB.
NB.                             N is INTEGER
NB.     [in]     N              The order of the matrix A.  N >= 0.
NB.
NB.                             NRHS is INTEGER
NB.     [in]     NRHS           The number of right hand sides, i.e., the number of columns
NB.                             of the matrix B.  NRHS >= 0.
NB.
NB.                             AP is DOUBLE PRECISION array, dimension (N*(N+1)/2)
NB.                             The triangular factor U or L from the Cholesky factorization
NB.                             A = U**T*U or A = L*L**T, packed columnwise in a linear
NB.     [in]     AP             array.  The j-th column of U or L is stored in the array AP
NB.                             as follows:
NB.                             if UPLO = 'U', AP(i + (j-1)*j/2) = U(i,j) for 1<=i<=j;
NB.                             if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = L(i,j) for j<=i<=n.
NB.
NB.                             B is DOUBLE PRECISION array, dimension (LDB,NRHS)
NB.     [in,out] B              On entry, the right hand side matrix B.
NB.                             On exit, the solution matrix X.
NB.
NB.                             LDB is INTEGER
NB.     [in]     LDB            The leading dimension of the array B.  LDB >= max(1,N).
NB.
NB.                             INFO is INTEGER
NB.     [out]    INFO           = 0:  successful exit
NB.                             < 0:  if INFO = -i, the i-th argument had an illegal value

require 'math/lapack2'

NB. cholesky factorization sample data
NB.
NB. solve AX = B where
NB. A is 4 x 4 symmetric matric
NB.   4.16
NB.  _3.12   5.03
NB.   0.56  _0.83   0.76
NB.  _0.10   1.18   0.34   1.18
NB.
NB. B is 2 x 4 matrix
NB.   8.70   8.30
NB. _13.35   2.13
NB.   1.89   1.61
NB.  _4.14   5.00

do_dpptrf=: 3 : 0
'uplo a b'=. y
assert. 1=#$a      NB. packed format
n=. , <:>.%:+:#a   NB. order of unpacked format matrix , dimension (N*(N+1)/2)
assert. (-:n*n+1)=#a
assert. 2=#$b
assert. n={.$b
nrhs=. ,{:$b

assert. 0= LASTINFO=: _1{::cdrc=. dpptrf_jlapack2_ (,uplo);n;a;,_1

ap=. 3{::cdrc         NB. packed U or L

assert. 0= LASTINFO=: _1{::cdrc=. dpptrs_jlapack2_ (,uplo);n;nrhs;ap;(|:b);(1>.n);,_1

|: 5{::cdrc
)

NB. symmetric positive definite matrix
a=: ".;._2[0 : 0
 4.16 _3.12  0.56 _0.1
_3.12  5.03 _0.83 1.18
 0.56 _0.83  0.76 0.34
 _0.1  1.18  0.34 1.18
)

NB. lower triagular packed format
ap=: 4.16 _3.12 0.56 _0.10 5.03 _0.83 1.18 0.76 0.34 1.18

NB. rhs
b=: 4 2 $ 8.70 8.30 _13.35 2.13 1.89 1.61 _4.14 5.00

x=: do_dpptrf 'L';ap;b
echo a ; b ; x ; a (+/ .*) x

