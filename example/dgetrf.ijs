NB. subroutine dgetrf ( integer                                M,
NB.                     integer                                N,
NB.                     double precision, dimension( lda, * )  A,
NB.                     integer                                LDA,
NB.                     integer, dimension( * )                IPIV,
NB.                     integer                                INFO
NB.                   )
NB.
NB. DGETRF
NB.
NB. DGETRF VARIANT: iterative version of Sivan Toledo's recursive LU algorithm
NB.
NB. DGETRF VARIANT: left-looking Level 3 BLAS version of the algorithm.
NB.
NB. Purpose:
NB.
NB.      DGETRF computes an LU factorization of a general M-by-N matrix A
NB.      using partial pivoting with row interchanges.
NB.
NB.      The factorization has the form
NB.         A = P * L * U
NB.      where P is a permutation matrix, L is lower triangular with unit
NB.      diagonal elements (lower trapezoidal if m > n), and U is upper
NB.      triangular (upper trapezoidal if m < n).
NB.
NB.      This is the right-looking Level 3 BLAS version of the algorithm.
NB.
NB. Parameters
NB.
NB.                             M is INTEGER
NB.     [in]     M              The number of rows of the matrix A.  M >= 0.
NB.
NB.                             N is INTEGER
NB.     [in]     N              The number of columns of the matrix A.  N >= 0.
NB.
NB.                             A is DOUBLE PRECISION array, dimension (LDA,N)
NB.                             On entry, the M-by-N matrix to be factored.
NB.     [in,out] A              On exit, the factors L and U from the factorization
NB.                             A = P*L*U; the unit diagonal elements of L are not stored.
NB.
NB.                             LDA is INTEGER
NB.     [in]     LDA            The leading dimension of the array A.  LDA >= max(1,M).
NB.
NB.                             IPIV is INTEGER array, dimension (min(M,N))
NB.     [out]    IPIV           The pivot indices; for 1 <= i <= min(M,N), row i of the
NB.                             matrix was interchanged with row IPIV(i).
NB.
NB.                             INFO is INTEGER
NB.                             = 0:  successful exit
NB.                             < 0:  if INFO = -i, the i-th argument had an illegal value
NB.     [out]    INFO           > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
NB.                                   has been completed, but the factor U is exactly
NB.                                   singular, and division by zero will occur if it is used
NB.                                   to solve a system of equations.

NB. subroutine dgetri ( integer                                N,
NB.                     double precision, dimension( lda, * )  A,
NB.                     integer                                LDA,
NB.                     integer, dimension( * )                IPIV,
NB.                     double precision, dimension( * )       WORK,
NB.                     integer                                LWORK,
NB.                     integer                                INFO
NB.                   )
NB.
NB. DGETRI
NB.
NB. Purpose:
NB.
NB.      DGETRI computes the inverse of a matrix using the LU factorization
NB.      computed by DGETRF.
NB.
NB.      This method inverts U and then computes inv(A) by solving the system
NB.      inv(A)*L = inv(U) for inv(A).
NB.
NB. Parameters
NB.
NB.                              N is INTEGER
NB.     [in]     N               The order of the matrix A.  N >= 0.
NB.
NB.                              A is DOUBLE PRECISION array, dimension (LDA,N)
NB.                              On entry, the factors L and U from the factorization
NB.     [in,out] A               A = P*L*U as computed by DGETRF.
NB.                              On exit, if INFO = 0, the inverse of the original matrix A.
NB.
NB.                              LDA is INTEGER
NB.     [in]     LDA             The leading dimension of the array A.  LDA >= max(1,N).
NB.
NB.                              IPIV is INTEGER array, dimension (N)
NB.     [in]     IPIV            The pivot indices from DGETRF; for 1<=i<=N, row i of the
NB.                              matrix was interchanged with row IPIV(i).
NB.
NB.                              WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
NB.     [out]    WORK            On exit, if INFO=0, then WORK(1) returns the optimal LWORK.
NB.
NB.                              LWORK is INTEGER
NB.                              The dimension of the array WORK.  LWORK >= max(1,N).
NB.                              For optimal performance LWORK >= N*NB, where NB is
NB.                              the optimal blocksize returned by ILAENV.
NB.     [in]     LWORK
NB.                              If LWORK = -1, then a workspace query is assumed; the routine
NB.                              only calculates the optimal size of the WORK array, returns
NB.                              this value as the first entry of the WORK array, and no error
NB.                              message related to LWORK is issued by XERBLA.
NB.
NB.                              INFO is INTEGER
NB.                              = 0:  successful exit
NB.     [out]    INFO            < 0:  if INFO = -i, the i-th argument had an illegal value
NB.                              > 0:  if INFO = i, U(i,i) is exactly zero; the matrix is
NB.                                    singular and its inverse could not be computed.

require 'math/lapack2'

NB. inverse of matrix using LU factorization sample data
NB.
NB. 1 2 3
NB. 6 5 4
NB. 8 9 7

do_dgetrf=: 3 : 0
assert. 2=#$y
assert. =/$y   NB. square matric
a=. y
'm n'=. ,"0 $a
mn=. m<.n

NB. lapack expect column major order |:a
assert. 0= LASTINFO=: _1{::cdrc=. dgetrf_jlapack2_ m;n;(|:a);(1>.m);(mn$_1);,_1

val=. 3{::cdrc         NB. val is in column major order
ipiv=. 5{::cdrc        NB. pivot indices

NB. call with lwork = _1 to query optimal workspace size
NB. val already in column major order, no need to transpose
assert. 0= LASTINFO=: _1{::cdrc=. dgetri_jlapack2_ n;val;(1>.n);ipiv;(1$0.0);(,_1);,_1

lwork=. <. _3{::cdrc

NB. call again with lwork
assert. 0= LASTINFO=: _1{::cdrc=. dgetri_jlapack2_ (_3}.}.cdrc),(lwork$0.0);lwork;,_1

|: 2{::cdrc
)

r=: do_dgetrf a=: 3 3$1.0 2 3 6 5 4 8 9 7
echo a;r; a (+/ .*) r
