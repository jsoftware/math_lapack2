NB. subroutine dgeqrf ( integer                                M,
NB.                     integer                                N,
NB.                     double precision, dimension( lda, * )  A,
NB.                     integer                                LDA,
NB.                     double precision, dimension( * )       TAU,
NB.                     double precision, dimension( * )       WORK,
NB.                     integer                                LWORK,
NB.                     integer                                INFO
NB.                   )
NB.
NB. DGEQRF VARIANT: left-looking Level 3 BLAS version of the algorithm.
NB.
NB. Purpose:
NB.
NB.  DGEQRF computes a QR factorization of a real M-by-N matrix A:
NB.  A = Q * R.
NB.
NB.  This is the left-looking Level 3 BLAS version of the algorithm.
NB.
NB. Parameters
NB.
NB.                              M is INTEGER
NB.     [in]     M               The number of rows of the matrix A.  M >= 0.
NB.
NB.                              N is INTEGER
NB.     [in]     N               The number of columns of the matrix A.  N >= 0.
NB.
NB.                              A is DOUBLE PRECISION array, dimension (LDA,N)
NB.                              On entry, the M-by-N matrix A.
NB.                              On exit, the elements on and above the diagonal of the array
NB.                              contain the min(M,N)-by-N upper trapezoidal matrix R (R is
NB.     [in,out] A               upper triangular if m >= n); the elements below the diagonal,
NB.                              with the array TAU, represent the orthogonal matrix Q as a
NB.                              product of min(m,n) elementary reflectors (see Further
NB.                              Details).
NB.
NB.                              LDA is INTEGER
NB.     [in]     LDA             The leading dimension of the array A.  LDA >= max(1,M).
NB.
NB.                              TAU is DOUBLE PRECISION array, dimension (min(M,N))
NB.     [out]    TAU             The scalar factors of the elementary reflectors (see Further
NB.                              Details).
NB.
NB.                              WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
NB.     [out]    WORK            On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
NB.
NB.                              LWORK is INTEGER
NB.
NB.                              The dimension of the array WORK. The dimension can be divided into three parts.
NB.
NB.                              1) The part for the triangular factor T. If the very last T is not bigger
NB.                                 than any of the rest, then this part is NB x ceiling(K/NB), otherwise,
NB.                                 NB x (K-NT), where K = min(M,N) and NT is the dimension of the very last T
NB.
NB.                              2) The part for the very last T when T is bigger than any of the rest T.
NB.                                 The size of this part is NT x NT, where NT = K - ceiling ((K-NX)/NB) x NB,
NB.     [in]     LWORK              where K = min(M,N), NX is calculated by
NB.                                       NX = MAX( 0, ILAENV( 3, 'DGEQRF', ' ', M, N, -1, -1 ) )
NB.
NB.                              3) The part for dlarfb is of size max((N-M)*K, (N-M)*NB, K*NB, NB*NB)
NB.
NB.                              So LWORK = part1 + part2 + part3
NB.
NB.                              If LWORK = -1, then a workspace query is assumed; the routine
NB.                              only calculates the optimal size of the WORK array, returns
NB.                              this value as the first entry of the WORK array, and no error
NB.                              message related to LWORK is issued by XERBLA.
NB.
NB.                              INFO is INTEGER
NB.     [out]    INFO            = 0:  successful exit
NB.                              < 0:  if INFO = -i, the i-th argument had an illegal value

NB. subroutine dorgqr ( integer                                M,
NB.                     integer                                N,
NB.                     integer                                K,
NB.                     double precision, dimension( lda, * )  A,
NB.                     integer                                LDA,
NB.                     double precision, dimension( * )       TAU,
NB.                     double precision, dimension( * )       WORK,
NB.                     integer                                LWORK,
NB.                     integer                                INFO
NB.                   )
NB.
NB. DORGQR
NB.
NB. Purpose:
NB.
NB.      DORGQR generates an M-by-N real matrix Q with orthonormal columns,
NB.      which is defined as the first N columns of a product of K elementary
NB.      reflectors of order M
NB.
NB.            Q  =  H(1) H(2) . . . H(k)
NB.
NB.      as returned by DGEQRF.
NB.
NB. Parameters
NB.
NB.                              M is INTEGER
NB.     [in]     M               The number of rows of the matrix Q. M >= 0.
NB.
NB.                              N is INTEGER
NB.     [in]     N               The number of columns of the matrix Q. M >= N >= 0.
NB.
NB.                              K is INTEGER
NB.     [in]     K               The number of elementary reflectors whose product defines the
NB.                              matrix Q. N >= K >= 0.
NB.
NB.                              A is DOUBLE PRECISION array, dimension (LDA,N)
NB.                              On entry, the i-th column must contain the vector which
NB.                              defines the elementary reflector H(i), for i = 1,2,...,k, as
NB.     [in,out] A               returned by DGEQRF in the first k columns of its array
NB.                              argument A.
NB.                              On exit, the M-by-N matrix Q.
NB.
NB.                              LDA is INTEGER
NB.     [in]     LDA             The first dimension of the array A. LDA >= max(1,M).
NB.
NB.                              TAU is DOUBLE PRECISION array, dimension (K)
NB.     [in]     TAU             TAU(i) must contain the scalar factor of the elementary
NB.                              reflector H(i), as returned by DGEQRF.
NB.
NB.                              WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
NB.     [out]    WORK            On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
NB.
NB.                              LWORK is INTEGER
NB.                              The dimension of the array WORK. LWORK >= max(1,N).
NB.                              For optimum performance LWORK >= N*NB, where NB is the
NB.                              optimal blocksize.
NB.     [in]     LWORK
NB.                              If LWORK = -1, then a workspace query is assumed; the routine
NB.                              only calculates the optimal size of the WORK array, returns
NB.                              this value as the first entry of the WORK array, and no error
NB.                              message related to LWORK is issued by XERBLA.
NB.
NB.                              INFO is INTEGER
NB.     [out]    INFO            = 0:  successful exit
NB.                              < 0:  if INFO = -i, the i-th argument has an illegal value

require 'math/lapack2'

NB. QR factorization sample data
NB.
NB. 12 _51   4
NB.  6 167 _68
NB. _4  24 _41

NB. return q r
do_dgeqrf=: 3 : 0
assert. 2=#$y
a=. y
'm n'=. ,"0 $a
mn=. m<.n

NB. call with lwork = _1 to query optimal workspace size
NB. lapack expect column major order |:a
assert. 0= LASTINFO=: _1{::cdrc=. dgeqrf_jlapack2_ m;n;(|:a);(1>.m);(mn$0.0);(,0.0);(,_1);,_1

lwork=. <. _3{::cdrc

NB. call again with lwork
assert. 0= LASTINFO=: _1{::cdrc=. dgeqrf_jlapack2_ (_3}.}.cdrc),(lwork$0.0);lwork;,_1

val=. 3{::cdrc         NB. val is in column major order
tau=. 5{::cdrc
r=. mn {. utri_jlapack2_ |:val   NB. upper triangular of row major

NB. create orthogonal matrix q from temp val
NB. call with lwork = _1 to query optimal workspace size
NB. val already in column major order, no need to transpose
assert. 0= LASTINFO=: _1{::cdrc=. dorgqr_jlapack2_ m;n;n;val;(1>.m);tau;(,0.0);(,_1);,_1

lwork=. <. _3{::cdrc

NB. call again with lwork
assert. 0= LASTINFO=: _1{::cdrc=. dorgqr_jlapack2_ (_3}.}.cdrc),(lwork$0.0);lwork;,_1

(|:4{::cdrc);r
)

'q r'=: do_dgeqrf a=: 3 3 $12.0 _51 4 6 167 _68 _4 24 _41
echo a ; q ; r ; (q (+/ .*) r)
