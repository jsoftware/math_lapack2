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

require 'math/lapack2'
cocurrent 'base'
coinsert 'jlapack2'
match=: matchclean;;
matchf=: matchcleanf;;

NB. =========================================================
tdgeqrf=: 4 : 0
assert. ismatrix y
zero=. (2|x){::0.0;0j0
'm n'=. $y
mn=. m <. n
lwork=. , 1 >. 10 * m >. n
assert. 0= _1{::cdrc=. dgeqrf`zgeqrf`sgeqrf`cgeqrf@.x (,m);(,n);(|:y);(,1>.m);(mn$zero);(lwork$zero);lwork;,_1
'val tau'=. 3 5{cdrc
val=. |: val
h=. mn {."1 (idmat m,n) + sltri val
q=. mn {."1 mp/ (idmat m) -"2 tau * (* +)"0/~"1 |: h
r=. mn {. utri val

echo q;r
echo r=. y match`matchf@.(x>1) q mp r
0{::r
)

NB. =========================================================
testdgeqrf=: 3 : 0
m0=. 0 0$0.0
m1=. ?.4 6$10
m2=. ?.6 4$10
m3=. ?.6 6$10
m4=. 0 0$0j0
m5=. j./ ?. 2 4 6$10
m6=. j./ ?.2 6 4$10
m7=. j./ ?.2 6 6$10
assert. 0&tdgeqrf &> m1;m2;m3
assert. 1&tdgeqrf &> m5;m6;m7
assert. 2&tdgeqrf &> m1;m2;m3
assert. 3&tdgeqrf &> m5;m6;m7
EMPTY
)
