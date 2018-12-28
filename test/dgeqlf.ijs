NB. subroutine dgeqlf ( integer                                M,
NB.                     integer                                N,
NB.                     double precision, dimension( lda, * )  A,
NB.                     integer                                LDA,
NB.                     double precision, dimension( * )       TAU,
NB.                     double precision, dimension( * )       WORK,
NB.                     integer                                LWORK,
NB.                     integer                                INFO
NB.                   )
NB.
NB. DGEQLF
NB.
NB. Purpose:
NB.
NB.      DGEQLF computes a QL factorization of a real M-by-N matrix A:
NB.      A = Q * L.
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
NB.                              On exit,
NB.                              if m >= n, the lower triangle of the subarray
NB.                              A(m-n+1:m,1:n) contains the N-by-N lower triangular matrix L;
NB.     [in,out] A               if m <= n, the elements on and below the (n-m)-th
NB.                              superdiagonal contain the M-by-N lower trapezoidal matrix L;
NB.                              the remaining elements, with the array TAU, represent the
NB.                              orthogonal matrix Q as a product of elementary reflectors
NB.                              (see Further Details).
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
NB.                              The dimension of the array WORK.  LWORK >= max(1,N).
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
NB.                              < 0:  if INFO = -i, the i-th argument had an illegal value

require 'math/lapack2'
cocurrent 'base'
coinsert 'jlapack2'
match=: matchclean;;
matchf=: matchcleanf;;

NB. =========================================================
tdgeqlf=: 4 : 0
if. (3=x) *. 807>:0".}.({.~ i.&'/')9!:14'' do. 1 return. end.
zero=. (2|x){::dzero;zzero

a=. zero + y
'm n'=. $a
mn=. m <. n
d=. n-m
lda=. 1 >. m
tau=. mn $ zero
lwork=. 1 >. 10 * m >. n

assert. 0= _1{::cdrc=. dgeqlf`zgeqlf`sgeqlf`cgeqlf@.x (,m);(,n);(|:a);(,1>.m);(tau=. mn$zero);(lwork$zero);(,lwork);,_1
'val tau'=. 3 5{cdrc
val=. |: val

Q=. H=. L=. 0

H=. (m,(-mn)) {. (d idmat m,n) + (d sutri val)
L=. ((-mn),n) {. d ltri val
Q=. (m,(-mn)) {. mp/ (idmat m) -"2 |. tau * (* +)"0/~"1 |: H

echo Q;L
echo r=. a match`matchf@.(x>1) Q mp L
0{::r
)

NB. =========================================================
testdgeqlf=: 3 : 0
m0=. 0 0$0
m1=. ?.4 6$10
m2=. ?.6 4$10
m3=. ?.6 6$10
m4=. 0 0$zzero
m5=. j./ ?. 2 4 6$10
m6=. j./ ?.2 6 4$10
m7=. j./ ?.2 6 6$10
assert. 0&tdgeqlf &> m0;m1;m2;m3
assert. 1&tdgeqlf &> m4;m5;m6;m7
assert. 2&tdgeqlf &> m0;m1;m2;m3
assert. 3&tdgeqlf &> m4;m5;m6;m7
EMPTY
)
