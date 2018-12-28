NB. subroutine dgelqf ( integer                                M,
NB.                     integer                                N,
NB.                     double precision, dimension( lda, * )  A,
NB.                     integer                                LDA,
NB.                     double precision, dimension( * )       TAU,
NB.                     double precision, dimension( * )       WORK,
NB.                     integer                                LWORK,
NB.                     integer                                INFO
NB.                   )
NB.
NB. DGELQF
NB.
NB. Purpose:
NB.
NB.      DGELQF computes an LQ factorization of a real M-by-N matrix A:
NB.      A = L * Q.
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
NB.                              On exit, the elements on and below the diagonal of the array
NB.     [in,out] A               contain the m-by-min(m,n) lower trapezoidal matrix L (L is
NB.                              lower triangular if m <= n); the elements above the diagonal,
NB.                              with the array TAU, represent the orthogonal matrix Q as a
NB.                              product of elementary reflectors (see Further Details).
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
NB.                              The dimension of the array WORK.  LWORK >= max(1,M).
NB.                              For optimum performance LWORK >= M*NB, where NB is the
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
tdgelqf=: 4 : 0
if. (3=x) *. 807>:0".}.({.~ i.&'/')9!:14'' do. 1 return. end.
zero=. (2|x){::dzero;zzero

a0=. a:{ a=. zero + y
'm n'=. $a
mn=. m <. n

assert. 0= _1{::cdrc=. dgelqf`zgelqf`sgelqf`cgelqf@.x (,m);(,n);(|:a);(,1>.m);(tau=. mn$zero);(lwork$zero);(,lwork=. 1 >. 10 * m >. n);,_1
'val tau'=. 3 5{cdrc
val=. |: val

Q=. H=. L=. 0

L=. (m,mn) {. ltri val
H=. (mn,n) {. (idmat m,n) + sutri val
Q=. mn {. mp/ (idmat n) -"2 |. (+ tau) * (* +)"0/~"1 + H

echo L;Q
echo r=. a match`matchf@.(x>1) L mp Q
0{::r
)

NB. =========================================================
testdgelqf=: 3 : 0
m0=. 0 0$0
m1=. ?.4 6$10
m2=. ?.6 4$10
m3=. ?.6 6$10
m4=. 0 0$zzero
m5=. j./ ?. 2 4 6$10
m6=. j./ ?.2 6 4$10
m7=. j./ ?.2 6 6$10
assert. 0&tdgelqf &> m0;m1;m2;m3
assert. 1&tdgelqf &> m4;m5;m6;m7
assert. 2&tdgelqf &> m0;m1;m2;m3
assert. 3&tdgelqf &> m4;m5;m6;m7
EMPTY
)
