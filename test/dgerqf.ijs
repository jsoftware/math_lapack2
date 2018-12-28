NB. subroutine dgerqf ( integer                                M,
NB.                     integer                                N,
NB.                     double precision, dimension( lda, * )  A,
NB.                     integer                                LDA,
NB.                     double precision, dimension( * )       TAU,
NB.                     double precision, dimension( * )       WORK,
NB.                     integer                                LWORK,
NB.                     integer                                INFO
NB.                   )
NB.
NB. DGERQF
NB.
NB. Purpose:
NB.
NB.      DGERQF computes an RQ factorization of a real M-by-N matrix A:
NB.      A = R * Q.
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
NB.                              if m <= n, the upper triangle of the subarray
NB.                              A(1:m,n-m+1:n) contains the M-by-M upper triangular matrix R;
NB.     [in,out] A               if m >= n, the elements on and above the (m-n)-th subdiagonal
NB.                              contain the M-by-N upper trapezoidal matrix R;
NB.                              the remaining elements, with the array TAU, represent the
NB.                              orthogonal matrix Q as a product of min(m,n) elementary
NB.                              reflectors (see Further Details).
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
NB.                              For optimum performance LWORK >= M*NB, where NB is
NB.                              the optimal blocksize.
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
tdgerqf=: 4 : 0
if. (3=x) *. 807>:0".}.({.~ i.&'/')9!:14'' do. 1 return. end.
zero=. (2|x){::dzero;zzero
a=. zero + y
'm n'=. $a
mn=. m<.n
d=. n-m
assert. 0= _1{::cdrc=. dgerqf`zgerqf`sgerqf`cgerqf@.x (,m);(,n);(|:a);(,1>.m);(tau=. mn$zero);(lwork$zero);(,lwork=. 1 >. 10 * m >. n);,_1

r=. h=. q=. 0

val=. |: 3{::cdrc
tau=. 5{::cdrc
r=. (m,(-mn)) {. d utri val
h=. (-mn) {. (d idmat m,n) + (d sltri val)
q=. (-mn) {. mp/ (idmat n) -"2 (+ tau) * (* +)"0/~"1 + h

echo r;q
echo r=. a match`matchf@.(x>1) r mp q
0{::r
)

NB. =========================================================
testdgerqf=: 3 : 0
m0=. 0 0$0
m1=. ?.4 6$10
m2=. ?.6 4$10
m3=. ?.6 6$10
m4=. 0 0$zzero
m5=. j./ ?. 2 4 6$10
m6=. j./ ?.2 6 4$10
m7=. j./ ?.2 6 6$10
assert. 0&tdgerqf &> m0;m1;m2;m3
assert. 1&tdgerqf &> m4;m5;m6;m7
assert. 2&tdgerqf &> m0;m1;m2;m3
assert. 3&tdgerqf &> m4;m5;m6;m7
EMPTY
)
