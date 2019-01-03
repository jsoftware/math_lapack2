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

require 'math/lapack2'
cocurrent 'base'
coinsert 'jlapack2'
match=: matchclean;;
matchf=: matchcleanf;;

NB. =========================================================
NB. DGETRF ZGETRF computes an LU factorization of a general M-by-N matrix A using partial pivoting with row interchanges.

tdgetrf=: 4 : 0
if. (3=x) *. 807>:0".}.({.~ i.&'/')9!:14'' do. 1 return. end.
zero=. (2|x){::dzero;zzero
a=. zero + y
'm n'=. $a
mn=. m<.n
assert. 0= _1{::cdrc=. dgetrf`zgetrf`sgetrf`cgetrf@.x (,m);(,n);(|:a);(,1>.m);(mn$izero);,_1

'r ipiv'=. 3 5{cdrc
r=. |:r
u=. l=. izero
l=. (idmat m,n) + sltri r
if. m < n do.
  l=. (_,mn) {. l
end.
u=. utri r
if. m > n do.
  u=. (mn,_) {. u
end.

echo l;u;ipiv
echo r=. a match`matchf@.(x>1) ipiv invperm~ l mp u
0{::r
)

NB. =========================================================
testdgetrf=: 3 : 0
m0=. 0 0$0
m1=. ?.4 6$10
m2=. ?.6 4$10
m3=. ?.6 6$10
m4=. 0 0$zzero
m5=. j./ ?. 2 4 6$10
m6=. j./ ?.2 6 4$10
m7=. j./ ?.2 6 6$10
assert. 0&tdgetrf &> m0;m1;m2;m3
assert. 1&tdgetrf &> m4;m5;m6;m7
assert. 2&tdgetrf &> m0;m1;m2;m3
assert. 3&tdgetrf &> m4;m5;m6;m7
EMPTY
)

