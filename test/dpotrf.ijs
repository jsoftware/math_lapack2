NB. subroutine dpotrf ( character                              UPLO,
NB.                     integer                                N,
NB.                     double precision, dimension( lda, * )  A,
NB.                     integer                                LDA,
NB.                     integer                                INFO
NB.                   )
NB.
NB. DPOTRF
NB.
NB. DPOTRF VARIANT: top-looking block version of the algorithm, calling Level 3 BLAS.
NB.
NB. Purpose:
NB.
NB.      DPOTRF computes the Cholesky factorization of a real symmetric
NB.      positive definite matrix A.
NB.
NB.      The factorization has the form
NB.         A = U**T * U,  if UPLO = 'U', or
NB.         A = L  * L**T,  if UPLO = 'L',
NB.      where U is an upper triangular matrix and L is lower triangular.
NB.
NB.      This is the block version of the algorithm, calling Level 3 BLAS.
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
NB.                             A is DOUBLE PRECISION array, dimension (LDA,N)
NB.                             On entry, the symmetric matrix A.  If UPLO = 'U', the leading
NB.                             N-by-N upper triangular part of A contains the upper
NB.                             triangular part of the matrix A, and the strictly lower
NB.                             triangular part of A is not referenced.  If UPLO = 'L', the
NB.     [in,out] A              leading N-by-N lower triangular part of A contains the lower
NB.                             triangular part of the matrix A, and the strictly upper
NB.                             triangular part of A is not referenced.
NB.
NB.                             On exit, if INFO = 0, the factor U or L from the Cholesky
NB.                             factorization A = U**T*U or A = L*L**T.
NB.
NB.                             LDA is INTEGER
NB.     [in]     LDA            The leading dimension of the array A.  LDA >= max(1,N).
NB.
NB.                             INFO is INTEGER
NB.                             = 0:  successful exit
NB.     [out]    INFO           < 0:  if INFO = -i, the i-th argument had an illegal value
NB.                             > 0:  if INFO = i, the leading minor of order i is not
NB.                                   positive definite, and the factorization could not be
NB.                                   completed.

require 'math/lapack2'
cocurrent 'base'
coinsert 'jlapack2'
match=: matchclean;;
matchf=: matchcleanf;;

NB. =========================================================
tdpotrf=: 4 : 0
if. (3=x) *. 807>:0".}.({.~ i.&'/')9!:14'' do. 1 return. end.
zero=. (2|x){::dzero;zzero
a=. zero + y
'm n'=. $a
assert. 0= _1{::cdrc=. dpotrf`zpotrf`spotrf`cpotrf@.x (,'L');(,n);(|:a);(,1>.m);,_1
echo L=. ltri |: 3{::cdrc
echo r=. a match`matchf@.(x>1) L mp +|:L
0{::r
)

NB. =========================================================
testdpotrf=: 3 : 0
m0=. 0 0$0
m1=. (mp |:) ?.4 4$10
m2=. (mp |:) _25 + ?.10 10$100
m3=. 0 0$zzero
m4=. (mp (+ @ |:)) j./ ?.2 4 4$10
m5=. (mp (+ @ |:)) _25 + j./ ?.2 10 10$100
assert. 0&tdpotrf &> m0;m1;m2;m3
assert. 1&tdpotrf &> m4;m5
assert. 2&tdpotrf &> m0;m1;m2;m3
assert. 3&tdpotrf &> m4;m5
EMPTY
)

