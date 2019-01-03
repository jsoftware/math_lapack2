NB. subroutine dgesv ( integer                                N,
NB.                    integer                                NRHS,
NB.                    double precision, dimension( lda, * )  A,
NB.                    integer                                LDA,
NB.                    integer, dimension( * )                IPIV,
NB.                    double precision, dimension( ldb, * )  B,
NB.                    integer                                LDB,
NB.                    integer                                INFO
NB.                  )
NB.
NB. DGESV computes the solution to system of linear equations A * X = B for GE matrices
NB.
NB. Purpose:
NB.
NB.      DGESV computes the solution to a real system of linear equations
NB.         A * X = B,
NB.      where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
NB.
NB.      The LU decomposition with partial pivoting and row interchanges is
NB.      used to factor A as
NB.         A = P * L * U,
NB.      where P is a permutation matrix, L is unit lower triangular, and U is
NB.      upper triangular.  The factored form of A is then used to solve the
NB.      system of equations A * X = B.
NB.
NB. Parameters
NB.
NB.                             N is INTEGER
NB.     [in]     N              The number of linear equations, i.e., the order of the
NB.                             matrix A.  N >= 0.
NB.
NB.                             NRHS is INTEGER
NB.     [in]     NRHS           The number of right hand sides, i.e., the number of columns
NB.                             of the matrix B.  NRHS >= 0.
NB.
NB.                             A is DOUBLE PRECISION array, dimension (LDA,N)
NB.                             On entry, the N-by-N coefficient matrix A.
NB.     [in,out] A              On exit, the factors L and U from the factorization
NB.                             A = P*L*U; the unit diagonal elements of L are not stored.
NB.
NB.                             LDA is INTEGER
NB.     [in]     LDA            The leading dimension of the array A.  LDA >= max(1,N).
NB.
NB.                             IPIV is INTEGER array, dimension (N)
NB.     [out]    IPIV           The pivot indices that define the permutation matrix P;
NB.                             row i of the matrix was interchanged with row IPIV(i).
NB.
NB.                             B is DOUBLE PRECISION array, dimension (LDB,NRHS)
NB.     [in,out] B              On entry, the N-by-NRHS matrix of right hand side matrix B.
NB.                             On exit, if INFO = 0, the N-by-NRHS solution matrix X.
NB.
NB.                             LDB is INTEGER
NB.     [in]     LDB            The leading dimension of the array B.  LDB >= max(1,N).
NB.
NB.                             INFO is INTEGER
NB.                             = 0:  successful exit
NB.     [out]    INFO           < 0:  if INFO = -i, the i-th argument had an illegal value
NB.                             > 0:  if INFO = i, U(i,i) is exactly zero.  The factorization
NB.                                   has been completed, but the factor U is exactly
NB.                                   singular, so the solution could not be computed.

require 'math/lapack2'
cocurrent 'base'
coinsert 'jlapack2'
match=: matchclean;;
matchf=: matchcleanf;;

NB. =========================================================
tdgesv=: 4 : 0
if. (3=x) *. 807>:0".}.({.~ i.&'/')9!:14'' do. 1 return. end.
zero=. (2|x){::dzero;zzero
'ma mvb'=. y
ma=. zero + ma
'm n'=. $ma
mvb=. zero + ,.^:(2>#@$)mvb
nrhs=. {:@$mvb
assert. 0= _1{::cdrc=. dgesv`zgesv`sgesv`cgesv@.x (,n);(,nrhs);(|:ma);(,1>.m);(ipv=. n$2-2);(|:ldb{.mvb);(,ldb=. 1>.n);,_1
R=. n{. |: 6{::cdrc
echo R
echo r=. mvb match`matchf@.(x>1) clean`cleanf@.(x>1) ma mp R
0{::r
)

NB. =========================================================
testdgesv=: 3 : 0
ma0=. 0 0$0
mb0=. 0 0$0
ma1=. ?. 10 10$100
mb1=. ?. 10 5$50
ma2=. 0 0$zzero
mb2=. 0 0$zzero
ma3=. j./ ?. 2 10 10$100
mb3=. j./ ?. 2 10 5$50
ma4=. 0 0$0
vb4=. 0$0
ma5=. ?. 10 10$100
vb5=. ?. 10$50
ma6=. 0 0$zzero
vb6=. 0$zzero
ma7=. j./ ?. 2 10 10$100
vb7=. j./ ?. 2 10$50
0&tdgesv &> (< ma0;mb0) , (< ma1;mb1) , (< ma2;mb2)
1&tdgesv &> (< ma3;mb3) , (< ma4;vb4) , (< ma5;vb5) , (< ma6;vb6) , (< ma7;vb7)
2&tdgesv &> (< ma0;mb0) , (< ma1;mb1) , (< ma2;mb2)
3&tdgesv &> (< ma3;mb3) , (< ma4;vb4) , (< ma5;vb5) , (< ma6;vb6) , (< ma7;vb7)
EMPTY
)
