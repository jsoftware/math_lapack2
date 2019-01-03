NB. subroutine dtrtrs ( character                              UPLO,
NB.                     character                              TRANS,
NB.                     character                              DIAG,
NB.                     integer                                N,
NB.                     integer                                NRHS,
NB.                     double precision, dimension( lda, * )  A,
NB.                     integer                                LDA,
NB.                     double precision, dimension( ldb, * )  B,
NB.                     integer                                LDB,
NB.                     integer                                INFO
NB.                   )
NB.
NB. DTRTRS
NB.
NB. Purpose:
NB.
NB.      DTRTRS solves a triangular system of the form
NB.
NB.         A * X = B  or  A**T * X = B,
NB.
NB.      where A is a triangular matrix of order N, and B is an N-by-NRHS
NB.      matrix.  A check is made to verify that A is nonsingular.
NB.
NB. Parameters
NB.
NB.                              UPLO is CHARACTER*1
NB.     [in]     UPLO            = 'U':  A is upper triangular;
NB.                              = 'L':  A is lower triangular.
NB.
NB.                              TRANS is CHARACTER*1
NB.                              Specifies the form of the system of equations:
NB.     [in]     TRANS           = 'N':  A * X = B  (No transpose)
NB.                              = 'T':  A**T * X = B  (Transpose)
NB.                              = 'C':  A**H * X = B  (Conjugate transpose = Transpose)
NB.
NB.                              DIAG is CHARACTER*1
NB.     [in]     DIAG            = 'N':  A is non-unit triangular;
NB.                              = 'U':  A is unit triangular.
NB.
NB.                              N is INTEGER
NB.     [in]     N               The order of the matrix A.  N >= 0.
NB.
NB.                              NRHS is INTEGER
NB.     [in]     NRHS            The number of right hand sides, i.e., the number of columns
NB.                              of the matrix B.  NRHS >= 0.
NB.
NB.                              A is DOUBLE PRECISION array, dimension (LDA,N)
NB.                              The triangular matrix A.  If UPLO = 'U', the leading N-by-N
NB.                              upper triangular part of the array A contains the upper
NB.                              triangular matrix, and the strictly lower triangular part of
NB.     [in]     A               A is not referenced.  If UPLO = 'L', the leading N-by-N lower
NB.                              triangular part of the array A contains the lower triangular
NB.                              matrix, and the strictly upper triangular part of A is not
NB.                              referenced.  If DIAG = 'U', the diagonal elements of A are
NB.                              also not referenced and are assumed to be 1.
NB.
NB.                              LDA is INTEGER
NB.     [in]     LDA             The leading dimension of the array A.  LDA >= max(1,N).
NB.
NB.                              B is DOUBLE PRECISION array, dimension (LDB,NRHS)
NB.     [in,out] B               On entry, the right hand side matrix B.
NB.                              On exit, if INFO = 0, the solution matrix X.
NB.
NB.                              LDB is INTEGER
NB.     [in]     LDB             The leading dimension of the array B.  LDB >= max(1,N).
NB.
NB.                              INFO is INTEGER
NB.                              = 0:  successful exit
NB.     [out]    INFO            < 0: if INFO = -i, the i-th argument had an illegal value
NB.                              > 0: if INFO = i, the i-th diagonal element of A is zero,
NB.                                   indicating that the matrix is singular and the solutions
NB.                                   X have not been computed.

require 'math/lapack2'
cocurrent 'base'
coinsert 'jlapack2'
match=: matchclean;;
matchf=: matchcleanf;;

NB. =========================================================
tdtrtrs=: 4 : 0
if. (3=x) *. 807>:0".}.({.~ i.&'/')9!:14'' do. 1 return. end.
zero=. (2|x){::dzero;zzero
'ma mvb'=. y
ma=. zero + ma
'm n'=. $ma
mvb=. zero + ,.^:(2>#@$)mvb
nrhs=. {:@$mvb
assert. 0= _1{::cdrc=. dtrtrs`ztrtrs`strtrs`ctrtrs@.x (,'U');(,'N');(,'N');(,n);(,nrhs);(|:ma);(,1>.m);(|:ldb{.mvb);(,ldb=. 1>.n);,_1
R=. n{. |: 8{::cdrc
echo R
echo r=. mvb match`matchf@.(x>1) clean`cleanf@.(x>1) ma mp R
0{::r
)

NB. =========================================================
testdtrtrs=: 3 : 0
ma0=. 0 0$0
mb0=. 0 0$0
ma1=. utri 1+ ?. 10 10$100
vb1=. ?. 10 5$50
ma2=. 0 0$zzero
mb2=. 0 0$zzero
ma3=. utri j./ ?. 2 10 10$100
vb3=. j./ ?. 2 10 5$50
ma4=. 0 0$0
vb4=. 0$0
ma5=. utri 1+ ?. 1 >. 10 10$100
vb5=. ?. 10$50
ma6=. 0 0$zzero
vb6=. 0$zzero
ma7=. utri j./ ?. 2 10 10$100
vb7=. j./ ?. 2 10$50
assert. 0&tdtrtrs &> (< ma0;mb0) , (< ma1;vb1) , (< ma2;mb2)
assert. 1&tdtrtrs &> (< ma3;vb3) , (< ma4;vb4) , (< ma5;vb5) , (< ma6;vb6) , (< ma7;vb7)
assert. 2&tdtrtrs &> (< ma0;mb0) , (< ma1;vb1) , (< ma2;mb2)
assert. 3&tdtrtrs &> (< ma3;vb3) , (< ma4;vb4) , (< ma5;vb5) , (< ma6;vb6) , (< ma7;vb7)
EMPTY
)
