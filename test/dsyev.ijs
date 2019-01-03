NB. subroutine dsyev ( character                              JOBZ,
NB.                    character                              UPLO,
NB.                    integer                                N,
NB.                    double precision, dimension( lda, * )  A,
NB.                    integer                                LDA,
NB.                    double precision, dimension( * )       W,
NB.                    double precision, dimension( * )       WORK,
NB.                    integer                                LWORK,
NB.                    integer                                INFO
NB.                  )
NB.
NB. DSYEV computes the eigenvalues and, optionally, the left and/or right eigenvectors for SY
NB. matrices
NB.
NB. Purpose:
NB.
NB.      DSYEV computes all eigenvalues and, optionally, eigenvectors of a
NB.      real symmetric matrix A.
NB.
NB. Parameters
NB.
NB.                              JOBZ is CHARACTER*1
NB.     [in]     JOBZ            = 'N':  Compute eigenvalues only;
NB.                              = 'V':  Compute eigenvalues and eigenvectors.
NB.
NB.                              UPLO is CHARACTER*1
NB.     [in]     UPLO            = 'U':  Upper triangle of A is stored;
NB.                              = 'L':  Lower triangle of A is stored.
NB.
NB.                              N is INTEGER
NB.     [in]     N               The order of the matrix A.  N >= 0.
NB.
NB.                              A is DOUBLE PRECISION array, dimension (LDA, N)
NB.                              On entry, the symmetric matrix A.  If UPLO = 'U', the
NB.                              leading N-by-N upper triangular part of A contains the
NB.                              upper triangular part of the matrix A.  If UPLO = 'L',
NB.                              the leading N-by-N lower triangular part of A contains
NB.     [in,out] A               the lower triangular part of the matrix A.
NB.                              On exit, if JOBZ = 'V', then if INFO = 0, A contains the
NB.                              orthonormal eigenvectors of the matrix A.
NB.                              If JOBZ = 'N', then on exit the lower triangle (if UPLO='L')
NB.                              or the upper triangle (if UPLO='U') of A, including the
NB.                              diagonal, is destroyed.
NB.
NB.                              LDA is INTEGER
NB.     [in]     LDA             The leading dimension of the array A.  LDA >= max(1,N).
NB.
NB.                              W is DOUBLE PRECISION array, dimension (N)
NB.     [out]    W               If INFO = 0, the eigenvalues in ascending order.
NB.
NB.                              WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
NB.     [out]    WORK            On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
NB.
NB.                              LWORK is INTEGER
NB.                              The length of the array WORK.  LWORK >= max(1,3*N-1).
NB.                              For optimal efficiency, LWORK >= (NB+2)*N,
NB.                              where NB is the blocksize for DSYTRD returned by ILAENV.
NB.     [in]     LWORK
NB.                              If LWORK = -1, then a workspace query is assumed; the routine
NB.                              only calculates the optimal size of the WORK array, returns
NB.                              this value as the first entry of the WORK array, and no error
NB.                              message related to LWORK is issued by XERBLA.
NB.
NB.                              INFO is INTEGER
NB.                              = 0:  successful exit
NB.     [out]    INFO            < 0:  if INFO = -i, the i-th argument had an illegal value
NB.                              > 0:  if INFO = i, the algorithm failed to converge; i
NB.                                    off-diagonal elements of an intermediate tridiagonal
NB.                                    form did not converge to zero.

require 'math/lapack2'
cocurrent 'base'
coinsert 'jlapack2'
match=: matchclean;;
matchf=: matchcleanf;;

NB. =========================================================
tdsyev=: 4 : 0
if. (3=x) *. 807>:0".}.({.~ i.&'/')9!:14'' do. 1 return. end.
a=. dzero + y
'm n'=. $a
assert. 0= _1{::cdrc=. dsyev`0:`ssyev`0:@.x (,'V');(,'U');(,n);(|:a);(,1>.m);(V=. n$dzero);(lwork$dzero);(,lwork=. 1>._1+3*n);,_1
'R V'=. 4 6{cdrc
R=. |:R
echo V;R
echo r=. (clean`cleanf@.(x>1) a mp R) match`matchf@.(x>1) (clean`cleanf@.(x>1) V *"1 R)
0{::r
)

NB. =========================================================
testdsyev=: 3 : 0
m0=. 0 0$dzero
m1=. (+ |:) ?.6 6$10
m2=. 0 0$dzero
assert. 0&tdsyev &> m0;m1;m2
assert. 2&tdsyev &> m0;m1;m2
EMPTY
)

