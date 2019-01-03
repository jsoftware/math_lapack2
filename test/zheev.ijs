NB. subroutine zheev ( character                         JOBZ,
NB.                    character                         UPLO,
NB.                    integer                           N,
NB.                    complex*16, dimension( lda, * )   A,
NB.                    integer                           LDA,
NB.                    double precision, dimension( * )  W,
NB.                    complex*16, dimension( * )        WORK,
NB.                    integer                           LWORK,
NB.                    double precision, dimension( * )  RWORK,
NB.                    integer                           INFO
NB.                  )
NB.
NB. ZHEEV computes the eigenvalues and, optionally, the left and/or right eigenvectors for HE
NB. matrices
NB.
NB. Purpose:
NB.
NB.      ZHEEV computes all eigenvalues and, optionally, eigenvectors of a
NB.      complex Hermitian matrix A.
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
NB.                              A is COMPLEX*16 array, dimension (LDA, N)
NB.                              On entry, the Hermitian matrix A.  If UPLO = 'U', the
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
NB.                              WORK is COMPLEX*16 array, dimension (MAX(1,LWORK))
NB.     [out]    WORK            On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
NB.
NB.                              LWORK is INTEGER
NB.                              The length of the array WORK.  LWORK >= max(1,2*N-1).
NB.                              For optimal efficiency, LWORK >= (NB+1)*N,
NB.                              where NB is the blocksize for ZHETRD returned by ILAENV.
NB.     [in]     LWORK
NB.                              If LWORK = -1, then a workspace query is assumed; the routine
NB.                              only calculates the optimal size of the WORK array, returns
NB.                              this value as the first entry of the WORK array, and no error
NB.                              message related to LWORK is issued by XERBLA.
NB.
NB.     [out]    RWORK           RWORK is DOUBLE PRECISION array, dimension (max(1, 3*N-2))
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
tzheev=: 4 : 0
if. (3=x) *. 807>:0".}.({.~ i.&'/')9!:14'' do. 1 return. end.
a=. zzero + y
'm n'=. $a
assert. 0= _1{::cdrc=. 0:`zheev`0:`cheev@.x (,'V');(,'U');(,n);(|:a);(,1>.m);(V=. n$dzero);(lwork$zzero);(,lwork=. 1>._1+2*n);(dzero$~1>._2+3*n);,_1
'R V'=. 4 6{cdrc
R=. |:R
echo V;R
echo r=. (clean`cleanf@.(x>1) a mp R) match`matchf@.(x>1) (clean`cleanf@.(x>1) V *"1 R)
0{::r
)

NB. =========================================================
testzheev=: 3 : 0
m3=. (+ (+@|:)) j./ ?.2 6 6$10
assert 1&tzheev m3
assert 3&tzheev m3
EMPTY
)

