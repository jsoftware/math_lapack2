NB. subroutine dgesvx ( character                               FACT,
NB.                     character                               TRANS,
NB.                     integer                                 N,
NB.                     integer                                 NRHS,
NB.                     double precision, dimension( lda, * )   A,
NB.                     integer                                 LDA,
NB.                     double precision, dimension( ldaf, * )  AF,
NB.                     integer                                 LDAF,
NB.                     integer, dimension( * )                 IPIV,
NB.                     character                               EQUED,
NB.                     double precision, dimension( * )        R,
NB.                     double precision, dimension( * )        C,
NB.                     double precision, dimension( ldb, * )   B,
NB.                     integer                                 LDB,
NB.                     double precision, dimension( ldx, * )   X,
NB.                     integer                                 LDX,
NB.                     double precision                        RCOND,
NB.                     double precision, dimension( * )        FERR,
NB.                     double precision, dimension( * )        BERR,
NB.                     double precision, dimension( * )        WORK,
NB.                     integer, dimension( * )                 IWORK,
NB.                     integer                                 INFO
NB.                   )
NB.
NB. DGESVX computes the solution to system of linear equations A * X = B for GE matrices
NB.
NB. Purpose:
NB.
NB.      DGESVX uses the LU factorization to compute the solution to a real
NB.      system of linear equations
NB.         A * X = B,
NB.      where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
NB.
NB.      Error bounds on the solution and a condition estimate are also
NB.      provided.
NB.
NB. Description:
NB.
NB.      The following steps are performed:
NB.
NB.      1. If FACT = 'E', real scaling factors are computed to equilibrate
NB.         the system:
NB.            TRANS = 'N':  diag(R)*A*diag(C)     *inv(diag(C))*X = diag(R)*B
NB.            TRANS = 'T': (diag(R)*A*diag(C))**T *inv(diag(R))*X = diag(C)*B
NB.            TRANS = 'C': (diag(R)*A*diag(C))**H *inv(diag(R))*X = diag(C)*B
NB.         Whether or not the system will be equilibrated depends on the
NB.         scaling of the matrix A, but if equilibration is used, A is
NB.         overwritten by diag(R)*A*diag(C) and B by diag(R)*B (if TRANS='N')
NB.         or diag(C)*B (if TRANS = 'T' or 'C').
NB.
NB.      2. If FACT = 'N' or 'E', the LU decomposition is used to factor the
NB.         matrix A (after equilibration if FACT = 'E') as
NB.            A = P * L * U,
NB.         where P is a permutation matrix, L is a unit lower triangular
NB.         matrix, and U is upper triangular.
NB.
NB.      3. If some U(i,i)=0, so that U is exactly singular, then the routine
NB.         returns with INFO = i. Otherwise, the factored form of A is used
NB.         to estimate the condition number of the matrix A.  If the
NB.         reciprocal of the condition number is less than machine precision,
NB.         INFO = N+1 is returned as a warning, but the routine still goes on
NB.         to solve for X and compute error bounds as described below.
NB.
NB.      4. The system of equations is solved for X using the factored form
NB.         of A.
NB.
NB.      5. Iterative refinement is applied to improve the computed solution
NB.         matrix and calculate error bounds and backward error estimates
NB.         for it.
NB.
NB.      6. If equilibration was used, the matrix X is premultiplied by
NB.         diag(C) (if TRANS = 'N') or diag(R) (if TRANS = 'T' or 'C') so
NB.         that it solves the original system before equilibration.
NB.
NB. Parameters
NB.
NB.                              FACT is CHARACTER*1
NB.                              Specifies whether or not the factored form of the matrix A is
NB.                              supplied on entry, and if not, whether the matrix A should be
NB.                              equilibrated before it is factored.
NB.                              = 'F':  On entry, AF and IPIV contain the factored form of A.
NB.     [in]     FACT                    If EQUED is not 'N', the matrix A has been
NB.                                      equilibrated with scaling factors given by R and C.
NB.                                      A, AF, and IPIV are not modified.
NB.                              = 'N':  The matrix A will be copied to AF and factored.
NB.                              = 'E':  The matrix A will be equilibrated if necessary, then
NB.                                      copied to AF and factored.
NB.
NB.                              TRANS is CHARACTER*1
NB.                              Specifies the form of the system of equations:
NB.     [in]     TRANS           = 'N':  A * X = B     (No transpose)
NB.                              = 'T':  A**T * X = B  (Transpose)
NB.                              = 'C':  A**H * X = B  (Transpose)
NB.
NB.                              N is INTEGER
NB.     [in]     N               The number of linear equations, i.e., the order of the
NB.                              matrix A.  N >= 0.
NB.
NB.                              NRHS is INTEGER
NB.     [in]     NRHS            The number of right hand sides, i.e., the number of columns
NB.                              of the matrices B and X.  NRHS >= 0.
NB.
NB.                              A is DOUBLE PRECISION array, dimension (LDA,N)
NB.                              On entry, the N-by-N matrix A.  If FACT = 'F' and EQUED is
NB.                              not 'N', then A must have been equilibrated by the scaling
NB.                              factors in R and/or C.  A is not modified if FACT = 'F' or
NB.                              'N', or if FACT = 'E' and EQUED = 'N' on exit.
NB.     [in,out] A
NB.                              On exit, if EQUED .ne. 'N', A is scaled as follows:
NB.                              EQUED = 'R':  A := diag(R) * A
NB.                              EQUED = 'C':  A := A * diag(C)
NB.                              EQUED = 'B':  A := diag(R) * A * diag(C).
NB.
NB.                              LDA is INTEGER
NB.     [in]     LDA             The leading dimension of the array A.  LDA >= max(1,N).
NB.
NB.                              AF is DOUBLE PRECISION array, dimension (LDAF,N)
NB.                              If FACT = 'F', then AF is an input argument and on entry
NB.                              contains the factors L and U from the factorization
NB.                              A = P*L*U as computed by DGETRF.  If EQUED .ne. 'N', then
NB.                              AF is the factored form of the equilibrated matrix A.
NB.
NB.                              If FACT = 'N', then AF is an output argument and on exit
NB.     [in,out] AF              returns the factors L and U from the factorization A = P*L*U
NB.                              of the original matrix A.
NB.
NB.                              If FACT = 'E', then AF is an output argument and on exit
NB.                              returns the factors L and U from the factorization A = P*L*U
NB.                              of the equilibrated matrix A (see the description of A for
NB.                              the form of the equilibrated matrix).
NB.
NB.                              LDAF is INTEGER
NB.     [in]     LDAF            The leading dimension of the array AF.  LDAF >= max(1,N).
NB.
NB.                              IPIV is INTEGER array, dimension (N)
NB.                              If FACT = 'F', then IPIV is an input argument and on entry
NB.                              contains the pivot indices from the factorization A = P*L*U
NB.                              as computed by DGETRF; row i of the matrix was interchanged
NB.                              with row IPIV(i).
NB.
NB.     [in,out] IPIV            If FACT = 'N', then IPIV is an output argument and on exit
NB.                              contains the pivot indices from the factorization A = P*L*U
NB.                              of the original matrix A.
NB.
NB.                              If FACT = 'E', then IPIV is an output argument and on exit
NB.                              contains the pivot indices from the factorization A = P*L*U
NB.                              of the equilibrated matrix A.
NB.
NB.                              EQUED is CHARACTER*1
NB.                              Specifies the form of equilibration that was done.
NB.                              = 'N':  No equilibration (always true if FACT = 'N').
NB.                              = 'R':  Row equilibration, i.e., A has been premultiplied by
NB.                                      diag(R).
NB.     [in,out] EQUED           = 'C':  Column equilibration, i.e., A has been postmultiplied
NB.                                      by diag(C).
NB.                              = 'B':  Both row and column equilibration, i.e., A has been
NB.                                      replaced by diag(R) * A * diag(C).
NB.                              EQUED is an input argument if FACT = 'F'; otherwise, it is an
NB.                              output argument.
NB.
NB.                              R is DOUBLE PRECISION array, dimension (N)
NB.                              The row scale factors for A.  If EQUED = 'R' or 'B', A is
NB.                              multiplied on the left by diag(R); if EQUED = 'N' or 'C', R
NB.     [in,out] R               is not accessed.  R is an input argument if FACT = 'F';
NB.                              otherwise, R is an output argument.  If FACT = 'F' and
NB.                              EQUED = 'R' or 'B', each element of R must be positive.
NB.
NB.                              C is DOUBLE PRECISION array, dimension (N)
NB.                              The column scale factors for A.  If EQUED = 'C' or 'B', A is
NB.                              multiplied on the right by diag(C); if EQUED = 'N' or 'R', C
NB.     [in,out] C               is not accessed.  C is an input argument if FACT = 'F';
NB.                              otherwise, C is an output argument.  If FACT = 'F' and
NB.                              EQUED = 'C' or 'B', each element of C must be positive.
NB.
NB.                              B is DOUBLE PRECISION array, dimension (LDB,NRHS)
NB.                              On entry, the N-by-NRHS right hand side matrix B.
NB.                              On exit,
NB.                              if EQUED = 'N', B is not modified;
NB.     [in,out] B               if TRANS = 'N' and EQUED = 'R' or 'B', B is overwritten by
NB.                              diag(R)*B;
NB.                              if TRANS = 'T' or 'C' and EQUED = 'C' or 'B', B is
NB.                              overwritten by diag(C)*B.
NB.
NB.                              LDB is INTEGER
NB.     [in]     LDB             The leading dimension of the array B.  LDB >= max(1,N).
NB.
NB.                              X is DOUBLE PRECISION array, dimension (LDX,NRHS)
NB.                              If INFO = 0 or INFO = N+1, the N-by-NRHS solution matrix X
NB.                              to the original system of equations.  Note that A and B are
NB.     [out]    X               modified on exit if EQUED .ne. 'N', and the solution to the
NB.                              equilibrated system is inv(diag(C))*X if TRANS = 'N' and
NB.                              EQUED = 'C' or 'B', or inv(diag(R))*X if TRANS = 'T' or 'C'
NB.                              and EQUED = 'R' or 'B'.
NB.
NB.                              LDX is INTEGER
NB.     [in]     LDX             The leading dimension of the array X.  LDX >= max(1,N).
NB.
NB.                              RCOND is DOUBLE PRECISION
NB.                              The estimate of the reciprocal condition number of the matrix
NB.                              A after equilibration (if done).  If RCOND is less than the
NB.     [out]    RCOND           machine precision (in particular, if RCOND = 0), the matrix
NB.                              is singular to working precision.  This condition is
NB.                              indicated by a return code of INFO > 0.
NB.
NB.                              FERR is DOUBLE PRECISION array, dimension (NRHS)
NB.                              The estimated forward error bound for each solution vector
NB.                              X(j) (the j-th column of the solution matrix X).
NB.                              If XTRUE is the true solution corresponding to X(j), FERR(j)
NB.     [out]    FERR            is an estimated upper bound for the magnitude of the largest
NB.                              element in (X(j) - XTRUE) divided by the magnitude of the
NB.                              largest element in X(j).  The estimate is as reliable as
NB.                              the estimate for RCOND, and is almost always a slight
NB.                              overestimate of the true error.
NB.
NB.                              BERR is DOUBLE PRECISION array, dimension (NRHS)
NB.                              The componentwise relative backward error of each solution
NB.     [out]    BERR            vector X(j) (i.e., the smallest relative change in
NB.                              any element of A or B that makes X(j) an exact solution).
NB.
NB.                              WORK is DOUBLE PRECISION array, dimension (4*N)
NB.                              On exit, WORK(1) contains the reciprocal pivot growth
NB.                              factor norm(A)/norm(U). The "max absolute element" norm is
NB.                              used. If WORK(1) is much less than 1, then the stability
NB.                              of the LU factorization of the (equilibrated) matrix A
NB.     [out]    WORK            could be poor. This also means that the solution X, condition
NB.                              estimator RCOND, and forward error bound FERR could be
NB.                              unreliable. If factorization fails with 0<INFO<=N, then
NB.                              WORK(1) contains the reciprocal pivot growth factor for the
NB.                              leading INFO columns of A.
NB.
NB.     [out]    IWORK           IWORK is INTEGER array, dimension (N)
NB.
NB.                              INFO is INTEGER
NB.                              = 0:  successful exit
NB.                              < 0:  if INFO = -i, the i-th argument had an illegal value
NB.                              > 0:  if INFO = i, and i is
NB.                                    <= N:  U(i,i) is exactly zero.  The factorization has
NB.                                           been completed, but the factor U is exactly
NB.                                           singular, so the solution and error bounds
NB.     [out]    INFO                         could not be computed. RCOND = 0 is returned.
NB.                                    = N+1: U is nonsingular, but RCOND is less than machine
NB.                                           precision, meaning that the matrix is singular
NB.                                           to working precision.  Nevertheless, the
NB.                                           solution and error bounds are computed because
NB.                                           there are a number of situations where the
NB.                                           computed solution can be more accurate than the
NB.                                           value of RCOND would suggest.

require 'math/lapack2'
cocurrent 'base'
coinsert 'jlapack2'
match=: matchclean;;
matchf=: matchcleanf;;

NB. =========================================================
tdgesvx=: 4 : 0
if. (3=x) *. 807>:0".}.({.~ i.&'/')9!:14'' do. 1 return. end.
zero=. (2|x){::dzero;zzero

'ma mvb'=. y
ma=. zero + ma
'm n'=. $ma
mvb=. zero + ,.^:(2>#@$)mvb
nrhs=. {:@$mvb

if. 0=2|x do.
  assert. 0= _1{::cdrc=. dgesvx`0:`sgesvx`0:@.x (,'N');(,'N');(,n);(,nrhs);(|:ma);(,1>.m);(af=. (n,1>.m)$zero);(,1>.n);(ipiv=. n$izero);(, equal=. ' ');(R=. n$dzero);(C=. n$dzero);(|:ldb{.mvb);(,ldb=. 1>.n);(X=. (nrhs,ldx)$zero);(,ldx=. 1>.n);(,dzero);(ferr=. nrhs$dzero);(berr=. nrhs$dzero);((1>.4*n)$zero);((1>.n)$izero);,_1
  'a af ipiv b X abyu'=. 5 7 9 13 15 20{cdrc
else.
  assert. 0= _1{::cdrc=. 0:`zgesvx`0:`cgesvx@.x (,'N');(,'N');(,n);(,nrhs);(|:ma);(,1>.m);(af=. (n,1>.m)$zero);(,1>.n);(ipiv=. n$izero);(, equal=. ' ');(R=. n$dzero);(C=. n$dzero);(|:ldb{.mvb);(,ldb=. 1>.n);(X=. (nrhs,ldx)$zero);(,ldx=. 1>.n);(,dzero);(ferr=. nrhs$dzero);(berr=. nrhs$dzero);((1>.2*n)$zero);((1>.2*n)$dzero);,_1
  'a af ipiv b X abyu'=. 5 7 9 13 15 21{cdrc
end.
a=. |: a
af=. m{. |:af
b=. n{. |: b
X=. n{. |: X
abyu=. {. abyu
u=. l=. 0
l=. (idmat n,m) + sltri af
u=. utri af

echo r=. mvb match`matchf@.(x>1) clean`cleanf@.(x>1) ma mp X
0{::r
)

NB. =========================================================
testdgesvx=: 3 : 0
ma0=. 0 0$0
mb0=. 0 0$0
ma1=. ? 10 10$100
mb1=. ? 10 5$50
ma2=. 0 0$zzero
mb2=. 0 0$zzero
ma3=. j./ ? 2 10 10$100
mb3=. j./ ? 2 10 5$50
ma4=. 0 0$0
vb4=. 0$0
ma5=. ? 10 10$100
vb5=. ? 10$50
ma6=. 0 0$zzero
vb6=. 0$zzero
ma7=. j./ ? 2 10 10$100
vb7=. j./ ? 2 10$50
assert. 0&tdgesvx &> (< ma0;mb0) , (< ma1;mb1) , (< ma4;vb4) , (< ma5;vb5)
assert. 1&tdgesvx &> (< ma2;mb2) , (< ma3;mb3) , (< ma6;vb6) , (< ma7;vb7)
assert. 2&tdgesvx &> (< ma0;mb0) , (< ma1;mb1) , (< ma4;vb4) , (< ma5;vb5)
assert. 3&tdgesvx &> (< ma2;mb2) , (< ma3;mb3) , (< ma6;vb6) , (< ma7;vb7)
EMPTY
)
