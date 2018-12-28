NB. subroutine dgelss ( integer                                M,
NB.                     integer                                N,
NB.                     integer                                NRHS,
NB.                     double precision, dimension( lda, * )  A,
NB.                     integer                                LDA,
NB.                     double precision, dimension( ldb, * )  B,
NB.                     integer                                LDB,
NB.                     double precision, dimension( * )       S,
NB.                     double precision                       RCOND,
NB.                     integer                                RANK,
NB.                     double precision, dimension( * )       WORK,
NB.                     integer                                LWORK,
NB.                     integer                                INFO
NB.                   )
NB.
NB. DGELSS solves overdetermined or underdetermined systems for GE matrices
NB.
NB. Purpose:
NB.
NB.      DGELSS computes the minimum norm solution to a real linear least
NB.      squares problem:
NB.
NB.      Minimize 2-norm(| b - A*x |).
NB.
NB.      using the singular value decomposition (SVD) of A. A is an M-by-N
NB.      matrix which may be rank-deficient.
NB.
NB.      Several right hand side vectors b and solution vectors x can be
NB.      handled in a single call; they are stored as the columns of the
NB.      M-by-NRHS right hand side matrix B and the N-by-NRHS solution matrix
NB.      X.
NB.
NB.      The effective rank of A is determined by treating as zero those
NB.      singular values which are less than RCOND times the largest singular
NB.      value.
NB.
NB. Parameters
NB.
NB.                              M is INTEGER
NB.     [in]     M               The number of rows of the matrix A. M >= 0.
NB.
NB.                              N is INTEGER
NB.     [in]     N               The number of columns of the matrix A. N >= 0.
NB.
NB.                              NRHS is INTEGER
NB.     [in]     NRHS            The number of right hand sides, i.e., the number of columns
NB.                              of the matrices B and X. NRHS >= 0.
NB.
NB.                              A is DOUBLE PRECISION array, dimension (LDA,N)
NB.                              On entry, the M-by-N matrix A.
NB.     [in,out] A               On exit, the first min(m,n) rows of A are overwritten with
NB.                              its right singular vectors, stored rowwise.
NB.
NB.                              LDA is INTEGER
NB.     [in]     LDA             The leading dimension of the array A.  LDA >= max(1,M).
NB.
NB.                              B is DOUBLE PRECISION array, dimension (LDB,NRHS)
NB.                              On entry, the M-by-NRHS right hand side matrix B.
NB.                              On exit, B is overwritten by the N-by-NRHS solution
NB.     [in,out] B               matrix X.  If m >= n and RANK = n, the residual
NB.                              sum-of-squares for the solution in the i-th column is given
NB.                              by the sum of squares of elements n+1:m in that column.
NB.
NB.                              LDB is INTEGER
NB.     [in]     LDB             The leading dimension of the array B. LDB >= max(1,max(M,N)).
NB.
NB.                              S is DOUBLE PRECISION array, dimension (min(M,N))
NB.     [out]    S               The singular values of A in decreasing order.
NB.                              The condition number of A in the 2-norm = S(1)/S(min(m,n)).
NB.
NB.                              RCOND is DOUBLE PRECISION
NB.                              RCOND is used to determine the effective rank of A.
NB.     [in]     RCOND           Singular values S(i) <= RCOND*S(1) are treated as zero.
NB.                              If RCOND < 0, machine precision is used instead.
NB.
NB.                              RANK is INTEGER
NB.     [out]    RANK            The effective rank of A, i.e., the number of singular values
NB.                              which are greater than RCOND*S(1).
NB.
NB.                              WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
NB.     [out]    WORK            On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
NB.
NB.                              LWORK is INTEGER
NB.                              The dimension of the array WORK. LWORK >= 1, and also:
NB.                              LWORK >= 3*min(M,N) + max( 2*min(M,N), max(M,N), NRHS )
NB.                              For good performance, LWORK should generally be larger.
NB.     [in]     LWORK
NB.                              If LWORK = -1, then a workspace query is assumed; the routine
NB.                              only calculates the optimal size of the WORK array, returns
NB.                              this value as the first entry of the WORK array, and no error
NB.                              message related to LWORK is issued by XERBLA.
NB.
NB.                              INFO is INTEGER
NB.                              = 0:  successful exit
NB.     [out]    INFO            < 0:  if INFO = -i, the i-th argument had an illegal value.
NB.                              > 0:  the algorithm for computing the SVD failed to converge;
NB.                                    if INFO = i, i off-diagonal elements of an intermediate
NB.                                    bidiagonal form did not converge to zero.

require 'math/lapack2'
cocurrent 'base'
coinsert 'jlapack2'
match=: matchclean;;
matchf=: matchcleanf;;

NB. =========================================================
tdgelss=: 4 : 0
if. (3=x) *. 807>:0".}.({.~ i.&'/')9!:14'' do. 1 return. end.
zero=. (2|x){::dzero;zzero

'ma mvb'=. y
ma=. zero + ma
'm n'=. $ma
mvb=. zero + ,.^:(2>#@$)mvb
nrhs=. {:@$mvb
mn=. m<.n

if. 0=2|x do.
  lwork=. 1 >. ((3*mn)+((2*mn)>.(m>.n)>.nrhs))
  assert. 0= _1{::cdrc=. dgelss`0:`sgelss`0:@.x (,m);(,n);(,nrhs);(|:ma);(,1>.m);(|:ldb{.mvb);(,ldb=. 1>.m>.n);(s=. mn$dzero);(,rcond=. 1e_10);(,rank=. 0);(lwork$zero);(,lwork);,_1
else.
  lwork=. 1 >. ((2*mn)+(m>.n>.nrhs))
  rwork=. (5*mn)$dzero
  assert. 0= _1{::cdrc=. 0:`zgelss`0:`cgelss@.x (,m);(,n);(,nrhs);(|:ma);(,1>.m);(|:ldb{.mvb);(,ldb=. 1>.m>.n);(s=. mn$dzero);(,rcond=. 1e_10);(,rank=. 0);(lwork$zero);(,lwork);rwork;,_1
end.
'A B S'=. 4 6 8{cdrc
B=. n {. |: B
A=. mn ({."1) |:A

echo r=. mvb match`matchf@.(x>1) ma mp B
0{::r
)

NB. =========================================================
testdgelss=: 3 : 0
dma0=. 0 0$0
dmb0=. 0 0$0
dma1=. ? 10 5$100          NB. match fails for this pair since solution is least squares
dmb1=. ? 10 3$50
dma2=. ? 5 10$100
dmb2=. ? 5 3$50
dma3=. 0 0$0
dvb3=. 0$0
dma4=. ? 10 5$100          NB. match fails for this pair since solution is least squares
dvb4=. ? 10$50
dma5=. ? 5 10$100
dvb5=. ? 5$50
zma0=. 0 0$zzero
zmb0=. 0 0$zzero
zma1=. j./ ? 2 10 5$100    NB. match fails for this pair since solution is least squares
zmb1=. j./ ? 2 10 3$50
zma2=. j./ ? 2 5 10$100
zmb2=. j./ ? 2 5 3$50
zma3=. 0 0$zzero
zvb3=. 0$zzero
zma4=. j./ ? 2 10 5$100    NB. match fails for this pair since solution is least squares
zvb4=. j./ ? 2 10$50
zma5=. j./ ? 2 5 10$100
zvb5=. j./ ? 2 5$50
assert. 0 1 0 0 1 0 +. 0&tdgelss &> (< dma0;dmb0) , (< dma1;dmb1) , (< dma2;dmb2) , (< dma3;dvb3) , (< dma4;dvb4) , (< dma5;dvb5)
assert. 0 1 0 0 1 0 +. 1&tdgelss &> (< zma0;zmb0) , (< zma1;zmb1) , (< zma2;zmb2) , (< zma3;zvb3) , (< zma4;zvb4) , (< zma5;zvb5)
assert. 0 1 0 0 1 0 +. 2&tdgelss &> (< dma0;dmb0) , (< dma1;dmb1) , (< dma2;dmb2) , (< dma3;dvb3) , (< dma4;dvb4) , (< dma5;dvb5)
assert. 0 1 0 0 1 0 +. 3&tdgelss &> (< zma0;zmb0) , (< zma1;zmb1) , (< zma2;zmb2) , (< zma3;zvb3) , (< zma4;zvb4) , (< zma5;zvb5)
EMPTY
)

NB. least squares tests from http://www.nag.co.uk/lapack-ex/lapack-ex.html
NB.
NB.    dma=. 6 4 $ _0.57 _1.28 _0.39 0.25 _1.93 1.08 _0.31 _2.14 2.30 0.24 0.40 _0.35 _1.93 0.64 _0.66 0.08 0.15 0.30 0.15 _2.13 _0.02 1.03 _1.43 0.50
NB.    dvb=. _2.67 _0.55 3.34 _0.77 0.48 4.10
NB.    dvx=. 1.5339 1.8707 _1.5241 0.0392
NB.    tgelss dma;dvb
NB. 1.53387 1.87075 _1.52407 0.039183
NB.
NB.    zma=. 6 4 $ 0.96j_0.81 _0.03j0.96 _0.91j2.06 _0.05j0.41 _0.98j1.98 _1.20j0.19 _0.66j0.42 _0.81j0.56 0.62j_0.46 1.01j0.02 0.63j_0.17 _1.11j0.60 _0.37j0.38 0.19j_0.54 _0.98j_0.36 0.22j_0.20 0.83j0.51 0.20j0.01 _0.17j_0.46 1.47j1.59 1.08j_0.28 0.20j_0.12 _0.07j1.23 0.26j0.26
NB.    zvb=. _2.09j1.93 3.34j_3.53 _4.94j_2.04 0.17j4.23 _5.19j3.63 0.98j2.53
NB.    zvx=. _0.5044j_1.2179 _2.4281j2.8574 1.4872j_2.1955 0.4537j2.6904
NB.    tgelss zma;zvb
NB. _0.504365j_1.21788 _2.42815j2.85742 1.48717j_2.19548 0.453708j2.69041
