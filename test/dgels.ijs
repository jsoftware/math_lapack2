NB. subroutine dgels ( character                              TRANS,
NB.                    integer                                M,
NB.                    integer                                N,
NB.                    integer                                NRHS,
NB.                    double precision, dimension( lda, * )  A,
NB.                    integer                                LDA,
NB.                    double precision, dimension( ldb, * )  B,
NB.                    integer                                LDB,
NB.                    double precision, dimension( * )       WORK,
NB.                    integer                                LWORK,
NB.                    integer                                INFO
NB.                  )
NB.
NB. DGELS solves overdetermined or underdetermined systems for GE matrices
NB.
NB. Purpose:
NB.
NB.      DGELS solves overdetermined or underdetermined real linear systems
NB.      involving an M-by-N matrix A, or its transpose, using a QR or LQ
NB.      factorization of A.  It is assumed that A has full rank.
NB.
NB.      The following options are provided:
NB.
NB.      1. If TRANS = 'N' and m >= n:  find the least squares solution of
NB.         an overdetermined system, i.e., solve the least squares problem
NB.                      minimize || B - A*X ||.
NB.
NB.      2. If TRANS = 'N' and m < n:  find the minimum norm solution of
NB.         an underdetermined system A * X = B.
NB.
NB.      3. If TRANS = 'T' and m >= n:  find the minimum norm solution of
NB.         an underdetermined system A**T * X = B.
NB.
NB.      4. If TRANS = 'T' and m < n:  find the least squares solution of
NB.         an overdetermined system, i.e., solve the least squares problem
NB.                      minimize || B - A**T * X ||.
NB.
NB.      Several right hand side vectors b and solution vectors x can be
NB.      handled in a single call; they are stored as the columns of the
NB.      M-by-NRHS right hand side matrix B and the N-by-NRHS solution
NB.      matrix X.
NB.
NB. Parameters
NB.
NB.                              TRANS is CHARACTER*1
NB.     [in]     TRANS           = 'N': the linear system involves A;
NB.                              = 'T': the linear system involves A**T.
NB.
NB.                              M is INTEGER
NB.     [in]     M               The number of rows of the matrix A.  M >= 0.
NB.
NB.                              N is INTEGER
NB.     [in]     N               The number of columns of the matrix A.  N >= 0.
NB.
NB.                              NRHS is INTEGER
NB.     [in]     NRHS            The number of right hand sides, i.e., the number of
NB.                              columns of the matrices B and X. NRHS >=0.
NB.
NB.                              A is DOUBLE PRECISION array, dimension (LDA,N)
NB.                              On entry, the M-by-N matrix A.
NB.                              On exit,
NB.     [in,out] A                 if M >= N, A is overwritten by details of its QR
NB.                                           factorization as returned by DGEQRF;
NB.                                if M <  N, A is overwritten by details of its LQ
NB.                                           factorization as returned by DGELQF.
NB.
NB.                              LDA is INTEGER
NB.     [in]     LDA             The leading dimension of the array A.  LDA >= max(1,M).
NB.
NB.                              B is DOUBLE PRECISION array, dimension (LDB,NRHS)
NB.                              On entry, the matrix B of right hand side vectors, stored
NB.                              columnwise; B is M-by-NRHS if TRANS = 'N', or N-by-NRHS
NB.                              if TRANS = 'T'.
NB.                              On exit, if INFO = 0, B is overwritten by the solution
NB.                              vectors, stored columnwise:
NB.                              if TRANS = 'N' and m >= n, rows 1 to n of B contain the least
NB.                              squares solution vectors; the residual sum of squares for the
NB.                              solution in each column is given by the sum of squares of
NB.     [in,out] B               elements N+1 to M in that column;
NB.                              if TRANS = 'N' and m < n, rows 1 to N of B contain the
NB.                              minimum norm solution vectors;
NB.                              if TRANS = 'T' and m >= n, rows 1 to M of B contain the
NB.                              minimum norm solution vectors;
NB.                              if TRANS = 'T' and m < n, rows 1 to M of B contain the
NB.                              least squares solution vectors; the residual sum of squares
NB.                              for the solution in each column is given by the sum of
NB.                              squares of elements M+1 to N in that column.
NB.
NB.                              LDB is INTEGER
NB.     [in]     LDB             The leading dimension of the array B. LDB >= MAX(1,M,N).
NB.
NB.                              WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
NB.     [out]    WORK            On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
NB.
NB.                              LWORK is INTEGER
NB.                              The dimension of the array WORK.
NB.                              LWORK >= max( 1, MN + max( MN, NRHS ) ).
NB.                              For optimal performance,
NB.                              LWORK >= max( 1, MN + max( MN, NRHS )*NB ).
NB.     [in]     LWORK           where MN = min(M,N) and NB is the optimum block size.
NB.
NB.                              If LWORK = -1, then a workspace query is assumed; the routine
NB.                              only calculates the optimal size of the WORK array, returns
NB.                              this value as the first entry of the WORK array, and no error
NB.                              message related to LWORK is issued by XERBLA.
NB.
NB.                              INFO is INTEGER
NB.                              = 0:  successful exit
NB.                              < 0:  if INFO = -i, the i-th argument had an illegal value
NB.     [out]    INFO            > 0:  if INFO =  i, the i-th diagonal element of the
NB.                                    triangular factor of A is zero, so that A does not have
NB.                                    full rank; the least squares solution could not be
NB.                                    computed.

require 'math/lapack2'
cocurrent 'base'
coinsert 'jlapack2'
match=: matchclean;;
matchf=: matchcleanf;;

NB. =========================================================
tdgels=: 4 : 0
if. (3=x) *. 807>:0".}.({.~ i.&'/')9!:14'' do. 1 return. end.
zero=. (2|x){::dzero;zzero
'ma mvb'=. y
ma=. zero + ma
'm n'=. $ma
mvb=. zero + ,.^:(2>#@$)mvb
nrhs=. {:@$mvb
lwork=. 1 >. (m<.n) + ((m<.n) >. nrhs)*4
assert. 0= _1{::cdrc=. dgels`zgels`sgels`cgels@.x (,'N');(,m);(,n);(,nrhs);(|:ma);(,1>.m);(|:ldb{.mvb);(,ldb=. 1>.m>.n);(lwork$zero);(,lwork);,_1
R=. n{. |: 7{::cdrc
echo R
echo r=. mvb match`matchf@.(x>1) clean`cleanf@.(x>1) ma mp R
0{::r
)

NB. =========================================================
testdgels=: 3 : 0
dma0=. 0 0$0
dmb0=. 0 0$0
dma1=. ?. 10 5$100          NB. match fails for this pair since solution is least squares
dmb1=. ?. 10 3$50
dma2=. ?. 5 10$100
dmb2=. ?. 5 3$50
dma3=. 0 0$0
dvb3=. 0$0
dma4=. ?. 10 5$100          NB. match fails for this pair since solution is least squares
dvb4=. ?. 10$50
dma5=. ?. 5 10$100
dvb5=. ?. 5$50
zma0=. 0 0$zzero
zmb0=. 0 0$zzero
zma1=. j./ ?. 2 10 5$100    NB. match fails for this pair since solution is least squares
zmb1=. j./ ?. 2 10 3$50
zma2=. j./ ?. 2 5 10$100
zmb2=. j./ ?. 2 5 3$50
zma3=. 0 0$zzero
zvb3=. 0$zzero
zma4=. j./ ?. 2 10 5$100    NB. match fails for this pair since solution is least squares
zvb4=. j./ ?. 2 10$50
zma5=. j./ ?. 2 5 10$100
zvb5=. j./ ?. 2 5$50
assert. 0 1 0 0 1 0 +. 0&tdgels &> (< dma0;dmb0) , (< dma1;dmb1) , (< dma2;dmb2) , (< dma3;dvb3) , (< dma4;dvb4) , (< dma5;dvb5)
assert. 0 1 0 0 1 0 +. 1&tdgels &> (< zma0;zmb0) , (< zma1;zmb1) , (< zma2;zmb2) , (< zma3;zvb3) , (< zma4;zvb4) , (< zma5;zvb5)
assert. 0 1 0 0 1 0 +. 2&tdgels &> (< dma0;dmb0) , (< dma1;dmb1) , (< dma2;dmb2) , (< dma3;dvb3) , (< dma4;dvb4) , (< dma5;dvb5)
assert. 0 1 0 0 1 0 +. 3&tdgels &> (< zma0;zmb0) , (< zma1;zmb1) , (< zma2;zmb2) , (< zma3;zvb3) , (< zma4;zvb4) , (< zma5;zvb5)
EMPTY
)
