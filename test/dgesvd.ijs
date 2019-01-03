NB. subroutine dgesvd ( character                               JOBU,
NB.                     character                               JOBVT,
NB.                     integer                                 M,
NB.                     integer                                 N,
NB.                     double precision, dimension( lda, * )   A,
NB.                     integer                                 LDA,
NB.                     double precision, dimension( * )        S,
NB.                     double precision, dimension( ldu, * )   U,
NB.                     integer                                 LDU,
NB.                     double precision, dimension( ldvt, * )  VT,
NB.                     integer                                 LDVT,
NB.                     double precision, dimension( * )        WORK,
NB.                     integer                                 LWORK,
NB.                     integer                                 INFO
NB.                   )
NB.
NB. DGESVD computes the singular value decomposition (SVD) for GE matrices
NB.
NB. Purpose:
NB.
NB.      DGESVD computes the singular value decomposition (SVD) of a real
NB.      M-by-N matrix A, optionally computing the left and/or right singular
NB.      vectors. The SVD is written
NB.
NB.           A = U * SIGMA * transpose(V)
NB.
NB.      where SIGMA is an M-by-N matrix which is zero except for its
NB.      min(m,n) diagonal elements, U is an M-by-M orthogonal matrix, and
NB.      V is an N-by-N orthogonal matrix.  The diagonal elements of SIGMA
NB.      are the singular values of A; they are real and non-negative, and
NB.      are returned in descending order.  The first min(m,n) columns of
NB.      U and V are the left and right singular vectors of A.
NB.
NB.      Note that the routine returns V**T, not V.
NB.
NB. Parameters
NB.
NB.                              JOBU is CHARACTER*1
NB.                              Specifies options for computing all or part of the matrix U:
NB.                              = 'A':  all M columns of U are returned in array U:
NB.                              = 'S':  the first min(m,n) columns of U (the left singular
NB.     [in]     JOBU                    vectors) are returned in the array U;
NB.                              = 'O':  the first min(m,n) columns of U (the left singular
NB.                                      vectors) are overwritten on the array A;
NB.                              = 'N':  no columns of U (no left singular vectors) are
NB.                                      computed.
NB.
NB.                              JOBVT is CHARACTER*1
NB.                              Specifies options for computing all or part of the matrix
NB.                              V**T:
NB.                              = 'A':  all N rows of V**T are returned in the array VT;
NB.                              = 'S':  the first min(m,n) rows of V**T (the right singular
NB.                                      vectors) are returned in the array VT;
NB.     [in]     JOBVT           = 'O':  the first min(m,n) rows of V**T (the right singular
NB.                                      vectors) are overwritten on the array A;
NB.                              = 'N':  no rows of V**T (no right singular vectors) are
NB.                                      computed.
NB.
NB.                              JOBVT and JOBU cannot both be 'O'.
NB.
NB.                              M is INTEGER
NB.     [in]     M               The number of rows of the input matrix A.  M >= 0.
NB.
NB.                              N is INTEGER
NB.     [in]     N               The number of columns of the input matrix A.  N >= 0.
NB.
NB.                              A is DOUBLE PRECISION array, dimension (LDA,N)
NB.                              On entry, the M-by-N matrix A.
NB.                              On exit,
NB.                              if JOBU = 'O',  A is overwritten with the first min(m,n)
NB.                                              columns of U (the left singular vectors,
NB.     [in,out] A                               stored columnwise);
NB.                              if JOBVT = 'O', A is overwritten with the first min(m,n)
NB.                                              rows of V**T (the right singular vectors,
NB.                                              stored rowwise);
NB.                              if JOBU .ne. 'O' and JOBVT .ne. 'O', the contents of A
NB.                                              are destroyed.
NB.
NB.                              LDA is INTEGER
NB.     [in]     LDA             The leading dimension of the array A.  LDA >= max(1,M).
NB.
NB.                              S is DOUBLE PRECISION array, dimension (min(M,N))
NB.     [out]    S               The singular values of A, sorted so that S(i) >= S(i+1).
NB.
NB.                              U is DOUBLE PRECISION array, dimension (LDU,UCOL)
NB.                              (LDU,M) if JOBU = 'A' or (LDU,min(M,N)) if JOBU = 'S'.
NB.                              If JOBU = 'A', U contains the M-by-M orthogonal matrix U;
NB.     [out]    U               if JOBU = 'S', U contains the first min(m,n) columns of U
NB.                              (the left singular vectors, stored columnwise);
NB.                              if JOBU = 'N' or 'O', U is not referenced.
NB.
NB.                              LDU is INTEGER
NB.     [in]     LDU             The leading dimension of the array U.  LDU >= 1; if
NB.                              JOBU = 'S' or 'A', LDU >= M.
NB.
NB.                              VT is DOUBLE PRECISION array, dimension (LDVT,N)
NB.                              If JOBVT = 'A', VT contains the N-by-N orthogonal matrix
NB.                              V**T;
NB.     [out]    VT              if JOBVT = 'S', VT contains the first min(m,n) rows of
NB.                              V**T (the right singular vectors, stored rowwise);
NB.                              if JOBVT = 'N' or 'O', VT is not referenced.
NB.
NB.                              LDVT is INTEGER
NB.     [in]     LDVT            The leading dimension of the array VT.  LDVT >= 1; if
NB.                              JOBVT = 'A', LDVT >= N; if JOBVT = 'S', LDVT >= min(M,N).
NB.
NB.                              WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
NB.                              On exit, if INFO = 0, WORK(1) returns the optimal LWORK;
NB.                              if INFO > 0, WORK(2:MIN(M,N)) contains the unconverged
NB.     [out]    WORK            superdiagonal elements of an upper bidiagonal matrix B
NB.                              whose diagonal is in S (not necessarily sorted). B
NB.                              satisfies A = U * B * VT, so it has the same singular values
NB.                              as A, and singular vectors related by U and VT.
NB.
NB.                              LWORK is INTEGER
NB.                              The dimension of the array WORK.
NB.                              LWORK >= MAX(1,5*MIN(M,N)) for the paths (see comments inside code):
NB.                                 - PATH 1  (M much larger than N, JOBU='N')
NB.                                 - PATH 1t (N much larger than M, JOBVT='N')
NB.                              LWORK >= MAX(1,3*MIN(M,N) + MAX(M,N),5*MIN(M,N)) for the other paths
NB.     [in]     LWORK           For good performance, LWORK should generally be larger.
NB.
NB.                              If LWORK = -1, then a workspace query is assumed; the routine
NB.                              only calculates the optimal size of the WORK array, returns
NB.                              this value as the first entry of the WORK array, and no error
NB.                              message related to LWORK is issued by XERBLA.
NB.
NB.                              INFO is INTEGER
NB.                              = 0:  successful exit.
NB.                              < 0:  if INFO = -i, the i-th argument had an illegal value.
NB.     [out]    INFO            > 0:  if DBDSQR did not converge, INFO specifies how many
NB.                                    superdiagonals of an intermediate bidiagonal form B
NB.                                    did not converge to zero. See the description of WORK
NB.                                    above for details.

require 'math/lapack2'
cocurrent 'base'
coinsert 'jlapack2'
match=: matchclean;;
matchf=: matchcleanf;;

NB. =========================================================
tdgesvd=: 4 : 0
if. (3=x) *. 807>:0".}.({.~ i.&'/')9!:14'' do. 1 return. end.
zero=. (2|x){::dzero;zzero
a=. zero + y
sa=. |. 'm n'=. $a
mn=. m<.n
sui=. m,m
svt=. n,n
lda=. ldu=. 1>.m
S=. mn$dzero
U=. sui$zero
VT=. svt$zero
lwork=. 1 >. (2|x) { (((3*mn)+(m>.n))>.(5*mn)) , ((2*mn)+(m>.n))
work=. lwork$zero
rwork=. ((3*mn)>.(_4+5*mn))$dzero
if. 0=2|x do.
  assert. 0= _1{::cdrc=. dgesvd`0:`sgesvd`0:@.x (,'A');(,'A');(,m);(,n);(|:a);(,1>.m);S;U;(,1>.m);VT;(,1>.n);(lwork$zero);(,lwork);,_1
else.
  assert. 0= _1{::cdrc=. 0:`zgesvd`0:`cgesvd@.x (,'A');(,'A');(,m);(,n);(|:a);(,1>.m);S;U;(,1>.m);VT;(,1>.n);(lwork$zero);(,lwork);rwork;,_1
end.
'S U VT'=. 7 8 10{cdrc
U=. |:U
VT=. +VT
S=. (m-n) diagmat S
echo U;S;VT
echo r=. a match`matchf@.(x>1) clean`cleanf@.(x>1) U mp S mp +|:VT
0{::r
)

NB. =========================================================
testdgesvd=: 3 : 0
m0=. 0 0$0
m1=. ?.4 4$10
m2=. ?.4 6$10
m3=. ?.6 4$10
m4=. 0 0$zzero
m5=. j./ ?.2 4 4$10
m6=. j./ ?.2 4 6$10
m7=. j./ ?.2 6 4$10
assert. 0&tdgesvd &> m0;m1;m2;m3;m4
assert. 1&tdgesvd &> m5;m6;m7
assert. 2&tdgesvd &> m0;m1;m2;m3;m4
assert. 3&tdgesvd &> m5;m6;m7
EMPTY
)

