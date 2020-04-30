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

NB. Singular value decomposition (svd) sample data
NB.
NB.      1    0.5 0.3333   0.25
NB.    0.5 0.3333   0.25    0.2
NB. 0.3333   0.25    0.2 0.1667

NB. return s u vt
do_dgesvd=: 3 : 0
assert. 2=#$y
a=. y
'm n'=. ,"0 $a
mn=. m<.n

NB. call with lwork = _1 to query optimal workspace size
NB. lapack expect column major order |:a
assert. 0= LASTINFO=: _1{::cdrc=. dgesvd_jlapack2_ (,'A');(,'A');m;n;(|:a);(1>.m);(mn$0.0);(0.0$~ldu,m);(ldu=. 1>.m);(0.0$~ldvt,n);(ldvt=. 1>.n);(,0.0);(,_1);,_1

lwork=. <. _3{::cdrc

NB. call again with lwork
assert. 0= LASTINFO=: _1{::cdrc=. dgesvd_jlapack2_ (_3}.}.cdrc),(lwork$0.0);lwork;,_1

's u vt'=. 7 8 10{cdrc
s;(m{.|:u);n{.|:vt
)

a=: ".;._2[0 : 0
     1    0.5 0.3333   0.25
   0.5 0.3333   0.25    0.2
0.3333   0.25    0.2 0.1667
)

's u vt'=: do_dgesvd a
echo (,.s);u;vt
echo sigma=:  ($a)$,s,("0 1) 0#~{:$a   NB. diagonal is s
echo a;u (+/ .*) sigma (+/ .*) vt
