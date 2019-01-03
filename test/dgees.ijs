NB. subroutine dgees ( character                               JOBVS,
NB.                    character                               SORT,
NB.                    external                                SELECT,
NB.                    integer                                 N,
NB.                    double precision, dimension( lda, * )   A,
NB.                    integer                                 LDA,
NB.                    integer                                 SDIM,
NB.                    double precision, dimension( * )        WR,
NB.                    double precision, dimension( * )        WI,
NB.                    double precision, dimension( ldvs, * )  VS,
NB.                    integer                                 LDVS,
NB.                    double precision, dimension( * )        WORK,
NB.                    integer                                 LWORK,
NB.                    logical, dimension( * )                 BWORK,
NB.                    integer                                 INFO
NB.                  )
NB.
NB. DGEES computes the eigenvalues, the Schur form, and, optionally, the matrix of Schur
NB. vectors for GE matrices
NB.
NB. Purpose:
NB.
NB.      DGEES computes for an N-by-N real nonsymmetric matrix A, the
NB.      eigenvalues, the real Schur form T, and, optionally, the matrix of
NB.      Schur vectors Z.  This gives the Schur factorization A = Z*T*(Z**T).
NB.
NB.      Optionally, it also orders the eigenvalues on the diagonal of the
NB.      real Schur form so that selected eigenvalues are at the top left.
NB.      The leading columns of Z then form an orthonormal basis for the
NB.      invariant subspace corresponding to the selected eigenvalues.
NB.
NB.      A matrix is in real Schur form if it is upper quasi-triangular with
NB.      1-by-1 and 2-by-2 blocks. 2-by-2 blocks will be standardized in the
NB.      form
NB.              [  a  b  ]
NB.              [  c  a  ]
NB.
NB.      where b*c < 0. The eigenvalues of such a block are a +- sqrt(bc).
NB.
NB. Parameters
NB.
NB.                               JOBVS is CHARACTER*1
NB.     [in]     JOBVS            = 'N': Schur vectors are not computed;
NB.                               = 'V': Schur vectors are computed.
NB.
NB.                               SORT is CHARACTER*1
NB.                               Specifies whether or not to order the eigenvalues on the
NB.     [in]     SORT             diagonal of the Schur form.
NB.                               = 'N': Eigenvalues are not ordered;
NB.                               = 'S': Eigenvalues are ordered (see SELECT).
NB.
NB.                               SELECT is a LOGICAL FUNCTION of two DOUBLE PRECISION arguments
NB.                               SELECT must be declared EXTERNAL in the calling subroutine.
NB.                               If SORT = 'S', SELECT is used to select eigenvalues to sort
NB.                               to the top left of the Schur form.
NB.                               If SORT = 'N', SELECT is not referenced.
NB.                               An eigenvalue WR(j)+sqrt(-1)*WI(j) is selected if
NB.                               SELECT(WR(j),WI(j)) is true; i.e., if either one of a complex
NB.     [in]     SELECT           conjugate pair of eigenvalues is selected, then both complex
NB.                               eigenvalues are selected.
NB.                               Note that a selected complex eigenvalue may no longer
NB.                               satisfy SELECT(WR(j),WI(j)) = .TRUE. after ordering, since
NB.                               ordering may change the value of complex eigenvalues
NB.                               (especially if the eigenvalue is ill-conditioned); in this
NB.                               case INFO is set to N+2 (see INFO below).
NB.
NB.                               N is INTEGER
NB.     [in]     N                The order of the matrix A. N >= 0.
NB.
NB.                               A is DOUBLE PRECISION array, dimension (LDA,N)
NB.     [in,out] A                On entry, the N-by-N matrix A.
NB.                               On exit, A has been overwritten by its real Schur form T.
NB.
NB.                               LDA is INTEGER
NB.     [in]     LDA              The leading dimension of the array A.  LDA >= max(1,N).
NB.
NB.                               SDIM is INTEGER
NB.                               If SORT = 'N', SDIM = 0.
NB.                               If SORT = 'S', SDIM = number of eigenvalues (after sorting)
NB.     [out]    SDIM                            for which SELECT is true. (Complex conjugate
NB.                                              pairs for which SELECT is true for either
NB.                                              eigenvalue count as 2.)
NB.
NB.     [out]    WR               WR is DOUBLE PRECISION array, dimension (N)
NB.
NB.                               WI is DOUBLE PRECISION array, dimension (N)
NB.                               WR and WI contain the real and imaginary parts,
NB.                               respectively, of the computed eigenvalues in the same order
NB.     [out]    WI               that they appear on the diagonal of the output Schur form T.
NB.                               Complex conjugate pairs of eigenvalues will appear
NB.                               consecutively with the eigenvalue having the positive
NB.                               imaginary part first.
NB.
NB.                               VS is DOUBLE PRECISION array, dimension (LDVS,N)
NB.                               If JOBVS = 'V', VS contains the orthogonal matrix Z of Schur
NB.     [out]    VS               vectors.
NB.                               If JOBVS = 'N', VS is not referenced.
NB.
NB.                               LDVS is INTEGER
NB.     [in]     LDVS             The leading dimension of the array VS.  LDVS >= 1; if
NB.                               JOBVS = 'V', LDVS >= N.
NB.
NB.                               WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
NB.     [out]    WORK             On exit, if INFO = 0, WORK(1) contains the optimal LWORK.
NB.
NB.                               LWORK is INTEGER
NB.                               The dimension of the array WORK.  LWORK >= max(1,3*N).
NB.                               For good performance, LWORK must generally be larger.
NB.
NB.     [in]     LWORK            If LWORK = -1, then a workspace query is assumed; the routine
NB.                               only calculates the optimal size of the WORK array, returns
NB.                               this value as the first entry of the WORK array, and no error
NB.                               message related to LWORK is issued by XERBLA.
NB.
NB.                               BWORK is LOGICAL array, dimension (N)
NB.     [out]    BWORK            Not referenced if SORT = 'N'.
NB.
NB.                               INFO is INTEGER
NB.                               = 0: successful exit
NB.                               < 0: if INFO = -i, the i-th argument had an illegal value.
NB.                               > 0: if INFO = i, and i is
NB.                                  <= N: the QR algorithm failed to compute all the
NB.                                        eigenvalues; elements 1:ILO-1 and i+1:N of WR and WI
NB.                                        contain those eigenvalues which have converged; if
NB.     [out]    INFO                      JOBVS = 'V', VS contains the matrix which reduces A
NB.                                        to its partially converged Schur form.
NB.                                  = N+1: the eigenvalues could not be reordered because some
NB.                                        eigenvalues were too close to separate (the problem
NB.                                        is very ill-conditioned);
NB.                                  = N+2: after reordering, roundoff changed values of some
NB.                                        complex eigenvalues so that leading eigenvalues in
NB.                                        the Schur form no longer satisfy SELECT=.TRUE.  This
NB.                                        could also be caused by underflow due to scaling.

require 'math/lapack2'
cocurrent 'base'
coinsert 'jlapack2'
match=: matchclean;;
matchf=: matchcleanf;;

NB. =========================================================
tdgees=: 4 : 0
if. (3=x) *. 807>:0".}.({.~ i.&'/')9!:14'' do. 1 return. end.
zero=. (2|x){::dzero;zzero
a=. zero + y
'm n'=. $a
if. 0=2|x do.
  assert. 0= _1{::cdrc=. dgees`0:`sgees`0:@.x (,'V');(,'N');(<0);(,n);(|:a);(,1>.m);(,0);(WR=. n$dzero);(WI=. n$dzero);(Z=. (,~n)$dzero);(,1>.n);(lwork$zero);(,lwork=. 1>.3*n);(,n$0);,_1
  'T WR WI Z'=. 5 8 9 10{cdrc
  T=. |:T
  Z=. |:Z
  W=. WR j. WI
else.
  assert. 0= _1{::cdrc=. 0:`zgees`0:`cgees@.x (,'V');(,'N');(<0);(,n);(|:a);(,1>.m);(,0);(W=. n$zzero);(Z=. (,~n)$zzero);(,1>.n);(lwork$zzero);(,lwork=. 1>.2*n);(n$dzero);(,n$0);,_1
  'T W Z'=. 5 8 9{cdrc
  T=. |:T
  W=. |:W
  Z=. |:Z
end.
echo Z;T;W
echo r=. a match`matchf@.(x>1) clean`cleanf@.(x>1) Z mp T mp +|:Z
0{::r
)

NB. =========================================================
testdgees=: 3 : 0
m0=. 0 0$0
m1=. ?.6 6$10
m2=. 0 0$zzero
m3=. j./ ?.2 6 6$10
assert. 0&tdgees &> m0;m1;m2
assert. 1&tdgees m3
assert. 2&tdgees &> m0;m1;m2
assert. 3&tdgees m3
EMPTY
)

