NB. subroutine dgeev ( character                               JOBVL,
NB.                    character                               JOBVR,
NB.                    integer                                 N,
NB.                    double precision, dimension( lda, * )   A,
NB.                    integer                                 LDA,
NB.                    double precision, dimension( * )        WR,
NB.                    double precision, dimension( * )        WI,
NB.                    double precision, dimension( ldvl, * )  VL,
NB.                    integer                                 LDVL,
NB.                    double precision, dimension( ldvr, * )  VR,
NB.                    integer                                 LDVR,
NB.                    double precision, dimension( * )        WORK,
NB.                    integer                                 LWORK,
NB.                    integer                                 INFO
NB.                  )
NB.
NB. DGEEV computes the eigenvalues and, optionally, the left and/or right eigenvectors for GE
NB. matrices
NB.
NB. Purpose:
NB.
NB.      DGEEV computes for an N-by-N real nonsymmetric matrix A, the
NB.      eigenvalues and, optionally, the left and/or right eigenvectors.
NB.
NB.      The right eigenvector v(j) of A satisfies
NB.                       A * v(j) = lambda(j) * v(j)
NB.      where lambda(j) is its eigenvalue.
NB.      The left eigenvector u(j) of A satisfies
NB.                    u(j)**H * A = lambda(j) * u(j)**H
NB.      where u(j)**H denotes the conjugate-transpose of u(j).
NB.
NB.      The computed eigenvectors are normalized to have Euclidean norm
NB.      equal to 1 and largest component real.
NB.
NB. Parameters
NB.
NB.                              JOBVL is CHARACTER*1
NB.     [in]     JOBVL           = 'N': left eigenvectors of A are not computed;
NB.                              = 'V': left eigenvectors of A are computed.
NB.
NB.                              JOBVR is CHARACTER*1
NB.     [in]     JOBVR           = 'N': right eigenvectors of A are not computed;
NB.                              = 'V': right eigenvectors of A are computed.
NB.
NB.                              N is INTEGER
NB.     [in]     N               The order of the matrix A. N >= 0.
NB.
NB.                              A is DOUBLE PRECISION array, dimension (LDA,N)
NB.     [in,out] A               On entry, the N-by-N matrix A.
NB.                              On exit, A has been overwritten.
NB.
NB.                              LDA is INTEGER
NB.     [in]     LDA             The leading dimension of the array A.  LDA >= max(1,N).
NB.
NB.     [out]    WR              WR is DOUBLE PRECISION array, dimension (N)
NB.
NB.                              WI is DOUBLE PRECISION array, dimension (N)
NB.                              WR and WI contain the real and imaginary parts,
NB.                              respectively, of the computed eigenvalues.  Complex
NB.     [out]    WI              conjugate pairs of eigenvalues appear consecutively
NB.                              with the eigenvalue having the positive imaginary part
NB.                              first.
NB.
NB.                              VL is DOUBLE PRECISION array, dimension (LDVL,N)
NB.                              If JOBVL = 'V', the left eigenvectors u(j) are stored one
NB.                              after another in the columns of VL, in the same order
NB.                              as their eigenvalues.
NB.                              If JOBVL = 'N', VL is not referenced.
NB.     [out]    VL              If the j-th eigenvalue is real, then u(j) = VL(:,j),
NB.                              the j-th column of VL.
NB.                              If the j-th and (j+1)-st eigenvalues form a complex
NB.                              conjugate pair, then u(j) = VL(:,j) + i*VL(:,j+1) and
NB.                              u(j+1) = VL(:,j) - i*VL(:,j+1).
NB.
NB.                              LDVL is INTEGER
NB.     [in]     LDVL            The leading dimension of the array VL.  LDVL >= 1; if
NB.                              JOBVL = 'V', LDVL >= N.
NB.
NB.                              VR is DOUBLE PRECISION array, dimension (LDVR,N)
NB.                              If JOBVR = 'V', the right eigenvectors v(j) are stored one
NB.                              after another in the columns of VR, in the same order
NB.                              as their eigenvalues.
NB.                              If JOBVR = 'N', VR is not referenced.
NB.     [out]    VR              If the j-th eigenvalue is real, then v(j) = VR(:,j),
NB.                              the j-th column of VR.
NB.                              If the j-th and (j+1)-st eigenvalues form a complex
NB.                              conjugate pair, then v(j) = VR(:,j) + i*VR(:,j+1) and
NB.                              v(j+1) = VR(:,j) - i*VR(:,j+1).
NB.
NB.                              LDVR is INTEGER
NB.     [in]     LDVR            The leading dimension of the array VR.  LDVR >= 1; if
NB.                              JOBVR = 'V', LDVR >= N.
NB.
NB.                              WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
NB.     [out]    WORK            On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
NB.
NB.                              LWORK is INTEGER
NB.                              The dimension of the array WORK.  LWORK >= max(1,3*N), and
NB.                              if JOBVL = 'V' or JOBVR = 'V', LWORK >= 4*N.  For good
NB.                              performance, LWORK must generally be larger.
NB.     [in]     LWORK
NB.                              If LWORK = -1, then a workspace query is assumed; the routine
NB.                              only calculates the optimal size of the WORK array, returns
NB.                              this value as the first entry of the WORK array, and no error
NB.                              message related to LWORK is issued by XERBLA.
NB.
NB.                              INFO is INTEGER
NB.                              = 0:  successful exit
NB.                              < 0:  if INFO = -i, the i-th argument had an illegal value.
NB.     [out]    INFO            > 0:  if INFO = i, the QR algorithm failed to compute all the
NB.                                    eigenvalues, and no eigenvectors have been computed;
NB.                                    elements i+1:N of WR and WI contain eigenvalues which
NB.                                    have converged.

require 'math/lapack2'
cocurrent 'base'
coinsert 'jlapack2'
match=: matchclean;;
matchf=: matchcleanf;;

NB. =========================================================
tdgeev=: 4 : 0
if. (3=x) *. 807>:0".}.({.~ i.&'/')9!:14'' do. 1 return. end.
zero=. (2|x){::dzero;zzero
a=. zero + y
'm n'=. $a

if. 0=2|x do.
  assert. 0= _1{::cdrc=. dgeev`0:`sgeev`0:@.x (,'V');(,'V');(,n);(|:a);(,1>.m);(wr=. n$zero);(wi=. n$zero);(vl=. (n,n)$zero);(,1>.n);(vr=. (n,n)$zero);(,1>.n);(lwork$zero);(,lwork=. 1>.4*n);,_1
  'wr wi vl vr'=. 6 7 8 10{cdrc
  w=. wr j. wi
else.
  assert. 0= _1{::cdrc=. 0:`zgeev`0:`cgeev@.x (,'V');(,'V');(,n);(|:a);(,1>.m);(w=. n$zero);(vl=. (n,n)$zero);(,1>.n);(vr=. (n,n)$zero);(,1>.n);(lwork$zero);(,lwork=. 1>.4*n);((+:n)$dzero);,_1
  'w vl vr'=. 6 7 9{cdrc
end.
V=. w
L=. n{. |:vl
R=. n{. |:vr

if. 0=2|x do.
  if. #cx=. I. wi ~: 0 do.
    L=. cx cxpair L
    R=. cx cxpair R
  end.
end.

echo L;V;R
echo a=. (y mp R) match`matchf@.(x>1) (V *"1 R)
echo b=. ((+|:L) mp y) match`matchf@.(x>1) (V * +|:L)
(0 pick a) *. 0 pick b
)

NB. =========================================================
testdgeev=: 3 : 0
m0=. 0 0$0
m1=. ?.6 6$10
m2=. 0 0$zzero
m3=. j./ ?.2 6 6$10
assert. 0&tdgeev &> m0;m1
assert. 1&tdgeev &> m1;m2;m3
assert. 2&tdgeev &> m0;m1
assert. 3&tdgeev &> m1;m2;m3
EMPTY
)
