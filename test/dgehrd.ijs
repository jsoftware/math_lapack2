NB. subroutine dgehrd ( integer                                N,
NB.                     integer                                ILO,
NB.                     integer                                IHI,
NB.                     double precision, dimension( lda, * )  A,
NB.                     integer                                LDA,
NB.                     double precision, dimension( * )       TAU,
NB.                     double precision, dimension( * )       WORK,
NB.                     integer                                LWORK,
NB.                     integer                                INFO
NB.                   )
NB.
NB. DGEHRD
NB.
NB. Purpose:
NB.
NB.      DGEHRD reduces a real general matrix A to upper Hessenberg form H by
NB.      an orthogonal similarity transformation:  Q**T * A * Q = H .
NB.
NB. Parameters
NB.
NB.                              N is INTEGER
NB.     [in]     N               The order of the matrix A.  N >= 0.
NB.
NB.     [in]     ILO             ILO is INTEGER
NB.
NB.                              IHI is INTEGER
NB.
NB.                              It is assumed that A is already upper triangular in rows
NB.     [in]     IHI             and columns 1:ILO-1 and IHI+1:N. ILO and IHI are normally
NB.                              set by a previous call to DGEBAL; otherwise they should be
NB.                              set to 1 and N respectively. See Further Details.
NB.                              1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.
NB.
NB.                              A is DOUBLE PRECISION array, dimension (LDA,N)
NB.                              On entry, the N-by-N general matrix to be reduced.
NB.                              On exit, the upper triangle and the first subdiagonal of A
NB.     [in,out] A               are overwritten with the upper Hessenberg matrix H, and the
NB.                              elements below the first subdiagonal, with the array TAU,
NB.                              represent the orthogonal matrix Q as a product of elementary
NB.                              reflectors. See Further Details.
NB.
NB.                              LDA is INTEGER
NB.     [in]     LDA             The leading dimension of the array A.  LDA >= max(1,N).
NB.
NB.                              TAU is DOUBLE PRECISION array, dimension (N-1)
NB.                              The scalar factors of the elementary reflectors (see Further
NB.     [out]    TAU             Details). Elements 1:ILO-1 and IHI:N-1 of TAU are set to
NB.                              zero.
NB.
NB.                              WORK is DOUBLE PRECISION array, dimension (LWORK)
NB.     [out]    WORK            On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
NB.
NB.                              LWORK is INTEGER
NB.                              The length of the array WORK.  LWORK >= max(1,N).
NB.                              For good performance, LWORK should generally be larger.
NB.
NB.     [in]     LWORK           If LWORK = -1, then a workspace query is assumed; the routine
NB.                              only calculates the optimal size of the WORK array, returns
NB.                              this value as the first entry of the WORK array, and no error
NB.                              message related to LWORK is issued by XERBLA.
NB.
NB.                              INFO is INTEGER
NB.     [out]    INFO            = 0:  successful exit
NB.                              < 0:  if INFO = -i, the i-th argument had an illegal value.

require 'math/lapack2'
cocurrent 'base'
coinsert 'jlapack2'
match=: matchclean;;
matchf=: matchcleanf;;

NB. =========================================================
tdgehrd=: 4 : 0
if. (3=x) *. 807>:0".}.({.~ i.&'/')9!:14'' do. 1 return. end.
zero=. (2|x){::dzero;zzero
'a ilo ihi'=. y

a0=. a:{ a=. zero + a
'm n'=. $a

if. (1 > ilo) +. ((ilo > ihi) *. ((0 ~: n) +. (1 ~: ilo) +. (0 ~: ihi))) +. (ihi > n) do.
  sminfo 'following should hold: 1 <: ILO <: IHI <: N'
  0 return.
end.

assert. 0= _1{::cdrc=. dgehrd`zgehrd`sgehrd`cgehrd@.x (,m);(,ilo);(,ihi);(|:a);(,1>.m);(,tau=. (0 >. n-1)$zero);(lwork$zero);(,lwork=. 1>.n);,_1

'ilo ihi a tau'=. 2 3 4 6{cdrc
a=. |: a
ilo=. {. ilo
iho=. {. iho
h=. _1 utri a

ldiff=. 0 >. ihi-ilo  NB. '>.' to fix case 0=n
qsize=. n , ldiff
nilo=. - {. ilo
qvec=. (nilo idmat qsize) + (((0 , <: ilo) ,: qsize) (nilo & sltri) ;. 0 a)
q=. mp/ (idmat n) -"2 (((0 >. <: ilo) ,: ldiff) ] ;. 0 tau) * (* +)"0/~"1 |:qvec

echo h;q
echo r=. h match`matchf@.(x>1) clean`cleanf@.(x>1) (+|:q) mp a0 mp q
0{::r
)

NB. =========================================================
testdgehrd=: 3 : 0
ios=. 1 0;2 0;3 0;4 0;5 0;6 0;7 0;8 0;9 0;9 1;9 2;9 3;9 4;9 5;9 6;9 7;9 8
ma0=. 0 0$0
ma1=. ? 10 10$100
ma2=. 0 ios } ma1
ma3=. 0 0$zzero
ma4=. j./ ? 2 10 10$100
ma5=. 0 ios } ma4
assert. 0&tdgehrd &> (< ma0;1;(#ma0)) , (< ma1;1;(#ma1)) , (< ma2;2;9)
assert. 1&tdgehrd &> (< ma3;1;(#ma3)) , (< ma4;1;(#ma4)) , (< ma5;2;9)
assert. 2&tdgehrd &> (< ma0;1;(#ma0)) , (< ma1;1;(#ma1)) , (< ma2;2;9)
assert. 3&tdgehrd &> (< ma3;1;(#ma3)) , (< ma4;1;(#ma4)) , (< ma5;2;9)
EMPTY
)
