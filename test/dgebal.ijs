NB. subroutine dgebal ( character                              JOB,
NB.                     integer                                N,
NB.                     double precision, dimension( lda, * )  A,
NB.                     integer                                LDA,
NB.                     integer                                ILO,
NB.                     integer                                IHI,
NB.                     double precision, dimension( * )       SCALE,
NB.                     integer                                INFO
NB.                   )
NB.
NB. DGEBAL
NB.
NB. Purpose:
NB.
NB.      DGEBAL balances a general real matrix A.  This involves, first,
NB.      permuting A by a similarity transformation to isolate eigenvalues
NB.      in the first 1 to ILO-1 and last IHI+1 to N elements on the
NB.      diagonal; and second, applying a diagonal similarity transformation
NB.      to rows and columns ILO to IHI to make the rows and columns as
NB.      close in norm as possible.  Both steps are optional.
NB.
NB.      Balancing may reduce the 1-norm of the matrix, and improve the
NB.      accuracy of the computed eigenvalues and/or eigenvectors.
NB.
NB. Parameters
NB.
NB.                              JOB is CHARACTER*1
NB.                              Specifies the operations to be performed on A:
NB.                              = 'N':  none:  simply set ILO = 1, IHI = N, SCALE(I) = 1.0
NB.     [in]     JOB                     for i = 1,...,N;
NB.                              = 'P':  permute only;
NB.                              = 'S':  scale only;
NB.                              = 'B':  both permute and scale.
NB.
NB.                              N is INTEGER
NB.     [in]     N               The order of the matrix A.  N >= 0.
NB.
NB.                              A is DOUBLE PRECISION array, dimension (LDA,N)
NB.                              On entry, the input matrix A.
NB.     [in,out] A               On exit,  A is overwritten by the balanced matrix.
NB.                              If JOB = 'N', A is not referenced.
NB.                              See Further Details.
NB.
NB.                              LDA is INTEGER
NB.     [in]     LDA             The leading dimension of the array A.  LDA >= max(1,N).
NB.
NB.     [out]    ILO             ILO is INTEGER
NB.
NB.                              IHI is INTEGER
NB.                              ILO and IHI are set to integers such that on exit
NB.     [out]    IHI             A(i,j) = 0 if i > j and j = 1,...,ILO-1 or I = IHI+1,...,N.
NB.                              If JOB = 'N' or 'S', ILO = 1 and IHI = N.
NB.
NB.                              SCALE is DOUBLE PRECISION array, dimension (N)
NB.                              Details of the permutations and scaling factors applied to
NB.                              A.  If P(j) is the index of the row and column interchanged
NB.                              with row and column j and D(j) is the scaling factor
NB.                              applied to row and column j, then
NB.     [out]    SCALE           SCALE(j) = P(j)    for j = 1,...,ILO-1
NB.                                       = D(j)    for j = ILO,...,IHI
NB.                                       = P(j)    for j = IHI+1,...,N.
NB.                              The order in which the interchanges are made is N to IHI+1,
NB.                              then 1 to ILO-1.
NB.
NB.                              INFO is INTEGER
NB.     [out]    INFO            = 0:  successful exit.
NB.                              < 0:  if INFO = -i, the i-th argument had an illegal value.

require 'math/lapack2'
cocurrent 'base'
coinsert 'jlapack2'
match=: matchclean;;
matchf=: matchcleanf;;

NB. =========================================================
tdgebal=: 4 : 0
if. (3=x) *. 807>:0".}.({.~ i.&'/')9!:14'' do. 1 return. end.
zero=. (2|x){::dzero;zzero

a=. zero + y
'm n'=. $a

assert. 0= _1{::cdrc=. dgebal`zgebal`sgebal`cgebal@.x (,'B');(,m);(|:a);(,1>.m);(,0);(,0);(n$dzero);,_1
'AB ILO IHI SCALE'=. 3 5 6 7{cdrc
AB=. |:AB
ILO=. {. ILO
IHI=. {. IHI

echo AB;ILO;IHI;SCALE

p=. (ILO , IHI) makeper SCALE
d=. ((>:&ILO *. <:&IHI) #\ i. # y)} 1 ,: SCALE
echo r=. AB match`matchf@.(x>1) d (*"1 % [) p ([ C."1 C.) y  NB. compare AB with D^_1 * P * A * P^_1 * D
0{::r
)

NB. =========================================================
testdgebal=: 3 : 0
m0=. 0 0$0
m1=. 7 7 $ 6 0 0 0 0 1 0 0 4 0 0.00025 0.0125 0.02 0.125 1 128 64 0 0 _2 16 0 16384 0 1 _400 256 _4000 _2 _256 0 0.0125 2 2 32 0 0 0 0 0 0 0 0 8 0 0.004 0.125 _0.2 3
m2=. 0 0$zzero
m3=. 6 6 $ 1j1 1j1 0 1j1 1j1 1j1 1j1 1j1 0 1j1 1j1 1j1 1j1 1j1 1j1 1j1 1j1 1j1 0 0 0 1j1 0 0 1j1 1j1 0 1j1 1j1 1j1 1j1 1j1 0 1j1 1j1 1j1
assert. 0&tdgebal &> m0;m1
assert. 1&tdgebal &> m2;m3
assert. 2&tdgebal &> m0;m1
assert. 3&tdgebal &> m2;m3
EMPTY
)
