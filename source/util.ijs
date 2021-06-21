NB. lapack utils
NB.
NB. z2d           convert complex to float datatype if zero imagine part
NB.
NB. matchclean    if clean x-y is all 0
NB.
NB. diagmat       rectangular diagonal matrix
NB. idmat         rectangular identity matrix with shifted diagonal
NB. ltmat         lower triangular (trapezoidal) matrix
NB. utmat         upper triangular (trapezoidal) matrix
NB.
NB. ltri          return only lower triangular (trapezoidal) matrix
NB. utri          return only upper triangular (trapezoidal) matrix
NB. sltri         return only strictly lower triangular (trapezoidal) matrix
NB. sutri         return only strictly upper triangular (trapezoidal) matrix
NB.
NB. cxpair        reconstruct complex columns
NB. xtoken        exclude tokens with indices in x from list y
NB. invperm       inverse permutation of x by pivot indices
NB.               from y
NB. makepermat    generate inverse permutation matrix P from
NB.               pivot indices y

izero=: 23-23
ione=: 23-22
dzero=: 1.1-1.1
done=: 2.1-1.1
zzero=: 1j1-1j1
zone=: 2j1-1j1

z2d=: [ ^: (-. @ -:) (9 & o.)

mp=: +/ . *

NB. from general/numeric
NB. =========================================================
NB.*clean v clean y to tolerance of x (default 1e_10)
NB. form: tolerance (default 1e_10) clean numbers
NB. sets values less than tolerance to 0
clean=: 1e_10&$: : (4 : 0)
if. L. y do.
  x clean each y
else.
  if. (3!:0 y) e. 16 16384 do.
    j./"1 y * (2*x) <: | y=. +.y
  else.
    y * x <: |y
  end.
end.
)

NB. =========================================================
matchclean=: 0: *./ . = clean @ , @: -

NB. single precision
neareq=: (2^_16) > [ |@:% -
matchcleanf=: 0: *./ . neareq 5e_5&clean @ , @: -
cleanf=: 5e_5&clean

NB. =========================================================
NB. diagmat   rectangular diagonal matrix with y on diagonal
NB. x=rows-columns , x=0 is default
NB. e.g.
NB.    diagmat 3 5 7
NB. 3 0 0
NB. 0 5 0
NB. 0 0 7
NB.    1 diagmat 3 5 7
NB. 3 0 0
NB. 0 5 0
NB. 0 0 7
NB. 0 0 0
NB.    _1 diagmat 3 5 7
NB. 3 0 0 0
NB. 0 5 0 0
NB. 0 0 7 0

diagmat=: (0 $: ]) :(((0 (>. , -@<.) [) + #@]) {. (* =@i.@#)@])

NB. =========================================================
NB. idmat   rectangular identity matrix with shifted diagonal
NB. e.g.
NB.    idmat 3
NB. 1 0 0
NB. 0 1 0
NB. 0 0 1
NB.    idmat 3 4
NB. 1 0 0 0
NB. 0 1 0 0
NB. 0 0 1 0
NB.    1 idmat 3 4
NB. 0 1 0 0
NB. 0 0 1 0
NB. 0 0 0 1

idmat=: (0 $: ]) :(= ({. -~/&i. {:))

NB. =========================================================
NB. ltmat   lower triangular (trapezoidal) boolean matrix
NB. e.g.
NB.    ltmat 3
NB. 1 0 0
NB. 1 1 0
NB. 1 1 1
NB.    ltmat 3 5
NB. 1 0 0 0 0
NB. 1 1 0 0 0
NB. 1 1 1 0 0
NB.    1 ltmat 3 5
NB. 1 1 0 0 0
NB. 1 1 1 0 0
NB. 1 1 1 1 0

ltmat=: (0 $: ]) :(>: ({. -~/&i. {:))

NB. =========================================================
NB. utmat   upper triangular (trapezoidal) boolean matrix
NB. e.g.
NB.    utmat 3
NB. 1 1 1
NB. 0 1 1
NB. 0 0 1
NB.    utmat 3 5
NB. 1 1 1 1 1
NB. 0 1 1 1 1
NB. 0 0 1 1 1
NB.    1 utmat 3 5
NB. 0 1 1 1 1
NB. 0 0 1 1 1
NB. 0 0 0 1 1

utmat=: (0 $: ]) :(<: ({. -~/&i. {:))

NB. =========================================================
NB. lhmat   lower Hessenberg boolean matrix
NB. H=. (ilo , ihi) lhmat size
NB. e.g.
NB.    2 4 lhmat 5
NB. 1 0 0 0 0
NB. 1 1 1 0 0
NB. 1 1 1 1 0
NB. 1 1 1 1 0
NB. 1 1 1 1 1
NB.    2 4 lhmat 5 5
NB. 1 0 0 0 0
NB. 1 1 1 0 0
NB. 1 1 1 1 0
NB. 1 1 1 1 0
NB. 1 1 1 1 1

lhmat=: 4 : '({. (>: +. (_1 = -) *. (ilo <: ]) *. (ihi >: >:@]))"0/&i. {:) y [ ''ilo ihi''=. x'

NB. =========================================================
NB. uhmat   upper Hessenberg boolean matrix
NB. H=. (ilo , ihi) uhmat size
NB. e.g.
NB.    2 4 uhmat 5
NB. 1 1 1 1 1
NB. 0 1 1 1 1
NB. 0 1 1 1 1
NB. 0 0 1 1 1
NB. 0 0 0 0 1
NB.    2 4 uhmat 5 5
NB. 1 1 1 1 1
NB. 0 1 1 1 1
NB. 0 1 1 1 1
NB. 0 0 1 1 1
NB. 0 0 0 0 1

uhmat=: 4 : '({. (<: +. ( 1 = -) *. (ilo <: [) *. (ihi >: >:@[))"0/&i. {:) y [ ''ilo ihi''=. x'

NB. =========================================================
NB. ltri   return only lower triangular (trapezoidal) matrix
NB. e.g.
NB.   ltri 3 5 $ 2
NB. 2 0 0 0 0
NB. 2 2 0 0 0
NB. 2 2 2 0 0
NB.    1 ltri 3 5 $ 2
NB. 2 2 0 0 0
NB. 2 2 2 0 0
NB. 2 2 2 2 0

ltri=: (0 $: ]) : (] * (>: ({. -~/&i. {:)@$))

NB. =========================================================
NB. utri   return only upper triangular (trapezoidal) matrix
NB. e.g.
NB.    utri 3 5 $ 2
NB. 2 2 2 2 2
NB. 0 2 2 2 2
NB. 0 0 2 2 2
NB.    1 utri 3 5 $ 2
NB. 0 2 2 2 2
NB. 0 0 2 2 2
NB. 0 0 0 2 2

utri=: (0 $: ]) : (] * (<: ({. -~/&i. {:)@$))

NB. =========================================================
NB. sltri   return only strictly lower triangular (trapezoidal) matrix
NB. e.g.
NB.    sltri 3 5 $ 2
NB. 0 0 0 0 0
NB. 2 0 0 0 0
NB. 2 2 0 0 0
NB.    1 sltri 3 5 $ 2
NB. 2 0 0 0 0
NB. 2 2 0 0 0
NB. 2 2 2 0 0

sltri=: (0 $: ]) : (] * (> ({. -~/&i. {:)@$))

NB. =========================================================
NB. sutri   return only strictly upper triangular (trapezoidal) matrix
NB. e.g.
NB.    sutri 3 5 $ 2
NB. 0 2 2 2 2
NB. 0 0 2 2 2
NB. 0 0 0 2 2
NB.    1 sutri 3 5 $ 2
NB. 0 0 2 2 2
NB. 0 0 0 2 2
NB. 0 0 0 0 2

sutri=: (0 $: ]) : (] * (< ({. -~/&i. {:)@$))

NB. =========================================================
NB. cxpair - reconstruct complex columns

cxpair=: 4 : 0
'i j'=: |: _2 [\ x
r=. i {"1 y
z=. j. j {"1 y
n=. (r+z) ,. r-z
n (i,j)}"1 y
)

NB. =========================================================
NB. pivot indices utilities
NB.
NB. Pivot indices come from xGEBAL, xGGBAL, xGEEVX, xLASYF,
NB. xGESV, xGESVX, xGETRF, xGETRI and xGETRS subroutines.
NB. This are 1-based indices of rows and/or columns
NB. interchanged in matrix.
NB. Pivot indices came from xGEBAL, xGGBAL and xGEEVX
NB. subroutines need to be separated from scaling factors
NB. and reordered.

NB. transform pivot indices to standard cycle representation of the permutation
ipiv2scrp=: ((}:^:({. -: {:))&.>)@:(<"1)@((i.@# ,. <:) : (((0 (1 i.@- {) [) ([ ,. <:@{) ]) , (,.~@(] + i.@-) <:)~/@[ , (1 { [) ((([ + i.@-~) #) ,. (- #) <:@{. ]) ]))

NB. do inverse permutation of x by pivot indices from y
invperm=: C.~ ipiv2scrp

NB. ---------------------------------------------------------
NB. Generate permutation from pivot indices.
NB.
NB. Syntax:
NB.   p=.             makeper    ipiv
NB.   p=. (ilo , ihi) makeper    scale
NB.   P=.             makepermat ipiv
NB.   P=. (ilo , ihi) makepermat scale
NB. where
NB.   ipiv  - n-vector, integers in range [1,n], incremented
NB.           pivot indices
NB.   ilo   > 0, integer, incremented IO first column of
NB.           scaled part of balanced matrix
NB.   ihi   ≥ ilo, integer, incremented IO last column of
NB.           scaled part of balanced matrix
NB.   scale - n-vector, float:
NB.              scale(j) = ipiv(j)    for j = 0    ,...,ilo-2
NB.                       = d(j)       for j = ilo-1,...,ihi-1
NB.                       = ipiv(j)    for j = ihi  ,...,n-1,
NB.           the order in which the interchanges are made is
NB.           n-1 to ihi, then 0 to ilo-2
NB.   p     - n-vector, non-negative integers, permutation
NB.   P     - n-by-n-matrix, non-negative integers,
NB.           permutation matrix
NB.
NB. Storage layout:
NB. - dyadic case forms the following pairs:
NB.     part    length       IO         localIO    0-based pair (j , scale(j)-1)
NB.     ----    ---------    -------    -------    -------------------------
NB.     T1      ilo-1        0          0          ilo-2 , scale(ilo-2)-1
NB.                          1          1          ilo-3 , scale(ilo-3)-1
NB.                          ...        ...        ...
NB.                          ilo-2      ilo-2      0     , scale(0)-1
NB.     D       ihi-ilo+1    ilo-1      0          ilo-1 , ilo-1
NB.                          ilo        1          ilo   , ilo
NB.                          ...        ...        ...
NB.                          ihi-1      ihi-ilo    ihi-1 , ihi-1
NB.     T2      n-ihi        ihi        0          ihi   , scale(ihi)-1
NB.                          ihi+1      1          ihi+1 , scale(ihi+1)-1
NB.                          ...        ...        ...
NB.                          n-1        n-ihi-1    n-1   , scale(n-1)-1

makeper=: C.@ipiv2scrp
makepermat=: ({ =)@makeper
