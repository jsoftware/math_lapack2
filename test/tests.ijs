cocurrent 'base'

routine=: 0 : 0
dgebal
dgees
dgeev
dgehrd
dgelqf
dgels
dgelss
dgeqlf
dgeqrf
dgerqf
dgesvd
dgesv
dgesvx
dgetrf
dpotrf
dsyev
dtrtrs
zheev
)

test=: 3 : 0
for_r. <;._2 routine do.
  echo >r
  load '~addons/math/lapack2/test/',(>r),'.ijs'
  ('test',(>r),'_base_')~''
end.
''
)

test''
