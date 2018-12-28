NB. build.ijs

mkdir_j_ jpath '~Addons/math/lapack2/test'
mkdir_j_ jpath '~Addons/math/lapack2/lib'
mkdir_j_ jpath '~addons/math/lapack2/test'
mkdir_j_ jpath '~addons/math/lapack2/lib'

writesourcex_jp_ '~Addons/math/lapack2/source';'~Addons/math/lapack2/lapack2.ijs'
('checklibrary$0',LF,'cocurrent ''base''',LF) fappend '~Addons/math/lapack2/lapack2.ijs'

(jpath '~addons/math/lapack2/lapack2.ijs') (fcopynew ::0:) jpath '~Addons/math/lapack2/lapack2.ijs'

f=. 3 : 0
(jpath '~Addons/math/lapack2/',y) fcopynew jpath '~Addons/math/lapack2/',y
(jpath '~addons/math/lapack2/',y) (fcopynew ::0:) jpath '~Addons/math/lapack2/',y
)

f 'manifest.ijs'
f 'lib/readme.txt'
f 'test/dgebal.ijs'
f 'test/dgees.ijs'
f 'test/dgeev.ijs'
f 'test/dgehrd.ijs'
f 'test/dgelqf.ijs'
f 'test/dgels.ijs'
f 'test/dgelss.ijs'
f 'test/dgeqlf.ijs'
f 'test/dgeqrf.ijs'
f 'test/dgerqf.ijs'
f 'test/dgesvd.ijs'
f 'test/dgesv.ijs'
f 'test/dgesvx.ijs'
f 'test/dgetrf.ijs'
f 'test/dpotrf.ijs'
f 'test/dsyev.ijs'
f 'test/dtrtrs.ijs'
f 'test/tests.ijs'
f 'test/zheev.ijs'
