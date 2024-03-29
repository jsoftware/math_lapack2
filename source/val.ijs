NB. lapack validation
NB.
NB. validation routines that check argument is a matrix:
NB.    vmatrix
NB.    vhermitian
NB.    vorthogonal
NB.    vsquare
NB.    vsymposdef
NB.
NB. validation routines that check argument is either
NB. a matrix or vector:
NB.    vmatrixorvector

isscalar=: 0 = #@$
isvector=: 1 = #@$
ismatrix=: 2 = #@$
ismatrixorvector=: 1 2 e.~ #@$
issquare=: =/ @ $
issymmetric=: -: |:
ishermitian=: -: +@|:

NB. =========================================================
NB. ishessenberg   dyadic adv. to check for lower/upper Hessenberg matrix
NB. syntax:
NB.   isH=. (ilo , ihi) (mkHmat ishessenberg) H
NB. where mkHmat is either lhmat or uhmat
NB. e.g.
NB.   isH=. 2 4 uhmat_jlapack2_ ishessenberg_jlapack2_ H

ishessenberg=: 1 : (':' ; '(-: (x u $)`(0&,:)}) y')

NB. =========================================================
isorthogonal=: 3 : 0
q=. y mp |: y
*./ 0 = clean ,q - idmat $q
)

NB. =========================================================
isunitary=: 3 : 0
q=. y mp + |: y
*./ 0 = clean ,q - idmat $q
)

NB. =========================================================
NB. f=. 2 : 'm&(13!:8)@(#&12)@(0 e. v)'
NB.
NB. vmatrix=: 'argument should be a matrix' f ismatrix
NB. vmatrixorvector=: 'argument should be either a matrix or vector' f ismatrixorvector
NB. vhermitian=: 'argument should be a hermitian matrix' f ishermitian [ vmatrix
NB. vorthogonal=: 'argument should be an orthogonal matrix' f isorthogonal [ vmatrix
NB. vsquare=: 'argument should be a square matrix' f issquare [ vmatrix
