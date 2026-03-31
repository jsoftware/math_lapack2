NB. lib version

NB. =========================================================
NB. library:
NB. First part of this script is making sure the library is loaded
3 : 0''
if. 0=4!:0<'liblapack' do. '' return. end.
if. 0: ~: 4!:0 @ <'IFWA64' do. IFWA64=. 0 end.
if. (<UNAME)e.'Linux';'FreeBSD';'OpenBSD' do.
  liblapack=: 'liblapack.so.3'
elseif. UNAME-:'Darwin' do.
  if. IFIOS do.
    liblapack=: '/System/Library/Frameworks/Accelerate.framework/Frameworks/vecLib.framework/vecLib'
  else.
    liblapack=: '/System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/vecLib'
    if. -.fexist liblapack do.
      liblapack=: '/System/Library/Frameworks/vecLib.framework/vecLib'
    end.
  end.
elseif. IFWIN do. NB. Win
  liblapack=: jpath '~addons/math/lapack2/lib/libopenblas',((-.IFWA64){::'_arm64';(-.IF64)#'_32'),'.dll'
  if. -.fexist liblapack do.
    liblapack=: 'libopenblas',((-.IFWA64){::'_arm64';(-.IF64)#'_32'),'.dll'
  end.
elseif. do.
  liblapack=: 'liblapack.so'
end.
)

NB. =========================================================
checklibrary=: 3 : 0
if. ((dquote liblapack) ,' dummyfunction n')&cd :: (1={.@cder) '' do.
  if. (<UNAME)e.'Linux';'FreeBSD';'OpenBSD' do.
    sminfo 'The binary needed for math/lapack2 has not yet been installed.',LF2,'Install libopenblas0 (or similar) package from linux distro'
  else.
    getbinmsg 'The binary needed for math/lapack2 has not yet been installed.',LF2,'To install, ' return.
  end.
end.
)

NB. =========================================================
NB. get lapack2 binary
NB. uses routines from pacman
getbin=: 3 : 0
if. -.IFWIN do. return. end.
require 'pacman'
path=. 'http://www.jsoftware.com/download/lapackbin/'
arg=. HTTPCMD_jpacman_
tm=. TIMEOUT_jpacman_
dq=. dquote_jpacman_ f.
to=. liblapack_jlapack2_
fm=. path,1 pick fpathname to
lg=. jpath '~temp/getbin.log'
cmd=. arg rplc '%O';(dquote to);'%L';(dquote lg);'%t';'3';'%T';(":tm);'%U';fm
res=. ''
fail=. 0
try.
  fail=. _1-: res=. shellcmd cmd
catch. fail=. 1 end.
if. fail +. 0 >: fsize to do.
  if. _1-:msg=. freads lg do.
    if. (_1-:msg) +. 0=#msg=. res do. msg=. 'Unexpected error' end. end.
  ferase to,lg
  smoutput 'Connection failed: ',msg
else.
  ferase lg
  smoutput 'LAPACK binary installed.'
end.
)

NB. =========================================================
getbinmsg=: 3 : 0
msg=. y,' run the getbin_jlapack2_'''' line written to the session.'
smoutput '   getbin_jlapack2_'''''
sminfo 'LAPACK';msg
)

NB. =========================================================
shellcmd=: 3 : 0
if. IFUNIX do.
  hostcmd_j_ y
else.
  spawn_jtask_ y
end.
)
