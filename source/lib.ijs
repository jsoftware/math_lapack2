NB. lib version

NB. =========================================================
NB. library:
NB. First part of this script is making sure the library is loaded
3 : 0''
if. 0=4!:0<'liblapack' do. '' return. end.
if. UNAME-:'Linux' do.
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
elseif. UNAME-:'Android' do.
  arch=. LF-.~ 2!:0'getprop ro.product.cpu.abi'
  if. IF64 < arch-:'arm64-v8a' do.
    arch=. 'armeabi-v7a'
  elseif. IF64 < arch-:'x86_64' do.
    arch=. 'x86'
  end.
  liblapack=: (jpath'~bin/../libexec/',arch,'/liblapack.so')
  if. -.fexist liblapack do.
    liblapack=: (({.~ i:&'/') LIBFILE),'/liblapack.so'
  end.
elseif. do.
  liblapack=: jpath '~addons/math/lapack2/lib/liblapack3',((-.IF64)#'_32'),'.dll'
end.
)

NB. =========================================================
checklibrary=: 3 : 0
if. ((dquote liblapack) ,' dummyfunction n')&cd :: (1={.@cder) '' do.
  if. (UNAME-:'Linux') do.
    sminfo 'The binary needed for math/lapack2 has not yet been installed.',LF2,'Install liblapack3 (or similar) package from linux distro'
  else.
    getbinmsg 'The binary needed for math/lapack2 has not yet been installed.',LF2,'To install, ' return.
  end.
end.
)

NB. =========================================================
NB. get lapack2 binary
NB. uses routines from pacman
getbin=: 3 : 0
if. +./ (UNAME-:'Darwin'),(UNAME-:'Linux'),(UNAME-:'Android') do. return. end.
require 'pacman'
path=. 'http://www.jsoftware.com/download/lapackbin/'
arg=. HTTPCMD_jpacman_
tm=. TIMEOUT_jpacman_
dq=. dquote_jpacman_ f.
to=. liblapack_jlapack2_
if. UNAME-:'Android' do.
  path=. 'http://www.jsoftware.com/download/'
  arch=. LF-.~ 2!:0'getprop ro.product.cpu.abi'
  if. IF64 < arch-:'arm64-v8a' do.
    arch=. 'armeabi-v7a'
  elseif. IF64 < arch-:'x86_64' do.
    arch=. 'x86'
  end.
  fm=. path,'android/libs/',z=. arch,'/liblapack.so'
  'res p'=. httpget_jpacman_ fm
  if. res do.
    smoutput 'Connection failed: ',z return.
  end.
  (<to) 1!:2~ 1!:1 <p
  2!:0 ::0: 'chmod 644 ', dquote to
  1!:55 ::0: <p
  smoutput 'LAPACK binary installed.'
  return.
end.
fm=. path,1 pick fpathname to
lg=. jpath '~temp/getbin.log'
cmd=. arg rplc '%O';(dquote to);'%L';(dquote lg);'%t';'3';'%T';(":tm);'%U';fm
res=. ''
fail=. 0
try.
  fail=. _1-: res=. shellcmd cmd
  2!:0 ::0:^:(UNAME-:'Linux') 'chmod 644 ', dquote to
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
