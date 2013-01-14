function base_path
  return, '/Users/bryna/Documents/Physics/bryna_svn/idl_working/'
end

PRO undefine, varname  
   tempvar = SIZE(TEMPORARY(varname))
END

pro upgrade
 ssw_upgrade,/hessi,/spex,/xray,/gen,/spawn,/loud
end

pro dhelp
print,''
print,'(fu_gauss,fu_lin,fu_quad,fu_poly,fu_polyoff)  x,a,f,pder'
print,'purgefit, x,y,w,a,sigmaa,frac,reps,function_name=function_name'
print,'datin,fname,ncol,d'
print,'set,x,xcol,d'
print,'setxy,x,y,xcol,ycol,d'
print,'cscale,x0,c'
print,'bin2,x,y,z,xd,yd,zd,xoff,yoff'
print,'bin2plot,xd,yd,zd,xoff,yoff,xtitle=xtitle,ytitle=ytitle,title=title'
print,''
print,'batboth,fname,num,col,pre,lf,lerr,bins,corr'
print,'batgraf,fname,num,col,pre,lf,lerr,corr'
print,'bathist,fname,num,col,pre,lf,lerr,bins,corr'
print,'batspec,fname,yr'
print,'batmap,fname'
print,'colorscatter,fname,cols,xcol,ycol,colorcol'
print,''
end

pro fu_gauss, x,a,f,pder
ex = exp(-(x-a(2))^2./(2.*a(1)^2))
f = a(0)*ex
if n_params(0) LE 3 then return
pder = dblarr(n_elements(x),n_elements(a))
pder(*,0) = ex
pder(*,1) = a(0)*ex*(x-a(2))^2./a(1)^3.
pder(*,2) = a(0)*ex*(x-a(2))/a(1)^2.
end


pro fu_lin, x,a,f,pder
f = a(1)*x + a(0)
if n_params(0) LE 3 then return
pder = dblarr(n_elements(x),n_elements(a))
pder(*,0) = 1.
pder(*,1) = x
end

pro fu_quad, x,a,f,pder
f = a(2)*x^2. + a(1)*x + a(0)
if n_params(0) LE 3 then return
pder = dblarr(n_elements(x),n_elements(a))
pder(*,0) = 1.
pder(*,1) = x
pder(*,2) = x^2.
end

pro fu_poly, x,a,f,pder
p=dindgen(n_elements(a))
f=dblarr(n_elements(x))
f(*)=0.
for i=0,n_elements(a)-1 do  f(*) = f(*) + a(i)*x(*)^p(i)
if n_params(0) LE 3 then return
pder = dblarr(n_elements(x),n_elements(a))
for i=0,n_elements(p)-1 do pder(*,i) = x(*)^p(i)
end

pro fu_polyoff, x,a,f,pder
p=dindgen(n_elements(a)-1)
f=dblarr(n_elements(x))
f(*)=0.
for i=1,n_elements(a)-1 do  f(*) = f(*) + a(i)*(x(*)-a(0))^p(i-1)
if n_params(0) LE 3 then return
pder = dblarr(n_elements(x),n_elements(a))
pder(*,*) = 0.
for i=2,n_elements(a)-1 do pder(*,0) = -p(i-1)*a(i)*(x(*)-a(0))^p(i-2)
for i=1,n_elements(a)-1 do pder(*,i) = (x(*)-a(0))^p(i-1)
end

pro purgefitrat, x,y,w,a,sigmaa,frac,reps,function_name=function_name
r = curvefit(x,y,w,a,sigmaa,function_name=function_name)
for i=1,reps do begin
  s = size(x)
  call_procedure, function_name,x,a,f
  for j=1,frac*(s(1)-1) do begin
    d = abs(y-f)/f
    x = x(where( d LT max(d)))
    y = y(where( d LT max(d)))
    w = w(where( d LT max(d)))
    f = f(where( d LT max(d)))
  end
  r = curvefit(x,y,w,a,sigmaa,function_name=function_name)
end
end

pro purgefitdist, x,y,w,a,sigmaa,frac,reps,function_name=function_name
r = curvefit(x,y,w,a,sigmaa,function_name=function_name)
for i=1,reps do begin
  s = size(x)
  call_procedure, function_name,x,a,f
  for j=1,frac*(s(1)-1) do begin
    d = abs(y-f)
    x = x(where( d LT max(d)))
    y = y(where( d LT max(d)))
    w = w(where( d LT max(d)))
    f = f(where( d LT max(d)))
  end
  r = curvefit(x,y,w,a,sigmaa,function_name=function_name)
end
end

pro purgefitn, x,y,w,a,sigmaa,n,function_name=function_name
r = curvefit(x,y,w,a,sigmaa,function_name=function_name)
call_procedure, function_name,x,a,f
for j=1,n do begin
    d = abs(y-f)
    x = x(where( d LT max(d)))
    y = y(where( d LT max(d)))
    w = w(where( d LT max(d)))
    f = f(where( d LT max(d)))
end
r = curvefit(x,y,w,a,sigmaa,function_name=function_name)
end

function datin,fname,ncol,d
openr,1,fname
s=''
first=-1l
last=-1l
i=0l
while not eof(1) and last EQ -1l do begin
  readf,1,s
  ss = strmid(strcompress(s,/remove_all),0,1)
  linetype=strpos('0123456789.-',ss)
  if strlen(ss) EQ 0 then linetype=-1
  if first GE  0l and linetype EQ -1 then last=i-1l
  if first EQ -1l and linetype GE  0l then first=i
  i=i+1
end
if last EQ -1l then last=i-1l
close,1
openr,1,fname
d=dblarr(ncol,last-first+1)
for i=0,first-1 do readf,1,s 
readf,1,d
close,1
return,last-first+1l
end

pro set,x,xcol,d
x = transpose(d(xcol,*))
end

pro setxy,x,y,xcol,ycol,d
x = transpose(d(xcol,*))
y = transpose(d(ycol,*))
end


pro cscale,x0,c
s0=size(x0)
s=size(x0)
x=x0
while s(1) GT s0(1)*.99 do begin
  x = x(where( x LT max(x)))
  s = size(x)
end
while s(1) GT s0(1)*.98 do begin
  x = x(where( x GT min(x)))
  s = size(x)
end
x1=min(x)
x2=max(x)
b=(!d.n_colors-20)/(x2-x1)
a= !d.n_colors-20 - b*x2
c=b*x0+a+4
end

pro cscale0,x0,c
x1=min(x0)
x2=max(x0)
b=(!d.n_colors-20)/(x2-x1)
a= !d.n_colors-20 - b*x2
c=b*x0+a+4
end

pro cscale2,x0,c
s0=size(x0)
s=size(x0)
x=x0
while s(1) GT s0(1)*.99 do begin
  x = x(where( x LT max(x)))
  s = size(x)
end
while s(1) GT s0(1)*.98 do begin
  x = x(where( x GT min(x)))
  s = size(x)
end
x1=0.
x2=max(x)
b=(!d.n_colors-20)/(x2-x1)
a= !d.n_colors-20 - b*x2
c=b*x0+a+4
end

pro cscale3,x0,c
x1=0.
x2=max(x0)
b=(!d.n_colors-20)/(x2-x1)
a= !d.n_colors-20 - b*x2
c= (!d.n_colors-20) - (b*x0+a-4)
end

pro powerlaw,x,a,f,pder
f=a[0]*x^a[1]
pder=[ [f/a[0]], [f*alog(x)]]
end


pro bin2plot,xd,yd,zd,xoff,yoff,xtitle=xtitle,ytitle=ytitle,title=title
s = size(zd)
sz = max(s)
sz1 = s(1)
x=[xoff-xd/10.,xoff+xd*(sz1+.1)]
y=[yoff-yd/10.,yoff+yd*(sz/sz1+.1)]
if not (keyword_set(xtitle)) then xtitle =''
if not (keyword_set(ytitle)) then ytitle =''
if not (keyword_set(title)) then title =''
plot,x,y,psym=3,xtitle=xtitle,ytitle=ytitle,title=title
cv = reform(zd,sz)
cscale,cv,c
cl = reform(c,sz1,sz/sz1)
for i=0,sz1-1 do begin
 for j=0,sz/sz1-1 do begin
   x=[xoff+xd*float(i),xoff+xd*float(i),xoff+xd*float(i+1),xoff+xd*float(i+1)]
   y=[yoff+yd*float(j),yoff+yd*float(j+1),yoff+yd*float(j+1),yoff+yd*float(j)]
   polyfill,x,y,color=cl(i,j)
 end
end
end
