pro prob,seed,n

below=0l

x = .22 * .0034 * .28 * .23d ;21n
;x = .33 * .084 * .79 * .74 ;15n
;x = .33 * .12 *.83 *.79 ; 13n
;x = .43 * .11 * 1 * .78 ; 21w
;x = .81 * .19 * 1 * 1 ;15w
;x = .74 * .28 * 1 * 1 ;13w



for i=0l,1000000l-1l do begin
a=randomu(seed,n)
pr=a[0]
for j=1,n-1 do pr=pr*a[j]
if pr LT x then below=below+1
end
print,below/1000000.d

end
