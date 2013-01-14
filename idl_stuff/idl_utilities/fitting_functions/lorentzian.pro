pro lorentzian, xx, a, f, pder
bx = (xx - a[1])^2 + a[2]^2
f = a[0]*a[2]^2/bx
pder=[[a[2]^2/bx], [2*a[0]*a[2]^2*(xx-a[1])/bx^2], $
      [2*a[0]*a[2]/bx(1-a[2]^2/bx)]]
end
