pro lorentzian6, xx, a, f, pder
bx = (xx - a[1])^2 + a[2]^2
f = a[0]*a[2]^2/bx + a[3] + a[4]*xx + a[5]*xx^2
pder=[[a[2]^2/bx], [2*a[0]*a[2]^2*(xx-a[1])/bx^2], $
      [2*a[0]*a[2]/bx(1-a[2]^2/bx)], [replicate(1.0, n_elements(xx))], [xx], $
     [xx^2]]
end
