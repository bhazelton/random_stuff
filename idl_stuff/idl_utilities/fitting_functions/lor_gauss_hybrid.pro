pro lor_gauss_hybrid, xx, a, f, pder
bx = (xx - a[4])^2 + a[5]^2
z=(xx-a[1])/a[2]
f = a[0]*exp(-1*z^2/2) + a[3]*a[5]^2/bx + a[6]
pder=[[exp(-1*z^2/2)], [a[0]*z/a[2]*exp(-1*z^2/2)], $
      [a[0]*z^2/a[2]*exp(-1*z^2/2)], [a[5]^2/bx], $
      [2*a[3]*a[5]^2*(xx-a[4])/bx^2], [2*a[3]*a[5]/bx(1-a[5]^2/bx)], $
      [replicate(1.0, n_elements(xx))]]
end
