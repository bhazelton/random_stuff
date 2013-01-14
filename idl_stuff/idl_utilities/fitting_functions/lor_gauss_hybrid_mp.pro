function lor_gauss_hybrid_mp, xx, a
  bx = (xx - a[4])^2 + a[5]^2
  z=(xx-a[1])/a[2]
  return, a[0]*exp(-1*z^2/2) + a[3]*a[5]^2/bx + a[6]
end
