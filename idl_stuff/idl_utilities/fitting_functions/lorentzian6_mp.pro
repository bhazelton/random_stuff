function lorentzian6_mp, xx, a
  bx = (xx - a[1])^2 + a[2]^2
  return, a[0]*a[2]^2/bx + a[3] + a[4]*xx + a[5]*xx^2
end
