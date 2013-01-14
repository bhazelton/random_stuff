function lorentzian_mp, xx, a
  bx = (xx - a[1])^2 + a[2]^2
  return, a[0]*a[2]^2/bx
end
