function exponential5_mp, xx, a
  return, a[0] * exp(a[1] * xx) + a[2] + a[3] * xx + a[4] * xx^2
end
