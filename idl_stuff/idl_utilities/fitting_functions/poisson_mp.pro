function poisson, xx, a
  return, exp(-a[0]) * a[0]^xx / factorial(xx)
end
