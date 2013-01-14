pro poisson, xx, a, f, pder
f = exp(-a[0]) * a[0]^xx / factorial(xx)
pder = [(xx/a[0] - 1) * f]
end
