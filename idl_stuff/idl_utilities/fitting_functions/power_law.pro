pro power_law, xx, a, f, pder
f = a[0] * xx ^ a[1]
pder = [[xx ^ a[1]], [a[0] * alog(xx) * xx ^ a[1]]]
end
