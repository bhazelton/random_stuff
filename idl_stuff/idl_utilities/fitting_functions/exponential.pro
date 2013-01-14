pro exponential, xx, a, f, pder
bx = exp(a[1] * xx) 
f = a[0] * bx
pder = [[bx], [a[0] * xx * bx]]
end
