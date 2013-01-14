pro exponential4, xx, a, f, pder
bx = exp(a[1] * xx) 
f = a[0] * bx + a[2] + a[3] * xx
pder = [[bx], [a[0] * xx * bx], [replicate(1.0, n_elements(xx))], [xx]]
end
