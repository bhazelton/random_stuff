


n_pts = 100
variances = 1/sin(findgen(n_pts)*30./n_pts.)^2.
variances[0]=1.

variances = fltarr(n_pts) + 1.
variances[indgen(n_pts/5.)*5] = 100.



n_trials = 100000
data = randomn(seed, n_pts, n_trials) * rebin(sqrt(variances), n_pts, n_trials,/sample)

data_ft = fft(data, dimension=1)

var_ft = total(abs(data_ft)^2., 2)

n_val = findgen(n_pts)-(n_pts/2. -1)
data_a0 = data_ft[where(n_val eq 0),*]
data_an = (data_ft[where(n_val gt 0),*] + data_ft[reverse(where(n_val lt 0)),*])/2.
data_bn = complex(0,1) * (data_ft[where(n_val gt 0),*] - data_ft[reverse(where(n_val lt 0)),*])/2.

var_cos = [total(abs(data_a0)^2.), total(abs(data_an)^2.,2)]/n_trials
var_sin = [0, total(abs(data_bn)^2.,2)]/n_trials

