pro poisson_errors, array, limits, symmetric = symmetric, feldcous = feldcous

num = n_elements(array)
if num eq 0 then message, 'array must have at least one element'
if num eq 1 then begin
   dimensions = size(array, /dimensions)
   if dimensions[0] eq 0 then array = [array]
endif

type = size(array, /type)

up_lim = array * 0.0d
low_lim = array * 0.0d

if keyword_set(feldcous) then begin
   ;; use limits from Feldman and Cousins for first 20, then switch to Geherls up to 50
   psn_upper = [1.29, 2.75, 4.25, 5.30, 6.78, 7.81, 9.28, 10.30, 11.32, 12.79, 13.81, 14.82, 16.29, 17.30, 18.32, 19.32, 20.80, $
                21.81, 22.82, 23.82, 25.30]

   psn_lower = [0.00, 0.37, 0.74, 1.10, 2.34, 2.75, 3.82, 4.25, 5.30, 6.33, 6.78, 7.81, 8.83, 9.28, 10.30, 11.32, 12.33, 12.79, $
                13.81, 14.82, 15.83]

endif else begin
   ;; use Gehrels numbers
   psn_upper = [1.814, 3.300, 4.638, 5.918, 7.163, 8.382, 9.584, 10.77, 11.95, 13.11, 14.27, 15.42, 16.56, 17.70, 18.83, 19.96, $
                21.08, 22.20, 23.32, 24.44, 25.55] 

   psn_lower = [0.000, 0.173, 0.708, 1.367, 2.086, 2.840, 3.620, 4.419, 5.232, 6.057, 6.891, 7.734, 8.585, 9.441, 10.30, 11.17, $
                12.04, 12.92, 13.80, 14.68, 15.57]
endelse
;; Use Geherls numbers for 21-50
psn_upper = [psn_upper, 26.66, 27.76, 28.87, 29.97, 31.07, 32.16, 33.26, 34.35, 35.45, 36.54, 37.63, 38.72, 39.80, 40.89, 41.97, $
             43.06, 44.14, 45.22, 46.30, 47.38, 48.46, 49.53, 50.61, 51.68, 52.76, 53.83, 54.90, 55.98, 57.05, 58.12]

psn_lower = [psn_lower, 16.45, 17.35, 18.24, 19.14, 20.03, 20.93, 21.84, 22.74, 23.65, 24.55, 25.46, 26.37, 27.28, 28.20, 29.11, $
             30.03, 30.94, 31.86, 32.78, 33.70, 34.62, 35.55, 36.47, 37.39, 38.32, 39.24, 40.17, 41.10, 42.02, 42.95]

if type eq 2 or type eq 3 then begin
   ;; make a mask for numbers le 50 in the array
   binsize = 51
   hist = histogram(array, binsize = binsize, min = 0, reverse_indices = ri)
   if hist[0] ne 0 then begin
      mask = ri[ri[0]:ri[1] - 1]

      ;;  set upper/lower errors
      up_lim[mask] = psn_upper[array[mask]]
      low_lim[mask] = psn_lower[array[mask]]
   endif

   ;; now deal with numbers > 50 using Geherls eq. 9 & 12 (same as 14 with beta = 0)
   if n_elements(hist) gt 1 then begin
      inv_mask = ri[ri[1]:ri[n_elements(hist)] - 1]
      up_lim[inv_mask] = (array[inv_mask] + 1.) * (1. - 1. / (9. * (array[inv_mask] + 1.)) + $
                                                   1. / (3. * sqrt(array[inv_mask] + 1.))) ^ 3.
      low_lim[inv_mask] = array[inv_mask] * (1. - 1. / (9. * array[inv_mask]) - 1. / (3. * sqrt(array[inv_mask]))) ^ 3
   endif
endif else begin
   ;; non-integer numbers, use Geherls eq. 9 & 12 (same as 14 with beta = 0)
   up_lim = (array + 1.) * (1. - 1. / (9. * (array + 1.)) + 1. / (3. * sqrt(array + 1.))) ^ 3.
   low_lim = array * (1. - 1. / (9. * array) - 1. / (3. * sqrt(array))) ^ 3
   ;; limit to 0 at low end
   if min(low_lim) lt 0 then low_lim[where(low_lim lt 0)] = 0d
endelse

limits = {upper:up_lim, lower: low_lim}

if arg_present(symmetric) then begin
   low_val = array - low_lim
   up_val = up_lim - array

   symmetric = (up_val + low_val) / 2
endif

end
