function ellipse_point, params

  n_pts = 100
  thetas = findgen(n_pts) * 2.*!pi / (n_pts-1)

  points = fltarr(2, n_pts)
  points[0,*] = params[2] + params[0]*cos(thetas)*cos(params[4]) - params[1]*sin(thetas)*sin(params[4])
  points[1,*] = params[3] + params[0]*cos(thetas)*sin(params[4]) + params[1]*sin(thetas)*cos(params[4])

  return, points
end

function inside_ellipse, x, y, params, count_inside=count_inside, count_outside=count_outside

  ;; ellipse equation: ax^2 + 2bxy + cy^2 + 2dx + 2fy + g = 0
  a = params[0]^2*sin(params[4])^2. + params[1]^2*cos(params[4])^2.
  b = 2*(params[1]^2-params[0]^2)*sin(params[4])*cos(params[4])
  c = params[0]^2*cos(params[4])^2. + params[1]^2*sin(params[4])^2.
  d = -2*a*params[2] - b*params[3]
  f = -b*params[2] - 2*c*params[3]
  g = a*params[2]^2. + b*params[2]*params[3] + c*params[3]^2. - params[0]^2*params[1]^2
  temp = a*x^2 + b*x*y + c*y^2 + d*x + f*y
  wh_inside = where(temp le -1*g, count_inside, ncomplement = count_outside)

  mask = intarr(n_elements(x))
  mask[wh_inside] = 1
  return, mask
end

function ellipse_cost, params
  common ellipse_xy, x, y

  area = !pi * params[0] * params[1]

  wh_inside = inside_ellipse(x,y,params, count_inside=count_inside, count_outside=count_outside)

  box_area = (max(x)-min(x))^2 + (max(y)-min(y))^2
  cost =  area*n_elements(x)*.2/box_area + count_outside
  ;e_data = ellipse_point(params)
  ;cgplot, x, y, psym=4, aspect = aspect
  ;cgplot, e_data[0,*], e_data[1,*],/over, psym=-3, color = 'red'

  return, cost
end

function fit_ellipse, xvals, yvals, start_params = start_params
  ;; params: [axis1, axis2, x0, y0, rot_angle]
  common ellipse_xy, x, y
  x = xvals
  y = yvals

  if n_elements(start_params) eq 0 then begin
    xlen = (max(x)-min(x))/2
    ylen = (max(y)-min(y))/2
    xmid = min(x)+xlen
    ymid = min(y)+ylen
    ;start_params = [ylen^2., 0, xlen^2., -2*ylen^2*xmid, -2*xlen^2.*ymid, ylen^2*xmid^2 + xlen^2*ymid^2 - xlen^2.*ylen^2. ]
    start_params = [xlen, ylen, xmid, ymid, 0]
  endif


  params = amoeba(1.0e-5, SCALE=[start_params[0:3],!pi] , p0 = start_params, function_name = 'ellipse_cost', function_value=fval)
  ; Check for convergence:
  IF n_elements(params) EQ 1 THEN message, 'amoeba failed to converge'

  return, params
end

pro crelox_test, file_name

  textfast, data, header, file_path = base_path('data') + 'crelox/'+file_name, first_line=1, column_list = indgen(3), /read

  data_dims = size(data, /dimension)
  if n_elements(data_dims) eq 2 then n_pts = data_dims[1]

  xrange = minmax(data[0,*])
  yrange = minmax(data[1,*])
  zrange = minmax(data[2,*])
  color_use = round((data[2,*] - zrange[0])*255/(zrange[1]-zrange[0]))

  aspect = (yrange[1]-yrange[0])/(xrange[1] - xrange[0])
  window, 1
  cgplot, data[0,*], data[1,*], psym=4, aspect = aspect
  ;cgplots, data[0,*], data[1,*], psym=4, color = color_use

  params = fit_ellipse(data[0,*], data[1,*], start_params = start_params)

  e_data_start = ellipse_point(start_params)
  e_data = ellipse_point(params)

  cgplot, e_data_start[0,*], e_data_start[1,*],/over, psym=-3, color = 'black'
  cgplot, e_data[0,*], e_data[1,*],/over, psym=-3, color = 'red'

  n_rand = n_pts*100
  rand_locs = randomu(seed, [2, n_rand]) * rebin([max(e_data[0,*])-min(e_data[0,*]), max(e_data[1,*])-min(e_data[1,*])], 2, n_rand) + $
    rebin([min(e_data[0,*]), min(e_data[1,*])], 2, n_rand)

  wh_ellipse = inside_ellipse(rand_locs[0,*],rand_locs[1,*],params, count_inside=count_inside, count_outside=count_outside)
  rand_locs = rand_locs[*,where(wh_ellipse)]
  n_rand = count_inside

  window,2
  cgplot, data[0,*], data[1,*], psym=4, aspect = aspect
  cgplot, rand_locs[0,*], rand_locs[1,*], psym=3, color='red',/over


  d=(rebin(transpose(data[0,*]),n_pts,n_pts,/SAMPLE)-rebin(data[0,*],n_pts,n_pts,/SAMPLE))^2 + $
    (rebin(transpose(data[1,*]),n_pts,n_pts,/SAMPLE)-rebin(data[1,*],n_pts,n_pts,/SAMPLE))^2 ;+ $
  ;(rebin(transpose(data[2,*]),n_pts,n_pts,/SAMPLE)-rebin(data[2,*],n_pts,n_pts,/SAMPLE))^2
  distances = sqrt(temporary(d))

  d=(rebin(transpose(rand_locs[0,*]),n_rand,n_rand,/SAMPLE)-rebin(rand_locs[0,*],n_rand,n_rand,/SAMPLE))^2 + $
    (rebin(transpose(rand_locs[1,*]),n_rand,n_rand,/SAMPLE)-rebin(rand_locs[1,*],n_rand,n_rand,/SAMPLE))^2 ;+ $
  ;(rebin(transpose(rand_locs[2,*]),n_rand,n_rand,/SAMPLE)-rebin(rand_locs[2,*],n_rand,n_rand,/SAMPLE))^2
  rand_distances = sqrt(temporary(d))

  d=(rebin(transpose(data[0,*]),n_pts,n_rand,/SAMPLE)-rebin(rand_locs[0,*],n_pts,n_rand,/SAMPLE))^2 + $
    (rebin(transpose(data[1,*]),n_pts,n_rand,/SAMPLE)-rebin(rand_locs[1,*],n_pts,n_rand,/SAMPLE))^2 ;+ $
  ;(rebin(transpose(data[2,*]),n_pts,n_rand,/SAMPLE)-rebin(rand_locs[2,*],n_pts,n_rand,/SAMPLE))^2
  cross_distances = sqrt(temporary(d))

  ;quick_histplot, distances, xtitle = 'distances between stomata pairs'


  binsize = 20.
  max_val = max([max(distances), max(rand_distances), max(cross_distances)])
  dist_hist = histogram(distances, binsize = binsize, locations = locs, max = max_val, min = binsize)
  rand_hist = histogram(rand_distances, binsize = binsize, locations = locs, max = max_val, min = binsize)
  cross_hist = histogram(cross_distances, binsize = binsize, locations = locs, max = max_val, min = binsize)
  undefine, distances, rand_distances, cross_distances

  window, 3
  cgplot, locs, dist_hist/n_pts^2.
  cgplot, locs, rand_hist/n_rand^2., /over, color='red'
  cgplot, locs, cross_hist/float(n_pts*n_rand), /over, color='blue'

  est1 = float(n_rand * dist_hist) / float(n_pts * cross_hist) -1
  est2 = float(dist_hist * rand_hist) / float(cross_hist^2) -1
  est3 = (dist_hist*(n_rand/float(n_pts))^2. - 2*cross_hist*n_rand/float(n_pts) + rand_hist)/float(rand_hist)

  window,4
  yrange=[-10,10]
  cgplot, locs, est1, yrange = yrange
  cgplot, locs, est2, /over, color='red'
  cgplot, locs, est3, /over, color='blue'



end