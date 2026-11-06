sensitivity_data = read_csv("Wayth2018_MWAPhs2_Sens128.csv")
k = sensitivity_data.field1
delta2_sens = sensitivity_data.field2
power_sens = delta2_sens / k^3

input_hours = 1000

; change this to change integration depth
desired_hours = 600

power_sens *= (float(input_hours)/float(desired_hours))

if desired_hours eq input_hours then begin
  plotfile = "Wayth2018_MWAPhs2_Sens128.pdf"
endif else begin
  plotfile = "Wayth2018_MWAPhs2_Sens128_" + number_formatter(desired_hours) + "hrs.pdf"
endelse

x_factor=2.3
y_factor=1.

base_size = 600
xsize = round(base_size * x_factor)
ysize = round(base_size * y_factor)

ps_aspect = y_factor / x_factor

if ps_aspect lt 1 then landscape = 1 else landscape = 0

charthick = 3
thick = 3
xthick = 3
ythick = 3
charsize = 2
font = 1

margin = [0.15, 0.2, 0.05, 0.1]
plot_pos = [margin[0], margin[1], (1-margin[2]), (1-margin[3])]

plot_aspect = (plot_pos[3] - plot_pos[1]) / (plot_pos[2] - plot_pos[0])

full_aspect = ps_aspect / plot_aspect

sizes = cgpswindow(LANDSCAPE=landscape, aspectRatio = full_aspect, /sane_offsets)

cgps_open, plotfile, /font, /nomatch, inches=sizes.inches, xsize=sizes.xsize, ysize=sizes.ysize, $
  xoffset=sizes.xoffset, yoffset=sizes.yoffset, landscape = landscape

DEVICE, /ISOLATIN1

xtitle = textoidl('k (Mpc^{-1})', font = font)
ytitle = textoidl('P_k (mK^2 Mpc^3)', font = font)

cgplot, k, power_sens, aspect = ps_aspect, /ylog, xrange = [0.017, 1.2], yrange = [1, 1e13], color="TG2", $
  xstyle=1, ystyle=1, axiscolor = 'black', xtitle = xtitle, ytitle = ytitle, psym=-3, $
  ytickformat = 'exponent', thick = thick, charthick = charthick, xthick = xthick, ythick = ythick, charsize = charsize, $
  font = font

cgps_close, /pdf, /delete_ps, density=600
