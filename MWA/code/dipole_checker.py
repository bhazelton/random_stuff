import struct
import math
import os
import sys
import getopt

# define global variables
global CH_MIN
global CH_MAX
global LOUD

def help_message():
	stars = '****************************************************************'
	print '\n%s' % (stars)
	print '  Usage:  dipole_checker.py evaluates the data stored in *_avg files to determine which tiles and polarizations have bad dipoles.  It both outputs results to screen and stores the results in a file in the data directory.\n\n  Mandatory argument:\n    -t YYYYMMDD.(decimal_day):  The starting time for the scan.  The first file from each dipole and receiver directly after this time will be evaluated.  (Specify receivers and dipoles by the -r and -d commands; see below)\n\n  Optional options:\n    -h: prints this help message\n    -f file1 file2, etc: analyze the data in the filenames listed after this flag.  Only include the file names here: paths may be specified with the -p flag, if necessary (see below).\n      Special options:\n        -f 0: prompt the user for file names\n        -f 1: use all files in the specified data directory (This is the default) (see -p below for info on setting the data directory)\n      Caveat:  If either the \'-r\' or \'-d\' flag is used,!
  then potentially only a subset of the files listed will be used.\n    -r recv: only analyze data taken by receiver number recv (default: use all receivers)\n    -T Tile: only output data for the tile Tile\n    -d dip: only analyze data taken with dipole dip (default: use all dipoles)\n    -y delay: only analyze data taken with this delay.  If y==-1, then use all delays (default: use only delay 0)\n    -p pathname: analyze data files located in directory pathname (default: current working directory)\n    -o outpathname: final results are written to a file located in outpathname (default: current working directory)\n    -i boolean: If boolean equals 1, then only consider the first file for each receiver/dipole pair after the start time (set by -t).  If boolean equals 0, then consider all times.  (Default: 1).\n    -X : Sets the expedition number to xnum.  (This is necessary for certain naming conventions). (Default: 12; ie, X12, the Feb 2010 site trip).\n    -l LOUD: LOUD==0!
  means don\'t print status messages to screen.  LOUD==1 means !
 print st
atus messages to screen\n    -n name: name of the ascii file containing the saved results.  (Default: bad_dipoles_dipX_recY.txt, where\'X\' and \'Y\' are the names of the dipoles and receivers to be considered (see the -d and -r flags), or the word \'All\' if none are specified.  Only include the file name here: use the -o flag to set the output directory).\n    -N: Don\'t save the results (good for debugging, results will still be printed to screen)\n\nError codes:\nL#: Gain for dipole is lower than typical by # decibals.\nH#: Gain for dipole is higher than typical by # decibals.\nD: Dipole appears to be dead (ie, no signal).\nF#: Data for dipole is garbage; data appears to be a flat noise floor at # decibals.\nR #1 #2: Data has bad RFI spikes.  Worst spike is at frequency #1 and is a #2 decibal spike above typical.\nE#: Data is fit too well at the highest frequencies by a parabola and/or \'effective potential\'.  # is related to the rms of the fit.  (This is sometimes caus!
 ed because the spectrum is too flat).\nU: Error is detected, but it doesn\'t fit into any of the above categories.\n'
	print '%s\n' % (stars)
	sys.exit(0)

def fit_line( d, nmpts, get_xs ): # Performs a chi-squared fit of the data in "d" (of length "nmpts") to a line, y = mx + b.  The function get_xs() inputs an element number ( 0 to nmpts-1 ) and outputs the value of x of that element.  Returns the tuple (m,b).

	# Calculate the necessary sums.
	dx_sum = d_sum = x2_sum = x_sum = 0.0;
	one_sum = float(nmpts);
	for i in range(0,nmpts):
		x_temp = get_xs(i)
		#print 'i = ', i, 'x = ', x_temp, 'data = ', d[i]
		dx_sum += d[i]*x_temp;
		d_sum += d[i];
		x2_sum += (x_temp*x_temp)
		x_sum += x_temp
	# Calculate the slope and intercept using the analytic solutions to the chi-sq fit
	denom = (one_sum*x2_sum-x_sum*x_sum)
	m = (one_sum*dx_sum - x_sum*d_sum)/denom
	b = (x2_sum*d_sum - x_sum*dx_sum)/denom
	#print 'm = ', m, 'b = ', b
	return ( m, b )

def fit_line_with_gaps( d, nmpts, indices2use, get_xs ): # Same as above, but now d is a vector of data with length greater than nmpts.  Only nmpts of the elements of d are used in the fit, those indices corresponding to the values held in indices2use (a vector of ints of length nmpts).
	# Calculate the necessary sums.
	dx_sum = d_sum = x2_sum = x_sum = 0.0;
	one_sum = float(nmpts);
	for i in range(0,nmpts):
		index = indices2use[i]
		x_temp = get_xs(index)
		#print 'index = ', index, 'x = ', x_temp, 'data = ', d[index]
		dx_sum += d[index]*x_temp;
		d_sum += d[index];
		x2_sum += (x_temp*x_temp)
		x_sum += x_temp
	# Calculate the slope and intercept using the analytic solutions to the chi-sq fit
	denom = (one_sum*x2_sum-x_sum*x_sum)
	m = (one_sum*dx_sum - x_sum*d_sum)/denom
	b = (x2_sum*d_sum - x_sum*dx_sum)/denom
	#print 'm = ', m, 'b = ', b
	return ( m, b )

def fit_effpot_with_gaps( d, nmpts, indices2use, get_xs, low_cut ): # Fit the coefficients of an 'effective potential' curve (a/f-b/f^2+c). The input d is a vector of data with length greater than nmpts.  Only nmpts of the elements of d are used in the fit, those indices corresponding to the values held in indices2use (a vector of ints of length nmpts).
	return fit_effpot_or_parab( d, nmpts, indices2use, get_xs, low_cut, eff_pot_func )

def fit_parab_with_gaps( d, nmpts, indices2use, get_xs, low_cut ): # Fit for a parabola (a*f-b*f^2+c). The input d is a vector of data with length greater than nmpts.  Only nmpts of the elements of d are used in the fit, those indices corresponding to the values held in indices2use (a vector of ints of length nmpts).
	return fit_effpot_or_parab( d, nmpts, indices2use, get_xs, low_cut, parab_func )

def fit_effpot_or_parab( d, nmpts, indices2use, get_xs, low_cut, func ): #option 0: fit effective potential, option 1: fit parabola

	freq0 = chan2freq( low_cut-1 ) # We need to avoid 1/0 in the below...
	#print freq0, low_cut
	# Calculate the necessary sums.
	F1 = F2 = F3 = F4 = d0 = d1 = d2 = 0.0
	F0 = float(nmpts);
	for i in range(0,nmpts):
		index = indices2use[i]
		f_temp = get_xs(index) # get the frequency
		x_temp =  func(f_temp,freq0)
		data = d[index]
		x2_temp = x_temp*x_temp
		#print 'index = ', index, 'x = ', x_temp, 'x2 = ', x2_temp, 'data = ', d[index], 'freq0 = ', freq0, 'f = ', get_xs(index)
		F1 += x_temp
		F2 += x2_temp
		F3 += x_temp*x2_temp
		F4 += x2_temp*x2_temp
		d0 += data
		d1 += data*x_temp
		d2 += data*x2_temp
	# Calculate the fit parameters using the analytic solutions to the chi-sq fit
	if( low_cut==99 ): print d0, d1, d2, d[0:3]
	G1 = (F3*F3-F4*F2)
	G2 = (F1*F3-F2*F2)
	G3 = (d2*F2-d1*F3)
	G4 = (F1*F3-F2*F2)
	G5 = (F1*F1-F0*F2)
	G6 = (d0*F2-d1*F1)
	H1 = (F2*G1)
	H2 = (F1*G1-G2*F3)
	H3 = (-G3*F3-d1*G1)
	H4 = (G5*G1-G2*G4)
	H5 = (G6*G1-G3*G4)
	a = (H5*H2-H3*H4)/(H1*H4)
	b = (H5*G2-G3*H4)/(G1*H4)
	c = -1.0*(H5)/(H4)
	# Calculate the rms error
	rms = 0.0
	for i in range(0,nmpts):
		index = indices2use[i]
		f_temp = get_xs(index) # get the frequency
		x_temp =  func(f_temp,freq0)
		x2_temp = x_temp*x_temp
		data = d[index]
		err = data - (a*x_temp + b*x2_temp+c)
		rms += err*err
	rms = math.sqrt( rms/float(nmpts) )

	return [ a, b, c, rms ]

def eff_pot_func( freq, freq0 ):
	return 1.0/(freq-freq0)

def parab_func( freq, freq0 ):
	return freq

def sub_line( d, nmpts, m, b, get_xs ): # Subtract the line defined by y = mx + b from the data d, a vector of nmpts data points.  get_xs() is a function which inputs an element number (0 to nmpts-1) and outputs the x-value of the data for that indice.  Returns the subtracted data d

	for i in range(0,nmpts):
		x_temp = get_xs(i)
		y_temp = m*x_temp + b
		d[i] -= y_temp

	return d

def chan2freq( chn ): # Inputs a channel number ( 0 to 255 ) and outputs the (minimum) frequency of that channel in MHz.
	return 1.28*float(chn)

def avg_data( d, nmpts, power ): # Inputs a data vector d of nmpts elements and outputs the average value of ( the data raised to the 'power' power )  Ex, power=2 outputs <x^2>
	avg = 0.0
	for i in range(0,nmpts):
		avg += math.pow(d[i],power)
		#print 'data = ', math.pow(d[i],power)

	return avg/float(nmpts)

def chan2freq_shift( chn ): # Same as above, but the channel number is shifted by CH_MIN from the input chn
	return 1.28*float(chn + CH_MIN)

def freqrange2chns( freq_min, freq_max ): # Inputs a minimum and maximum frequency in MHz and outputs a tuple (ch_min1, ch_max1) containing the minimum and maximum channel numbers within that range.
	ch_min1 = int(math.ceil(freq_min/1.28))
	ch_max1 = int(math.floor(freq_max/1.28))
	return (ch_min1, ch_max1)

def pol2XY( pol ): # Inputs a polarization (0 or 1) and outputs 'X' for 0 and 'Y' for 1 
	if( pol == 0 ):
		return 'X'
	else:
		return 'Y'

def closestchan2freq( freq ): # returns the channel number with frequency closest to the input frequency of freq
	ch1 = int(math.ceil(freq/1.28))
	ch2 = int(math.floor(freq/1.28))
	diff1 = math.fabs(chan2freq(ch1)-freq)
	diff2 = math.fabs(chan2freq(ch2)-freq)
	if( diff1 > diff2 ):
		return int(ch2)
	else:

		return int(ch1)

def slot2tile( slot_or_tile, recv, Xnum, mode ): # mode=0: inputs a slot number, and returns the tile number (which depends upon the expedition's setup, thus the 'Xnum'). mode=1: inputs a tile number, and returns a slot number and receiver number 
	if( mode==0 ): # slot input with slot_or_tile, return the tile number
		return int((slot_or_tile-1)*4+recv)
	if( mode==1 ): # tile input with slot_or_tile, return the slot number and receiver number
		rx_num = int(slot_or_tile) % 4
		if( rx_num==0 ):
			rx_num=4
		slot = (int(slot_or_tile) - int(rx_num))/4+1
		return [slot, rx_num]

def cutfiles( start_time, end_time, recv2use, dip2use, delays2use, times, dips, recvs, delays, filenames ): # Inputs a list of filenames and eliminates all files with time stamps before start_time.
	num_files = len(filenames)
	new_filenames = []
	new_dips = []
	new_recvs = []
	new_delays = []
	new_times = []
	# Cycle through and eliminate files that were created before start_time
	new_num_files = 0
	for i in range(num_files):
		time = times[i]
		receiver = recvs[i]
		dipole = dips[i]
		delay = delays[i]
		if( time_good( start_time, end_time, time ) and recv_good( recv2use, receiver ) and dip_good( dips2use, dipole ) and delay_good( delays2use, delay ) ):
			new_num_files = new_num_files + 1
			new_filenames.append(filenames[i])
			new_dips.append(dips[i])
			new_recvs.append(recvs[i])
			new_delays.append(delays[i])
			new_times.append(times[i])
	return [ new_times, new_dips, new_recvs, new_delays, new_filenames ]

def time_good( start_time, end_time, time ): # Returns 0 if the time is before the starting time, 1 otherwise
	if( start_time <= time and end_time >= time ):
		return 1
	else:
		return 0

def recv_good( recv2use, receiver ): # Returns 0 if this receiver is not to be used, 1 otherwise
	if( recv2use==-1 or recv2use==receiver ):
		return 1
	else:
		return 0

def dip_good( dip2use, dipole ): # Returns 0 if this dipole is not to be used, 1 otherwise
	if( dip2use==-1 or dip2use==dipole ):
		return 1
	else:
		return 0

def delay_good( delays2use, delay ): # Returns 0 if this delay is not to be used, 1 otherwise
	if( delays2use==-1 or delays2use==delay ):
		return 1
	else:
		return 0

def isolate_one_set( num_dips, num_recvs, num_delays, times, dips, recvs, delays, filenames ): # Inputs (sorted) lists, and returns only the earliest created file for each combination of dipole, receiver, and delay.
	num_files = len(filenames)
	new_filenames = []
	new_dips = []
	new_recvs = []
	new_delays = []
	new_times = []
	already_used = [  [ [ 0 for h in range(num_delays) ] for i in range(num_dips) ] for j in range(num_recvs) ]
	for i in range(0,num_files):
		receiver = recvs[i]
		dipole = dips[i]
		delay = delays[i]
		if( already_used[receiver-1][dipole-1][delay-1] == 0 ):
			new_filenames.append(filenames[i])
			new_dips.append(dips[i])
			new_recvs.append(recvs[i])
			new_delays.append(delays[i])
			new_times.append(times[i])
			already_used[receiver-1][dipole-1][delay-1] = 1
	return [ new_times, new_dips, new_recvs, new_delays, new_filenames ]

def slow_sort( times, dips, recvs, delays, filenames ): # Sorting algorithm: not the most graceful bit of code (hence, the 'slow').  No files with time stamps earlier than start_time are considered; these are all left out of the final solution.

	global LOUD

	num_files = len(filenames)
	new_filenames = []
	new_dips = []
	new_recvs = []
	new_delays = []
	new_times = []
	already_used = [ 0 for i in range(num_files) ]
	too_early = [ 0 for i in range(num_files) ]
	universal_min_time = 9.99e99
	universal_max_time = -1.0
	for i in range(num_files):
		time = times[i]
		if( universal_min_time > time ): universal_min_time = time
		if( universal_max_time < time ): universal_max_time = time
		new_min_time = 9.99e99
		new_min_dip = 9999
		new_min_recv = 9999
		new_min_delay = 9999
		min_spot = -1
		min_spot_old = -1
		for j in range(num_files):
			if( min_spot_old != min_spot ): # new spot was found last cycle, update the minimum values.
				new_min_time = times[min_spot]
				new_min_dip = dips[min_spot]
				new_min_recv = recvs[min_spot]
				new_min_delay = delays[min_spot]
				min_spot_old = min_spot
			if( already_used[j]==0 ):
				time = times[j]
				dip = int(dips[j])
				recv = recvs[j]
				delay = delays[j]
				if( recv < new_min_recv ):
					min_spot = j
				elif( recv == new_min_recv ):
					if( dip < new_min_dip ):
						min_spot = j
					elif( dip == new_min_dip ):
						if( delay < new_min_delay ):
							min_spot = j
						elif( delay == new_min_delay ):
							if( time < new_min_time ):
								min_spot = j
		already_used[min_spot] = 1
		new_filenames.append(filenames[min_spot])
		new_dips.append(dips[min_spot])
		new_recvs.append(recvs[min_spot])
		new_delays.append(delays[min_spot])
		new_times.append(times[min_spot])

	if LOUD: print 'slow_sort: min time = %14.5f, max time = %14.5f' % (universal_min_time, universal_max_time)

	return [ new_times, new_dips, new_recvs, new_delays, new_filenames ]

def valid_file( filename ): # Outputs 1 is the input filename is valid, 0 otherwise
	match_str1 = 'Dipole'
	match_str2 = '_Rx'
	match_str3 = '_avg'
	if( (match_str1 in filename) and (match_str2 in filename) and (match_str3 in filename) ):
		return 1
	else:
		return 0

def calc_avg_std( num_ref_freqs, num_pols, num_z, fixed_quant, option, data, bad_plots, use_mods ): # Calc the avg and std dev of the data located in data
	# option==0 --> file number is fixed, average over slots
	# option==1 --> slot number is fixed, average over files

	avg_val = [ [ 0 for j in range(num_ref_freqs) ] for k in range(num_pols) ]
	std_val = [ [ 0 for j in range(num_ref_freqs) ] for k in range(num_pols) ]
	#print 'calc_avg num_pols = %d' % (num_pols)
	# Cycle through reference frequencies
	for i in range(0,num_ref_freqs):
		for pol in range(0,num_pols):
			# Calculate the avg values for each polarizations
			tot = tot2 = 0.0 # tot = sum total, tot2 = sum of square values
			nmpts = 0
			for z in range(0,num_z): # cycle through zs
				if( option==0 ):
					bad_val = bad_plots[z][pol][fixed_quant]
				else:
					bad_val = bad_plots[fixed_quant][pol][z]
					#if (pol==1): print 'pol = 1, fixed = %d, z = %d, bad_val = %d' % (fixed_quant, z, bad_val)
				if( bad_val==0 ):
					temp = data[z][pol][i]
					tot += temp
					tot2 += temp*temp
					nmpts = nmpts+1
			#print 'not quite there i = %d, pol = %d, nmpts = %d' % (i, pol, nmpts)
			if( nmpts!=0 ):
				#print 'here! i = %d, pol = %d, nmpts = %d' % (i, pol, nmpts)
				N = float(nmpts)
				avg = tot/N
				std = math.sqrt( (tot2/N - avg*avg) )
				if( use_mods and nmpts>2 ):
					# Calculate the 'modified' averages and standard deviations.  These are the average and standard deviation of the set that includes all the data points, except for two data points.  The two data points excluded are those that lead to the lowest standard deviation when left out.
					min_std = 999e99
					for j in range(0,num_z):
						if( option==0 ):
							bad_val = bad_plots[j][pol][fixed_quant]
						else:
							bad_val = bad_plots[fixed_quant][pol][j]
						if( bad_val==0 ):
							for k in range(j+1,num_z):
								if( option==0 ):
									bad_val = bad_plots[k][pol][fixed_quant]
								else:
									bad_val = bad_plots[fixed_quant][pol][k]
								if( bad_val==0 ):
									temp1 = data[j][pol][i]
									temp2 = data[k][pol][i]
									temp_val1 = tot - temp1 - temp2
									temp_val2 = tot2 - temp1*temp1 - temp2*temp2
									temp_val1 /= float(N-2)
									#print temp_val2, temp_val1, N
									#print temp_val2/float(N-2), temp_val1*temp_val1
									temp_val2 = math.sqrt( math.fabs( (temp_val2/float(N-2) - temp_val1*temp_val1) ) )
									if( temp_val2 < min_std ):
										min_std = temp_val2
										mod_avg = temp_val1
										mod_std = temp_val2
										left_out1 = j+1
										left_out2 = k+1
										#if LOUD: print 'looking at (%d,%d) for pol %s, avg = %2.3f, std = %2.3f' % (j+1,k+1,pol2XY(pol),temp_val1, temp_val2)
					avg_val[pol][i] = mod_avg
					std_val[pol][i] = mod_std
					if LOUD: print '%s:  Modified vals calc\'d: freq = %d, pol = %s, left out zs (%d,%d), avg = %2.3f, std = %2.3f, mod avg = %2.3f, mod std = %2.3f' % (prog_name, i, pol2XY(pol), left_out1, left_out2, avg, std, avg_val[pol][i], std_val[pol][i])
				else:
					avg_val[pol][i] = avg
					std_val[pol][i] = math.sqrt( (tot2/N - avg*avg) ) # this is a stddev of the data set, not an unbiased estimate of a sqrt variance (ie, divided by sqrt(N), not sqrt(N-1) )
					if LOUD: print '%s:  Non-modified vals calc\'d: freq = %d, pol = %s, nmpts = %d, avg = %2.3f, std = %2.3f' % (prog_name, i, pol2XY(pol), nmpts, avg_val[pol][i], std_val[pol][i])

	return [ avg_val, std_val ]

def calc_avg_std2( slots, pols, freqs, files, allpows, bad_plots, use_mods ): # Calc the avg and std dev of the data located in data

	num_slots = len(slots)
	num_pols = len(pols)
	num_ref_freqs = len(freqs)
	num_files = len(files)

	avg_val = [ [ 0 for j in range(num_ref_freqs) ] for k in range(num_pols) ]
	std_val = [ [ 0 for j in range(num_ref_freqs) ] for k in range(num_pols) ]

	if( num_files==1 ):
		option = 0
		num_zs = num_slots
	elif( num_slots==1 ):
		option = 1
		num_zs = num_files
	else:
		print 'ERROR: Input to calc_std_avg() must have either len(slots)==1 or len(files)==1.  Aborting...'

	for i, fitem in enumerate(freqs):  # Cycle through reference frequencies
		for pol, pitem in enumerate(pols):  # Cycle through polarizations
			# First, determine how many dipoles are still potentially good.  Calculate the average value and value squared of these good dipoles.
			nmpts = 0
			tot = tot2 = 0.0
			for z in range(0,num_zs): # cycle through slots or files
				[slot,filenum] = find_slot_filenum( option, z, slots, files )
				if( bad_plots[slot][pol][filenum]==0 ):
					temp = allpows[slot][pol][i][filenum]
					tot += temp
					tot2 += temp*temp
					nmpts = nmpts+1
			if( nmpts!=0 ):
				N = float(nmpts)
				avg = tot/N
				std = math.sqrt( math.fabs( (tot2/N - avg*avg) ) )
				if( use_mods and nmpts>2 ):
					# Calculate the 'modified' averages and standard deviations.  These are the average and standard deviation of the set that includes all the data points, except for two data points.  The two data points excluded are those that lead to the lowest standard deviation when left out.
					min_std = 999e99
					for j in range(0,num_zs): # cycle through slots or files
						[slot1,filenum1] = find_slot_filenum( option, j, slots, files )
						#print 'j: slot, file, num_zs', slot1, filenum1, bad_plots[slot1][pol][filenum1], num_zs
						if( bad_plots[slot1][pol][filenum1]==0 ):
							for k in range(j+1,num_zs):
								[slot2,filenum2] = find_slot_filenum( option, k, slots, files )
								#print 'k: slot, file, bad', slot2, filenum2, bad_plots[slot2][pol][filenum2]
								if( bad_plots[slot2][pol][filenum2]==0 ):
									temp1 = allpows[slot1][pol][i][filenum1]
									temp2 = allpows[slot2][pol][i][filenum2]
									temp_val1 = tot - temp1 - temp2
									temp_val2 = tot2 - temp1*temp1 - temp2*temp2
									temp_val1 /= float(N-2)
									#print temp_val2, temp_val1, N
									#print temp_val2/float(N-2), temp_val1*temp_val1
									temp_val2 = math.sqrt( math.fabs( (temp_val2/float(N-2) - temp_val1*temp_val1) ) )
									if( temp_val2 < min_std ):
										min_std = temp_val2
										mod_avg = temp_val1
										mod_std = temp_val2
										if( option==0 ):
											left_out1 = slot1+1
											left_out2 = slot2+1
										else:
											left_out1 = filenum1+1
											left_out2 = filenum2+1
										#if LOUD: print 'looking at (%d,%d) for pol %s, avg = %2.3f, std = %2.3f' % (j+1,k+1,pol2XY(pol),temp_val1, temp_val2)
					avg_val[pol][i] = mod_avg
					std_val[pol][i] = mod_std
					if LOUD: print '%s:  Modified vals calc\'d: freq = %d, pol = %s, left out zs (%d,%d), avg = %2.3f, std = %2.3f, mod avg = %2.3f, mod std = %2.3f' % (prog_name, i, pol2XY(pol), left_out1, left_out2, avg, std, avg_val[pol][i], std_val[pol][i])
				else:
					avg_val[pol][i] = avg
					std_val[pol][i] = std # this is a stddev of the data set, not an unbiased estimate of a sqrt variance (ie, divided by sqrt(N), not sqrt(N-1) )
					if LOUD: print '%s:  Non-modified vals calc\'d: freq = %d, pol = %s, nmpts = %d, avg = %2.3f, std = %2.3f' % (prog_name, i, pol2XY(pol), nmpts, avg_val[pol][i], std_val[pol][i])

	return [ avg_val, std_val ]

def find_slot_filenum( option, z, slots, files ):
	if( option == 0 ):
		filenum = files[0]
		slot = slots[z]
	elif( option == 1 ):
		filenum = files[z]
		slot = slots[0]
	else:
		print 'ERROR: invalid input to find_slot_filenum.  option must be 0 or 1.  Aborting...'
		sys.exit(1)

	return [ slot, filenum ]

if __name__ == '__main__':

	global CH_MIN
	global CH_MAX
	global LOUD

	# Set important variables
	prog_name = "dipole_checker.py"
	num_chs = 256 # total number of coarse channels
	spr = 8 # slots per receiver
	num_pols = 2 # number of polarizations, X and Y
	num_dips = 16 # number of dipoles per tile
	num_rxs = 4 # number of receivers
	num_delays = 32 # max possible number of delay lines
	num_tiles = num_rxs*spr
	compare_option = 1


	print '\n%s: Commencing program\n' % (prog_name)


	# Set the defaults
	LOUD = 0
	dipole = -1
	receiver = -1
	get_files_from_dir = 1
	prompt_user = 0
	one_set = 1
	Xnum = 12
	use_r = use_d = use_p = use_f = use_o = use_n = use_t = use_y = nowrite = use_T = 0

	# Read in the command line arguments, if any
	try:
		options, extra_args = getopt.getopt(sys.argv[1:], 'hf:r:d:l:p:o:n:t:i:y:X:NT:')
	except getopt.error:
		if LOUD: print '%s:  ERROR: Unknown commandline argument entered. Printing function usage and then exiting...\n' % (prog_name)
		help_message()
		sys.exit(0)
	for input_opt, arg in options[:]:
		if input_opt == '-h':
			help_message()
		elif input_opt == '-t':
			use_t = 1
			start_time = float(arg)
		elif input_opt == '-f':
			use_f = 1
			if( arg == '0' ):
				if LOUD: print '%s: As per specified on the command line via -f 0, the user will be prompted for files.' % (prog_name)
				prompt_user = 1
			elif( arg=='1' ):
				if LOUD: print '%s: As per specified on the command line via -f 1, all files in the data directory will be considered.' % (prog_name)
				get_files_from_dir = 1
				prompt_user = 0
			else:
				if LOUD: print '%s: File names entered via the command line' % (prog_name)
				filenames = []
				filenames.append(arg)
				for item in extra_args:
					filenames.append(item)
				prompt_user = 0
		elif input_opt == '-r':
			use_r = 1
			receiver = int(arg)
			if LOUD: print '%s: Only files involving receiver %d will be examined.' % (prog_name, receiver)
		elif input_opt == '-d':
			use_d = 1
			dipole = int(arg)
			if LOUD: print '%s: Only files involving dipole %d will be examined.' % (prog_name, dipole)
		elif input_opt == '-y':
			use_y = 1
			delay = int(arg)
			if LOUD: print '%s: Only files involving delay %d will be examined.' % (prog_name, delay)
		elif input_opt == '-l':
			loudt = int(arg)
			if( loudt==0 ): print '%s: LOUD entered via commandline as 0.  Suppressing output...' % (prog_name)
			LOUD = loudt
		elif input_opt == '-p':
			use_p = 1
			data_dir = arg
		elif input_opt == '-o':
			use_o = 1
			out_dir = arg
		elif input_opt == '-n':
			use_n = 1
			outfilename = arg
		elif input_opt == '-i':
			one_set = int(arg)
		elif input_opt == '-T':
			use_T = 1
			use_tile_only = int(arg)
		elif input_opt == '-X':
			Xnum = int(arg)
		elif input_opt == '-N':
			nowrite = 1
		else: # Actually, this line shouldn't be reached because of the 'except getopt.error' line above
			print '%s: ERROR: Unknown input option %s entered on the command line.  Printing function usage and then exiting...' % (prog_name, input_opt)
			help_message()
	# Exit if no start time is entered
	if( use_t == 0 ):
		print '%s: ERROR: beginning time must be entered on the command line.  Printing function usage and then exiting...' % (prog_name)
		help_message()		
	# Determine the directories to use
	# (i) the data file directory...
	if( use_p == 0 ):
		data_dir = os.getcwd()
		if LOUD: print '%s: Data files will be pulled from the current working directory, %s (use -p dirname on the command line to change this)' % (prog_name, data_dir)
	else:
		if LOUD: print '%s: As specified via command line argument -p, files will be pulled from the directory, %s' % (prog_name, data_dir)
	# (ii) the output file directory...
	if( use_o == 0 ):
		out_dir = os.getcwd()
		if LOUD: print '%s: Final results will be saved to the current working directory, %s (use -o dirname on the command line to change this)' % (prog_name, out_dir)
	else:
		if LOUD: print '%s: As specified via command line argument -o, files will be saved to the directory, %s' % (prog_name, out_dir)
	# Determine the name of the output file.
	if( use_n == 0 ):
		if( use_d ):
			dipole_str = '%02d' % (dipole)
		else:
			dipole_str = 'All'
		if( use_r ):
			receiver_str = '%1d' % (receiver)
		else:
			receiver_str = 'All'
		if( use_y ):
			if( delay == -1 ):
				delay_str = 'All'
			else:
				delay_str = '%1d' % (delay)
		else:
			delay_str = '0'
		if( use_T ):
			outfilename = 'tile_%d_dip%s_rec%s_delay%s_t%014.5f.txt' % (use_tile_only, dipole_str, receiver_str, delay_str, start_time)
		else:
			outfilename = 'dip%s_rec%s_delay%s_t%014.5f.txt' % (dipole_str, receiver_str, delay_str, start_time)
		g_name = 'vals_%s' % (outfilename)
		b_name = 'errs_%s' % (outfilename)
		if LOUD: print '%s: Output file name will be %s and %s (use -n outfilename on the command line to change the back end here)' % (prog_name, g_name, b_name)
	else:
		g_name = 'vals_%s' % (outfilename)
		b_name = 'errs_%s' % (outfilename)
		if LOUD: print '%s: As specified via command line argument -n, the output file will be named %s and %s' % (prog_name, g_name, b_name)

	# If a user prompt was not requested, get all legitimate file names from the data directory, data_dir
	if( get_files_from_dir ):
		filenames = os.listdir(data_dir)
		filenames2 = []
		# Make sure that the file name is valid
		for item in filenames:
			if( valid_file( item ) ):
				filenames2.append(item)
		filenames = filenames2

	# If file names not entered via command line, then either prompt the user (if prompt_user==1) or use the default (if prompt_user==-1; this option is good for debugging)
	if( prompt_user == 1):
		filenames = []
		stay_in_loop = 1
		while(stay_in_loop==1):
			input_str = '%s: Please input a filename (Use the -f flag to avoid this prompt.  Input \'END\' if no more file names are to be entered):\n' % (prog_name)
			filename = raw_input(input_str)
			if( filename != 'END' ): # If a name is entered...
				# ... first check to make sure that the name is valid
				if( valid_file( filename )==0 ):
					print '%s: ERROR: Input file name  ( = %s ) is not valid.  Ignoring this input...' % (prog_name, filename)
				# ... and if it is, then add it to the list of files
				else:
					filenames.append(filename)
			else: # Otherwise, break the loop
				stay_in_loop = 0
	elif( prompt_user == -1): # else, just use the file listed below (good for debugging only)
		filenames.append('Dipole10_Rx1_20091123.46998_avg')

	# Determine the dipole number, receiver number, and time of files in the list
	num_files = len(filenames)
	times = []
	dips = []
	recvs = []
	delays = []
	if( Xnum==9 ):
		shift = 0
	else:
		shift = 1
	for item in filenames:
		parts = item.split('_')
		times.append(float(parts[2+shift]))
		temp_str = parts[0].split('e')
		dips.append(int(temp_str[1]))
		temp_str = parts[1+shift].split('x')
		recvs.append(int(temp_str[1]))
		if( Xnum != 9 ):
			delays.append(int(parts[1]))
		else:
			delays.append(-1)
			
	# Cut all files that are created before the determined starting time (-t), are not of the correct receiver number (-r), or are not of the correct dipole number (-d)
	if( use_r==1 ):
		recvs2use = receiver
	else:
		recvs2use = -1
	if( use_d==1 ):
		dips2use = dipole
	else:
		dips2use = -1
	if( use_y==1 ):
		delays2use = delay
	else:
		delays2use = 0
	time_length = 5.0 # only look for files created within time_length minutes of the start time
	end_time = start_time + time_length/(60.0*24.0)
	[ times, dips, recvs, delays, filenames ] = cutfiles( start_time, end_time, recvs2use, dips2use, delays2use, times, dips, recvs, delays, filenames )

	# Sort the files.
	[ times, dips, recvs, delays, filenames ] =  slow_sort( times, dips, recvs, delays, filenames )

	# If desired, isolate one set of 16 dipoles for all possible receivers
	if( one_set!=0 ):
		[ times, dips, recvs, delays, filenames ] = isolate_one_set( num_dips, num_rxs, num_delays, times, dips, recvs, delays, filenames )

	num_files = len(filenames)

	# Make sure that at least one file is chosen.  Abort if not.
	if(num_files==0):
		if LOUD: print '%s: No valid names entered.  Aborting...' % (prog_name)
		sys.exit(1)

	if LOUD:
		print '%s: Files to be examined:' % (prog_name)
		print filenames

	# Define the frequencies to consider
	num_ref_freqs = 3 # number of reference frequencies.  Here, =3 (80, 120, and 160 MHz)
	reffreq = [ 80.0, 120.0, 160.0 ] # reference frequencies in MHz

	# Create the array that holds all power readings
	allpows = [ [ [ [ 0 for h in range(len(filenames)) ] for i in range(num_ref_freqs) ] for j in range(num_pols) ] for k in range(spr) ]
	bad2 = [ [ [ 0 for h in range(len(filenames)) ] for j in range(num_pols) ] for k in range(spr) ]
	bad_props2 = [ [ [ 0 for h in range(len(filenames)) ] for j in range(num_pols) ] for k in range(spr) ]
	bad_props3 = [ [ [ 0 for h in range(len(filenames)) ] for j in range(num_pols) ] for k in range(spr) ]
	# The extra dimension 'i' in allms and allbs below is so that calc_avg_std() may be used on allms and allbs (this function was tailored for allpows above.
	allms = [ [ [ [ 0 for h in range(len(filenames)) ] for i in range(0,1) ] for j in range(num_pols) ] for k in range(spr) ]
	allbs = [ [ [ [ 0 for h in range(len(filenames)) ] for i in range(0,1) ] for j in range(num_pols) ] for k in range(spr) ]

	# Cycle through files
	
	for index, item in enumerate(filenames):

		# SECTION I: Read in the data

		filename = item
		filenumber = int(index)
		if LOUD: print '\n%s: Reading in the data from file %s' % (prog_name, filename)

		read_in_data = [ [ [ 0 for i in range(num_chs) ] for j in range(num_pols) ] for k in range(spr) ]
		read_from_file = '%s/%s' % (data_dir, filename)
		try:
			data_file = open(read_from_file, "rb") # open the file for reading
		except IOError:
			print '%s: Can\'t open file %s for reading.  Aborting...' % (prog_name, read_from_file)
	       		sys.exit(2)

		for slot in range(0,spr): # cycle through slots
			for pol in range (0,num_pols): # cycle through polarizations
				for ch in range(0,num_chs): # cycle through coarse channels
					try:
						s = data_file.read(4) # read in a float
					except IOError:
						print '%s: Can\'t read in float from file %s.  Aborting...' % (prog_name, filename)
	       					sys.exit(0)
					value = struct.unpack("f", s)
					value = value[0]
					if( value > 0.0 ):
						read_in_data[slot][pol][ch] = 10.0*math.log10(value) # 1.0/(chan2freq(ch)-chan2freq(38)) - 2.0/(chan2freq(ch)-chan2freq(38))/(chan2freq(ch)-chan2freq(38)) + 5.0 #effpot
						#if( slot==0 and pol==0 ): print 'ch, data = ', ch, read_in_data[slot][pol][ch]
					else:
						read_in_data[slot][pol][ch] = -1.0

		bad = [ [ 0 for j in range(num_pols) ] for k in range(spr) ] # Boolean array for bad slots/pols: 0 = good, 1 = bad
		bad_props = [ [ 0.0 for j in range(num_pols) ] for k in range(spr) ] # float array containing information about bad slots/pols

		# SECTION II: Fit and subtract a line from the data, using only freq<50MHz and freq>300MHz (ie, fit and subtract the noise floor)

		sub_noise_floor = 1
		if( sub_noise_floor==1 ):
			temp_str = 'and subtracting '
		else:
			temp_str = ''
		if LOUD: print '\n%s: Fitting %sthe out-of-freq-range noise floor' % (prog_name, temp_str)
		# The extra unnecessary dimension is so that I can use the function calc_avg_std(), created for another part of the code
		# Determine the frequency bins to use in the calculation
		low_only = 0 # if low_only==1, then only fit a line to freqs less than 50MHz.  if low_only==0, then include freqs greater than 300MHz as well
		ch50 = closestchan2freq( 50.0 ) # channel number closest to 50MHz
		if( low_only==1 ):
			nchans2fit = ch50
		else:
			ch300 = closestchan2freq( 300.0 ) # channel number closest to 300MHz
			nchans2fit = ch50 + (num_chs - ch300) # Notice: we're excluding the DC component here (bin # 0)
		# Create indices2use[], a list of all the elements of the data to be used in the fit (needed for fit_line_with_gaps)
		indices2use =  [ 0 for i in range(0,nchans2fit) ]
		for i in range(1,ch50+1):
			indices2use[i-1] = i
		if( low_only==0 ):
			for i in range(0,num_chs-ch300):
				indices2use[i+ch50] = ch300+i
		# Cycle through slots and polarizations, fitting lines to the noise floor, and then subsequently subtracting that line
		for slot in range(0,spr): # cycle through slots
			for pol in range (0,num_pols): # cycle through polarizations
				# Fit a line to the extremes of this data
				(m,b) = fit_line_with_gaps( read_in_data[slot][pol][0:num_chs], nchans2fit, indices2use, chan2freq )
				allms[slot][pol][0][index] = m
				allbs[slot][pol][0][index] = b
				if LOUD: print '%s: slot = %d%s, fit noise line: slope=%3.3g, intercept=%3.3g' % (prog_name,slot+1,pol2XY(pol),m,b)
				if( sub_noise_floor ):
					# Subtract the fit line from the data
					read_in_data[slot][pol][0:num_chs] = sub_line( read_in_data[slot][pol][0:num_chs], num_chs, m, b, chan2freq )
		# Search for bad files based upon these values for the slope and y intercepts
		use_mods = 1
		avg_slope = [ [ 0 for j in range(1) ] for k in range(num_pols) ]
		std_slope = [ [ 0 for j in range(1) ] for k in range(num_pols) ]
		avg_intcp = [ [ 0 for j in range(1) ] for k in range(num_pols) ]
		std_intcp = [ [ 0 for j in range(1) ] for k in range(num_pols) ]
		[ avg_intcp, std_intcp ] = calc_avg_std2( range(0,spr), range(0,num_pols), range(0,1), range(index,index+1), allbs, bad2, use_mods )
		crit1 = 0 # if 1, then compare the fit intercept to the other fit values.
		crit2 = 1 # if 1, then look to see if the noise floor is greater than 40dB
		for pol in range(0,num_pols):
			# Cycle through all slots, looking for problems
			for slot in range(0,spr): # cycle through slots
				temp = allbs[slot][pol][0][index]
				if( crit1==1 ):
					avg = avg_intcp[pol][0]
					std = std_intcp[pol][0]
					diff = math.fabs( temp - avg )
					if( (diff>std) and (diff>5.0)  ):
						if LOUD: print '%s: slot = %d%s has funny noise floor.  fit intcp = %2.3f, avg = %2.3f, diff = %2.3f, std = %2.3f' % (prog_name,slot+1,pol2XY(pol),temp,avg,diff,std)
						bad2[slot][pol][index] = 5
						bad_props2[slot][pol][index] = temp - avg
				if( crit2==1 ):
					if( temp >= 40.0  ):
						if LOUD: print '%s: slot = %d%s has funny noise floor.  fit intcp = %2.3f: seems too high' % (prog_name,slot+1,pol2XY(pol),temp)
						bad2[slot][pol][index] = 5
						bad_props2[slot][pol][index] = temp


		# SECTION IIb: Search for bad RFI


		# Determine the frequency bins to use in the calculation
		sat_min = closestchan2freq( 240.0 ) # min channel number for satellite band
		sat_max = closestchan2freq( 275.0 ) # max channel number for satellite band
		num_sats = sat_max - sat_min + 1
		orb_min = closestchan2freq( 137.0 ) # min channel number for orbcom band
		orb_max = closestchan2freq( 139.0 ) # max channel number for orbcom band
		num_orbs = orb_max - orb_min + 1
		nchans2use = num_chs - (num_sats + num_orbs) - 2 # The final -1 is to exclude the DC component and final freq
		#print num_chs, nchans2use
		#print sat_min, sat_max, orb_min, orb_max
		# Create indices2use[], a list of all the elements of the data to be used when searching for RFI
		indices2use =  [ 0 for i in range(0,nchans2use) ]
		i = 0
		for j in range(0,num_chs):
			if( (j!=0) and (j<orb_min or j>orb_max) and (j<sat_min or j>sat_max) and (j!=num_chs-1) ):
				indices2use[i] = j
				i = i+1

		bad_spike = 5.0 # decimal jump of what is considered a bad RFI spike
		real_bad_spike = 20.0
		toss_threshold = 50.0
		for slot in range(0,spr): # cycle through slots
			for pol in range (0,num_pols): # cycle through polarizations
				if(bad2[slot][pol][index]==0):
					max_jump = -999999.9
					num_spikes = 0
					toss = 0.0
					for i in range(0,nchans2use):
						spot = indices2use[i]
						jump = read_in_data[slot][pol][spot] - (read_in_data[slot][pol][spot-1] + read_in_data[slot][pol][spot+1])/2.0
						if( jump > max_jump ):
							max_jump = jump
							bad_chan = spot
						if( jump > bad_spike ):
							num_spikes = num_spikes+1
							toss = toss + jump
					if( (max_jump>real_bad_spike) or (toss>toss_threshold) ):
						if LOUD: print '%s: Bad spike(s) found in file %s!  pol = %s, slot = %d, freq = %3.2f, jump = %2.2f, num_spikes = %d, toss = %2.2f' % (prog_name, filename, pol2XY(pol), slot+1, chan2freq(bad_chan), max_jump, num_spikes, toss)
						bad2[slot][pol][index] = 6
						bad_props2[slot][pol][index] = chan2freq(bad_chan)
						bad_props3[slot][pol][index] = max_jump


		fit_effpot = 0
		if( fit_effpot==1 ):
			effpot_str = 'an \'effective potential-esque\' function'
			allrms_eff = [ [ 99999.99e99 for j in range(spr) ] for k in range(num_pols) ]
			allrms_chs_eff = [ [ 0.0 for j in range(spr) ] for k in range(num_pols) ]
			all_a1s = [ [ 0.0 for j in range(spr) ] for k in range(num_pols) ]
			all_b1s = [ [ 0.0 for j in range(spr) ] for k in range(num_pols) ]
			all_c1s = [ [ 0.0 for j in range(spr) ] for k in range(num_pols) ]
		else:
			effpot_str = ''
		fit_parab = 1
		if( fit_parab==1 ):
			parab_str = 'a parabolic function'
			allrms_par = [ [ 99999.99e99 for j in range(spr) ] for k in range(num_pols) ]
			allrms_chs_par = [ [ 0.0 for j in range(spr) ] for k in range(num_pols) ]
			all_a2s = [ [ 0.0 for j in range(spr) ] for k in range(num_pols) ]
			all_b2s = [ [ 0.0 for j in range(spr) ] for k in range(num_pols) ]
			all_c2s = [ [ 0.0 for j in range(spr) ] for k in range(num_pols) ]
		else:
			parab_str = ''
		if( fit_effpot and fit_parab ):
			and_str = ' and '
		else:
			and_str = ''
		if LOUD: print '\n%s: Fitting for %s%s%s' % (prog_name,effpot_str,and_str,fit_parab)

		if( fit_effpot==1 and fit_parab==1 ):
			rms_cut = 5.0 # acceptable rms lower cut.  If rms < rms_cut, then we've found the bad shape: mark this file as bad
			kmin = 70
			kmax = 80
		elif( fit_effpot==1 ):
			rms_cut = 5.0
			kmin = 70
			kmax = 80
		elif( fit_parab==1 ):
			rms_cut = 1.0
			kmin = 80
			kmax = 80

		for k in range(kmin,kmax+1):
			low_cut = k #closestchan2freq( low_freq ) # lowest frequency channel to use
			low_only = 1
			if( low_only==1 ):
				nchans2fit = num_chs - low_cut
				indices2use =  [ i for i in range(low_cut,num_chs) ]
			for slot in range(0,spr): # cycle through slots
				for pol in range (0,num_pols): # cycle through polarizations
					if(bad2[slot][pol][index]==0):
						# Fit for an effective potential
						if( fit_effpot ):
							[ a1, b1, c1, rms1 ] = fit_effpot_with_gaps( read_in_data[slot][pol][0:num_chs], nchans2fit, indices2use, chan2freq, low_cut )
							if( rms1 < allrms_eff[pol][slot] ):
								allrms_eff[pol][slot] = rms1
								allrms_chs_eff[pol][slot] = k
								all_a1s[pol][slot] = a1
								all_b1s[pol][slot] = b1
								all_c1s[pol][slot] = c1
						# Fit for a parabola
						if( fit_parab ):
							[ a2, b2, c2, rms2 ] = fit_parab_with_gaps( read_in_data[slot][pol][0:num_chs], nchans2fit, indices2use, chan2freq, low_cut )
							if( rms2 < allrms_par[pol][slot] ):
								allrms_par[pol][slot] = rms2
								allrms_chs_par[pol][slot] = k
								all_a2s[pol][slot] = a2
								all_b2s[pol][slot] = b2
								all_c2s[pol][slot] = c2
		for slot in range(0,spr): # cycle through slots
			for pol in range (0,num_pols): # cycle through polarizations
				if(bad2[slot][pol][index]==0):
					if( fit_effpot==1 and fit_parab==1 ):
						value = (allrms_eff[pol][slot]*allrms_par[pol][slot])
					elif( fit_effpot==1 ):
						value = (allrms_eff[pol][slot])
					elif( fit_parab==1 ):
						value = (allrms_par[pol][slot])
					if( value<rms_cut ):
						if LOUD:
							print '%s: Bad form at slot %d%s for file %s!  ' % (prog_name,slot+1,pol2XY(pol),filename)
							if( fit_effpot ):
								chan1 = allrms_chs_eff[pol][slot]
								print '    Eff pot: rms = %1.2f, vals=(%3.4f, %3.4f, %3.4f), freq=ch %d, %03.2fMHz ' % (allrms_eff[pol][slot],all_a1s[pol][slot],all_b1s[pol][slot],all_c1s[pol][slot],chan1,chan2freq(chan1) )
							if( fit_parab ):
								chan2 = allrms_chs_par[pol][slot]
								print '    Parabola: rms = %1.2f, vals=(%3.4f, %3.4f, %3.4f), freq=ch %d, %03.2fMHz' % (allrms_par[pol][slot],all_a2s[pol][slot],all_b2s[pol][slot],all_c2s[pol][slot],chan2,chan2freq(chan2) )
							print '    Overall: prod = %1.2f, cut = %1.2f' % ( value,rms_cut )
						bad2[slot][pol][index] = 7
						bad_props2[slot][pol][index] = value

		# SECTION III: Calculate the power at various frequencies

		if LOUD: print '\n%s: Calculating the power at various frequencies' % (prog_name)
		#Determine the frequency bins associated with these three frequencies
		channels = [ closestchan2freq( reffreq[i] ) for i in range(num_ref_freqs) ]
		#Set the number of bins around this bin to average over
		bins2avg = 2 # number of bins of each side, total bins in average is 2*bins2avg+1
		totbins = 2*bins2avg+1
		# Create the array to hold these power levels
		refpow = [ [ [ 0 for i in range(num_ref_freqs) ] for j in range(num_pols) ] for k in range(spr) ]
		# Cycle through slots and polarizations, and estimate the relevent quantities
		for slot in range(0,spr): # cycle through slots
			for pol in range (0,num_pols): # cycle through polarizations
				for i in range(0,num_ref_freqs):
					chan_temp = channels[i]
					refpow[slot][pol][i] = avg_data( read_in_data[slot][pol][chan_temp-bins2avg:chan_temp+bins2avg+1], totbins, 1 )
					# Add these power readings to the full data array
					allpows[slot][pol][i][index] = refpow[slot][pol][i]
		# Print the calculated values
		for slot in range(0,spr): # cycle through slots
			for pol in range (0,num_pols): # cycle through polarizations
				if LOUD: print '%s: slot = %d%s, power levels = (%2.3f,%2.3f,%2.3f) at MHz = (80,120,160), respectively' % (prog_name,slot+1,pol2XY(pol),refpow[slot][pol][0],refpow[slot][pol][1],refpow[slot][pol][2])

		if LOUD: print '\n%s: Looking for dead dipoles across single receivers' % (prog_name)

		# Find dead dipoles
		use_mod = 1
		avg_val = [ [ 0 for j in range(num_ref_freqs) ] for k in range(num_pols) ]
		std_val = [ [ 0 for j in range(num_ref_freqs) ] for k in range(num_pols) ]
		[ avg_val, std_val ] = calc_avg_std2( range(0,spr), range(0,num_pols), reffreq, range(index,index+1), allpows, bad2, use_mods )
		dead_cut = 3.0 # plots of dead_cut db or less in all three frequency channels are mostly likely dead
		for pol in range(0,num_pols):
			# Cycle through all slots, looking for outliers
			for slot in range(0,spr): # cycle through slots
				if( bad2[slot][pol][index]==0 ):
					dead = 0
					for i in range(0,num_ref_freqs):
						temp = refpow[slot][pol][i]
						diff = math.fabs( temp - avg_val[pol][i] )
						if( (temp<0.2*avg_val[pol][i]) or temp<dead_cut ):
							dead = dead + 1
							if LOUD: print '%s: slot = %d%s may be dead.  At freq = %3.3f MHz., power = %2.3f avg = %2.3f' % (prog_name,slot+1,pol2XY(pol),reffreq[i],temp,avg_val[pol][i])
					if( dead == num_ref_freqs ):
						bad2[slot][pol][index] = 2





		if( compare_option==0 ): # In this case, we're comparing different slots in the same receiver

			# SECTION IV: Make a quick recommendation based upon these values of the power

			# Look for other peculiarities
			[ avg_val, std_val ] = calc_avg_std2( range(0,spr), range(0,num_pols), reffreq, range(index,index+1), allpows, bad2, use_mods )
			for pol in range(0,num_pols):
				# Cycle through all slots, looking for outliers
				for slot in range(0,spr): # cycle through slots
					if( bad2[slot][pol][index] == 0 ):
						low = high = some_bad = 0
						low_avg = high_avg = 0.0
						for i in range(0,num_ref_freqs):
							temp = refpow[slot][pol][i]
							avg = avg_val[pol][i]
							diff = math.fabs( temp - avg )
							if( (diff>std_val[pol][i]) and (diff>3.0)  ):
								some_bad = 1
								if( temp < avg ):
									low = low + 1
									low_avg = low_avg + diff
								elif( temp > avg ):
									high = high + 1
									high_avg = high_avg + diff
								if LOUD: print '%s: slot = %d%s may be bad.  At freq = %3.3f MHz., power = %2.3f avg = %2.3f, diff = %2.3f, std dev = %2.3f' % (prog_name,slot+1,pol2XY(pol),reffreq[i],refpow[slot][pol][i],avg_val[pol][i],diff,std_val[pol][i])
						if( some_bad==1 ):
							if( low==num_ref_freqs ):
								bad2[slot][pol][index] = 3
								bad_props2[slot][pol][index] = low_avg/float(num_ref_freqs)
							elif( high==num_ref_freqs ):
								bad2[slot][pol][index] = 4
								bad_props2[slot][pol][index] = high_avg/float(num_ref_freqs)
							else:
								bad2[slot][pol][index] = 1


			# SECTION V: Print results to screen and save the results.

			# Determine the dipole number, receiver number, delay and time
			dip_num = dips[index]
			recv_num = recvs[index]
			time = times[index]
			delay_num = delays[index]
			if LOUD: print '%s: File %s is for dipole %d, receiver %d, delay %d, time %8.5f' % (prog_name, filename, dip_num, recv_num, delay_num, time)

			# Print out the values
			stars = '************************************************************************'
			out_txt = ''
			write_to_file = '%s/%s' % (out_dir, g_name)
			if LOUD: print'\n%s' % (stars)
			if LOUD: print '\n%s: Output added to %s from file %s\n' % (prog_name, g_name, filename)
			out_txt1 = 'time %9.5f Rx %d Dip %d' % (time, recv_num, dip_num)
			for pol in range(0,num_pols):
				out_txt2 = '%s%s Delay %d slot ' % (out_txt1, pol2XY(pol), delay_num )
				for slot in range(0,spr):
					out_txt = '%s%s%d tile %d ' % (out_txt, out_txt2, slot+1, slot2tile( slot+1, recv_num, Xnum, 0 ) )
					for i in range(0,num_ref_freqs):
						out_txt = '%s%2.2f ' % (out_txt, refpow[slot][pol][i] )
					out_txt = '%s\n' % (out_txt)
			if LOUD: print '%s' % (out_txt)
			if LOUD: print'%s' % (stars)

			# Save the results
			if( nowrite == 0 ):
				#Open the file
				if(index == 0): # If this is the first file, then the outfile needs to be created
					try:
			       			outfile = open(write_to_file, 'w')
					except IOError:
						print '%s: Can\'t open file %s for writing.  Aborting...' % (prog_name, write_to_file)
			       			sys.exit(0)
				else: # Otherwise, the outfile needs to be opened, and the new results appended to the end
					try:
			       			outfile = open(write_to_file, 'a')
					except IOError:
						print '%s: Can\'t open file %s for writing.  Aborting...' % (prog_name, write_to_file)
			       			sys.exit(0)
				# Write the data
				outfile.write(out_txt)
				# Close the file
				outfile.close

			out_txt = ''
			write_to_file = '%s/%s' % (out_dir, b_name)
			if LOUD: print'\n%s' % (stars)
			if LOUD: print '\n%s: Output added to %s from file %s\n' % (prog_name, b_name, filename)
			out_txt1 = 'time %9.5f Rx %d Dip %d' % (time, recv_num, dip_num)
			for slot in range(0,spr):
				for pol in range(0,num_pols):
					out_txt2 = '%s%s Delay %d slot ' % (out_txt1, pol2XY(pol), delay_num )
					value = bad2[slot][pol][index]
					if( value!=0 ):
						out_txt = '%s%s%d tile %d ' % (out_txt, out_txt2, slot+1, slot2tile( slot+1, recv_num, Xnum, 0 ) )
						if( value==1 ): # 'U' for Unknown (or, Unusual)
							out_txt = '%sU' % (out_txt)
						elif( value==2 ): # 'D' for Dead dipole
							out_txt = '%sD' % (out_txt)
						elif( value==3 ): # 'L' for Low signal
							out_txt = '%sL%2.1f' % (out_txt, bad_props2[slot][pol][index])
						elif( value==4 ): # 'H' for High signal
							out_txt = '%sH%2.1f' % (out_txt, bad_props2[slot][pol][index])
						elif( value==5 ): # 'F' for unusual noise Floor
							out_txt = '%sF%2.1f' % (out_txt, bad_props2[slot][pol][index])
						elif( value==6 ): # 'R' for Rfi spikes
							out_txt = '%sR %3.2f %3.2f' % (out_txt, bad_props2[slot][pol][index], bad_props3[slot][pol][index])
						out_txt = '%s\n' % (out_txt)
			if LOUD: print '%s' % (out_txt)
			if LOUD: print'%s' % (stars)

			# Save the results
			if( nowrite == 0 ):
				#Open the file
				if(index == 0): # If this is the first file, then the outfile needs to be created
					try:
			       			outfile = open(write_to_file, 'w')
					except IOError:
						print '%s: Can\'t open file %s for writing.  Aborting...' % (prog_name, write_to_file)
			       			sys.exit(0)
				else: # Otherwise, the outfile needs to be opened, and the new results appended to the end
					try:
			       			outfile = open(write_to_file, 'a')
					except IOError:
						print '%s: Can\'t open file %s for writing.  Aborting...' % (prog_name, write_to_file)
			       			sys.exit(0)
				# Write the data
				outfile.write(out_txt)
				# Close the file
				outfile.close
		# End, if( compare_option==0 )

	# End, cycle through files

	if( compare_option==1):

		first_time_through = 1
		if( use_T ):
			start_tile = use_tile_only
			stop_tile = use_tile_only+1
		else:
			start_tile = 1
			stop_tile = num_tiles
		for tilenum in range(start_tile,stop_tile):

			# Determine the slot number and receiver of this tile
			[ slot, recv ] = slot2tile( tilenum, -1, Xnum, 1 )
			# As a check, make sure that slot2tile outputs the tile number for this slot number and receiver
			tile = slot2tile( slot, recv, Xnum, 0 )
			if( tile != tilenum ):
				if LOUD: print '%s: ERROR: function slot2tile() does not give consistent results.  tilenum=%d.  output (mode 1): slot = %d, recv = %d.  output (mode 0) using these values: tile = %d.  Aborting...' % (prog_name, tilenum, slot, recv, tile)
				sys.exit(1)
			else:
				if LOUD: print '%s: Working on tile number %d, slot number %d, recv number %d...' % (prog_name, tilenum, slot, recv)

			rel_files = [] # Determine the relevent files for this recv
			for i in range(0,num_files):
				if( recvs[i]== recv ):
					rel_files.append(i)
			num_rel_files = len(rel_files)

			if( num_rel_files != 0 ): # If at least one file contains this receiver (and thus tile), then proceed...

				for delay in range(0,num_delays):
					if LOUD: print '\n%s: Working on delay %d...' % (prog_name, delay)

					rel_files2 = [] # Determine the relevent files for this recv
					for i in range(0,num_files):
						if( recvs[i]== recv and delays[i]==delay ):
							rel_files2.append(i)
					num_rel_files2 = len(rel_files2)

					if( num_rel_files2 != 0 ): # If at least one file contains this delay, then proceed...
		
						# SECTION IV: Make a quick recommendation based upon these values of the power

						if LOUD: print '\n%s: Making quick recommendations on this tile using %d files' % (prog_name, num_rel_files2)

						use_mod = 1
						avg_val = [ [ 0 for j in range(num_ref_freqs) ] for k in range(num_pols) ]
						std_val = [ [ 0 for j in range(num_ref_freqs) ] for k in range(num_pols) ]

						# Find dead dipoles (Search was also done across slots in the same receiver above: this helps in the case when the entire tile is dead.)
						[ avg_val, std_val ] = calc_avg_std2( range(slot-1,slot), range(0,num_pols), reffreq, rel_files2, allpows, bad2, use_mods )
						dead_cut = 3.0 # plots of dead_cut db or less in all three frequency channels are mostly likely dead
						for pol in range(0,num_pols):
							for files in range(0,num_rel_files2): # cycle through relevent files (which should be 16 files, one for each dipole)
								filenum = rel_files2[files]
								if( bad2[slot-1][pol][filenum]==0 ):
									dead = 0
									for i in range(0,num_ref_freqs):
										temp = allpows[slot-1][pol][i][filenum]
										avg = avg_val[pol][i]
										diff = math.fabs( temp - avg )
										if( (temp < 0.2*avg) or temp<dead_cut ):
											dead = dead + 1
											if LOUD: print '%s: dipole = %d%s for tile %d may be dead.  At freq = %3.3f MHz, power = %2.3f avg = %2.3f' % (prog_name,dips[filenum],pol2XY(pol),tilenum,reffreq[i],temp,avg)
									if( dead == num_ref_freqs ):
										bad2[slot-1][pol][filenum] = 2

						# Look for other peculiarities
						[ avg_val, std_val ] = calc_avg_std2( range(slot-1,slot), range(0,num_pols), reffreq, rel_files2, allpows, bad2, use_mods )
						for pol in range(0,num_pols):
							for files in range(0,num_rel_files2): # cycle through relevent files (which should be 16 files, one for each dipole)
								filenum = rel_files2[files]
								if( bad2[slot-1][pol][filenum] == 0 ):
									low = high = some_bad = 0
									low_avg = high_avg = 0.0
									for i in range(0,num_ref_freqs):
										temp = allpows[slot-1][pol][i][filenum]
										avg = avg_val[pol][i]
										#print 'avg = %2.2f, avg_val = %2.2f, pol = %d, files = %d, i = %d' % ( avg, avg_val[pol][i], pol, files, i )
										std = std_val[pol][i]
										diff = math.fabs( temp - avg )
										if( (diff>std) and (diff>3.0)  ):
											some_bad = 1
											if( temp < avg ):
												low = low + 1
												low_avg = low_avg + diff
											elif( temp > avg ):
												high = high + 1
												high_avg = high_avg + diff
											if LOUD: print '%s: dipole = %d%s for tile %d may be bad.  At freq = %3.3f MHz, power = %2.3f avg = %2.3f, diff = %2.3f, std dev = %2.3f' % (prog_name,dips[filenum], pol2XY(pol), tilenum, reffreq[i],temp,avg,diff,std)
									if( some_bad==1 ):
										if( low==num_ref_freqs ):
											bad2[slot-1][pol][filenum] = 3
											bad_props2[slot-1][pol][filenum] = low_avg/float(num_ref_freqs)
										elif( high==num_ref_freqs ):
											bad2[slot-1][pol][filenum] = 4
											bad_props2[slot-1][pol][filenum] = high_avg/float(num_ref_freqs)
										else:
											bad2[slot-1][pol][filenum] = 1


						# SECTION V: Print results to screen and save the results.

						# Print out the values
						stars = '************************************************************************'
						write_to_file = '%s/%s' % (out_dir, g_name)
						if LOUD: print'\n%s' % (stars)
						if LOUD: print '\n%s: Output added to %s for tile %d\n' % (prog_name, g_name, tilenum)
						out_txt = 'TILE %03d (RX %02d, SLOT %d) for DELAY %02d at TIME ~%9.3f\n' % (tilenum, recv, slot, delay, times[rel_files2[num_rel_files2/2]])
						for files in range(0,num_rel_files2): # cycle through relevent files (which should be 16 files, one for each dipole)
							filenum = rel_files2[files]
							dip_num = dips[filenum]
							#print filenames[filenum], dip_num
							for pol in range(0,num_pols):
								out_txt = '%sdipole %02d%s ' % (out_txt, dip_num, pol2XY(pol) )
								for i in range(0,num_ref_freqs):
									out_txt = '%s%2.2f ' % (out_txt, allpows[slot-1][pol][i][filenum] )
								out_txt = '%s\n' % (out_txt)
						out_txt = '%s\n' % (out_txt)
						if LOUD: print '%s' % (out_txt)
						if LOUD: print'%s' % (stars)

						# Save the results
						if( nowrite == 0 ):
							#Open the file
							if(first_time_through==1): # If this is the first file, then the outfile needs to be created
								try:
						       			outfile = open(write_to_file, 'w')
								except IOError:
									print '%s: Can\'t open file %s for writing.  Aborting...' % (prog_name, write_to_file)
						       			sys.exit(0)
							else: # Otherwise, the outfile needs to be opened, and the new results appended to the end
								try:
						       			outfile = open(write_to_file, 'a')
								except IOError:
									print '%s: Can\'t open file %s for writing.  Aborting...' % (prog_name, write_to_file)
						       			sys.exit(0)
							# Write the data
							outfile.write(out_txt)
							# Close the file
							outfile.close

						out_txt = ''
						write_to_file = '%s/%s' % (out_dir, b_name)
						if LOUD: print'\n%s' % (stars)
						if LOUD: print '\n%s: Output added to %s for tile %d\n' % (prog_name, b_name, tilenum)
						out_txt = 'TILE %03d (RX %02d, SLOT %d) for DELAY %02d at TIME ~%9.3f\n' % (tilenum, recv, slot, delay, times[rel_files2[num_rel_files2/2]])
						for files in range(0,num_rel_files2): # cycle through relevent files (which should be 16 files, one for each dipole)
							filenum = rel_files2[files]
							dip_num = dips[filenum]
							#print filenames[filenum], dip_num
							for pol in range(0,num_pols):
								value = bad2[slot-1][pol][filenum]
								if( value!=0 ):
									out_txt = '%sdipole %02d%s ' % (out_txt, dip_num, pol2XY(pol) )
									if( value==1 ): # 'U' for Unknown (or, Unusual)
										out_txt = '%sU' % (out_txt)
									elif( value==2 ): # 'D' for Dead dipole
										out_txt = '%sD' % (out_txt)
									elif( value==3 ): # 'L' for Low signal
										out_txt = '%sL%2.1f' % (out_txt, bad_props2[slot-1][pol][filenum])
									elif( value==4 ): # 'H' for High signal
										out_txt = '%sH%2.1f' % (out_txt, bad_props2[slot-1][pol][filenum])
									elif( value==5 ): # 'F' for unusual noise Floor
										out_txt = '%sF%2.1f' % (out_txt, bad_props2[slot-1][pol][filenum])
									elif( value==6 ): # 'R' for Rfi spikes
										out_txt = '%sR %3.2f %3.2f' % (out_txt, bad_props2[slot-1][pol][filenum], bad_props3[slot-1][pol][filenum])
									elif( value==7 ): # 'E' for Effective potential
										out_txt = '%sE %3.2f' % (out_txt, bad_props2[slot-1][pol][filenum])
										#print 'HERE! file = %s' % (filenames[filenum])
									out_txt = '%s\n' % (out_txt)
						out_txt = '%s\n' % (out_txt)
						if LOUD: print '%s' % (out_txt)
						if LOUD: print'%s' % (stars)

						# Save the results
						if( nowrite == 0 ):
							#Open the file
							if(first_time_through==1): # If this is the first file, then the outfile needs to be created
								try:
						       			outfile = open(write_to_file, 'w')
								except IOError:
									print '%s: Can\'t open file %s for writing.  Aborting...' % (prog_name, write_to_file)
						       			sys.exit(0)
								first_time_through=0
							else: # Otherwise, the outfile needs to be opened, and the new results appended to the end
								try:
						       			outfile = open(write_to_file, 'a')
								except IOError:
									print '%s: Can\'t open file %s for writing.  Aborting...' % (prog_name, write_to_file)
						       			sys.exit(0)
							# Write the data
							outfile.write(out_txt)
							# Close the file
							outfile.close



					else: # no files contain this tile number and delay
						if LOUD: print '%s: No input files contain this tile number and delay.' % (prog_name)

				# End, cycle through delays

			else: # no files contain this tile number
				if LOUD: print '%s: No files contain tile number %d' % (prog_name, tilenum)

		# End, cycle through tiles

	# End, if( compare_option==1 )



		

	print '\n%s: Program complete.\n' % (prog_name)


