
pro waterfall_script

  froot = '/Users/bryna/Documents/Physics/MWA_data/'
  cc_file = froot + 'alldas.lccspc'

  cc_arr  = read_cross(cc_file, nchans = 768)

  index13 = ccindex(13, cc_arr, labels = labels13)

  index = index13[where(labels13 eq 55)]

  wf_arr = reform(cc_arr[*, index, *])

  loadct, 0
  waterfall, atan(wf_arr,/phase)



  loadct, 39
end
