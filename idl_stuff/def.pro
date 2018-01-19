git_path = '/Users/bryna/Projects/Physics/hazelton_git/idl_stuff/'

path_dirs = git_path + ['', 'idl_utilities', 'idl_utilities/fitting_functions',$
                        'coyote/' + ['', 'public/'], 'astron/pro/' + ['','jhuapl/'], $
                        'mwa', 'power_spectrum', 'single_use', $
                        'Healpix_3.11/src/idl/' + ['', 'examples', 'fits', 'interfaces', $
                                                   'misc', 'toolkit', 'visu', 'ximview/' + $
                                                   ['', 'docs', 'gscroll', 'hpx', 'utilities'], $
                                                   'zzz_external/' + ['cgis', 'obsolete_astron']]]
ps_path_dirs = '/Users/bryna/Projects/Physics/eppsilon/' + $
  ['', 'fhdps_utils', 'ps_core', 'ps_compare', 'ps_plotting', $
  'ps_setup', 'ps_utils', 'ps_wrappers', 'textoidl']

;; add plus sign to get subdirectories. Don't want fhdps_utils added from here
;; so don't include the root directory in that part
fhd_path = '/Users/bryna/Projects/Physics/FHD/'
fhd_path_dirs = fhd_path + ['catalog_data', 'fhd_core', 'fhd_output', 'fhd_utils', $
                            'instrument_config', 'Observations',  'simulation']

temp = fhd_path
for i=0, n_elements(fhd_path_dirs)-1 do temp = [temp, expand_path('+' + fhd_path_dirs[i])]
fhd_path_dirs = temp

path_dirs = [ps_path_dirs, path_dirs, fhd_path_dirs]
path_string = strjoin(path_dirs, ':') + ':'

!path = path_string + !path

init_healpix
defsysv, '!HEALPIX', exist = hpx_is_defined
if (hpx_is_defined) then print, 'environment variable HEALPIX exists'

device, decomposed = 0
loadct,39
!p.background = 255
!p.color = 0

;; added this to prevent BadMatch errors on Lion,
;; see coyote tip: http://www.idlcoyote.com/misc_tips/badmatch.php
Device, RETAIN=2

;; added this to get imagemagick to be found by IDL
;; took it out in Jan 2018 because it broke spawn commands and doesn't seem to be needed for imagemagick
;; SetEnv, 'PATH=/opt/local/bin:$PATH'
