git_path = '/Users/bryna/Projects/Physics/random_stuff/idl_stuff/'

path_dirs = git_path + ['', 'idl_utilities', 'idl_utilities/fitting_functions',$
                        'coyote/' + ['', 'public/'], 'astron/pro/' + ['','jhuapl/'], $
                        'mwa', 'power_spectrum', 'single_use', $
                        'Healpix_3.11/src/idl/' + ['', 'examples', 'fits', 'interfaces', $
                                                   'misc', 'toolkit', 'visu', 'ximview/' + $
                                                   ['', 'docs', 'gscroll', 'hpx', 'utilities'], $
                                                   'zzz_external/' + ['cgis', 'obsolete_astron']]]
ps_path = '/Users/bryna/Projects/Physics/eppsilon/'
;; get subdirectories
ps_path_dirs = expand_path('+' + ps_path)

fhd_path = '/Users/bryna/Projects/Physics/FHD/'
;; get subdirectories
fhd_path_dirs = expand_path('+' + fhd_path)

pipeline_scripts_path = '/Users/bryna/Projects/Physics/pipeline_scripts'
pipeline_scripts_dirs = expand_path('+' + pipeline_scripts_path)

fhdps_utils_dirs = '/Users/bryna/Projects/Physics/fhdps_utils'

path_dirs = [ps_path_dirs, fhdps_utils_dirs, path_dirs, fhd_path_dirs, pipeline_scripts_dirs]
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
