;base_path = '/home/bryna/Documents/Physics/idl_working/'
;base_path = '/Users/bryna/Documents/Physics/bryna_svn/idl_working/'
git_path = '/Users/bryna/Documents/Physics/hazelton_git/idl_stuff/'

path_dirs = git_path + ['', 'idl_utilities', 'idl_utilities/fitting_functions','coyote/' + ['', 'public/'], 'astron/pro/' + ['','jhuapl/'], $
  'mwa', 'power_spectrum', $;;'fhd_sims', $
  'single_use', 'Healpix_3.11/src/idl/' + ['', 'examples', 'fits', 'interfaces', 'misc', 'toolkit', 'visu', 'ximview/' + $
  ['', 'docs', 'gscroll', 'hpx', 'utilities'], 'zzz_external/' + ['cgis', 'obsolete_astron']]] ;; , $
;; 'UCSC/' + ['', 'ADELE', 'A0535', 'crab', 'preflight', 'solarfss/' + ['', 'visibilities', 'sas_temp'], $
;;           'TGF_science/' + ['', 'land_sea/' + ['', 'lis'], 'wwlln', 'geant/' + ['', 'deadtime']]]]
  
;fhdps_utils_path = '/Users/bryna/Documents/Physics/fhdps_utils/'
;fhdps_utils_dirs = fhdps_utils_path
  
ps_path_dirs = '/Users/bryna/Documents/Physics/PS/' + ['', 'fhdps_utils','ps_utils', 'ps_core', 'ps_wrappers', 'textoidl']

fhd_path = '/Users/bryna/Documents/Physics/FHD/'
fhd_path_dirs = fhd_path + ['', 'fhd_utils', 'fhd_output', 'fhd_core', 'Observations', 'catalog_data', 'instrument_config', $
    'simulation', 'simulation/simulation_wrappers']

;path_dirs = [fhdps_utils_dirs, ps_path_dirs, path_dirs, fhd_path_dirs]
path_dirs = [ps_path_dirs, path_dirs, fhd_path_dirs]
path_string = strjoin(path_dirs, ':') + ':'

;!path = path_string + ':/home/dsmith/rhessi/code:' + !path
!path = path_string + !path
;search_network, /enable

init_healpix
defsysv, '!HEALPIX', exist = hpx_is_defined
if (hpx_is_defined) then print, 'environment variable HEALPIX exists'

device, decomposed = 0
loadct,39
!p.background = 255
!p.color = 0

;; added this to prevent BadMatch errors on Lion, see coyote tip: http://www.idlcoyote.com/misc_tips/badmatch.php
Device, RETAIN=2

