;base_path = '/home/bryna/Documents/Physics/idl_working/'
;base_path = '/Users/bryna/Documents/Physics/bryna_svn/idl_working/'
base_path = '/Users/bryna/Documents/Physics/bryna_git/idl_stuff/'

path_dirs = base_path + ['', 'idl_utilities', 'idl_utilities/fitting_functions','coyote', 'astron/pro', 'mwa', $
                         'power_spectrum', 'fhd_sims', 'fhd_sims/simulations', 'single_use', 'textoidl', $
                         'Healpix_2.15a/src/idl/' + ['', 'examples', 'fits', 'interfaces', 'misc', 'toolkit', 'visu', 'ximview/' + $
                                                     ['', 'docs', 'gscroll', 'hpx', 'utilities'], 'zzz_external/' + $
                                                     ['cgis', 'obsolete_astron']]] ;; , $
                         ;; 'UCSC/' + ['', 'ADELE', 'A0535', 'crab', 'preflight', 'solarfss/' + ['', 'visibilities', 'sas_temp'], $
                         ;;           'TGF_science/' + ['', 'land_sea/' + ['', 'lis'], 'wwlln', 'geant/' + ['', 'deadtime']]]]

ps_path = '/Users/bryna/Documents/Physics/PS/' + ['', 'ps_utils', 'power_spectrum']

ian_path = '/Users/bryna/Documents/Physics/RadCos/sullivan/'
ian_path_dirs = ian_path + ['FHD', 'fhd\ utils', 'fhd_output', 'mwa\ analysis', 'single_use', 'utilities']

path_dirs = [path_dirs, ian_path_dirs]
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
