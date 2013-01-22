PRO observation_setup,filename=filename,data_directory=data_directory
IF N_Elements(data_directory) EQ 0 THEN data_directory='DATA\r4\clmw\X14\PicA_121_20100924211938'
IF N_Elements(filename) EQ 0 THEN filename='PicA_121_20100924211938.cal'
;filename='PicA_121_20100924211938'
ext='.UVFITS'
header_filename=filename+'_header'
tile_gain_x_filename=filename+'_tile_gains_x'
tile_gain_y_filename=filename+'_tile_gains_y'
lat=-26.703319 ;degrees
lon=116.67081 ;degrees
alt=377.83 ;meters

ntiles=32.
nfreq_bin=24.
gain_array=fltarr(17,nfreq_bin*ntiles)+1.
gain_array[0,*]=Floor(indgen(nfreq_bin*ntiles)/nfreq_bin)+1

tile_gain_x_filepath=filepath(tile_gain_x_filename,root_dir=rootdir('mwa'),subdir=data_directory)
tile_gain_y_filepath=filepath(tile_gain_y_filename,root_dir=rootdir('mwa'),subdir=data_directory)

;do not overwrite a gain_array if one already exists (it's either real data, or the same default data as this!)
IF file_test(tile_gain_x_filepath) EQ 0 THEN $
    textfast,gain_array,/write,filename=tile_gain_x_filename,root=rootdir('mwa'),filepathfull=data_directory
IF file_test(tile_gain_y_filepath) EQ 0 THEN $
    textfast,gain_array,/write,filename=tile_gain_y_filename,root=rootdir('mwa'),filepathfull=data_directory

i_test=Strpos(filename,'.cal')
IF i_test EQ -1 THEN filename2=filename ELSE filename2=Strmid(filename,0,i_test)
file_path=filepath(filename2,root_dir=rootdir('mwa'),subdir=data_directory)
data_struct=mrdfits(file_path+'.cal'+ext,0,data_header,/silent)
gcount=sxpar(data_header,'gcount') 
n_grp_params=sxpar(data_header,'pcount')
;    uu_arr=reform(data_struct.params[0])
;    vv_arr=reform(data_struct.params[1])
;    ww_arr=reform(data_struct.params[2])
;    baseline_arr=reform(data_struct.params[3])
;    date_arr=reform(data_struct.params[4])
uvw_baseline=fltarr(n_grp_params,gcount)
FOR i=0,n_grp_params-1 DO uvw_baseline[i,*]=data_struct.params[i]
sxaddpar,data_header,'lat',lat,'degrees',after='INSTRUME'
sxaddpar,data_header,'lon',lon,'degrees',after='lat'
sxaddpar,data_header,'alt',alt,'meters',after='lon'

FitsFast,uvw_baseline,data_header,/write,filename=header_filename,root=rootdir('mwa'),filepathfull=data_directory
;FitsFast,gain_array,data_header,/write,filename=filename+'_tile_gains_y',root=rootdir('mwa'),filepathfull=data_directory

END
