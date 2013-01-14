import time
import os


prefix='PicA_121_20100924211938'

datapath = '/Users/bryna/Documents/Physics/MWA/data/'
savepath = datapath + prefix + '/'
fitsdata = datapath + prefix + '.UVFITS'

if not os.path.exists(savepath):
    os.mkdir(savepath)

# Set up some useful variables (these will also be set later on)
msfile = savepath + prefix + '.ms'
btable = savepath + prefix + '.bcal'
gtable = savepath + prefix + '.gcal'
ftable = savepath + prefix + '.fluxscale'
splitms = savepath + prefix + '.src.split.ms'
imname = savepath + prefix + '.cleanimg'

exportcalfits = savepath + prefix + '.cal.uvfits'


print '--Import--'

# Safest to start from task defaults
default('importuvfits')

# Set up the MS filename and save as new global variable
#msfile = prefix + '.ms'

# Use task importuvfits
fitsfile = fitsdata
vis = msfile

saveinputs('importuvfits',savepath + prefix + '.importuvfits.saved')

importuvfits()

listobs()

 ## To get the Az/El for the observation, used times out of listobs() with single observation:
 ## single_observation.py --starttime='2010-09-24,21:19:42' --stoptime='2010-09-24,21:24:34' --source='PicA' --nodatabase --creator='bryna' --frequencies=100,24 --useazel
 ## Result for PicA: (05:19:49.67, -45:46:44.40)
 ## Observation from 969398392 to 969398688
 ## Calculations for PicA from 969398392 to 969398688; shift every 296 s
 ## WARNING: No mode specified: using default
 ## GPSsec[s]	Az[deg]		El[deg]
 ## 969398392	181:36:43	+70:55:52
 ## WARNING: No gain_control_value specified: using default

## then used idl routine stokes_off_zenith to get Stokes parameters for unpolarized source at field center.
## az_deg = 181d + 36/60d + 43/3600d
## el_deg = 70d + 55/60d + 52/3600d
## stokes = stokes_off_zenith(az_deg, el_deg)
## No Stokes at zenith provided, assuming 1 Jy unpolarized.
## unpolarized fraction:        1
## Ex_mag, Ey_mag, Ez_mag, Emag: 
##       0.70707692      0.66833691      0.23101515       1.0000000
## Stokes I, Q, U, V
##       0.94663200     0.053283540       0.0000000      -0.0000000


print '--Setjy--'
default('setjy')

vis = msfile

fluxdensity = [0.94663200,     0.053283540,       0.0000000,      0.0000000]

saveinputs('setjy',savepath + prefix + '.setjy.saved')
setjy()


print '--Bandpass--'
default('bandpass')

vis = msfile
caltable = btable

refant='Tile01'
saveinputs('bandpass',savepath + prefix + '.bandpass.saved')
bandpass()



print '--Plotcal (bandpass)--'
default('plotcal')

caltable = btable
xaxis='chan'
yaxis='amp'
iteration='antenna'

plotcal()

yaxis='phase'
plotcal()


# time dependent calibration -- usually small for MWA, can skip
## print '--Gaincal--'
## default('gaincal')

## vis = msfile
## caltable = gtable
## gaintable = btable
## solint = '16s'

## gaincal()

## inp plotcal

## xaxis='time'
## plotcal()



print '--ApplyCal--'
default('applycal')

vis = msfile

gaintable = btable

# if did gaincal:
## gaintable.append(gtable)

applycal()



default('exportuvfits')

vis = msfile
fitsfile = exportcalfits

exportuvfits()
