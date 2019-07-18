import pdb
from obsplanning.pandexo_prep_dev import tepcat, hstsim
#import tepcat, jwstsim

"""
Copy this script to your working directory. Edit the import
statements to find the pandexo_prep routines on your system.

Currently, you also need to download the following TEPCat
catalogues and put them in the working directory:
  'For planning observations' --> save as 'tepcat1.txt'
  'Well-studied transiting planets' --> save as 'tepcat2.txt'

Then there are two steps to edit below:
  1. Specify instrument modes by setting 'inst_modes' variable.
  2. Specify planets by setting 'planets' variable.
"""

z = tepcat.load( download_latest=False ) # load the TEPCat catalogues

# 1. Specify a list of strings containing the JWST instrument
#    modes that are recognised by PandExo, 
#    e.g. [ 'NIRSpec G395H', 'NIRSpec G235H', 'NIRISS SOSS' ] 
#    or set to 'all' to run calculations for all instrument modes.
inst_modes = 'all'

# 2. Define a dictionary containing the properties of the star and planet:

planetdict = { 'name':'HD209458b', 'star':{}, 'planet':{} }

planetdict['star']['type'] = 'phoenix' # use the provided Phoenix spectra 
planetdict['star']['mag'] = 6.31 # stellar magnitude
planetdict['star']['ref_wave'] = 2.22 # K band central wavelength in micron
planetdict['star']['temp'] = 6080 # stellar effective temperature in Kelvin
planetdict['star']['metal'] = 0 # stellar metallicity as log10( [Fe/H] )
planetdict['star']['logg'] = 4.4 # stellar surface gravity as log10( g (c.g.s.) )
planetdict['star']['radius'] = 1.19
planetdict['star']['r_unit'] = 'R_sun' # stellar radius unit

planetdict['planet']['transit_duration'] = 0.127 # transit duration in days
planetdict['planet']['depth'] = ( 0.121 )**2.
planetdict['planet']['period'] = 3.5247 # planet orbital period in days
planetdict['planet']['i'] = 86.7
planetdict['planet']['ars'] = 8.8
planetdict['planet']['ecc'] = 0.
planetdict['planet']['w'] = 90.
planetdict['planet']['radius'] = 1.4
planetdict['planet']['r_unit'] = 'R_jup' # planet radius unit
planetdict['planet']['w_unit'] = 'um' # wavelength unit is micron; other options include 'Angs', secs" (for phase curves)
planetdict['planet']['f_unit'] = 'rp^2/r*^2' # options are 'rp^2/r*^2' or 'fp/f*'
planetdict['planet']['td_unit'] = 'd' # transit duration unit

scan = 'Round Trip' # 'Round Trip' or 'Forward'
nchan = 14
useFirstOrbit = False

##########################
# Can only do one planet at a time currently because for HST calculations
# PandExo requires various planet properties not provided by TEPCat; one
# possibly solution is to switch from TEPCat to NASA Exoplanet Archive:
hstsim.main( planetdict, scan=scan, nchan=nchan, inst_modes=inst_modes, useFirstOrbit=useFirstOrbit )

