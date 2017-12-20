import pdb
#from pandexo_prep_dev.jwst import tepcat, jwstsim
import tepcat, jwstsim

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

z = tepcat.load() # load the TEPCat catalogues

# 1. Specify a list of strings containing the JWST instrument
#    modes that are recognised by PandExo, 
#    e.g. [ 'NIRSpec G395H', 'NIRSpec G235H', 'NIRISS SOSS' ] 
#    or set to 'all' to run calculations for all instrument modes.
inst_modes = 'all'

# 2. Specify a list of strings containing the names of the 
#    systems to run calculations for, in TEPCat format, 
#    e.g. [ 'WASP-121', 'WASP-109' ]
planets = z['names'] # this will do all TEPCat planets

##########################
# Below here is automatic:
nplanets = len( planets )
for i in range( nplanets ):
    print( '\nPlanet {0} of {1}'.format( i+1, nplanets ) )
    jwstsim.main( planets[i], z, inst_modes=inst_modes )

