import pdb
#from pandexo_prep import tepcat, jwstsim
import tepcat, jwstsim

#inst_modes = 'all' # will run on all modes
inst_modes = ['NIRSpec G395H'] # just test one for now
z = tepcat.load()
allnames = z['names']
allnames = [ 'WASP-121', 'WASP-019' ] # temporary
nplanets = len( allnames )
for i in range( nplanets ):
    print( '\nPlanet {0} of {1}'.format( i+1, nplanets ) )
    jwstsim.main( allnames[i], z, inst_modes=inst_modes )

