import pdb
import tepcat, jwstsim

inst_modes = ['NIRSpec G395H'] # just test one for now

z = tepcat.load()
allnames = z['names']


planet_label = 'WASP-121'
jwstsim.main( planet_label, z, inst_modes=inst_modes )

