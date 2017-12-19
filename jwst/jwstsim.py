from __future__ import print_function
import os, pdb, sys, time
import numpy as np
import pandexo.engine.justdoit as jdi
import pandexo.engine.justplotit as jpi
import matplotlib.pyplot as plt



def main( planet_label, tepcat, sat_level=80, sat_unit='%', noise_floor_ppm=20, inst_modes='all', outdir='.' ):
    t1 = time.time()

    # Prepare the output directory:
    cwd = os.getcwd()
    if outdir=='.':
        outdir = cwd
    odirfull = os.path.join( outdir, planet_label )
    if os.path.isdir( odirfull )==False:
        os.makedirs( odirfull )

    # Identify the tepcat index:
    ix = ( tepcat['names']==planet_label )

    # Prepare the PandExo inputs:
    z = jdi.load_exo_dict()
    z['observation']['sat_level'] = 80 # saturation level
    z['observation']['sat_unit'] = '%' # saturation level unit
    z['observation']['noccultations'] = 1 # number of transits
    z['observation']['R'] = None # do not pre-bin the output spectra
    z['observation']['baseline'] = 1.0 # out-of-transit baseline quantity
    z['observation']['baseline_unit'] ='frac' # baseline quantity is fraction of time in transit versus out
    z['observation']['noise_floor'] = noise_floor_ppm # noise floor in p.p.m.
    z['star']['type'] = 'phoenix' # use the provided Phoenix spectra 
    z['star']['mag'] = float( tepcat['kmags'][ix] ) # stellar magnitude
    z['star']['ref_wave'] = 2.22 # K band central wavelength in micron
    z['star']['temp'] = float( tepcat['tstar'][ix] ) # stellar effective temperature in Kelvin
    z['star']['metal'] = float( tepcat['metalstar'][ix] ) # stellar metallicity as log10( [Fe/H] )
    z['star']['logg'] = float( tepcat['loggstar'][ix] ) # stellar surface gravity as log10( g (c.g.s.) )
    z['planet']['w_unit'] = 'um' # wavelength unit is micron; other options include 'Angs', secs" (for phase curves)
    z['planet']['f_unit'] = 'fp/f*'               #options are 'rp^2/r*^2' or 'fp/f*'
    z['planet']['transit_duration'] = float( tepcat['tdurs'][ix] )*24.*60.*60.   #transit duration in seconds
    # Use a null spectrum for the planet:
    z['planet']['exopath'] = get_nullpath()
    if os.path.isfile( z['planet']['exopath'] )==False:
        generate_nullspec()

    # Run PandExo over requested instrument modes:
    if inst_modes=='all':
        inst_modes = list( jdi.ALL.keys() )
    inst_modes.remove( 'WFC3 G141' ) # remove HST modes
    nmodes = len( inst_modes )

    if nmodes==1:
        modestr = '{0} instrument mode:\n'.format( nmodes )
    else:
        modestr = '{0} instrument modes:\n'.format( nmodes )
    for m in inst_modes: modestr += '{0}, '.format( m )
    print( '\n{0}\nRunning PandExo for {1}\n{2}\n{0}\n'.format( 50*'#', planet_label, modestr[:-2] ) )
    #print( '{0}\n'.format( modesstr[:-2] ) )
    for k in range( nmodes ):
        oname = '{0}.txt'.format( inst_modes[k].replace( ' ', '-' ) )
        opath = os.path.join( odirfull, oname )
        y = jdi.run_pandexo( z, [inst_modes[k]], save_file=False )
        wav = y['FinalSpectrum']['wave']
        err = y['FinalSpectrum']['error_w_floor']*( 1e6 )
        outp = np.column_stack( [ wav, err ] )
        np.savetxt( opath, outp )
        print( '\nSaved noise: {0}'.format( opath ) )
    t2 = time.time()
    print( 'Total time taken = {0:.2f} minutes'.format( (t2-t1)/60. ) )
    return None


def generate_nullspec():
    n = int( 1e4 )
    x = np.linspace( 0.2, 40, n )
    y = np.zeros( n )
    nullpath = get_nullpath()
    np.savetxt( nullpath, np.column_stack( [ x, y ] ) )
    return nullpath

def get_nullpath():
    return os.path.join( os.getcwd(), 'zeros.txt' )
