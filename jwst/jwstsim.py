from __future__ import print_function
import os, pdb, sys, time
import numpy as np
import pandexo.engine.justdoit as jdi
import pandexo.engine.justplotit as jpi
import matplotlib.pyplot as plt



def main( planet_label, tepcat, sat_level=80, sat_unit='%', noise_floor_ppm=20, inst_modes='all', outdir='.' ):
    """
    Routine called by the run_jwst.py script to run PandExo over specified modes for specified planet.
    """
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
    if ix.sum()==0:
        print( 'Could not match {0} to any TEPCat planets - skipping'.format( planet_label ) )
        return None

    # Prepare the PandExo inputs:
    z = jdi.load_exo_dict()
    z['observation']['sat_level'] = sat_level # default to 80% for basic run
    z['observation']['sat_unit'] = sat_unit
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
    z['star']['r_unit'] = 'R_sun' # stellar radius unit
    z['planet']['r_unit'] = 'R_jup' # planet radius unit
    z['planet']['w_unit'] = 'um' # wavelength unit is micron; other options include 'Angs', secs" (for phase curves)
    z['planet']['f_unit'] = 'fp/f*' # options are 'rp^2/r*^2' or 'fp/f*'
    z['planet']['transit_duration'] = float( tepcat['tdurs'][ix] )*24.*60.*60. # transit duration in seconds
    z['planet']['td_unit'] = 's' # transit duration unit
    
    # Use a null spectrum for the planet:
    z['planet']['exopath'] = get_nullpath()
    if os.path.isfile( z['planet']['exopath'] )==False:
        generate_nullspec()

    # Run PandExo over requested instrument modes:
    if inst_modes=='all':
        inst_modes = list( jdi.ALL.keys() )
        inst_modes.remove( 'WFC3 G141' ) # remove HST modes
    if 'WFC3 G141' in inst_modes:
        pdb.set_trace() # this module is for JWST only, i.e. not WFC3 which requires batman etc.
        # todo = write another module called wfc3sim.py to handle WFC3, or generalise this
        # module to handle JWST and WFC3 and rename it noisesim.py or something.
    nmodes = len( inst_modes )

    if nmodes==1:
        modestr = '{0} instrument mode:\n'.format( nmodes )
    else:
        modestr = '{0} instrument modes:\n'.format( nmodes )
    for m in inst_modes: modestr += '{0}, '.format( m )
    print( '\n{0}\nRunning PandExo for {1}\n{2}\n{0}\n'.format( 50*'#', planet_label, modestr[:-2] ) )

    for k in range( nmodes ):
        #G140 has two filters, need to run the second (non-default) manually
        if 'G140' in inst_modes[k]:
            #G140 non-default:
            g140 = jdi.load_mode_dict(inst_modes[k])
            g140['configuration']['instrument']['filter'] = 'f100lp'
            oname = '{0}-F100LP.txt'.format( inst_modes[k].replace( ' ', '-' ) )
            oname_obs = '{0}-F100LP.obspar.txt'.format( inst_modes[k].replace( ' ', '-' ) )
            opath = os.path.join( odirfull, oname )
            opath_obs = os.path.join( odirfull, oname_obs)
            y = jdi.run_pandexo( z, g140, save_file=False )
            wav = y['FinalSpectrum']['wave']
            err = y['FinalSpectrum']['error_w_floor']*( 1e6 )
            outp = np.column_stack( [ wav, err ] )
            np.savetxt( opath, outp )
            #print( '\nSaved noise: {0}\n{1}\n'.format( opath, 50*'#' ) )
            print( '\nSaved noise: {0}'.format( opath, 50*'#' ) )
            save_obspar( opath_obs, y )
            #G140 prepare default:
            oname = '{0}-F070LP.txt'.format( inst_modes[k].replace( ' ', '-' ) )
            oname_obs = '{0}-F070LP.obspar.txt'.format( inst_modes[k].replace( ' ', '-' ) )
        else:
            oname = '{0}.txt'.format( inst_modes[k].replace( ' ', '-' ) )
            oname_obs = '{0}.obspar.txt'.format( inst_modes[k].replace( ' ', '-' ) )
        opath = os.path.join( odirfull, oname )
        opath_obs = os.path.join( odirfull, oname_obs)
        y = jdi.run_pandexo( z, [inst_modes[k]], save_file=False )
        wav = y['FinalSpectrum']['wave']
        err = y['FinalSpectrum']['error_w_floor']*( 1e6 )
        outp = np.column_stack( [ wav, err ] )
        np.savetxt( opath, outp )
        print( '\nSaved noise: {0}'.format( opath, 50*'#' ) )
        save_obspar( opath_obs, y )

    t2 = time.time()
    print( 'Total time taken = {0:.2f} minutes'.format( (t2-t1)/60. ) )
    return None


def save_obspar( opath, y ):
    with open( opath, 'w' ) as f:
        if y['warning']['Saturated?']=='All good':
            f.write( 'NOT SATURATED\n' )
        else:
            f.write( 'SATURATED\n(see warning below)\n' )
        f.write( '\nWARNINGS\n{0}\n'.format( 50*'-' ) )
        for key, value in y['warning'].items():
            f.write( '*** {0}\n--> {1}\n'.format( key, value ) )
        f.write( '\nSETUP\n{0}'.format( 50*'-' ) )
        ikeys = [ 'Instrument', 'Mode', 'Aperture', 'Disperser', 'Subarray', 'Readmode', 'Filter' ]
        for key in ikeys:
            f.write( '{0}:  {1}\n'.format( key, y['input'][key] ) )
        f.write( '\nEXPOSURES\n{0}\n'.format( 50*'-' ) )
        for key, value in y['timing'].items():
            f.write('{0}:  {1}\n'.format(key, value))
    f.close()
    print( '\nSaved observation parameters: {0}\n{1}\n'.format( opath, 50*'#' ) )
    return None
    


def generate_nullspec():
    """
    Generates a null spectrum and saves it in the working directory.
    """
    n = int( 1e4 )
    x = np.linspace( 0.2, 40, n )
    y = np.zeros( n )
    nullpath = get_nullpath()
    np.savetxt( nullpath, np.column_stack( [ x, y ] ) )
    return nullpath

def get_nullpath():
    """
    Returns the path of the null spectrum.
    """
    return os.path.join( os.getcwd(), 'zeros.txt' )
