from __future__ import print_function
import os, pdb, sys, time
import numpy as np
import pandexo.engine.justdoit as jdi
import pandexo.engine.justplotit as jpi
import matplotlib.pyplot as plt



def main( planetdict, sat_level=80, sat_unit='%', noise_floor_ppm=20, useFirstOrbit=False, \
          inst_modes=['WFC3 G141'], scan='Round Trip', nchan=14, outdir='.' ):
    """
    Routine called by the run_jwst.py script to run PandExo over specified 
    modes for specified planet.
    """
    t1 = time.time()

    # Prepare the output directory:
    cwd = os.getcwd()
    if outdir=='.':
        outdir = cwd
    planet_label = planetdict['name']        
    odirfull = os.path.join( outdir, planet_label )
    if os.path.isdir( odirfull )==False:
        os.makedirs( odirfull )

    # Identify the tepcat index:
    #ix = ( tepcat['names']==planet_label )
    #if ix.sum()==0:
    #    print( 'Could not match {0} to any TEPCat planets - skipping'.format( planet_label ) )
    #    return None

    # Prepare the PandExo inputs:
    z = jdi.load_exo_dict()
    z['observation']['sat_level'] = sat_level # default to 80% for basic run
    z['observation']['sat_unit'] = sat_unit
    z['observation']['noccultations'] = 1 # number of transits
    z['observation']['R'] = None # do not pre-bin the output spectra
    z['observation']['baseline'] = 1.0 # out-of-transit baseline quantity
    z['observation']['baseline_unit'] ='frac' # baseline quantity is fraction of time in transit versus out
    z['observation']['noise_floor'] = noise_floor_ppm # noise floor in p.p.m.
    z['star'] = planetdict['star']
    z['planet'] = planetdict['planet']
    z['planet']['type'] = 'constant'
    HST_ORB_PERIOD_DAYS = 96./60./24.
    tobs_days = 2*z['planet']['transit_duration']
    norb = int( np.ceil( tobs_days/HST_ORB_PERIOD_DAYS ) )
    
    # Use a null spectrum for the planet:
    z['planet']['exopath'] = get_nullpath()
    if os.path.isfile( z['planet']['exopath'] )==False:
        generate_nullspec()

    # Run PandExo over requested instrument modes:
    if inst_modes=='all':
        inst_modes = [ 'WFC3 G141' ]  # G102 not implemented?
    nmodes = len( inst_modes )

    if nmodes==1:
        modestr = '{0} instrument mode:\n'.format( nmodes )
    else:
        modestr = '{0} instrument modes:\n'.format( nmodes )
    for m in inst_modes: modestr += '{0}, '.format( m )
    print( '\n{0}\nRunning PandExo for {1}\n{2}\n{0}\n'.format( 50*'#', planet_label, modestr[:-2] ) )

    for k in range( nmodes ):
        s1 = inst_modes[k].replace( ' ', '-' )
        if useFirstOrbit==True:
            s2 = 'keepFirstOrbit'
        else:
            s2 = 'dropFirstOrbit'
        if scan=='Round Trip':
            s3 = 'RTscan'
        elif scan=='Forward':
            s3 = 'Fscan'
        else:
            pdb.set_trace()
        oname = '{0}.{1}.{2}.txt'.format( s1, s2, s3 )
        oname_obs = oname.replace( '.txt', '.obspar.txt' )
        opath = os.path.join( odirfull, oname )
        opath_obs = os.path.join( odirfull, oname_obs)
        wfc3 = jdi.load_mode_dict( inst_modes[k] )
        wfc3['strategy']['calculateRamp'] = useFirstOrbit
        if useFirstOrbit==True:
            wfc3['strategy']['norbits'] = norb
        else:
            wfc3['strategy']['norbits'] = norb+1
        wfc3['strategy']['nchan'] = nchan
        wfc3['strategy']['schedulability'] = 100
        wfc3['strategy']['scanDirection'] = scan
        y = jdi.run_pandexo( z, wfc3, save_file=False )
        #pdb.set_trace()
        #y = jdi.run_pandexo( z, [inst_modes[k]], save_file=False )
        wav = y['planet_spec']['binwave']
        err = (1e6)*y['planet_spec']['error']*np.ones( nchan )
        outp = np.column_stack( [ wav, err ] )
        np.savetxt( opath, outp )
        print( '\nSaved noise: {0}'.format( opath, 50*'#' ) )
        save_obspar( opath_obs, y )

    t2 = time.time()
    print( 'Total time taken = {0:.2f} minutes'.format( (t2-t1)/60. ) )
    return None


def save_obspar( opath, y ):
    #print( y.keys() )
    #pdb.set_trace()
    with open( opath, 'w' ) as f:
        f.write( '\nSETUP\n{0}\n'.format( 50*'-' ) )
        f.write( 'Number of HST orbits = {0:.0f}\n'\
                 .format( int( y['wfc3_TExoNS']['info']['Number of HST orbits'] ) ) )
        f.write( 'Use first orbit = {0}\n'\
                 .format( y['wfc3_TExoNS']['info']['Use first orbit']  ) )
        f.write( 'NSAMP = {0}\n'\
                 .format( int( y['wfc3_TExoNS']['info']['WFC3 parameters: NSAMP'] ) ) )
        f.write( 'SAMP_SEQ = {0}\n'\
                 .format( y['wfc3_TExoNS']['info']['WFC3 parameters: SAMP_SEQ'] ) )
        f.write( 'exposure time = {0:.0f} (s)\n'\
                 .format( float( y['wfc3_TExoNS']['info']['exposure time'] ) ) )
        scanrate = y['wfc3_TExoNS']['info']['Recommended scan rate (arcsec/s)']
        f.write( 'Recommended scan rate = {0:.4f} (arcsec/s)\n'\
                 .format( float( scanrate ) ) )
        f.write( 'Scan height = {0:.0f} pixels\n'\
                 .format( float( y['wfc3_TExoNS']['info']['Scan height (pixels)'] ) ) )
        f.write( 'Maximum fluence = {0:.0f} (electrons)\n'\
                 .format( float( y['wfc3_TExoNS']['info']['Maximum pixel fluence (electrons)'] ) ) )
        f.write( 'Estimated duty cycle = {0:.2f} %\n'\
                 .format( float( y['wfc3_TExoNS']['info']['Estimated duty cycle (outside of Earth occultation)'] ) ) )
        f.write( 'Transit depth uncertainty = {0:.0f} ppm per channel\n'\
                 .format( float( y['wfc3_TExoNS']['info']['Transit depth uncertainty(ppm)'] ) ) )
        f.write( 'Number of channels = {0:.0f}\n'\
                 .format( int( y['wfc3_TExoNS']['info']['Number of channels'] ) ) )
        f.write( 'Number of transits = {0:.0f}\n'\
                 .format( int( y['wfc3_TExoNS']['info']['Number of Transits'] ) ) )
        f.write( 'Start observations between orbital phases = {0}\n'\
                 .format( y['wfc3_TExoNS']['info']['Start observations between orbital phases'] ) )
        #for key in ikeys:
        #    f.write( '{0}:  {1}\n'.format( key, y['input'][key] ) )
        #f.write( '\nEXPOSURES\n{0}\n'.format( 50*'-' ) )
        #for key, value in y['timing'].items():
        #    f.write('{0}:  {1}\n'.format(key, value))
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
