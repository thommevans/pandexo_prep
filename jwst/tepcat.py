from __future__ import print_function
import pdb, sys, os
import numpy as np
try:
    # For Python 3.0 and later:
    from urllib.request import urlopen
    from urllib.error import URLError
except ImportError:
    # For earlier Python:
    from urllib2 import urlopen, URLError

HPLANCK_SI = 6.62607e-34 # planck's constant in J*s
C_SI = 2.99792e8 # speed of light in vacuum in m s^-1
KB_SI = 1.3806488e-23 # boltzmann constant in J K^-1
G_SI = 6.67428e-11 # gravitational constant in m^3 kg^-1 s^-2
DAY_SI = 24*60*60 # day in seconds
RSUN_SI = 6.9551e8 # solar radius in m
MSUN_SI = 1.99e30 # solar mass in kg
RJUP_SI = 7.1492e7 # jupiter radius in m
MJUP_SI = 1.89852e27 # jupiter mass in kg
AU_SI = 1.49598e11 # au to metres conversion factor
RGAS_SI = 8.314 # gas constant in J mol^-1 K^-1
MUJUP_SI = 2.22e-3 # jupiter atmosphere mean molecular weight in kg mole^-1

def load( download_latest=True ):
    """
    Routine called by the run_jwst.py script to load the TEPCat catalogues into
    a format that can be used by the routines in the jwstsim.py module.
    """

    if internet_on() and download_latest==True:
        print("TEPCat connection successful: Downloading latest tables")
        download_tables()
    elif not internet_on():
        print("No internet connection to TEPCat: Using saved tables")
    else:
        print("Download not requested: Using saved tables")


    # Read contents of the first tepcat file into numpy arrays:
    ifile1 = open( 'tepcat1.txt', 'r' )
    header1 = ifile1.readline() # skip first header line
    names1 = []
    vmags1 = []
    kmags1 = []
    tdurs1 = []
    tdepths1 = []
    periods1 = []
    for line in ifile1:
        z = line.split()
        names1 += [ str( z[0] ) ]
        vmags1 += [ float( z[8] ) ]
        kmags1 += [ float( z[9] ) ]
        tdurs1 += [ float( z[10] ) ]
        tdepths1 += [ float( z[11] ) ]
        periods1 += [ float( z[14] ) ]
        print( names1[-1], periods1[-1] )
    ifile1.close()
    names1 = np.array( names1 )
    vmags1 = np.array( vmags1 )
    kmags1 = np.array( kmags1 )
    tdurs1 = np.array( tdurs1 )
    tdepths1 = np.array( tdepths1 )
    periods1 = np.array( periods1 )
    # Read contents of the second tepcat file into numpy arrays:
    ifile2 = open( 'tepcat2.txt', 'r' )
    ifile2.seek( 0 )
    header2 = ifile2.readline() # skip first header line
    names2 = []
    tstar = []
    metalstar = []
    mstar = []
    rstar = []
    loggstar = []
    a = []
    mplanet = []
    rplanet = []
    littleg = []    
    rhoplanet = []
    tplanet_tepcat = []
    for line in ifile2:
        z = line.split()
        mstar_i = float( z[7] )
        rstar_i = float( z[10] )
        loggstar_i = float( z[11] )
        mplanet_i = float( z[26] )
        rplanet_i = float( z[29] )
        littleg_i = float( z[32] )
        rhoplanet_i = float( z[35] )
        if ( mstar_i>0 )*( rstar_i>0 )*( mplanet_i>0 )*( rplanet_i>0 ):
            names2 += [ str( z[0] ) ]
            tstar += [ float( z[1] ) ]
            metalstar += [ float( z[2] ) ]
            mstar += [ mstar_i ]
            rstar += [ rstar_i ]
            loggstar += [ loggstar_i ]
            a += [ float( z[23] ) ]
            mplanet += [ mplanet_i ]
            rplanet += [ rplanet_i ]
            littleg += [ littleg_i ]
            rhoplanet += [ rhoplanet_i ]
            tplanet_tepcat += [ float( z[38] ) ]
    ifile2.close()
    names2 = np.array( names2 )
    tstar = np.array( tstar )
    metalstar = np.array( metalstar )
    mstar = np.array( mstar )
    rstar = np.array( rstar )
    loggstar = np.array( loggstar )
    a = np.array( a )
    mplanet = np.array( mplanet )
    rplanet = np.array( rplanet )
    littleg = np.array( littleg )
    rhoplanet = np.array( rhoplanet )
    tplanet_tepcat = np.array( tplanet_tepcat )
    # Merge the two catalogues:
    print( 'Merging both catalogues:' )
    n = len( names2 )
    m = len( names1 )
    names = names2
    vmags = np.zeros( n )
    kmags = np.zeros( n )
    tdurs = np.zeros( n )
    tdepths = np.zeros( n )
    periods = np.zeros( n )
    for i in range( n ):
        for j in range( m ):
            if names[i]==names1[j]:
                vmags[i] = vmags1[j]
                kmags[i] = kmags1[j]
                tdurs[i] = tdurs1[j]
                tdepths[i] = tdepths1[j]
                periods[i] = periods1[j]
                print( '... matched {0} --> {1}'.format( names[i], names1[j] ) )
                break
        if j==m-1:
            print( 'Could not match {0}!'.format( names[i] ) )
            pdb.set_trace()

    littleg = G_SI*mplanet*MJUP_SI/( ( rplanet*RJUP_SI )**2. )
    volplanet = ( (4*np.pi/3.)*( ( rplanet*RJUP_SI )**3. ) )
    rhoplanet = mplanet/volplanet
    term1 = ( periods*DAY_SI/2/np.pi )**2.
    term2 = G_SI*( mstar*MSUN_SI+mplanet*MJUP_SI )
    a = ( ( term1*term2 )**(1./3.) )/AU_SI
    aRs = ( a*AU_SI )/( rstar*RSUN_SI )
    tplanet = calc_teq( tstar, aRs, Ab=0, fprime=0.25 )
    RpRs = ( rplanet*RJUP_SI )/( rstar*RSUN_SI )

    ##############################################
    # Calculate and normalise the thermal signal assuming
    # planet radiates as blackbody at equilibrium temperature:
    nearir_um = 2.2
    ecdepth_nearir = calc_ecdepth( tplanet, tstar, RpRs, nearir_um )
    ix = ( names=='WASP-121' )
    tref = tstar[ix]
    kref = kmags[ix]
    fratios_nearir = calc_fratio( nearir_um, tstar, kmags, tref, kref )
    sn_em = ecdepth_nearir*np.sqrt( fratios_nearir )
    sn_em /= sn_em[ix]

    # Calculate and normalise the transmission signal
    # assuming a hydrogen-dominated atmosphere:
    hatm = RGAS_SI*tplanet/MUJUP_SI/littleg
    hdepth = 2*hatm*(rplanet*RJUP_SI)/( ( rstar*RSUN_SI )**2. )
    ix = ( names=='WASP-121' )
    dkmag = kmags-kmags[ix]
    fratio = 10**( -dkmag/2.5 )
    sn_tr = hdepth*np.sqrt( fratio )
    ##############################################

    ixs = ( vmags>0 )*( kmags>0 ) # restrict to those with reliable brightnesses
    tepcat = {}
    tepcat['names'] = names[ixs]
    tepcat['a'] = a[ixs]
    tepcat['periods'] = periods[ixs]
    tepcat['mplanet'] = mplanet[ixs]
    tepcat['RpRs'] = RpRs[ixs]
    tepcat['mstar'] = mstar[ixs]
    tepcat['rstar'] = rstar[ixs]
    tepcat['loggstar'] = loggstar[ixs]
    tepcat['rplanet'] = rplanet[ixs]
    tepcat['tstar'] = tstar[ixs]
    tepcat['metalstar'] = metalstar[ixs]
    tepcat['tplanet'] = tplanet[ixs]
    tepcat['rplanet'] = rplanet[ixs]
    tepcat['kmags'] = kmags[ixs]
    tepcat['vmags'] = vmags[ixs]
    tepcat['littleg'] = littleg[ixs]
    tepcat['rhoplanet'] = rhoplanet[ixs]
    tepcat['tdurs'] = tdurs[ixs]
    tepcat['tdepths'] = tdepths[ixs]
    tepcat['sn_tr'] = sn_tr[ixs]
    tepcat['sn_em'] = sn_em[ixs]

    # Sort in order of decreasing transmission signal:
    ixs = np.argsort( tepcat['sn_tr'] )
    keys = list( tepcat.keys() )
    for k in keys:
        tepcat[k] = tepcat[k][ixs][::-1]

    print( '\nFinished reading TEPCat.\n' )
    return tepcat


def calc_teq( tstar, aRs, Ab=0, fprime=0.25 ):
    redist = fprime*( 1-Ab )
    return tstar*( np.sqrt( 1./aRs ) )*( redist**0.25 )

def calc_ecdepth( tplanet, tstar, RpRs, wav_um ):
    wav_m = wav_um*(1e-6)
    bratio = planck( wav_m, tplanet )/planck( wav_m, tstar )
    ecdepth = bratio*( RpRs**2. )    
    return ecdepth

def planck( wav_m, temp ):
    """
    Evaluates the Planck function for given values of wavelength
    and temperature. Wavelength should be provided in metres and
    temperature should be provided in Kelvins.
    """
    term1 = 2 * HPLANCK_SI * ( C_SI**2. ) / ( wav_m**5. )
    term2 = np.exp( HPLANCK_SI * C_SI / KB_SI / wav_m / temp ) - 1
    bbflux = term1 / term2
    return bbflux

def calc_fratio( wav_um, t, kmag, tref, kmagref ):
    wav_m = wav_um*(1e-6)
    k_m = 2.2e-6
    bratio = planck( wav_m, t )/planck( k_m, t )
    bratioref = planck( wav_m, tref )/planck( k_m, tref )
    termA = bratio/bratioref
    delk = kmag-kmagref
    termB = 10**( -delk/2.5 )
    return termA*termB

def internet_on():
    try:
        urlopen('http://www.astro.keele.ac.uk/jkt/tepcat/allplanets-ascii.txt', timeout=1)
        return True
    except URLError as err:
        return False

def download_tables():
    table = urlopen('http://www.astro.keele.ac.uk/jkt/tepcat/observables.txt').read()
    with open('tepcat1.txt', 'wb') as f:
        f.write(table)

    table2 = urlopen('http://www.astro.keele.ac.uk/jkt/tepcat/allplanets-ascii.txt').read()
    with open('tepcat2.txt', 'wb') as g:
        g.write(table2)
    return 
