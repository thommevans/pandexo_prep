from __future__ import print_function
import os, pdb, sys, time
import numpy as np
import pandexo.engine.justdoit as jdi
import pandexo.engine.justplotit as jpi
import matplotlib.pyplot as plt



def main( planet_label, tepcat, sat_level=80, sat_unit='%', noise_floor_ppm=20, inst_modes='all' ):

    ix = ( tepcat['names']==planet_label )
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

    # Load a null spectrum for the planet:
    z['planet']['exopath'] = 'zeros.txt'

    # Run PandExo over requested instrument modes:
    if inst_modes=='all':
        inst_modes = jdi.print_instruments()
    y = jdi.run_pandexo( z, inst_modes )

    # Extract output (TODO properly)
    x1,y1,e1 = jpi.jwst_1d_spec(y, R=100)

    data1 = np.transpose([x1[0], y1[0], e1[0]])


    plt.figure()
    ax = plt.axes()

    data1 = data1[data1[:,2] < 0.01] #Masks values with very high errors
    ax.errorbar(data1[:,0], data1[:,1], yerr=data1[:,2], ms=3, linestyle='None', marker='s', capsize=0, markeredgewidth=0)

    plt.show()
    pdb.set_trace()
    return None


