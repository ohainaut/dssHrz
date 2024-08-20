#!/usr/bin/python3
# dssHrz
#
# Create a DSS finding chart, for a solar system object, 
# for one or more nights
#
# Call JPL Horizons for the ephemerides;
# Get DSS from ESO archive.
#
# replacement for dss_ephAuto and dss_ephAuto.prg
#
Oversion = "dssHrz 2024-08-19 ohainaut"

#.history
# = 2013-09-05 creation
#
#
import argparse
from astropy.time import Time
import astropy.units as u
from astroquery.jplhorizons import Horizons
from astropy.coordinates import Angle
from astropy.io import fits
from astropy.wcs import WCS

import matplotlib.pyplot as plt
import numpy as np
import sys
import os



def getDSS(RAdeg, DECdeg, FoV, survey='any', dssFile='DSS.fits'):
    '''Retreive a FITS from DSS'''

    # assemble web command line and URL
    dssGET   =  'wget -O ' + dssFile + ' '
    dssHTML  = 'http://archive.eso.org/dss/dss/image?mime-type=download-fits'
    dssHTML += f'&ra={RAdeg:.3f}&dec={DECdeg:.3f}'
    dssHTML += f'&x={FoV}&y={FoV}'

    if survey == 'any':
        surveys = ['&Sky-Survey=DSS2-red', '&Sky-Survey=DSS2-blue', '&Sky-Survey=DSS1']
    else:
        surveys = [f'&Sky-Survey={survey}']

    for dssMode in surveys:
        # call the archive
        print( "Calling Archive: ", dssHTML+dssMode )

        os.system(f'{dssGET} "{dssHTML}{dssMode}"')

        # check if it was returned
        test = os.system('grep ERROR ' + dssFile)
        if test > 0 :
            break   # found a DSS!
        else:
            raise RuntimeError( 'ERROR- no DSS found')

    return dssFile, dssMode


def makeEllipse( xy0, a, b, theta):
    '''Compute an ellipse,
    centre x0,y0
    axis a, b
    orientation theta'''
    tr = -np.radians(theta) # from the def of theta (from RA to N)

    th = np.radians(np.arange(0.,361.,15.))
    x1 = a * np.cos(th)
    y1 = b * np.sin(th)
    x =  x1 * np.cos(tr) + y1 * np.sin(tr)
    y = -x1 * np.sin(tr) + y1 * np.cos(tr)
    x = xy0[0] + x/3600.
    y = xy0[1] + y/3600.
    return x,y
#------------------------------------------------------------------------------


def plotIt(startTime, object,
           observatory='309', FoV=10., outFile="DSS", ESOprogId=""):


    Ts = Time(startTime)
    Te = Ts + 24.*u.h
    epochs = {'start': Ts.value,
            'stop' : Te.value,
            'step' : '1h' }
    # read the ephem.

    print(epochs, observatory, object)
    ephall = Horizons( id=object, location=observatory, epochs=epochs ).ephemerides()
    print('Ephemerides in')
    if len(ephall) == 0:
        raise RuntimeError('Ephemerides is empty for {object} for {startTime} - object not observable')
    else:
        print(f'Ephem: {len(ephall)} lines')

    # canonical name of the object:
    Oobject = ephall['targetname'][0]

    # filter ephemerides
    ephVis = ephall[  ephall['EL'] >= 27. ]
    if len(ephVis)== 0:
        raise RuntimeError('Ephemerides is empty for Elev>27- object not observable')
    else:
        print(f'{len(ephVis)} observable points')

    # coordinates of the centre:
    raCentre  = np.average(ephVis['RA']) # I know... dangerous
    decCentre = np.average(ephVis['DEC'])
    print('Centre',raCentre, decCentre)

    wra = Angle(raCentre/15., 'hourangle').hms
    wde = Angle(decCentre, 'deg').signed_dms
    if wde[0] >0:
        Ocenter = f'{int(wra[0])}:{int(wra[1])}:{wra[2]:.2f}  +{int(wde[1])}:{int(wde[2])}:{wde[3]:.2f}'
    else:
        Ocenter = f'{int(wra[0])}:{int(wra[1])}:{wra[2]:.2f}  -{int(wde[1])}:{int(wde[2])}:{wde[3]:.2f}'
    print( Ocenter)

    # check field-of-view
    wFoV = np.sqrt( (ephVis['RA'][-1]-ephVis['RA'][0])**2
                   + (ephVis['DEC'][-1]-ephVis['DEC'][0])**2)*60. #arcmin

    if wFoV > FoV:
        print(f'WARNING: FoV={FoV} is too small ({wFoV}arcmin needed)')
        if wFoV > 60.:
            FoV = 60.
            print(f'WARNING: Truncating to {FoV} arcmin')
        else:
            FoV = wFoV


    # time range
    Odate = f'{ephVis[0]["datetime_str"]} .. {ephVis[-1]["datetime_str"]}'
    print( Odate)

    # centre of ephem
    icent = int(len(ephVis)/2.)


    # DSS image
    dssFile, dssMode = getDSS(raCentre, decCentre, FoV, survey='any', dssFile=outFile+'.fits')

    hdu = fits.open(dssFile)
    image = hdu[0].data
    wcs = WCS(hdu[0].header)

    # cuts for the image
    im_av = np.average(image)
    im_std = np.std(image)

    plt.figure(figsize=(8.3,11))
    ax = plt.subplot(projection=wcs)

    ax.imshow(image, cmap='gray_r', vmin=im_av-im_std, vmax=im_av+ 10*im_std)
    ax.grid(color='b', ls='solid', alpha=0.4)
    ax.set_xlabel('RA')
    ax.set_ylabel('Dec')

    ax.scatter(ephVis['RA'], ephVis['DEC'], color='red',  marker='.', s=10, transform=ax.get_transform('world'))

    def rot(rar,decr):
        # returns label rotation for a dot where the rates are rar, decr
        rot = np.degrees(np.arctan2( rar, decr ))
        rot = rot if np.cos(np.radians(rot)) > 0 else rot + 180.
        return rot

    ax.text(ephVis['RA'][0], ephVis['DEC'][0],    " "+ephVis['datetime_str'][0],
            rotation=rot(ephVis['RA_rate'][0], ephVis['DEC_rate'][0]  ),
            size=8,  horizontalalignment='left', verticalalignment='center', rotation_mode='anchor',
            color='r', transform=ax.get_transform('world'))
    ax.text(ephVis['RA'][-1], ephVis['DEC'][-1],  " "+ephVis['datetime_str'][-1][12:],
            rotation=rot(ephVis['RA_rate'][0], ephVis['DEC_rate'][0]  ),
            size=8,  horizontalalignment='left', verticalalignment='center', rotation_mode='anchor',
            color='r', transform=ax.get_transform('world'))

    # label extra days

    for ephLine in ephVis[1:-1]:
        if ephLine['datetime_str'][12:] == "00:00":
            ax.text(ephLine['RA'], ephLine['DEC'],   " "+ephLine['datetime_str'],
                rotation=rot(ephLine['RA_rate'], ephLine['DEC_rate']  ),
                size=8,  horizontalalignment='left', verticalalignment='center', rotation_mode='anchor',
                color='r', transform=ax.get_transform('world'))

    # plot uncertainty
    if ephVis["SMAA_3sigma"][icent] > 5.:
        trans = ax.get_transform('world').inverted()
        x,y =  makeEllipse(trans.transform( (150., 200)),
                       ephVis["SMAA_3sigma"][icent],
                       ephVis["SMIA_3sigma"][icent],
                       ephVis["Theta_3sigma"][icent]
                       )
        ax.fill(  x, y, color='r', alpha=0.3, transform=ax.get_transform('world'))
        [x,y] = trans.transform( (150., 150))
        ax.text( x,y, "Orbit uncertainty", color='r', transform=ax.get_transform('world'))


    mag = ""
    if "V" in ephVis.columns:
        mag +=  f'{ephVis["V"][icent]:.2f} '
    elif 'Tmag' in ephVis.columns:
        mag += f'{ephVis["Tmag"][icent]:.2f} (Tmag) '
        if "Nmag" in ephVis.columns:
            mag += f'{ephVis["Nmag"][icent]:.2f} (Nmag)'


    # labels
    ax.text(0.2, 1.25, 'Object:',
                  size=20, horizontalalignment='right', transform=ax.transAxes)
    ax.text(0.22, 1.25,  Oobject,
                  size=20, horizontalalignment='left', transform=ax.transAxes)

    ax.text(0.2, 1.20, 'Date:',
                  size=15, horizontalalignment='right', transform=ax.transAxes)
    ax.text(0.22, 1.20, Odate,
                  size=15, horizontalalignment='left', transform=ax.transAxes)

    ax.text(0.2, 1.15, 'Field Center:',
                  size=15, horizontalalignment='right', transform=ax.transAxes)
    ax.text(0.22, 1.15, Ocenter,
                  size=15, horizontalalignment='left', transform=ax.transAxes)



    ax.text(0.2, 1.08, 'Field of View:',
                  size=12, horizontalalignment='right', transform=ax.transAxes)
    ax.text(0.22, 1.08, f'{FoV:.1f}\'',
                  size=12, horizontalalignment='left', transform=ax.transAxes)

    ax.text(0.6, 1.08, 'Observatory code:',
                  size=12, horizontalalignment='right', transform=ax.transAxes)
    ax.text(0.62, 1.08, observatory,
                  size=12, horizontalalignment='left', transform=ax.transAxes)
    
    ax.text(0.2, 1.05,'Uncertainty:',
                  size=12, horizontalalignment='right', transform=ax.transAxes)
    ax.text(0.22, 1.05, f'{ephVis["SMAA_3sigma"][icent]:.1f} \"',
                  size=12, horizontalalignment='left', transform=ax.transAxes)
    ax.text(0.6, 1.05, 'Apparent velocity:',
                  size=12, horizontalalignment='right', transform=ax.transAxes)
    ax.text(0.62, 1.05, f'{ephVis["RA_rate"][icent]:.2f}, {ephVis["DEC_rate"][icent]:.2f}"/h',
                  size=12, horizontalalignment='left', transform=ax.transAxes)

    ax.text(0.2, 1.02, 'Galactic lat:',
                  size=12, horizontalalignment='right', transform=ax.transAxes)
    ax.text(0.22, 1.02, f'{ephVis["GlxLat"][icent]:.1f} deg',
                  size=12, horizontalalignment='left', transform=ax.transAxes)

    ax.text(0.6, 1.02, 'Magnitude:',
                  size=12, horizontalalignment='right', transform=ax.transAxes)
    ax.text(0.62, 1.02, mag,
                  size=12, horizontalalignment='left', transform=ax.transAxes)



    ax.text(1., -0.1,  ESOprogId,
                  size=20, horizontalalignment='right', transform=ax.transAxes)

    ax.text(0.0, -0.1,  Oversion,
                  size=8, horizontalalignment='left', color='grey', transform=ax.transAxes)

    ax.text(0.5, -0.1,  f'Survey: {dssMode}',
                  size=8, horizontalalignment='center', transform=ax.transAxes)



    plt.savefig(outFile+'.pdf')
    print(f'>> Finder in {outFile}.pdf')


#############################################################################################

if __name__ == "__main__":
    #--- command line arguments

    parser = argparse.ArgumentParser(description='''Generate a DSS finding chart for
                                     a solar system object for 1 night''')
    parser.add_argument('-f','--outFile', default='DSS',
                            help='Root of the output files.')
    parser.add_argument('-o','--object',
                            help='''Designation of the object;
                            must be resolved by Horizon;
                            in case of doubt use the Unique JPL ID
                            (in the ephem header, 1st line Rec #:''')
    parser.add_argument('-n','--night',
                            help='''Night to be considered: YYYY-MM-DD[Thh:mm];
                            the ephem will go from t0 to t0+24h''')
    parser.add_argument('-e','--endNight',
                            help='''If provided, will iterate from night to endNight''')
    parser.add_argument('-t','--tel', default=309,
                        help="IAU Observatory code or name, must be resolved by JPL")
    parser.add_argument('--FoV', default=10,
                        help="Field of view, arcmin; will be extended if needed")
    parser.add_argument('-p','--ESOprogId', default=" ",
                        help="ESO program ID")

    parser.add_argument('-s','--survey', default="any",
                        help="survey, from DSS2-red DSS2-blue DSS1 or any")

    myargs = parser.parse_args()

    Ts = Time(myargs.night)

    if myargs.endNight:
        Te = Time(myargs.endNight) + 1.*u.day
    else:
        Te = Ts + 1.*u.day

    FoV = float(myargs.FoV)

    for jd in np.arange(Ts.jd, Te.jd):
        night = Time(jd, format='jd').iso[:10]
        outFile = f'{myargs.outFile}_{night}'
        plotIt(Time(jd, format='jd').iso, myargs.object, observatory=myargs.tel, FoV=FoV, outFile=outFile, ESOprogId=myargs.ESOprogId  )
