from   astropy import time
from   astropy.coordinates import AltAz, EarthLocation, SkyCoord
import astropy.units as u
from   astropy.coordinates import get_sun, get_moon
from   astropy import coordinates as coord
from   astropy import time
from   astroplan import moon
from   astroquery.vizier import Vizier
from   astroquery.utils.tap.core import TapPlus
from   misc import bcolors
import pandas as pd
from   plotsettings import *
from   standard_libraries import *


# Optical
obs_GemminiSouth  = EarthLocation.of_site('gemini_south')
obs_LaSilla       = EarthLocation.of_site('lasilla')
obs_LBT           = EarthLocation.of_site('lbt')
obs_NOT           = EarthLocation.of_site('lapalma')
obs_Palomar       = EarthLocation.of_site('Palomar')
obs_VLT           = EarthLocation.of_site('paranal')

# Radio
obs_ALMA          = EarthLocation.of_site('ALMA')
obs_ATCA          = EarthLocation.from_geodetic(lat=-30.3128 * u.deg, lon=149.5502 * u.deg, height=0 * u.m)
obs_GMRT          = EarthLocation.from_geodetic(lat=-19.06   * u.deg, lon=74.03 *    u.deg, height=0 * u.m)
obs_VLA           = EarthLocation.of_site('vla')

def get_atca_timeslots(DATE, URL):
       
    # process obs calendar
    tables = pd.read_csv(URL, sep=',', header=None)
    tables = tables.rename(columns={0: 'date_start', 1: 'UT_start', 2: 'date_stop', 3: 'UT_stop', 4: 'program', 5: 'configuration'})
    
    # Filter the free slots
    tables = tables[(tables.program == ' Directors time')]
    #tables = tables[('Reconfigure' not in tables.program)]
    
    # Get duration
    
    tstart = [tables.date_start.values.tolist()[ii] + 'T' + tables.UT_start.values.tolist()[ii] for ii in range(len(tables))]
    tstop  = [tables.date_stop.values.tolist()[ii]  + 'T' + tables.UT_stop.values.tolist()[ii]  for ii in range(len(tables))]

    tstart = [x.replace(' ', '') for x in tstart]
    tstop  = [x.replace(' ', '') for x in tstop]
    
    duration = [(time.Time(tstop[ii], format='isot', scale='utc') - time.Time(tstart[ii], format='isot', scale='utc')).value*24 for ii in range(len(tstart))]
    
    
    tables['duration'] = duration
    #tables['duration'].format = '.1f'
    
    return tables.query("@DATE <= date_start")

def check_moon(coords, obstime, observatory, avoid=30.*u.degree):

    if not isinstance(avoid, u.Quantity):
        avoid = float(avoid)*u.degree
    else:
        avoid = avoid.to(u.degree)

    _moon = get_moon(obstime, obs_NOT)

    target = coords
    
    moon_alt = _moon.transform_to(AltAz(obstime=obstime, location=observatory)).alt.to(u.deg)
    
    #if moon_alt < 0*u.degree:
        #print('Moon is down')
    #    return {'sep': 0, 'fli':np.nan}
    
    """
    else:
        sep = _moon.separation(target)
        fli = moon.moon_illumination(obstime)
        
        msg_sep = 'Moon separation = {sep:.0f} deg, fli = {fli:.0f}%'.format(sep=sep.to(u.degree).value, fli=fli*100)
        #msg_fli = 'Full moon illumination {:.0f}%'.format(fli*100)
        msg_fail = 'Object too close to the moon!'
       
        if sep > avoid:
            print(bcolors.OKGREEN + msg_sep + bcolors.ENDC)
            #print(bcolors.OKGREEN + msg_fli + bcolors.ENDC)
            #print(bcolors.OKGREEN + msg_fail + bcolors.ENDC)
        else:
            print(bcolors.FAIL + bcolors.BOLD + msg_sep + bcolors.ENDC)
            #print(bcolors.FAIL + bcolors.BOLD + msg_fli + bcolors.ENDC)

        return {'sep': sep.to(u.degree).value, 'fli':fli*100}
    """
    
def check_sun(coords, obstime, observatory, avoid=30.*u.degree):

    if not isinstance(avoid, u.Quantity):
        avoid = float(avoid)*u.degree
    else:
        avoid = avoid.to(u.degree)

    _sun = get_sun(obstime)

    target = coords
    
    sun_alt = _sun.transform_to(AltAz(obstime=obstime, location=observatory)).alt.to(u.deg)
    
    if sun_alt < 0*u.degree:
        print('Sun is down')
        return {'sep': 0}
    
    else:
        sep = _sun.separation(target)
        
        msg_sep = 'Sun separation = {sep:.0f} deg'.format(sep=sep.to(u.degree).value)
        msg_fail = 'Object too close to the moon!'
        
        """
        if sep > avoid:
            print(bcolors.OKGREEN + msg_sep + bcolors.ENDC)
            #print(bcolors.OKGREEN + msg_fli + bcolors.ENDC)
            #print(bcolors.OKGREEN + msg_fail + bcolors.ENDC)
        else:
            print(bcolors.FAIL + bcolors.BOLD + msg_sep + bcolors.ENDC)
            #print(bcolors.FAIL + bcolors.BOLD + msg_fli + bcolors.ENDC)
        """
        return {'sep': sep.to(u.degree).value}

def plot_objects(DATE, CAT, OBS=obs_NOT, SAVE=False):

    # Time array

    midnight_utc = time.Time(DATE, format='isot', scale='utc')#+1
    delta_midnight = np.linspace(0, 24, 1000)*u.hour

    # Plot

    plt.figure(figsize=(9*np.sqrt(2), 9))
    ax = plt.subplot(111)

    frame_time  = AltAz(obstime=midnight_utc + delta_midnight, location=OBS)

    for ii in range(len(CAT)):

        print(bcolors.BOLD + bcolors.HEADER + CAT['NAME'][ii] + bcolors.ENDC)
        
        coords      = SkyCoord(CAT['RA'][ii], CAT['DEC'][ii], unit=(u.hour, u.deg))
        obj_altazs  = coords.transform_to(frame_time)
        obj_airmass = obj_altazs.secz

        moon_distance = check_moon(coords, midnight_utc, OBS)        
        
        mask = (obj_airmass < 5) & (obj_airmass > 0)
        #label='{obj} (moon distance: {sep:.0f}°)'.format(obj=CAT['NAME'][ii], sep=moon_distance['sep'])
        label='{obj}'.format(obj=CAT['NAME'][ii])
        
        ax.plot(delta_midnight, obj_altazs.alt,
                label=label, lw=4, color='crimson', path_effects=[PathEffects.withStroke(linewidth=8, foreground="white")])

        sun_distance = check_sun(coords, midnight_utc, OBS)        


    # Night

    sunaltazs = get_sun(midnight_utc + delta_midnight).transform_to(frame_time)

    ax.fill_between(delta_midnight, 0*u.deg, 90*u.deg,
                    sunaltazs.alt < -0*u.deg, color='navy', zorder=0, alpha=0.2)

    ax.fill_between(delta_midnight, 0*u.deg, 90*u.deg,
                    sunaltazs.alt < -18*u.deg, color='navy', zorder=0)

    # Moon

    moonaltazs = get_moon(midnight_utc + delta_midnight).transform_to(frame_time)
    ax.plot(delta_midnight, moonaltazs.alt, lw=2, color='black', path_effects=[PathEffects.withStroke(linewidth=8, foreground="white")], label='Moon')

    # Sun

    sunaltazs = get_sun(midnight_utc + delta_midnight).transform_to(frame_time)
    ax.plot(delta_midnight, sunaltazs.alt, lw=2, color='yellow', path_effects=[PathEffects.withStroke(linewidth=8, foreground="white")], label='Sun')


    # Airmass limits

    ax.axhline(90 - np.arccos(1/2.)*180/np.pi, color=vigit_color_12)
    #ax.axhline(90 - np.arccos(1/3.)*180/np.pi, color=vigit_color_12, ls='--')

    # Midnight

    #ax.axvline(0., color='white', lw=2)

    # Legend

    ax.legend(loc='upper left', fontsize=legend_size-4, bbox_to_anchor=(1, 0.25, 0.5, 0.5))

    # Prettify

    xmin, xmax = 0, 24

    ax.set_xlim(xmin*u.hour, xmax*u.hour)
    labels = ['{:.0f} h'.format(x)for x in (np.arange(xmax-xmin)+xmin) if x % 2 == 0 ]
    ax.set_xticks((2*np.arange(12))*u.hour, labels=labels, rotation=45)

    ax.set_ylim(10*u.deg, 90*u.deg)

    ax.set_xlabel('Universal time', fontsize=label_size)
    ax.set_ylabel('Altitude', fontsize=label_size)

    try:
        ax.set_title(DATE + ' (full moon illumination = ' + str(int(moon_distance['fli'])) + '%)', fontsize=label_size)
    except:
        ax.set_title(DATE, fontsize=label_size)

    if SAVE:
        plt.savefig('visibility_targets_{}.pdf'.format(DATE))

def query_askap(COORD, RADIUS):
    
    # Reference: Hale et al. (2021), https://arxiv.org/abs/2109.00956
    # Frequency: 887.5 MHz
    
    casdatap = TapPlus(url="https://casda.csiro.au/casda_vo_tools/tap")

    query = "SELECT * FROM  AS110.racs_dr1_gaussians_galacticcut_v2021_08_v01 \
            where 1=CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', {ra:.4f}, {dec:.4f}, {radius:.4f}))".format(ra=COORD.ra.deg, dec=COORD.dec.deg, radius=RADIUS.to(u.deg).value)
    
    job = casdatap.launch_job_async(query)

    result = job.get_results()
    
    if len(result) > 0:
        
        _coords_askap = coord.SkyCoord(result['ra'], result['dec'])
        
        keys = result.keys()
        
        result['_r'] = _coords_askap.separation(COORD).to(u.arcmin).value
        result['_r'].format = '.4f'
        result.sort('_r')
        
        return result[['_r'] + keys]
    
    else:
        return None

def query_vlass(COORD, RADIUS):
    
    # Reference: Gordon et al. (2021), https://arxiv.org/abs/2102.11753
    # Frequency: 3 GHz
    
    v = Vizier(columns = ['all'], catalog = 'J/ApJS/255/30', row_limit = -1)
    result = v.query_region(COORD, radius = RADIUS, catalog='J/ApJS/255/30')

    if len(result) > 0:
        result = result['J/ApJS/255/30/comp']
        result.sort('_r')
        return result
    
    else:
        return None

def query_first(COORD, RADIUS):
    
    # Reference: White et al. (1997), https://ui.adsabs.harvard.edu/abs/1997ApJ...475..479W/abstract
    # Frequency: 1.4 GHz
    
    v = Vizier(columns = ['all'], catalog = 'VIII/90', row_limit = -1)
    result = v.query_region(COORD, radius = RADIUS, catalog='VIII/90')

    if len(result) > 0:
        result = result['VIII/90/first12']
        result.sort('_r')
        return result
    
    else:
        return None

def query_nvss(COORD, RADIUS):

    # Reference: White et al. (1997), https://ui.adsabs.harvard.edu/abs/1997ApJ...475..479W/abstract
    # Frequency: 1.4 GHz
    
    v = Vizier(columns = ['all'], catalog = 'VIII/65', row_limit = -1)
    result = v.query_region(COORD, radius = RADIUS, catalog='VIII/65')

    if len(result) > 0:
        result = result['VIII/65/nvss']
        result.sort('_r')
        return result
    
    else:
        return None
