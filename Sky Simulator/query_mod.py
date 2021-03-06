from astroquery.simbad import Simbad
from astroquery.gaia import Gaia
import astropy.coordinates as coord
from astropy import wcs
from astropy import units as u
import numpy as np
import json

from Logger import hist
from astropy.logger import log

#filename = "Antola_data.json"
filename = "San_Pedro_data.json"
with open(filename) as f:
    DATA = json.load(f)

    
def Coordinator(coord, center, CCD_structure, catalog):
    '''
    The function takes the coordinates of a star and uses the WCS keywords
    to give back the position on the CCD in pixel

    Parameters
    ----------

    coord : SkyCoord
        It is the coordinates of the star already elaborated by astropy

    center : array
        It must contain the position of the center of the CCD and it is used to obtain
        the distance of the star from it

    CCD_res : float
        It is the resolution of the CCD in arcsec/pixel

    Outputs
    -------

    x : float
        It is the position of the star on the grid of the CCD along the axis x

    y : float
        It is the position of the star on the grid of the CCD along the axis y
    '''

    log.info(hist())
    if catalog == 'Simbad':
        offset_x = CCD_structure[0]/CCD_structure[2]
        offset_y = CCD_structure[1]/CCD_structure[2]
        w = wcs.WCS(naxis=2)
        w.wcs.crpix = [1, 1]
        w.wcs.crval = [center[0], center[1]]
        w.wcs.ctype = ["RA", "DEC"]    
        x,y = wcs.utils.skycoord_to_pixel(coord, w)
        x = ((x*3600)/CCD_structure[3]) + offset_x/2 #(deg*arc/sec*deg)*arc/sec*pixel+offset
        y = ((y*3600)/CCD_structure[3]) + offset_y/2
    elif catalog == 'Gaia':
        offset_x = CCD_structure[0]/CCD_structure[2]
        offset_y = CCD_structure[1]/CCD_structure[2]
        w = wcs.WCS(naxis=2)
        w.wcs.crpix = [1, 1]
        w.wcs.crval = [center[0], center[1]]
        w.wcs.ctype = ["RA", "DEC"]    
        x,y = wcs.utils.skycoord_to_pixel(coord, w)
        x = ((x*3600)/CCD_structure[3]) + offset_x/2 #(deg*arc/sec*deg)*arc/sec*pixel+offset
        y = ((y*3600)/CCD_structure[3]) + offset_y/2

    return y, x


def VEGA_to_AB(magnitudo, photo_filter):
    '''
    This function takes the magnitudo of a star in the VEGA system and converts it in the AB system.
    It download the conversion table calling data_loader

    Parameters
    ----------

    magnitudo : array (float elements)
        An element of magnitudo is the magnitudo of a star in a specific range (photo_filter)

    photo_filter : list (string elements)
        It is a list where a element is a string that identifies a filter, e.g. "U","R","I"

    Output
    ------

    magnitudo : array
        The array contains the elements converted        
    '''

    log.debug(hist(photo_filter))

    sub = DATA["Convertion_table"]
    convertion_table = [sub[k] for k in photo_filter if k in sub]
    for i in range(len(magnitudo)-1):
        if magnitudo[i]== 0:
            magnitudo[i]= 40
        else:
            magnitudo[i] = magnitudo[i] - convertion_table[i]
        magnitudo[i] = magnitudo[i] - convertion_table[i]
    return magnitudo


def atmospheric_attenuation(magnitudo, photo_filter, AirMass):
    '''
    This function takes the magnitudo of a star and gives back the same magnitudo attenuated by the atmosphere.
    It download the extinction coefficients calling data_loader

    Parameters
    ----------

    magnitudo : array (float elements)
        An element of magnitudo is the magnitudo of a star in a specific range (photo_filter)

    photo_filter : list (string elements)
        It is a list where a element is a string that identifies a filter, e.g. "U","R","I"

    AirMass : float
        The AirMass provide the thickness of the atmosphere crossed by the light 

    Output
    ------

    magnitudo : array
        The array contains the elements attenuated        
    '''
    log.debug(hist())

    sub = DATA["Extinction_coefficient"]
    extinction_coefficient = [sub[k] for k in photo_filter if k in sub]
    
    for i in range(len(magnitudo)-1):
        #if magnitudo[i] == 0:
         #   magnitudo[i] = 0
        #else:
         #   magnitudo[i] = magnitudo[i]-AirMass*extinction_coefficient[i]
        magnitudo[i] = magnitudo[i]+AirMass*extinction_coefficient[i]
    return magnitudo


def magnitudo_to_photons(magnitudo, photo_filter):
    '''
    This function convertes the magnitudo of a star in number of photons/second.
    It download informations calling data_loader

    Parameters
    ----------

    magnitudo : array (float elements)
        An element of magnitudo is the magnitudo of a star in a specific range (photo_filter)

    photo_filter : list (string elements)
        It is a list where a element is a string that identifies a filter, e.g. "U","R","I"

    Output
    ------

    photons : array
        The elements of the array are the number of photons for a specific range        
    '''
    log.debug(hist())
    
    sub = DATA["Central_wavelenght"]
    central_wavelenght = [sub[k] for k in photo_filter if k in sub]
    
    number = magnitudo   
    for i in range(len(magnitudo)-1):
        if magnitudo[i] == 0:
            number[i] = 0
        else:
            exp = 6.74-0.4*magnitudo[i]
            number[i] = (10**exp)#/(central_wavelenght[i])
    return number


def photons_to_electrons(photons, photo_filter):
    '''
    This function takes the number of photons/sec and converts them in number of electron/sec on the CCD.
    It calls data_loader to download the quantum efficiency of the CCD

    Parameters
    ----------

    photons : array
        The elements of the array are the number of photons for a specific range

    photo_filter : list (string elements)
        It is a list where a element is a string that identifies a filter, e.g. "U","R","I"

    Output

    electron : array
        The elements of the array are the number of electrons for a specific range 
    '''
    log.debug(hist())

    sub = DATA["Quantum_efficiency"]
    quantum_efficiency = [sub[k] for k in photo_filter if k in sub]
    
    electrons = photons
    for i in range(len(photons)-1):
        electrons[i] = photons[i]*quantum_efficiency[i]

    return electrons


def Data_structure():
    '''
    This function simply creates the data structure used in query
    '''
    log.debug(hist())
    
    Coord_x = []
    Coord_y = []
    Flux_tot = []
    data = [Coord_x, Coord_y, Flux_tot]
    return data


def radius_sky_portion(CCD_structure):
    '''
    This function uses the information on the CCD to estimate the radius of the portion of the sky visualized

    Parameters
    ----------

    CCD_structure : array (float elements)
        This object contains the information about the measure of the CCD
        [0] : the length of the CCD on the x axis in mm
        [1] : the length of the CCD on the y axis in mm
        [2] : the length of a single pixel in mm

    Output
    ------
     radius : astropy.units.quantity.Quantity
         it is the radius in an astropy undestandable format
    '''
    log.debug(hist())
    
    pixel_on_x = CCD_structure[0]/CCD_structure[2]
    pixel_on_y = CCD_structure[1]/CCD_structure[2]
    diagonal_diameter = np.sqrt((pixel_on_x)**2+(pixel_on_y)**2)
    radius_on_pixel = diagonal_diameter/2
    radius = radius_on_pixel*CCD_structure[3]/60
    radius = (radius/60)*u.deg
    return radius


def magnitudo_to_electrons(magnitudo, photo_filters, AirMass, exposure_time, Controll):
    '''
    This functions calls VEGA_to_AB, atmospheric_attenuation, magnitudo_to_photons, photons_to_electrons
    to obtain the total electrons generated in a CCD by the light of a star
    
    Parameters
    ----------

    magnitudo : array (float elements)
        An element of magnitudo is the magnitudo of a star in a specific range (photo_filter)

    photo_filter : list (string elements)
        It is a list where a element is a string that identifies a filter, e.g. "U","R","I"

    AirMass : float
        The AirMass provide the thickness of the atmosphere crossed by the light

    exposure_time : int
        It is the exposure time used to obtain the image, it is in seconds

    Output
    ------

    tot : float
        It is the total number of electrons generated by the star on the CCD
    '''
    log.debug(hist())
    magnitudo = VEGA_to_AB(magnitudo, photo_filters)
    magnitudo = atmospheric_attenuation(magnitudo, photo_filters, 1)
    photons = magnitudo_to_photons(magnitudo, photo_filters)
    electrons = photons_to_electrons(photons, photo_filters)
    
    if Controll:
        for n, i in enumerate(photo_filters):
            if i == "V":
                if electrons[n] <= 3e-4:
                    electrons[n] = electrons[-1]
        tot = (sum(electrons)-electrons[-1])*exposure_time
    else:
        tot = sum(electrons)*exposure_time
    return tot


def query(coordi, photo_filters, CCD_structure, exposure_time, catalog):
    '''
    This funtion uses all the functions above to query Simbad and extract the information about the stars in a
    specific area of the sky

    Parameters
    ----------

    coordi : string
        It contains the center of the area of interest

    photo_filter : list (string elements)
        It is a list where a element is a string that identifies a filter, e.g. "U","R","I"

    CCD_structure : array (float elements)
        This object contains the information about the measure of the CCD
        [0] : the length of the CCD on the x axis in mm
        [1] : the length of the CCD on the y axis in mm
        [2] : the length of a single pixel in mm

    exposure_time : int
        It is the exposure time used to obtain the image, it is in seconds

    Outputs
    -------

    data : list
        It is a list of array that contains
         [0] : array, the position of the stars along the x axis
         [1] : array, the position of the stars along the y axis
         [2] : array, the number of electrons generates by the stars on the CCD

    center : list
        It contains the center of the sky portion in "RA" and "DEC"
    '''
    log.info(hist())
    
    radius = radius_sky_portion(CCD_structure)
    coordi = coord.SkyCoord(coordi, frame = 'icrs', unit=(u.hourangle, u.deg))

    center = [0,0]
    center[0] = coordi.ra.degree
    center[1] = coordi.dec.degree
    data = Data_structure()
    AirMass = 1 #simple Antola

    if catalog == 'Simbad':
        Controll_V = False

        for n, i in enumerate(photo_filters):
            if i == "g":
                photo_filters[n] = "V"

        photo_filters_temp = photo_filters
    
        if "V" in photo_filters_temp:
            photo_filters_temp.append("g")
            Controll_V = True

        f =['flux({0})'.format(x) for x in photo_filters_temp]
        [Simbad.add_votable_fields(g) for g in f]
        result_table = Simbad.query_region(coordi, radius)

        Stars_number = len(result_table)
        len_filter_list = len(photo_filters_temp)
        
        if Controll_V:
            log.warning('If filter V is not available, g filter will be used instead')
            log.info(hist(f"Creating {Stars_number} stars:"))
        for i in range(Stars_number):

            log.info(f"Star n° {i+1}/{Stars_number}")
            row = result_table[i]
            tot = 0
            magnitudo = []
            for j in range(11, 11+len_filter_list):  #controll
                if row[j] is np.ma.masked:
                    row[j] = 0
                magnitudo.append(row[j])
        
            tot = magnitudo_to_electrons(magnitudo, photo_filters_temp, AirMass, exposure_time, Controll_V)
    
            if tot != 0: #if total flux is 0 the object is not passed
                c = coord.SkyCoord(row[1], row[2], frame = 'icrs', unit=(u.hourangle, u.deg))##change here
                x, y = Coordinator(c, center, CCD_structure, catalog)
                data[0].append(x)
                data[1].append(y)
                data[2].append(tot)
    elif catalog == 'Gaia':
        Gaia.ROW_LIMIT = -1
        tables = Gaia.cone_search_async(coordi, radius)
        tables = tables.get_results()
        photonis = tables['phot_g_mean_flux'][:]
        mag = tables['phot_g_mean_mag'][:]
        pos_ra = tables['ra'][:]
        pos_dec = tables['dec'][:]
        number_of_stars = len(photonis)
        for i in range(number_of_stars):
            c = coord.SkyCoord(pos_ra[i], pos_dec[i], frame = 'icrs', unit=(u.degree, u.degree))
            x, y = Coordinator(c, center, CCD_structure, catalog)
            mag_atm_att = AirMass*0.2   #extinction_coefficient_mean
            mag[i] = mag[i] + mag_atm_att
            mag[i] = 10**(6.74-0.4*mag[i]) #here becomes photons
            quantum_efficiency = 0.3
            transmissivity = 0.5              
            multiplier = quantum_efficiency * transmissivity * exposure_time *0.5
            data[0].append(x)
            data[1].append(y)
            data[2].append(mag[i] * multiplier)
    return data, center

def sky_brightness(plate_scale, x_pix, y_pix, photo_filters, exposure_time, moon_phase=3):
    
    log.info(hist())
    '''
    moon_phase = int, optional
        moon_phase can go to 0 (new moon) to 3 (full moon) with intermedian phases 1 (7 day from new moon) and 2 (10 day from new moon)
    '''
    sub = DATA["Sky_brightness"]
    for n, i in enumerate(photo_filters):
        if i == "g":
            photo_filters[n] = "V"
    sky_brightness = [sub[k] for k in photo_filters if k in sub]
    x_arcsec = plate_scale * x_pix
    y_arcsec = plate_scale * y_pix
    Area = x_arcsec * y_arcsec
    magnitudo_ab_sky = []
    magnitudo_ab_sky.append(sky_brightness[0][moon_phase] - 2.5*np.log10(Area))
    photons = magnitudo_to_photons(magnitudo_ab_sky, photo_filters)
    electrons = photons_to_electrons(photons, photo_filters)
    tot = sum(electrons)*exposure_time
    return tot

def main():
    data = query()
    print(data)

if __name__ == '__main__':
    main()
