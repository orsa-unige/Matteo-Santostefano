from astroquery.simbad import Simbad
import astropy.coordinates as coord
from astropy import wcs
from astropy import units as u
import numpy as np
import json

from Logger import hist
from astropy.logger import log


def Coordinator(c, center, CCD_res):
    log.info(hist())
    
    #offset_x = CCD_structure[0]/CCD_structure[2]
    #offset_y = CCD_structure[1]/CCD_structure[2]
    w = wcs.WCS(naxis=2)
    w.wcs.crpix = [1, 1]
    w.wcs.crval = [center[0], center[1]]
    w.wcs.ctype = ["RA", "DEC"]    
    x,y = wcs.utils.skycoord_to_pixel(c, w)
    x = (x*3600)/CCD_res #+ offset_x/2 #(deg*arc/sec*deg)*arc/sec*pixel+offset
    y = (y*3600)/CCD_res #+ offset_y/2

    return x, y

def data_loader(key, photo_filter):
    log.info(hist())
    
    f = open("Antola_data.json", "r")
    data = json.load(f)
    key_data = data[key]
    datum = []
    for i in photo_filter:
          datum.append(key_data[str(i)])

    f.close()
    return datum

def VEGA_to_AB(magnitudo, photo_filter):
    log.info(hist())
    
    convertion_table = data_loader("Convertion_table", photo_filter)
    for i in range(len(magnitudo)):
        if magnitudo[i]== 0:
            magnitudo[i]=0
        else:
            magnitudo[i] = magnitudo[i] - convertion_table[i]
    return magnitudo

def atmospheric_attenuation(magnitudo, photo_filter, AirMass):
    log.info(hist())
    
    extinction_coefficient = data_loader("Extinction_coefficient", photo_filter)
    for i in range(len(magnitudo)):
        if magnitudo[i] == 0:
            magnitudo[i] = 0
        else:
            magnitudo[i] = magnitudo[i]-AirMass*extinction_coefficient[i]
    return magnitudo

def magnitudo_to_photons(magnitudo, photo_filter):
    log.info(hist())
    
    central_wavelenght = data_loader("Central_wavelenght", photo_filter) #in Armstrong
    number = magnitudo
    for i in range(len(magnitudo)):
        if magnitudo[i] == 0:
            number[i] = 0
        else:
            exp = 6.74-0.4*magnitudo[i]
            number[i] = (10**exp)/(central_wavelenght[i])
    return number

def photons_to_electrons(photons, photo_filter):
    log.info(hist())
    
    QE = data_loader("Quantum_efficiency", photo_filter) #in Armstrong
    electrons = photons
    for i in range(len(photons)):
        electrons[i] = photons[i]*QE[i]
    return electrons

def Data_structure():
    log.info(hist())
    
    Coord_x = []
    Coord_y = []
    Flux_tot = []
    data = [Coord_x, Coord_y, Flux_tot]
    return data

def radius_sky_portion(CCD_structure):
    log.info(hist())
    
    pixel_on_x = CCD_structure[0]/CCD_structure[2]
    pixel_on_y = CCD_structure[1]/CCD_structure[2]
    diagonal_diameter = np.sqrt((pixel_on_x)**2+(pixel_on_y)**2)
    radius_on_pixel = diagonal_diameter/2
    radius = radius_on_pixel*CCD_structure[3]/60
    return (radius/60)*u.deg

def magnitudo_to_electrons(magnitudo, photo_filters, AirMass, exposure_time):
    log.info(hist())
    
    magnitudo = VEGA_to_AB(magnitudo, photo_filters)
    magnitudo = atmospheric_attenuation(magnitudo, photo_filters, 1)
    photons = magnitudo_to_photons(magnitudo, photo_filters)
    electrons = photons_to_electrons(photons, photo_filters)
    tot = sum(electrons)*exposure_time
    return tot
    
def query(coordi, photo_filters, CCD_structure, exposure_time):
    log.info(hist())
    
    radius = radius_sky_portion(CCD_structure)
    coordi = coord.SkyCoord(coordi, frame = 'icrs', unit=(u.hourangle, u.deg))
    f =['flux({0})'.format(x) for x in photo_filters]
    [Simbad.add_votable_fields(g) for g in f]
    result_table = Simbad.query_region(coordi, radius)
    
    center = [0,0]
    center[0] = coordi.ra.degree
    center[1] = coordi.dec.degree

    data = Data_structure()

    Stars_number = len(result_table)
    len_filter_list = len(photo_filters)
    #exposure_time = 60
    AirMass = 1
    
    for i in range (Stars_number):
        row = result_table[i]
        tot = 0
        magnitudo = []
        for j in range(11, 11+len_filter_list):  #sum the fluxies
            if row[j] is np.ma.masked:
                row[j] = 0
            magnitudo.append(row[j])
        
        tot = magnitudo_to_electrons(magnitudo, photo_filters, AirMass, exposure_time)
    
        if tot != 0: #if total flux is 0 the object is not passed
            c = coord.SkyCoord(row[1], row[2], frame = 'icrs', unit=(u.hourangle, u.deg))##change here
            x, y = Coordinator(c, center, CCD_structure[3])
            data[0].append(x)
            data[1].append(y)
            data[2].append(tot)
    return data, center



def main():
    data = query()
    print(data)

if __name__ == '__main__':
    main()
