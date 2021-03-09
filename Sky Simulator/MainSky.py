import Query_modulus as Qm
import Telescope_modulus as Tm
import json
import matplotlib.pyplot as plt
import numpy as np
from scipy import signal
from astropy.io import fits

from Logger import hist
from astropy.logger import log

def load_measure():
    '''Load the basic data of the telescope'''
    log.info(hist())
    
    f = open("Antola_data.json", "r")
    data = json.load(f)
    Telescope = data['Telescope']
    f_l = Telescope['focal_lenght']
    ape = Telescope['aperture']
    obs = Telescope['obstruction']
    wav = Telescope['wavelenght']
    CCD = data['CCD']
    C_x = CCD['CCD_x']
    C_y = CCD['CCD_y']
    pix = CCD['pixels']
    f.close()
    return f_l, ape, obs, wav, C_x, C_y, pix

def simple_point(x, y, amp=5.0, x_0=0, y_0=0):
    '''Creates a point in the CCD where will be a star'''
    log.info(hist())
    
    KX = x-x_0
    KY = y-y_0
    r=np.sqrt(KX**2+KY**2)
    point = np.piecewise(r, [r, r<2], [0,1])
    return point*amp

def save_fits(data, defocus_distance, center):
    '''Save a fits with a simple header'''
    log.info(hist())
    
    hdul = fits.PrimaryHDU(data) #save the file
    hdr = hdul.header
    hdr['RA'] = center[0]
    hdr['DEC'] = center[1]
    hdr['IMGTYPE'] = 'object'
    hdr['IMAGETYP'] = 'object'
    hdr['BZERO'] = 32768
    hdr['DEFOCUS'] = defocus_distance
    hdul.scale('int16', bzero=32768)
    hdul.writeto(f'try.fits')


def main():
    #"07 59 05.8395618539 +15 23 29.240025256"
    log.info(hist())
    
    photo_filters = ['U', 'B', 'V', 'R', 'I']
    default_defocus = 2.0
    defocus_distance = float(input(f'The defocus distance? [mm] Default: {default_defocus} \n') or default_defocus) #mm

    default_coordinates = "07 59 05.83 +15 23 29.24"
    coordinates = input(f'Coordinates? Default: {default_coordinates} \n' or default_coordinates)
    #defocus_distance = 0 #mm
    #coordinates = "07 59 05.8395618539 +15 23 29.240025256"

    f_l, ape, obs, wav, C_x, C_y, pix = load_measure()  #load the measure
    C_r = Tm.plate_scale(pix, f_l)  #calculate the resolution of the CCD arcsec/pixel
    CCD_structure = (C_x, C_y, pix, C_r) #put the information of the CCD in a array
    CCD_sample = Tm.telescope_on_CCD(CCD_structure[3], f_l, ape, obs, defocus_distance, wav, True, True) #generates the sample of a star with the right measure
    CCD, x, y = Tm.physical_CCD(CCD_structure) #generate the CCD

    sky, center = Qm.query(coordinates, photo_filters, CCD_structure) #call a function that gives back positions, fluxs of the stars and a data for the header 
    p_x = sky[0]
    p_y = sky[1] 
    flux = sky[2]
    data = CCD
    for i in range(len(p_x)):
        data +=  simple_point(x, y, flux[i], p_y[i], p_x[i]) #positioning of the stars 

    data = signal.fftconvolve(data, CCD_sample, mode='same')  #convolve the position with the sample

    save_fits(data, defocus_distance, center) #save the fits

if __name__ == '__main__':
    main()
