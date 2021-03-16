import query_mod as Qm
import telescope_mod as Tm
import json
import matplotlib.pyplot as plt
import numpy as np
from scipy import signal
from astropy.io import fits
import ccd_mod 

from Logger import hist
from astropy.logger import log

"""
Prova push
"""

def load_measure():
    '''
    Load the basic data of the telescope and the CCD from a JSON file

    Outputs
    -------

    f_l : int
        It is the focal length of the telescope in mm

    ape : int
        It is the aperture of the telescope in mm

    obs : int
        It is the central obstruction of the telescope in mm

    wav : float
        It is the central wavelength of the sensible spectrum of the CCD in mm

    C_x : float
        It is the length of the CCD along the x axis in mm

    C_y : float
        It is the length of the CCD along the y axis in mm

    pix : float
        It is the lenght of a single pixel in mm
    '''
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

def save_fits(data, defocus_distance, center):
    '''
    Save a fits with a simple header.

    Parameters
    ----------

    data : numpy array
        Image to save

    defocus_distance : float
        Information about the defocus distance to insert in the header

    center : tuple, array
        Invormation about the coordinate to insert in the header
    
    '''
    log.info(hist())
    
    hdul = fits.PrimaryHDU(data) #save the file
    hdr = hdul.header
    hdr['RA'] = (center[0], "Right Ascension in decimal hours" )
    hdr['DEC'] = (center[1], "Declination in decimal degrees")
    hdr['IMGTYPE'] = 'object'
    hdr['IMAGETYP'] = 'object'
    hdr['BZERO'] = 32768
    hdr['DEFOCUS'] = (defocus_distance, "[mm] Distance from focal plane")
    hdul.scale('int16', bzero=32768)
    hdul.writeto(f'try.fits', overwrite = True)
 

def main():
    '''
    This is the main function of the system, it coordinates the Telescope_mod, the ccd_mod and the Query_mod.
    First it calls load_measure to obtain some data for Telescope_mod and Query_mod, the uses Telescope_mod
    with some additional parameters in input to create the PSF with the right measure for the telescope and the CCD.
    
    At this point the function calls again the Telescope_mod to creats an empty image of the CCD, and the Query_mod
    to obtain information about the position and the flux of the stars in a specific reagion of the sky,
    the center of that is an imput from the operator.

    The function uses the information of Query_mod to "light_up" some pixel on the CCD's image in position where shoud be the stars
    and then convolve this updated CCD with the calculated PSF from Telescope_mod obtainig the CCD with the stars drawn.

    If the value of realistic is True, the main calls the ccd_mod to creates the bias frame, the dark frame, the read out noise,
    the backround noise and the flat frame and then combines all the frames to obtain a more realistic image on the CCD.    
    '''
    log.info(hist())
    binning = 2
    photo_filters = ['U', 'B', 'V', 'R', 'I']
    gain = 2.0
    exposure_time = 600 #second
    default_defocus = 0.9 #mm
    default_coordinates = "07 59 05.83 +15 23 29.24"

    #Standard data for the noise formation 
    realistic = True
    dark = 0.1
    sky_counts = 20
    bias_level = 1100
    read_noise_electrons = 5

    defocus_distance = float(input(f'The defocus distance? [mm] Default: {default_defocus} \n') or default_defocus) #mm
    coordinates = input(f'Coordinates? Default: {default_coordinates} \n') or default_coordinates

    gain = gain**binning

    f_l, ape, obs, wav, C_x, C_y, pix = load_measure()  #load the measure
    pix = pix * binning
    C_r = Tm.plate_scale(pix, f_l)  #calculate the resolution of the CCD arcsec/pixel
    CCD_structure = (C_x, C_y, pix, C_r) #put the information of the CCD in a array
    CCD_sample = Tm.telescope_on_CCD(CCD_structure[3], f_l, ape, obs, defocus_distance, wav, True, True) #generates the sample of a star with the right measure
    CCD, x, y = Tm.physical_CCD(CCD_structure) #generate the CCD
    size = CCD.shape
    sky, center = Qm.query(coordinates, photo_filters, CCD_structure, exposure_time) #call a function that gives back positions, fluxs of the stars and a data for the header

    p_x = sky[0]
    p_y = sky[1] 
    flux = sky[2]

    for i in range(len(p_x)):
        if p_x[i] >=0 and p_x[i] <=size[0]:
            if p_y[i] >=0 and p_y[i] <=size[1]:
                CCD[int(p_x[i]), int(p_y[i])] = flux[i]*gain #bind the multiplier to the gain
                
    CCD = signal.fftconvolve(CCD, CCD_sample, mode='same')  #convolve the position with the sample

    if realistic: 


        flat = ccd_mod.sensitivity_variations(CCD, vignetting=True, dust=False)
        bias_only = ccd_mod.bias(CCD, bias_level, realistic=True)
        noise_only = ccd_mod.read_out_noise(CCD, read_noise_electrons, gain=gain) 
        dark_only = ccd_mod.dark_current(CCD, dark, exposure_time, gain=gain)
        sky_only = ccd_mod.sky_background(CCD, sky_counts, gain=gain)

        CCD = bias_only + noise_only + dark_only + flat * (sky_only + CCD)

        CCD = ccd_mod.saturation_controll(CCD)
    
    save_fits(CCD, defocus_distance, center) #save the fits

if __name__ == '__main__':
    main()
