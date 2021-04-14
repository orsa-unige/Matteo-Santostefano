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


def load_measure():
    '''
    Load the basic data of the telescope and the CCD from a JSON file

    Outputs
    -------
    telescope_structure : list
        a list that contains the information about the telescope
        telescope_structure[0] = f_l : int
            It is the focal length of the telescope in mm            
        telescope_structure[1] = ape : int
            It is the aperture of the telescope in mm
        telescope_structure[2] = obs : int
            It is the central obstruction of the telescope in mm
        telescope_structure[3] = wav : float
            It is the central wavelength of the sensible spectrum of the CCD in mm

    ccd_structure : array
        a array that contains the information about the CCD
        ccd_structure[0] = C_x : float
            It is the length of the CCD along the x axis in mm
        ccd_structure[0] = C_y : float
            It is the length of the CCD along the y axis in mm
        ccd_structure[0] = pix : float
            It is the lenght of a single pixel in mm

    ccd_data : list
        a list that cointains information about the errors generated in the CCD
        ccd_data[0] = gain : float
            the gain of the CCD
        ccd_data[1] = read_out_electrons : float
            the electrons that generate the read out noise
    '''
    log.info(hist())

    #filename = "Antola_data.json"
    filename = "San_Pedro_data.json"
    f = open(filename, "r")
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
    gain = CCD['gain']
    read_out_electrons = CCD["read_out_electrons"]
    f.close()
    telescope_structure = (f_l, ape, obs, wav)
    ccd_structure = [C_x, C_y, pix]
    ccd_data = (gain, read_out_electrons)
    return telescope_structure, ccd_structure, ccd_data


def save_fits(data, defocus_distance, center, choose):
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

    choose : int
        is a variable that determines the star used
    
    '''
    log.info(hist())
 
    name = ['Kolmogorov', 'Speckle_Sum', 'Simple_seeing']
    choose = choose-1
    hdul = fits.PrimaryHDU(data) #save the file
    hdr = hdul.header
    hdr['RA'] = (center[0], "Right Ascension in decimal hours" )
    hdr['DEC'] = (center[1], "Declination in decimal degrees")
    hdr['IMGTYPE'] = 'object'
    hdr['IMAGETYP'] = 'object'
    hdr['BZERO'] = 32768
    hdr['DEFOCUS'] = (defocus_distance, "[mm] Distance from focal plane")
    hdul.scale('int16', bzero=32768)
    hdul.writeto(f'try_{name[choose]}.fits', overwrite = True)

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
    #photo_filters = ['U', 'B', 'V', 'R', 'I']
    photo_filters = ['V']

    exposure_time = 60 #second
    default_defocus = 0.9 #mm
    default_coordinates = "07 59 08.445 +15 24 42.00"
    default_seeing = 3

    choose = int(input('do you want Kolmogorov (1), spekle method(2) or gaussian approximation(3)? \n'))
    defocus_distance = float(input(f'The defocus distance? [mm] Default: {default_defocus} \n') or default_defocus) #mm
    coordinates = input(f'Coordinates? Default: {default_coordinates} \n') or default_coordinates
    seeing = float(input(f'The seeing? [arcsec] Default: {default_seeing} \n') or default_seeing) 

    telescope_structure, ccd_structure, ccd_data = load_measure()

    #Standard data for the noise formation 
    realistic = True
    dark = 1.04
    bias_level = 2000
    
    gain = ccd_data[0]
    read_noise_electrons = ccd_data[1]

    ccd_structure[2] = ccd_structure[2] * binning
    ccd_structure.append(Tm.plate_scale(ccd_structure[2], telescope_structure[0])) #calculate the resolution of the CCD arcsec/pixel

    if choose == 1:
        CCD_sample = Tm.telescope_on_CCD(ccd_structure[3], binning, telescope_structure, defocus_distance, True, True)
    else:
        CCD_sample = Tm.telescope_on_CCD(ccd_structure[3], binning, telescope_structure, defocus_distance, True, False) #generates the sample of a star with the right measure

    CCD = Tm.physical_CCD(ccd_structure) #generate the CCD
    size = CCD.shape

    if choose == 2:
        CCD_seeing = seeing/(ccd_structure[3]*2)
        number_of_spekle = 100 * exposure_time
        seeing_image = Tm.main_spekle(number_of_spekle, CCD_seeing)
        seeing_image = seeing_image/number_of_spekle
    elif choose == 3:
        CCD_seeing = seeing/(ccd_structure[3]*2)
        seeing_image = Tm.seeing(CCD_seeing)

    sky, center = Qm.query(coordinates, photo_filters, ccd_structure, exposure_time) #call a function that gives back positions, fluxs of the stars and a data for the header
    sky_counts = Qm.sky_brightness(ccd_structure[3], size[0], size[1], photo_filters, exposure_time)

    photons_collection_area = (np.pi/400)*(telescope_structure[1]**2-telescope_structure[2]**2)

    if choose == 1:
        multiplier = 2 * gain * photons_collection_area *0.05
    else:
        multiplier = 1 * gain * photons_collection_area *0.05 #photons_collection_area * gain * binning**2 #* 200 #I don't know where 200 cames from I'm investigating
                                    #the flux is in ph cm^-2 so it has to be multiplied for the effective area of the aperure in cm^2

    for i in range(len(sky[0])):
        if sky[0][i] >=0 and sky[0][i] <=size[0]:
            if sky[1][i] >=0 and sky[1][i] <=size[1]:
                CCD[int(sky[0][i])][int(sky[1][i])] = sky[2][i] * multiplier

    if choose == 2 or choose == 3:
                CCD = signal.fftconvolve(CCD, seeing_image, mode='same')       

    CCD = signal.fftconvolve(CCD, CCD_sample, mode='full')  #convolve the position with the sample

    if realistic: 

        flat = ccd_mod.sensitivity_variations(CCD, vignetting=True, dust=False)
        bias_only = ccd_mod.bias(CCD, bias_level, realistic=True)
        noise_only = ccd_mod.read_out_noise(CCD, read_noise_electrons, gain=gain) 
        dark_only = ccd_mod.dark_current(CCD, dark, exposure_time, gain=gain)
        sky_only = ccd_mod.sky_background(CCD, sky_counts, gain=gain)
        #cosmic_rays = ccd_mod.make_cosmic_rays(CCD, np.random.randint(10,30))

        CCD = bias_only + noise_only + dark_only + flat * (sky_only + CCD)
        CCD = ccd_mod.saturation_controll(CCD)
    
    save_fits(CCD, defocus_distance, center, choose) #save the fits
    

if __name__ == '__main__':
    main()
