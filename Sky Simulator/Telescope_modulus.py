import numpy as np
import matplotlib.pyplot as plt
from scipy import fft, signal
from aotools import zernikeArray as ZA
from aotools.turbulence.infinitephasescreen import PhaseScreenKolmogorov

from Logger import hist
from astropy.logger import log

scale = 10 #data for the atmospheric turbolence
D = 8
r0 = 0.164 #(tipycally between 0.10 and 0.20)
L0 = 100
modes = 10

def phase_aberration(aperture, scale, D, r0, L0, pupil):
    '''
    Creates an image with the Kolmogorov algorithm
    with the phase atmospheric aberration.
    This function works but usually is not used
    because the telescope is limited by the seeing.

    Parameters
    ----------

    aperture : int
        The aperture of the telescope in mm

    scale : int
        The scale of the image generated

    D : int #inutile e rindondante con aperture, da cambiare e togliere

    r0 : float
        The Fried parameter of the 'seeing'

    L0 : float, int
        The outer scale of the 'seeing'

    pupil : numpy ndarray
        The image of the aperture


    Output
    ------

    phase_screen : numpy ndarray
        The image of the turbolence overlapped with the image of the aperture
    '''

    log.info(hist())
        
    nx_size= aperture//scale
    plx_scale = D/nx_size
    phase_screen = PhaseScreenKolmogorov(nx_size, plx_scale, r0, L0)
    phase_screen.add_row()
    phase_screen = phase_screen.scrn
    phase_screen = reshape_aber(phase_screen, scale, pupil)
    return phase_screen

def intensity_aberration(aperture, scale, modes, pupil):
    '''Creates the amplitude atmospheric aberration.

    Parameters
    ----------

    aperture : int
        The aperture of the telescope in mm

    scale : int
        The scale of the image generated

    modes : int
        Number of zernikle terms to use 

    pupil : numpy ndarray
        The image of the aperture


    Output
    ------

    zernike : numpy ndarray
        The image of the aberration overlapped with the image of the aperture
    '''

    log.info(hist())
        
    nx_size= aperture//scale
    zernike = []
    zernike_array = ZA(modes, nx_size)
    for i in range(modes):
        zernike.append(reshape_aber(zernike_array[i], scale, pupil))
    return zernike

def reshape_aber(image, scale, pupil):
    '''
    Adapt the aberration to the pupil size

    Parameters
    ----------

    image : int
        The aperture of the telescope in mm

    scale : int
        The scale of the image generated

    pupil : numpy ndarray
        The image of the aperture


    Output
    ------

    image : numpy ndarray
        The image of the aberration with the right shape to be combinated
        with the aperture    
    '''

    log.info(hist())
        
    image = np.repeat(np.repeat(image, scale, axis = 0), scale, axis = 1)
    units= int((np.shape(pupil)[0]-np.shape(image)[0])/2)
    image = np.pad(image,(units), mode = 'constant')
    return image

def plate_scale(pixel, focal_lenght):
    '''
    Convert to arcsec/pixel

    Parameters
    ----------

    pixel : float
        The length of the pixel in mm

    focal_lenght : int
        The focal length of the telescope in mm

    Output
    ------

    CCD_r : float
        The resolution of the CCD in arcsec/pixel 
    '''
    CCD_r = 206265*(pixel)/focal_lenght
    return CCD_r
    
def converter_to_pixel(CCD_resolution, focal_lenght, quantity_to_convert):  #TO REVISIONATE
    '''Convert a quantity on the aperture to pixels on CCD '''
    log.info(hist())
        
    converter_rad_arc = 180*3600/np.pi
    quantity_to_convert = quantity_to_convert
    quantity_angle = np.arctan(quantity_to_convert/focal_lenght) * converter_rad_arc
    Q_on_CCD = quantity_angle / CCD_resolution
    return Q_on_CCD

def telescope(focal_lenght, aperture, obstruction, trellis=True):
    '''
    Creates the figures of the telescope

    Parameters
    ----------

    focal_lenght : int
        The focal length of the telescope in mm

    aperture : int
        The aperture of the telescope in mm

    obstruction : int
        The central obstraction of the telescope in mm

    trellis : bool, optional
        If trellis == True the structure that hold on the obstruction is drawn,
        else the structure is not drawn

    Output
    ------

    pupil : numpy ndarray
        The image of the aperture

    r : numpy ndarray
        The radial variable on the pupil   
    '''
    log.info(hist())

    x,y = np.mgrid[-aperture:aperture, -aperture:aperture] # creates the 2D grid for the 2D function *4
    r = np.sqrt(x**2+y**2)
    print(type(r))
    pupil = np.piecewise(r, [r < aperture/2, r > aperture/2, r < obstruction/2], [1, 0, 0]) #creates the aperture
    if trellis:
        structure_x = np.piecewise(x, [x, x>5, x<-5], [0,1,1]) #creates the structure that keep the obstruction
        structure_y = np.piecewise(y, [y, y>5, y<-5], [0,1,1])
        pupil = pupil*(structure_x*structure_y)
    return pupil, r

def physical_CCD(CCD_structure):
    '''
    Generates a CCD with a given structure

    Parameters
    ----------

    CCD_structure : array
        The varies entries must have the measure on the CCD
        [0] : int, float; the x length of the CCD in mm
        [1] : int, float; the y length of the CCD in mm
        [2] : float     ; the length of a pixel in mm

    Outputs
    -------

    CCD : numpy ndarray
        The image of the CCD

    x : numpy ndarray
        The x variable on the CCD

    y : numpy ndarray
        The y variable on the CCD   
    '''
    log.info(hist())
        
    len_y = int(CCD_structure[0]/CCD_structure[2])//2
    len_x = int(CCD_structure[1]/CCD_structure[2])//2
    x, y = np.mgrid[-len_x:len_x,-len_y:len_y]
    r = np.sqrt(x**2+y**2)
    CCD = np.piecewise(r, [r], [0])
    return CCD, x, y    

def sky_background_aperture(focal_lenght, aperture, obstruction, trellis=True, atmosphere=False):
    '''
    Calls telescope, intensity_aberration and phase_aberration to create the optical
    figure at the aperure
    
    Parameters
    ----------

    focal_lenght : int
        The focal length of the telescope in mm

    aperture : int
        The aperture of the telescope in mm

    obstruction : int
        The central obstraction of the telescope in mm

    trellis : bool, optional
        If trellis == True the structure that hold on the obstruction is drawn,
        else the structure is not drawn

    atmosphere : bool, optional
        If atmosphere == True intensity_aberration and phase_aberration are called,
        else they are not

    Output
    ------

    pupil : numpy ndarray
        The image of the aperture with all the aberrations

    r : numpy ndarray
        The radial variable on the pupil
    '''
    log.info(hist())
        
    pupil, r = telescope(focal_lenght, aperture, obstruction, trellis)

    if atmosphere:
        zer = intensity_aberration(aperture, scale, modes, pupil)
        zernike = zer[0]

        for i in range(modes-1):
            a = np.random.random()
            zernike += (a/2)*zer[i]

        phase = phase_aberration(aperture, scale, D, r0, L0, pupil)
        kernel = pupil * zernike
        phase = np.exp(1j*phase*0.01+0j)
        pupil = kernel*phase
    return pupil, r

def defocus(pupil, defocus_distance, r, wavelenght):
    '''
    Takes the image and created the FT at some distance

    Parameters
    ----------

    pupil : numpy ndarray
        The image of the aperture to transform

    defocus_distance : float
        The distance of the defocus in mm in range(-2.5,2.5)

    r : numpy ndarray
        The radial variable on the pupil

    wavelenght : float
        The center wavelength in mm

    Output
    ------

    image : numpy ndarray
        The FT of the pupil at defocus_distance from the focus plane, AKA the PSF
    '''
    log.info(hist())
        
    phaseAngle = 1j*defocus_distance*np.sqrt((2*np.pi/wavelenght)**2-r**2+0j) #unnecessary 0j but keeping it for complex reasons
    kernel = np.exp(phaseAngle)
    defocusPupil = pupil * kernel
    defocusPSFA = fft.fftshift(fft.fft2(defocusPupil))
    image = np.abs(defocusPSFA)
    return image

def image_processing(image, units, aperture):
    '''
    This function takes an image and some parameters to obtain the same image cut and strechetd to fit the CCD.
    Also, this function takes only the essentioal information thus reducing the total weight of the final image 

    Parameters
    ----------

    image : numpy ndarray
        The image, usually the PSF

    units : float
        The convertion coefficient to a pixel on the aperture plane to the CCD

    aperture : int
        The aperture of the telescope in mm

    Output
    ------

    np.abs(image) : numpy ndarray
        The image cleaned with only the good parts and with the right measure
    '''
    log.info(hist())
        
    zoom = (aperture-8, aperture+8) #*4
    image = image[zoom[0]:zoom[1], zoom[0]:zoom[1]]  #takes only the good part
    image = np.repeat(np.repeat(image,units, axis=0), units, axis=1) #Strech the image to fit the CCD scale
    image = image**2  #the intensity is the amplitude squared
    image = image/(np.max(image))#TODO change to an integral normalization
    return np.abs(image)

def telescope_on_CCD(CCD_resolution, focal_lenght, aperture, obstruction, defocus_distance, wavelenght, trellis=True, atmosphere=False):
    '''
    This function calls sky_backgroud_aperture, defocus and image_processing to create the sample on the PSF on the CCD

    Parameters
    ----------

    CCD_resolution : float
        The resolution on the CCD in arcsec/pixel
    
    focal_lenght : int
        The focal length of the telescope in mm

    aperture : int
        The aperture of the telescope in mm

    obstruction : int
        The central obstraction of the telescope in mm

    defocus_distance : float
        The distance of the defocus in mm in range(-2.5,2.5)

    wavelenght : float
        The center wavelength in mm

    trellis : bool, optional
        If trellis == True the structure that hold on the obstruction is drawn,
        else the structure is not drawn

    atmosphere : bool, optional
        If atmosphere == True intensity_aberration and phase_aberration are called,
        else they are not

    Output
    ------

    image : numpy ndarray
        The image of the PSF
    '''
    log.info(hist())
        
    #units =  converter_to_pixel(CCD_resolution, focal_lenght, 1) #convert 1mm on the aperture in pixel on CCD
    units = CCD_resolution//0.12 #empiric value
    image, r = sky_background_aperture(focal_lenght, aperture, obstruction, trellis, atmosphere) #creates the aperture 
    image = defocus(image, defocus_distance, r, wavelenght) #creates the image on the screen
    image = image_processing(image, units, aperture) 
    return image
