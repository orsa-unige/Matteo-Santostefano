import numpy as np
import matplotlib.pyplot as plt
from scipy import fft, signal
from aotools import zernikeArray as ZA
from aotools.turbulence.infinitephasescreen import PhaseScreenKolmogorov

scale = 10 #data for the atmospheric turbolence
D = 8
r0 = 0.164 #(tipycally between 0.10 and 0.20)
L0 = 100
modes = 10

def phase_aberration(aperture, scale, D, r0, L0, pupil): #creates the phase atmospheric aberration
    nx_size= aperture//scale
    plx_scale = D/nx_size
    phase_screen = PhaseScreenKolmogorov(nx_size, plx_scale, r0, L0)
    phase_screen.add_row()
    phase_screen = phase_screen.scrn
    phase_screen = reshape_aber(phase_screen, scale, pupil)
    return phase_screen

def intensity_aberration(aperture, scale, modes, pupil): #creates the amplitude atmospheric aberration
    nx_size= aperture//scale
    zernike = []
    zernike_array = ZA(modes, nx_size)
    for i in range(modes):
        zernike.append(reshape_aber(zernike_array[i], scale, pupil))
    return zernike

def reshape_aber(zer_ar, scale, pupil): #adapt the aberration to the pupil size
    zer_ar = np.repeat(np.repeat(zer_ar, scale, axis = 0), scale, axis = 1)
    units= int((np.shape(pupil)[0]-np.shape(zer_ar)[0])/2)
    zer_ar = np.pad(zer_ar,(units), mode = 'constant')
    return zer_ar

def plate_scale(pixel, focal_lenght): #convert to arcsec/pixel
    CCD_r = 206265*(pixel)/focal_lenght
    return CCD_r
    
def converter_to_pixel(CCD_resolution, focal_lenght, quantity_to_convert): #convert a quantity on the aperture to pixels on CCD 
    converter_rad_arc = 180*3600/np.pi
    quantity_to_convert = quantity_to_convert/16
    quantity_angle = np.arctan(quantity_to_convert/focal_lenght) * converter_rad_arc
    Q_on_CCD = quantity_angle / CCD_resolution
    return Q_on_CCD

def telescope(focal_lenght, aperture, obstruction, trellis=True): #creates the figures of the telescope
    aperture_angle = np.arctan(aperture/focal_lenght)
    x,y = np.mgrid[-aperture*4:aperture*4, -aperture*4:aperture*4] # creates the 2D grid for the 2D function
    r = np.sqrt(x**2+y**2)
    pupil = np.piecewise(r, [r < aperture/2, r > aperture/2, r < obstruction/2], [1, 0, 0]) #creates the aperture
    if trellis:
        structure_x = np.piecewise(x, [x, x>5, x<-5], [0,1,1]) #creates the structure that keep the obstruction
        structure_y = np.piecewise(y, [y, y>5, y<-5], [0,1,1])
        pupil = pupil*(structure_x*structure_y)
    return pupil, r

def physical_CCD(CCD_structure):  #generates a CCD with a given structure
    len_x = int(CCD_structure[0]/CCD_structure[2])//2
    len_y = int(CCD_structure[1]/CCD_structure[2])//2
    x, y = np.mgrid[-len_x:len_x,-len_y:len_y]
    r = np.sqrt(x**2+y**2)
    CCD = np.piecewise(r, [r], [0])
    return CCD, x, y    

def sky_background_aperture(focal_lenght, aperture, obstruction, trellis, atmosphere): #unites the amperture with the optical aberrarion ##TODO phase aberration
    pupil, r = telescope(focal_lenght, aperture, obstruction, trellis)
    if atmosphere:
        zer = intensity_aberration(aperture, scale, modes, pupil)
        zernike = zer[0]
        for i in range(modes-1):
            a = np.random.random()
            zernike += (a/2)*zer[i]
        phase = phase_aberration(aperture, scale, D, r0, L0, pupil)
        kernel = pupil * zernike
        phase = np.exp(1j*phase*0.05+0j)
    pupil = kernel*phase
    return pupil, r

def defocus(pupil, defocus_distance, r, wavelenght): #takes the image and created the FT at some distance
    phaseAngle = 1j*defocus_distance*np.sqrt((2*np.pi/wavelenght)**2-r**2+0j) #unnecessary 0j but keeping it for complex reasons
    kernel = np.exp(phaseAngle)
    defocusPupil = pupil * kernel
    defocusPSFA = fft.fftshift(fft.fft2(defocusPupil))
    return np.abs(defocusPSFA)

def image_processing(image, units, aperture):
    zoom = (aperture*4-40, aperture*4+40)
    image = image/(np.max(image)/100) #normalize the image
    image = image[zoom[0]:zoom[1], zoom[0]:zoom[1]]  #takes only the good part
    image = np.repeat(np.repeat(image,units, axis=0), units, axis=1) #Strech the image to fit the CCD scale
    image = image**2  #the intensity is the amplitude squared    
    return np.abs(image)

def telescope_on_CCD(CCD_resolution, focal_lenght, aperture, obstruction, defocus_distance, wavelenght, trellis, atmosphere):
    units =  converter_to_pixel(CCD_resolution, focal_lenght, 1) #convert 1mm on the aperture in pixel on CCD
    units = units//2
    image, r = sky_background_aperture(focal_lenght, aperture, obstruction, trellis, atmosphere) #creates the aperture 
    image = defocus(image, defocus_distance, r, wavelenght) #creates the image on the screen
    image = image_processing(image, units, aperture) 
    return image







