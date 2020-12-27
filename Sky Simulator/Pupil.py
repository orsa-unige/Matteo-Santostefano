from scipy import fft
from scipy.special import j1 #Bessel function
import matplotlib.pyplot as plt
import numpy as np

#-------------Global parameter set-------------------#
plt.style.use('dark_background')
plt.rcParams['image.cmap'] = 'plasma' #'gray'
zoom = 40  # Visualized pixels
imgSize = 513 # Odd number to preserve rotational symmetry about the center

wavelength = 1 #"dimensionless"
NA = 0.8 # Numerical aperture of the objective works well in range(1 to 0.5)
pixelSize  = 0.1 # Arcsec per Pixel (the right value isn't so important)


def Take_profile(imgSize, psf): #Create the profile of a PSF
    Profile = np.real(psf[int(np.floor(imgSize/2)), int(np.floor(imgSize/2)-20) : int(np.floor(imgSize/2)+21)])
    return Profile

def Shower(defocusPSF ,imgSize, devc, vmax, Airy, bool=False):
    Name = f'Defocus: {devc} lambda'
    if bool:
        Name ='In-focus'
    fig, (ax, ax1) = plt.subplots(1,2)
    ax.imshow(defocusPSF, vmin =0, vmax = vmax, interpolation = 'nearest')
    Profile = Take_profile(imgSize, defocusPSF)
    ax.set_xlim((np.floor(imgSize/2)-zoom/2, np.floor(imgSize/2)+zoom/2))
    ax.set_ylim((np.floor(imgSize/2)-zoom/2, np.floor(imgSize/2)+zoom/2))
    ax.set_xlabel('x, pixels')
    ax.set_ylabel('y, pixels')
    ax.set_title(Name)
    ax1.plot(Airy, linewidth = 2,  label = 'Airy PSF')
    if not bool:
        ax1.plot(Profile, linewidth=2,  label = 'Defocus PSF')
    ax1.set_xlabel('Pixels')
    ax1.set_ylabel('Irradiance')
    ax1.legend()
    plt.tight_layout()
    plt.show()

def Defocusor(imgSize, wavelenght, structure):
    pupil = structure[0]
    KX = structure[1]
    KY = structure[2]
    dk = structure[3]
    defocusDistance = np.arange(-2.75, 0.1, 0.25) # Creates an arry of test distances
    defocusPSF = np.zeros((imgSize, imgSize, defocusDistance.size))
    
    # Classic Integration through Fourier Transformation
    for ctr, z in enumerate(defocusDistance): 
        phaseAngle = 1j* z*np.sqrt((2*np.pi/wavelength)**2-KX**2-KY**2+0j) #unnecessary 0j but keeping it for complex reasons
        kernel = np.exp(phaseAngle)
        defocusPupil = pupil * kernel
        defocusPSFA = fft.fftshift(fft.fft2(defocusPupil))*dk**2
        defocusPSF[:,:,ctr] = np.abs(defocusPSFA)

    vmax = np.max(defocusPSF[:,:,10])  #find the maximum in PSF (in focus) to scale the others (can be eliminated)

    return defocusDistance, defocusPSF, vmax

def PupilStructure(wavelength, NA, pixelSize, imgSize):
    # Working with FFT a new set of
    kMax = 2*np.pi/pixelSize # Value of k at the maximum extent of the pupil function
    kNA = 2*np.pi*NA/wavelength
    dk = 2*np.pi/(imgSize*pixelSize)

    kx = np.arange((-kMax+dk)/2, (kMax+dk)/2, dk)
    ky = np.arange((-kMax+dk)/2, (kMax+dk)/2, dk)
    KX, KY = np.meshgrid(kx, ky) # Coordinate system for the pupil function

    maskRadius = kNA/dk # Radius of amplitude mask for defining the pupil
    maskCenter = np.floor(imgSize/2)
    W, H = np.meshgrid(np.arange(0, imgSize), np.arange(0, imgSize))
    mask = np.sqrt((W-maskCenter)**2 + (H-maskCenter)**2) < maskRadius

    plt.imshow(mask, extent = (kx.min(), kx.max(), ky.min(), ky.max()))
    plt.xlabel('kx')
    plt.ylabel('ky')
    plt.colorbar()
    plt.show()

    amp = np.ones((imgSize, imgSize))*mask
    phase = 2j*np.pi*np.ones((imgSize, imgSize))
    pupil = amp*np.exp(phase)
    structure = [pupil, KX, KY, dk]
    psfA_Unaberrated = fft.fftshift(fft.fft2(pupil))*dk**2
    psf = np.abs(psfA_Unaberrated)
    return psf, structure

def main():
    psf, structure = PupilStructure(wavelength, NA, pixelSize, imgSize)
    vmax = np.max(psf)
    Airy_profile= Take_profile(imgSize, psf)
    Shower(psf, imgSize, 0, vmax, Airy_profile, True)
    defocusDistance, defocusPSF, maximum=Defocusor(imgSize, wavelength, structure)
    for i in range (defocusDistance.size):
        Shower(defocusPSF[:,:, i],imgSize,defocusDistance[i], maximum, Airy_profile)

if __name__ =="__main__":
    main()