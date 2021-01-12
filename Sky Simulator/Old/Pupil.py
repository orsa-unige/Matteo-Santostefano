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

def Shower(defocusPSF ,imgSize, devc, vmax, Airy, bool=False): #Shows a psf with its profile
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

def Defocusor(imgSize, wavelenght, structure): #creates the defocused psf
    pupil = structure[0] #it use the data of the pupil
    KX = structure[1]
    KY = structure[2]
    dk = structure[3]
    defocusDistance = np.arange(-2.75, 0.1, 0.25) # Creates an arry of test distances
    defocusPSF = np.zeros((imgSize, imgSize, defocusDistance.size)) #creates an array of defocused psf
    
    # Classic Integration through Fourier Transformation
    for ctr, z in enumerate(defocusDistance): #a control number and the actual distance
        phaseAngle = 1j* z*np.sqrt((2*np.pi/wavelength)**2-KX**2-KY**2+0j) #unnecessary 0j but keeps it to make sure sqrt knows is a complex
        kernel = np.exp(phaseAngle) 
        defocusPupil = pupil * kernel #bind the amplitude with the phase
        defocusPSFA = fft.fftshift(fft.fft2(defocusPupil))*dk**2 #creates the psf
        defocusPSF[:,:,ctr] = np.abs(defocusPSFA) #adds the psf at the array of psfs

    vmax = np.max(defocusPSF[:,:,10])  #find the maximum in PSF (in focus) to scale the others (can be eliminated)

    return defocusDistance, defocusPSF, vmax

def PupilStructure(wavelength, NA, pixelSize, imgSize):
    # Working with FFT a new set of variables is created
    kMax = 2*np.pi/pixelSize # Value of k at the maximum extent of the pupil function
    kNA = 2*np.pi*NA/wavelength
    dk = 2*np.pi/(imgSize*pixelSize)

    kx = np.arange((-kMax+dk)/2, (kMax+dk)/2, dk)
    ky = np.arange((-kMax+dk)/2, (kMax+dk)/2, dk)
    KX, KY = np.meshgrid(kx, ky) # Coordinate system for the pupil function

    maskRadius = kNA/dk # Radius of amplitude mask for defining the pupil
    maskCenter = np.floor(imgSize/2)
    W, H = np.meshgrid(np.arange(0, imgSize), np.arange(0, imgSize))
    mask = np.sqrt((W-maskCenter)**2 + (H-maskCenter)**2) < maskRadius #create the figure of a circular dot

    plt.imshow(mask, extent = (kx.min(), kx.max(), ky.min(), ky.max()))
    plt.xlabel('kx')
    plt.ylabel('ky')
    plt.colorbar()
    plt.show()

    amp = np.ones((imgSize, imgSize))*mask #creates the amplitude
    phase = 2j*np.pi*np.ones((imgSize, imgSize)) #creates the phase angle
    pupil = amp*np.exp(phase) #Bind the amplitude with the phase in order to obtain the pupil
    structure = [pupil, KX, KY, dk] #the structural parameter of the pupil
    psfA_Unaberrated = fft.fftshift(fft.fft2(pupil))*dk**2 #the psf of the pipul
    psf = np.abs(psfA_Unaberrated) #the absolute value of the psf
    return psf, structure

def main():
    psf, structure = PupilStructure(wavelength, NA, pixelSize, imgSize) #creates a psf and the structure of the pupil (dimensions)
    vmax = np.max(psf) #takes the maximum of the psf only for the comparison with the intensity of the defocused figure
    Airy_profile= Take_profile(imgSize, psf) #takes the profile of the psf for comparison's purpose
    Shower(psf, imgSize, 0, vmax, Airy_profile, True) #Shows the psf near its profile
    defocusDistance, defocusPSF, maximum=Defocusor(imgSize, wavelength, structure) #creates the defocused figured
    for i in range (defocusDistance.size):
        Shower(defocusPSF[:,:, i],imgSize,defocusDistance[i], maximum, Airy_profile) #shows the defocused psfs near their profile 

if __name__ =="__main__":
    main()
