from astropy.modeling.models import custom_model
import numpy as np
import matplotlib.pyplot as plt
from scipy import fft


@custom_model
def modellino(x, y, amp=5.0, x_0=0.0, y_0=0.0, defocus_distance = 0.06):
    
    KX = x-x_0
    KY = y-y_0
    wavelength = 0.005 #value to change
    r=np.sqrt(KX**2+KY**2)
    
    pupil = np.piecewise(r, [r < 500, r > 500, r < 5], [1, 0, 0])#value to chance

    phaseAngle = 1j*defocus_distance*np.sqrt((2*np.pi/wavelength)**2-KX**2-KY**2+0j) #unnecessary 0j but keeping it for complex reasons
    kernel = np.exp(phaseAngle)
    defocusPupil = pupil * kernel

    defocusPSFA = fft.fftshift(fft.fft2((defocusPupil)))
    
    return amp*np.abs(defocusPSFA)
    

def test(): #visualise a test sample star with 
    w,z = np.mgrid[0:200, 0:200]
    pest_0 = modellino(1,100,100,0)(w,z)
    pest_1 = modellino(1,100,100,0.005)(w,z) #values to change
    pest_2 = modellino(1,100,100,0.01)(w,z) #values to change
    pest_3 = modellino(1,100,100,0.02)(w,z) #values to change
    pest_4 = modellino(1,100,100,0.04)(w,z) #values to change
    pest_5 = modellino(1,100,100,0.06)(w,z) #values to change
    pest_6 = modellino(1,100,100,0.10)(w,z) #values to change
    pest_7 = modellino(1,100,100,0.20)(w,z) #values to change

    
    plt.subplot(241)
    plt.imshow(pest_0, origin = 'lower')
    plt.title('In-focus')

    plt.subplot(242)
    plt.imshow(pest_1, origin = 'lower')
    plt.title('Defocus 0.005')

    plt.subplot(243)
    plt.imshow(pest_2, origin = 'lower')
    plt.title('Defocus 0.01')

    plt.subplot(244)
    plt.imshow(pest_3, origin = 'lower')
    plt.title('Defocus 0.02')

    plt.subplot(245)
    plt.imshow(pest_4, origin = 'lower')
    plt.title('Defocus 0.04')

    plt.subplot(246)
    plt.imshow(pest_5, origin = 'lower')
    plt.title('Defocus 0.06')

    plt.subplot(247)
    plt.imshow(pest_6, origin = 'lower')
    plt.title('Defocus 0.1')

    
    plt.subplot(248)
    plt.imshow(pest_7, origin = 'lower')
    plt.title('Defocus 0.2')

    
    plt.show()



def main():
    test()

if __name__ == "__main__":
    main()

