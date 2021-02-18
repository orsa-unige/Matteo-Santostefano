from astropy.modeling.models import custom_model
import numpy as np
import matplotlib.pyplot as plt
import astropy.coordinates as coord
from astropy.io import fits
from datetime import datetime
from scipy.special import j1
import Query as Q
from astropy import units as u
from scipy import fft, signal

@custom_model
def modellino(x, y, amp=5.0, x_0=0.0, y_0=0.0, defocus_distance = 6):
    
    wavelength = 0.005
    KX = x_0-x
    KY = y_0-y
    r = np.sqrt(KX**2+KY**2)
    pupil = np.piecewise(r, [r < 50, r > 50, r < 2], [1, 0, 0])

    phaseAngle = 1j*defocus_distance*np.sqrt((2*np.pi/wavelength)**2-KX**2-KY**2+0j) #unnecessary 0j but keeping it for complex reasons
    kernel = np.exp(phaseAngle)
    defocusPupil = pupil * kernel
    defocusPSFA = fft.fftshift(fft.fft2(defocusPupil))

    return amp*np.abs(defocusPSFA)
  
@custom_model
def simple_point(x,y, amp=5.0, x_0=0, y_0=0):
    KX = x_0-x
    KY = y_0-y
    r=np.sqrt(KX**2+KY**2)
    point = np.piecewise(r, [r < 0.5, r > 0.5], [1,0])

    return point*amp

#the initial fase of the program it allows for different choice    
def intro(YES):
    a = input('Do you want a In-focus star? (Y/N)\n')
    if a.upper() in YES:
        choice = 0
    else:
        b = input('Do you want a sample star? (Y/N)\n ')
        if b.upper() in YES:
            choice = 1
        else:
            choice = 2
    return choice


#dependig of the choice in Intro, building creates a test star         
def building(YES): 
    choice = intro(YES)
    if choice == 0:
        test(0)
        return 0
    elif choice == 1:
        test(6)
        return 6
    else:
        defocus = float(input(print('How much defocus? between(0,10)\n')))
        
        test(defocus)
        return defocus

def test(defocus): #visualise a test sample star with 
    AMP=6000
    w,z =np.mgrid[0:200, 0:200]
    test = modellino(AMP, 100, 100, defocus)(w,z)
    plt.imshow(test, origin='lower')
    plt.show()
    
def filter_choice():

    photo_filters = ['U', 'B', 'V', 'R', 'I', 'G', 'J', 'H','K']
    list_filter = [item for item in input('choose the filters: ').split()]
    for filter in list_filter:
        if filter.upper() == 'ALL':
            return photo_filters
        elif filter.upper() in photo_filters:
           flag = True
        else:
           flag = False

    if flag:
        return list_filter
    else:
        return 'R'

def datum(coordinates, defocus=6, YES= ['Y', 'YES', 'S', 'SI', 'Sì']): #create and at will dispay the sky
    
    filters = filter_choice()
    print("Initialing the CCD\n")
    Dx, Dy = 1024, 1032
    x, y = np.mgrid[0:Dx, 0:Dy]
    data = modellino(0,0,0) (x, y) #inizialize a datum to 0
    AMP= 700
    ph = np.random.randint(0,360) #the initial fase is randomly choosen

    p_x,p_y,flux = Q.query(coordinates, filters)
    kernel = modellino(1,512,518,defocus)(x,y)
    tot = len(p_x)
    print(tot)
    for i in range(tot):
        data += simple_point(flux[i]*AMP, p_x[i], p_y[i])(x,y)
        print(f"{(i+1)*100//(tot)}% Sky Completed")

    data = signal.fftconvolve(data, kernel, mode= 'same')
    
    choice = input('Do you want do visualize? (Y/N)\n')

    if choice.upper() in YES:
        plt.imshow(data, origin='lower')
        plt.show()

    choice = input('Do you want to save it? (Y/N)\n')
    if choice.upper() in YES:
        saver(data)
    return data

def saver(data): #Save the simulated sky in a fits
    hdul = fits.PrimaryHDU(data)
    hdul.scale('int16', bzero=32768)
    hdul.writeto(f'{datetime.now().strftime("%Y%m%d-%H%M%S")}.fits')


def main():
    print('Build your Star!\n')
    YES = ['Y', 'YES', 'S', 'SI', 'Sì']
    a = False
    while not a: #this cycle reapet itself until the sample model is good for the user
        defocus = building(YES)
        b = input('Do you like your model? (Y/N)\n')
        if b.upper() in YES:
            a = True
        else:
            a = False
    coordinates = input('Sky Coordinates: ')

    print('Bulding your sky, it can take few seconds\n')
    
    data = datum(coordinates, defocus, YES)
  
if __name__ == "__main__":
    main()


