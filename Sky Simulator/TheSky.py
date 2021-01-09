from astropy.modeling.models import custom_model
import numpy as np
import matplotlib.pyplot as plt
from photutils.datasets import make_noise_image
from astropy.io import fits
from datetime import datetime

#Pocus (sorry for the name) create the model star with parameter set by the user°
@custom_model
def Pocus(x, y, amp=5.0, x_0=0.0, y_0=0.0, phase_1=0, phase_2=0, d=1, c=1, b=1, a=1, e=1):  
    r = np.sqrt((x-x_0)**2+(y-y_0)**2)
    cosp= (x-x_0)/(r+0.001) #the +0.001 is used to prevent the function from exploding at zero
    sinp= (y-y_0)/(r+0.001)
    t_cosp1= cosp*np.cos(phase_1)+sinp*np.sin(phase_1)#rotation of the system
    t_cosp2= cosp*np.cos(phase_2)+sinp*np.sin(phase_2)
    z_s = d*amp*(r**4-r**2) #spherical aberration
    z_c = c*amp*(r**3-r*2/3)*t_cosp1 #coma aberration
    z_a = a*amp*(t_cosp2**2-0.5)*r**2 #astigmatic aberration
    z_g = b*amp*np.exp(-r**2/(4+0.001)) #speudo-gaussian
    return (z_a+z_s+z_c+z_g)*z_g
#------------------------------------------------------------------------------------#
# It was choosen to use the forth order spherical aberration to simulate the defocus
# it isn't the most clean and matematical correct way to implemented the defocus
# but it is fast and simple to calibrate and the overall error seems neglettable
# Also to create the star it was choosen a pseudo-gaussian because the Airy function
# explode at zero saturing the 'pixel
#------------------------------------------------------------------------------------#

#the initial fase of the program it allows for different choice
def Intro(YES):
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
    choice = Intro(YES)
    if choice == 0:
        test(0, 0, 1, 0)
        return 0, 0, 1, 0, 0
    elif choice == 1:
        test(1, 1, 1, 1)
        return 1, 1, 1, 1, 0
    else:
        defocus = float(input('How much defocus? '))
        coma = float(input('How much coma? '))
        bessel = float(input('How much bessel? '))
        astig = float(input('How much astigmatism? '))
        phase = float(input('The phase between coma and stigmatism? (deg) '))
        test(defocus, coma, bessel, astig , phase)
        return defocus, coma, bessel, astig, phase

def test(d,c,b,a, ph=0): #visualise a test sample star with 
    w,z =np.mgrid[0:50, 0:50]
    test = Pocus(20, 25, 25, 0, ph, d, c, b, a) (w,z)
    plt.imshow(test, origin='lower')
    plt.show()

def datum(defo, com, bess, ast, ph, YES): #create and at will dispay the sky
    x, y = np.mgrid[0:500, 0:500] 

    data = Pocus(0) (x, y) #inizialize a datum to 0

    phase1 = np.random.randint(0,360) #the initial fase is randomly choosen

    #the amplitude and the position are randomy choosen
    #the choice from catalogs can be implemented
    n= np.random.randint(10,100)
    for N in range(n):
        amp= np.random.randint(1, 50)
        X= np.random.randint(0, 500)
        Y= np.random.randint(0, 500)
        data += Pocus(amp,X,Y, phase1, phase1+ph, defo, com, bess, ast)(x,y)

    data +=  make_noise_image(data.shape, distribution='gaussian', mean=10.,
                          stddev=5., seed=123) #add background noise
    c = input('Do you want do visualize? (Y/N)\n')
    if c.upper() in YES:
        plt.imshow(data, origin='lower')
        plt.show()
    return data

def saver(data): #Save the simulated sky in a fits
    hdul = fits.PrimaryHDU(data)
    hdul.writeto(f'{datetime.now().strftime("%Y%m%d-%H%M%S")}.fits')

def main():
    print('Build your Star!\n')
    YES = ['Y', 'YES', 'S', 'SI', 'Sì']
    a = False
    while not a: #this cycle reapet itself until the sample model is good for the user
        defocus, coma, bessel, astig, phase = building(YES)
        b = input('Do you like your model? (Y/N)\n')
        if b.upper() in YES:
            a = True
        else:
            a = False
    print('Bulding your sky, it can take few seconds\n')
    data = datum(defocus, coma, bessel, astig, phase, YES)
    c = input('Do you want to save it? (Y/N)\n')
    if c.upper() in YES:
        saver(data)

    
if __name__ == "__main__":
    main()


