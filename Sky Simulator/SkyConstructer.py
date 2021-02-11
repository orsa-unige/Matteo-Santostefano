from astropy.modeling.models import custom_model
import numpy as np
import matplotlib.pyplot as plt
from photutils.datasets import make_noise_image
import astropy.coordinates as coord
from astropy.io import fits
from datetime import datetime
from scipy.special import j1
import Query as Q
from astropy import units as u


#Pocus (sorry for the name) create the model star with parameter set by the user#
@custom_model
def Pocus(x, y, amp=6000.0, x_0=0.0, y_0=0.0, phase_1=0, phase_2=0, phase_3=0, z11=0, z20=0, z40=0, z31=0, z22=0):  
    rho = 5+z20*5 #the first minimum of the Airy disk
    r = (np.sqrt((x-x_0)**2+(y-y_0)**2)+0.001)/rho
    r1 = r*rho
    cosp= (x-x_0)/(r1) #the +0.001 is used to prevent the function from exploding at zero
    sinp= (y-y_0)/(r1)
    t_cosp1= cosp*np.cos(phase_1)+sinp*np.sin(phase_1)#rotation of the system
    t_cosp2= cosp*np.cos(phase_2)+sinp*np.sin(phase_2)
    t_cosp3= cosp*np.cos(phase_3)+sinp*np.sin(phase_3)
    z_t = z11*r*(t_cosp1+1) #tilt
    z_d = z20*(2*(r**2)-1) #defocus Field curvature
    z_s = z40*(r**4-r**2+1/6) #primary spherical aberration
    z_c = z31*(r**3-r*2/3)*(t_cosp2+1) #primary coma aberration
    z_a = z22*(r**2)*(t_cosp3**2) #astigmatic aberration
    R= np.pi*1.2197/rho
    AiryDisk = (2*j1(r1*R)/(r1*R))**2
    Box = np.piecewise(r1, [r1 < rho+5*z20, r1 > rho+5*z20], [1, 0])
    sphere = (-(r1**2)+225)/225

    return (amp*((z_t+z_d+z_s+z_c+z_a)*sphere+AiryDisk)*Box)

#------------------------------------------------------------------------------------#
# It was choosen to use the forth order spherical aberration to simulate the defocus
# it isn't the most clean and matematical correct way to implemented the defocus
# but it is fast and simple to calibrate and the overall error seems neglettable
# Also to create the star it was choosen a pseudo-gaussian because the Airy function
# explode at zero saturing the 'pixel'
#------------------------------------------------------------------------------------#


def Pass(AMP, z11_t, z20, z40, z31, z22): #gives them back multiplied for their Zernikle's coefficient
    tot = z11_t+z40*6*np.sqrt(5)+z31*np.sqrt(8)*3/2+z22*np.sqrt(6)+z20+1
    AMP = AMP/tot                         #normalize the amplitude
    return AMP, z11_t, z20, z40*6*np.sqrt(5), z31*np.sqrt(8)*3/2, z22*np.sqrt(6)

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
        test(0, 0, 0, 0, 0, 0, 0, 0)
        return 0, 0, 0, 0, 0, 0, 0, 0
    elif choice == 1:
        test(1, 1, 0, 1, 1, 0, 120, 240)
        return 1, 1, 1, 1, 1, 0, 120, 240
    else:
        print('Every aberration can have a value between 0 and 1, exept for angles (0,360)')
        tilt = float(input('How much tilt? '))
        ph1 = float(input('The angle of the tilt? '))
        defocus = float(input('How much defocus? '))
        spher = float(input('How much spheric? '))
        coma = float(input('How much coma? '))
        ph2 = float(input('The angle of the coma? '))
        astig = float(input('How much astigmatism? '))
        ph3 = float(input('The angle of the astigmatism? '))
        
        test(tilt, defocus, spher, coma , astig, ph1, ph2, ph3)
        return tilt, defocus, spher, coma , astig, ph1, ph2, ph3

def test(z11, z20, z40, z31, z22, ph1, ph2, ph3): #visualise a test sample star with 
    AMP=6000
    w,z =np.mgrid[0:200, 0:200]
    AMP, z11, z20, z40, z31, z22 = Pass(AMP, z11, z20, z40, z31, z22)
    test = Pocus(AMP, 100, 100, ph1, ph2, ph3, z11, z20, z40, z31, z22)(w,z)
    plt.imshow(test, origin='lower')
    plt.show()
    

def datum(coordinates,z11, z20, z40, z31, z22, ph1, ph2, ph3, YES= ['Y', 'YES', 'S', 'SI', 'Sì']): #create and at will dispay the sky
    print("Initialing the CCD\n")
    Dx, Dy = 1024, 1032
    x, y = np.mgrid[0:Dx, 0:Dy]
    data = Pocus(0) (x, y) #inizialize a datum to 0
    AMP=7000
    AMP, z11, z20, z40, z31, z22 = Pass(AMP, z11, z20, z40, z31, z22)
    ph = np.random.randint(0,360) #the initial fase is randomly choosen

    p_x,p_y,f_u,f_b,f_v,f_r,f_i,f_g,f_j,f_h,f_k = Q.query(coordinates)
    #TODO, pass the total flux weithed for the efficency spectrum

    tot = len(p_x)
    print(tot)
    for i in range(tot):
        f_tot=f_u[i]+f_b[i]+f_v[i]+f_r[i]+f_i[i]+f_g[i]+f_j[i]+f_h[i]+f_k[i]
        f_tot=f_tot*10  #TODO multiply for the gain of the CCD
        data += Pocus(f_tot,p_y[i],p_x[i], ph1+ph, ph2+ph, ph3+ph, z11, z20, z40, z31, z22)(x,y)
        print(f"{(i+1)*100/(tot)}% Sky Completed")

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
data_cache = {}

def main():
    print('Build your Star!\n')
    YES = ['Y', 'YES', 'S', 'SI', 'Sì']
    a = False
    while not a: #this cycle reapet itself until the sample model is good for the user
        tilt, defocus, sphere, coma, astig, ph1, ph2, ph3= building(YES)
        b = input('Do you like your model? (Y/N)\n')
        if b.upper() in YES:
            a = True
        else:
            a = False
    c = input('Sky Coordinates: ')
    coordinates = coord.SkyCoord(c, frame='icrs', unit=(u.hourangle, u.deg))
    print('Bulding your sky, it can take few seconds\n')
    data = datum(coordinates, tilt, defocus, sphere, coma, astig, ph1, ph2, ph3, YES)
    c = input('Do you want to save it? (Y/N)\n')
    if c.upper() in YES:
        saver(data)


    
if __name__ == "__main__":
    main()


