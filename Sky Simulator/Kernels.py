import numpy as np
from scipy import signal
from pathlib import Path
from PIL import Image

def star(klen): #Easter Egg
    ker = Image.open('star.png').convert('L')
    ker = ker.resize((klen,klen))
    return ker

def k1(keln): #cross, just to test manual kernel
    kernel=np.array([
        1,0,0,0,0,0,0,0,0,1,
        0,1,0,0,0,0,0,0,1,0,
        0,0,1,0,0,0,0,1,0,0,
        0,0,0,1,0,0,1,0,0,0,
        0,0,0,0,.5,.5,0,0,0,0,
        0,0,0,0,.5,.5,0,0,0,0,
        0,0,0,1,0,0,1,0,0,0,
        0,0,1,0,0,0,0,1,0,0,
        0,1,0,0,0,0,0,0,1,0,
        1,0,0,0,0,0,0,0,0,1,
        ])
    kernel = kernel.reshape(10,10)
    return kernel

def k2(klen=10, std=3, anis=1): #Naive gaussian kernel, test 
    gker1d= signal.gaussian(klen, std=std).reshape(klen).reshape(klen,1)
    gker1d_b= signal.gaussian(klen, std=std*anis).reshape(klen).reshape(klen,1)
    gker2d=np.outer(gker1d, gker1d_b)
    return gker2d

def k3(klen, anis=1): #Naive defocus, test from image
    ker = Image.open('defocus1.png').convert('L')
    ker = ker.resize((klen,klen*anis))
    return ker

def k4(klen, anis=1): #Naive astigmatism, test from image
    ker = Image.open('astigmatims.png').convert('L')
    ker = ker.resize((klen,klen*anis))
    return ker

def k5(klen, anis=1): #Naive coma, test from image
    ker = Image.open('coma.jpg').convert('L')
    ker = ker.resize((klen,klen*anis))
    return ker

def k6(klen,anis=1): #Naive spherical, test from image
    ker = Image.open('spherical.jpg').convert('L')
    ker= ker.resize((klen,klen*anis))
    return ker

