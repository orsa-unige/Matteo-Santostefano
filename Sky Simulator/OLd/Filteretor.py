from PIL import Image
from pathlib import Path
import os
from scipy import signal
import matplotlib.pyplot as plt
import Kernels
import numpy as np

def main():
    path = Path.cwd()
    path = os.path.join(path, 'test.png')
    new_im = Image.open(path).convert('L')
    dimker=40
    kernel = Kernels.k2(dimker)  #si può passare(dimensione, anisotropia) / (dim, std, anis) nel caso di k2
    arry = signal.fftconvolve(new_im, kernel, mode='same')

    fig, (normal, ker, convol)= plt.subplots(3,1)
    normal.imshow(new_im, cmap='gray')
    normal.set_title('Original')
    normal.set_axis_off()
    ker.imshow(kernel, cmap='gray')
    ker.set_title('Kernel')
    ker.set_axis_off()
    convol.imshow(arry, cmap='gray')
    convol.set_title('Convoluted')
    convol.set_axis_off()
    fig.show()
    return arry
   #plt.savefig('my.png')

if __name__ == '__main__':
    main()

