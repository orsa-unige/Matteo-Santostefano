import numpy as np
from astropy.modeling.models import Gaussian2D, RickerWavelet2D, Const2D
from photutils.aperture import EllipticalAperture
from Logger import hist
from astropy.logger import log

def bias(image, value, realistic=False):

    log.info(hist())

    bias_image= np.zeros_like(image)+value

    if realistic:
        shape = image.shape
        number_of_columns = 5  #to randomize
        columns = np.random.randint(0, shape[1], size=number_of_columns)
        col_pattern = np.random.randint(0, int(0.1*value), size=shape[0]) #add a little pseudo-random noise

        for c in columns:
            bias_image[:, c] = value + col_pattern

    return bias_image

def dark_current(image, current, exposure_time, gain=1.0, hot_pixels=False):

    log.info(hist())

    base_current = current*exposure_time/gain
    dark_bias = np.random.poisson(base_current, size = image.shape)

    if hot_pixels:
        '''set the probability of 0.01% of a pixel to be hot'''
        y_max, x_max = dark_bias.shape

        numb_hot = int(0.0001*x_max*y_max)
        hot_x = np.random.randint(0, x_max, size=numb_hot)
        hot_y = np.random.randint(0, y_max, size=numb_hot)

        hot_current = 10000*current
        dark_bias[[hot_y, hot_x]] = hot_current *exposure_time#/gain

    return dark_bias
    
def read_out_noise(image, amount, gain=1.0):

    log.info(hist())
    
    shape = image.shape
    noise = np.random.normal(scale=amount/gain, size=shape)

    return noise

def synthetic_ccd(shape):
    
    log.info(hist())

    ccd_image = np.zeros(shape)
    
    return ccd_image

def sky_background(image, sky_count, gain=1):

    sky_image = np.random.poisson(sky_count*gain, size= image.shape) / gain

    return sky_image

def make_cosmic_rays(image, number, strength=10000):

    log.info(hist())

    cosmic_image = np.zeros_like(image)

    # Yes, the order below is correct. The x axis is the column, which
    # is the second index.
    max_y, max_x = cosmic_image.shape

    # Get the smallest dimension to ensure the cosmic rays are within the image
    maximum_pos = np.min(cosmic_image.shape)
    # These will be center points of the cosmic rays, which we place away from
    # the edges to ensure they are visible.
    xy_cosmic = np.random.randint(0.1 * maximum_pos, 0.9 * maximum_pos,
                          size=[number, 2])

    cosmic_length = 5  # pixels, a little big
    cosmic_width = 2
    theta_cosmic = 2 * np.pi * np.random.rand()
    apertures = EllipticalAperture(xy_cosmic, cosmic_length, cosmic_width, theta_cosmic)
    masks = apertures.to_mask(method='center')
    for mask in masks:
        cosmic_image += strength * mask.to_image(shape=cosmic_image.shape)

    return cosmic_image

def make_one_donut(center, diameter=10, amplitude=0.25):
    '''dust simulation'''
    log.info(hist())
    
    sigma = diameter / 2
    mh = RickerWavelet2D(amplitude=amplitude,
                      x_0=center[0], y_0=center[1],
                      sigma=sigma)
    gauss = Gaussian2D(amplitude=amplitude,
                       x_mean=center[0], y_mean=center[1],
                       x_stddev=sigma, y_stddev=sigma)
    return Const2D(amplitude=1) + (mh - gauss)


def add_donuts(image, number=20):

    log.info(hist())

    y, x = np.indices(image.shape)

    rng = np.random.RandomState(43901) #TONOT
    shape = np.array(image.shape)
    border_padding = 50

    # dust specks range from 1% to 5% of the image size

    min_diam = int(0.02 * shape.max())
    max_diam = int(0.05 * shape.max())

    # Weight towards the smaller donuts because it looks more like real flats..
    diameters = np.random.choice([min_diam, min_diam, min_diam, max_diam],
                           size=number)

    # Add a little variation in amplitude
    amplitudes = np.random.normal(0.25, 0.05, size=number)
    center_x = np.random.randint(border_padding,
                           high=shape[1] - border_padding, size=number)
    center_y = np.random.randint(border_padding,
                           high=shape[0] - border_padding, size=number)
    centers = [[x, y] for x, y in zip(center_x, center_y)]

    donut_model = make_one_donut(centers[0], diameter=diameters[0],
                                 amplitude=amplitudes[0])
    donut_im = donut_model(x, y)
    idx = 1
    for center, diam, amplitude in zip(centers[1:],
                                       diameters[1:],
                                       amplitudes[1:]):
        idx += 1
        donut_model = make_one_donut(center, diameter=diam,
                                      amplitude=amplitude)
        donut_im += donut_model(x, y)

    donut_im /= number

    return donut_im


def sensitivity_variations(image, vignetting=True, dust=True):
    '''the sensivity isn't constant, this function provides for that'''
    log.info(hist())
    
    sensitivity = np.zeros_like(image) + 1.0
    shape = np.array(sensitivity.shape)

    if dust or vignetting:
        # I don't know why, but y,x not x,y
        y, x = np.indices(sensitivity.shape)

    if vignetting: #TODO, centro gaussiana da spostare
        # Generate very wide gaussian centered on the center of the image,
        # multiply the sensitivity by it.
        vign_model = Gaussian2D(amplitude=1,
                                x_mean=shape[0] / 2, y_mean=shape[1] / 2,
                                x_stddev=2 * shape.max(),
                                y_stddev=2 * shape.max())
        vign_im = vign_model(x, y)
        sensitivity *= vign_im

    if dust:
        dust_im = add_donuts(image, number=20)
        dust_im = dust_im / dust_im.max()
        sensitivity *= dust_im

    return sensitivity

def saturation_controll(image):
    log.info(hist())
    size = image.shape
    for i in range(size[0]):
        for j in range(size[1]):
            if image[i,j] >= 65535:
                image[i,j] = 65535
    return image


#image = synthetic_ccd([2000,2000])

#gain = 1.0
#exposure = 30.0
#dark = 0.1
#sky_counts = 20
#bias_level = 1100
#read_noise_electrons = 5

#flat = sensitivity_variations(image)
#bias_only = bias(image, bias_level, realistic=True)
#noise_only = read_out_noise(image, read_noise_electrons, gain=gain)
#dark_only = dark_current(image, dark, exposure, gain=gain)
#sky_only = sky_background(image, sky_counts, gain=gain)

#final_image = bias_only + noise_only + dark_only + flat * (sky_only)

#data = final_image
#hdul = fits.PrimaryHDU(data)
#hdul.writeto(f'sample.fits', overwrite = True)
#plt.imshow(dark_only, cmap='gray')
#plt.show()
