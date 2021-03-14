import numpy as np
from astropy.modeling.models import Gaussian2D, RickerWavelet2D, Const2D
from photutils.aperture import EllipticalAperture
from Logger import hist
from astropy.logger import log

def bias(image, value, realistic=False):
    '''
    This function creates a matrix with the shape of the CCD with the bias and
    some bad colums
    
    Parameters
    ----------

    image : numpy ndarray  #can be changed with the only shape 
        The image of the CCD, used by the function to extract the shape

    value: int
        It is the bias level to add at the CCD

    realistic : bool, optional
        If True there will be add some bad columns at the CCD

    Output
    ------

    bias_image : numpy ndarray
        The image of the CCD with the bais generated
    '''

    log.info(hist())

    bias_image= np.zeros_like(image)+value

    if realistic:
        shape = image.shape
        number_of_columns = np.random.randint(1,6)  
        columns = np.random.randint(0, shape[1], size=number_of_columns)
        col_pattern = np.random.randint(0, int(0.1*value), size=shape[0]) #add a little pseudo-random noise

        for c in columns:
            bias_image[:, c] = value + col_pattern

    return bias_image

def dark_current(image, current, exposure_time, gain=1.0, hot_pixels=False):
    '''
    This function creates a matrix with the shape of the CCD with the errors introduced
    by the dark current with a poissonian distribution. It also provide for the presence
    of hot pixels in the CCD

    Parameters
    ----------

    image : numpy ndarray  #can be changed with the only shape 
        The image of the CCD, used by the function to extract the shape

    current : float
        This is the value of the dark current for 1 second

    exposure_time : int
        The number of second used for obtainig the CCD's image

    gain : float, optional
        The value of the gain used in the CCD, more the gain, more the errors

    hot_pixel : bool, optional
        This flag allows to choose if there will by the hot pixels or not
        if True there will be hot pixels

    Output
    ------

    dark_bias : numpy ndarray
        The image of the CCD with the dark current noise generated
    '''

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
    '''
    This function provides the errors introduced with the readout of the CCD
    
    Parameters
    ----------

    image : numpy ndarray  #can be changed with the only shape 
        The image of the CCD, used by the function to extract the shape

    amount : float, int
        The central value of the distribution of the errors, it should be
        related with the background (the sqrt)

    gain : float, optional
        The value of the gain used by the CCD

    Output
    ------

    noise : numpy ndarray
        The image of the CCD with the readout noise generated
    '''

    log.info(hist())
    
    shape = image.shape
    noise = np.random.normal(scale=amount/gain, size=shape)

    return noise

def synthetic_ccd(shape):
    '''
    Given a shape this function will creates an empty image with that shape

    Parameters
    ----------

    shape: tuple
        shape[0] : int, the shape of the CCD along the x axis in number of pixel
        shape[1] : int, the shape of the CCD along the y axis in number of pixel

    Output
    ------

    ccd_image : numpy ndarray
        The empy image with the choosen shape
    '''
    
    log.info(hist())

    ccd_image = np.zeros(shape)
    
    return ccd_image

def sky_background(image, sky_count, gain=1):
    '''
    This function creates a image with the shape of the CCD and provides for
    the background light as the moonlight or the city lights.
    Possible improvment: add a gradient
    
    Parameters
    ----------

    image : numpy ndarray  #can be changed with the only shape 
        The image of the CCD, used by the function to extract the shape

    sky_count : int
        The central number of the electrons generated by the background

    gain : float, optional
        The value of the gain used by the CCD

    Output
    ------

    sky_image : numpy ndarray
        The image of the CCD with the background generated
    '''

    sky_image = np.random.poisson(sky_count*gain, size= image.shape) / gain

    return sky_image

def make_cosmic_rays(image, number, strength=10000):
    '''
    It can appens that during an acquisition of an image some pixel are "saturated"
    by the cosmic rays, this function provides for a CCD image with this effect

    Parameters
    ----------

    image : numpy ndarray  #can be changed with the only shape
        The image of the CCD, used by the function to extract the shape

    number : int
        This number is the number of cosmic ray within a single image

    strenght : int, optional
        This is the intensity of a cosmic ray on a pixel

    Output
    ------

    cosmic_image :numpy ndarray
        The image of the CCD with the pixel overflowed by the cosmic rays generated
    '''

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
    '''
    This fuction is pretty eavy, it uses 3 mathematical fuction creates the
    image of a single donut to by added to the flat image

    Parameters
    ----------

    center : numpy array
        center[0] : the position along the x axis of the center of the donut
        center[1] : the position along the y axis of the center of the donut

    diameter : int, float
        Can be interpreted as the diameter in pixel of the donut.
        In reality it is 2sigma of the Gaussian and RickerWavelet (ex MexicaHat)
        functions used to simulate the donut

    amplitude : float, optional
        The peak intensity of the aforementioned functions

    Output
    ------

    Const2D(amplitude=1) + (mh - gauss) : numpy array
        It is a image of the donut created 
    '''
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
    '''
    This function will creates an image with the same shape of the ccd with a specific
    number of dust-donuts drawn on it.

    Parameters
    ----------

    image : numpy ndarray  #can be changed with the only shape 
        The image of the CCD, used by the function to extract the shape

    number : int, optional
        The number of the dust-donut on the CCD

    Output
    ------

    donut_im : numpy ndarray
        The image with the donuts drawn and the shape of the ccd
    '''

    log.info(hist())

    y, x = np.indices(image.shape)

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
    '''
    The sensivity isn't constant, but can vary over the CCD, tipically
    with a gaussian trend and can be worsen by the presence of the dust.
    This function provides for this, creating the flat frame
    

    Parameters
    ----------

    image : numpy ndarray  #can be changed with the only shape nect
        The image of the CCD, used by the function to extract the shape

    vignetting : bool, optional
        If True, the gaussian figure is created on a image with the shape
        of the CCD

    dust : bool, optional
        If True, there will be added the donut generated by the dust on the
        CCD image


    Output
    ------

    sensitivity :numpy ndarray
        The image of the CCD with a big gaussuan curve
        to simulate the variability of the sensivity
    '''
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
    '''
    The value of a pixel can't overcome the maximum value of a 16bit memory (65535) and it is
    called saturated. This function will controll that no pixel would surpass this value and,
    if it does, this function will reset its value to 65535

    Parameter
    ---------

    image : numpy ndarray
        It is the image to controll

    Output
    ------

    image : numpy ndarray
        The image controlled and adjusted
    '''
    
    log.info(hist())
    
    size = image.shape
    
    for i in range(size[0]):
        for j in range(size[1]):
            if image[i,j] >= 65535:
                image[i,j] = 65535
    return image



