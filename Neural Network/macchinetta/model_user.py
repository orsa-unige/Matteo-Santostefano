# example of loading a pix2pix model and using it for one-off image translation
from keras.models import load_model
from keras.preprocessing.image import img_to_array
from keras.preprocessing.image import load_img
from numpy import load
from numpy import expand_dims
from numpy import zeros
from numpy import rot90, fliplr
from matplotlib import pyplot
from astropy.io import fits
from keras import backend as K

# load an image
def load_image(filename, size=(256,256)):
	# load image with the preferred size
	pixels = load_img((filename), color_mode = 'grayscale',target_size=size)
	# convert to numpy array
	pixels = img_to_array(pixels)
	bit8 = 127.5
	bit16 = 32767.5
	# scale from [0,255] to [-1,1]
	pixels = (pixels - bit16) / bit16
	# reshape to 1 sample
	X_1 = zeros((1,pixels.shape[0], pixels.shape[1],3))
	X_1[:,:,:] = pixels
	pixels = expand_dims(pixels, 0)
	return X_1

def save_fits(data):
    hdul = fits.PrimaryHDU(data) #save the file
    hdr = hdul.header

    hdul.scale('int16', bzero=32768)
    hdul.writeto('try.fits', overwrite = True)

def main():
        # load source image
        src_image = load_image('try.png')
        print('Loaded', src_image.shape)
        # load model
        model = load_model('model_005600.h5',compile=False)
        # generate image from source
        gen_image = model.predict(src_image)
        K.clear_session()
        # scale from [-1,1] to [0,1]
        bit16 = 32767.5
        gen_image = ((gen_image + 1) / 2.0)
        gen_image_1 = gen_image * 2 * bit16
        data = zeros((256,256))
        dato = zeros((256,256))
        for i in range (gen_image.shape[3]):
                data += gen_image_1[0,:,:,i]

        dato = rot90(data, 2)
        dato = fliplr(dato)
         # plot the image
        save_fits(dato)
        #pyplot.imshow(gen_image[0])
        #pyplot.axis('off')
        #pyplot.show()

if __name__ == '__main__':
        main()
