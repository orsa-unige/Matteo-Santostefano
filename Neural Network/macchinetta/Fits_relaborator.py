from keras.models import load_model
from numpy import load
from numpy import expand_dims
from numpy import zeros
from astropy.io import fits
from keras import backend as K

# load an image
def load_image(pixels, size=(256,256)):
	bit8 = 127.5
	bit16 = 32767.5
	# scale from [0,255] to [-1,1]
	pixels = expand_dims(pixels, 2)  #sembra funzionarw, poi però controlla
	pixels = (pixels - bit16) / bit16
	# reshape to 1 sample
	X_1 = zeros((1,pixels.shape[0], pixels.shape[1],3))
	X_1[:,:,:] = pixels
	return X_1

def save_fits(data):
    hdul = fits.PrimaryHDU(data) #save the file
    hdr = hdul.header
    hdul.scale('int16', bzero=32768)
    hdul.writeto('try.fits', overwrite = True)

def load_fits(filename):
        hdul = fits.open(filename)
        data = hdul[0].data
        return data

def image_slicer(image):
        image_array = []
        size = image.shape
        number_x = size[0]//256
        number_y = size[1]//256
        for i in range(number_x):
                for j in range(number_y):
                        x = image[(256*i):(256*i)+256,(256*j):(256*j)+256]
                        image_array.append(x)
        return image_array, number_x, number_y

def image_recomposer(image_array, number_x, number_y):
        new_image = zeros((number_x*256, number_y*256))
        a = 0
        for ctrl_x in range(number_x):
                for ctrl_y in range(number_y):
                        for i in range(256-1):
                                for j in range(256-1):
                                        new_image[i+256*ctrl_x][j+256*ctrl_y] = image_array[a][i][j]
                        a += 1
        return new_image              

def main():
        # load source image
        image_array = []
        src_image_array = []
        gen_image_array = []
        data_array= []

        filename = 'try_1.fits'
        image = load_fits(filename)
        image_array, num_x, num_y = image_slicer(image)
        for image in image_array:
                src_image_array.append(load_image(image))

        model = load_model('model_042800.h5', compile=False)

        # generate image from source
        for src_image in src_image_array:
                gen_image_array.append(model.predict(src_image))
        K.clear_session()
        # scale from [-1,1] to [0,1]
        bit16 = 32767.5
        for gen_image in gen_image_array:
                gen_image = ((gen_image + 1) * bit16) # from[-1,1] to [0,65535]
                data = zeros((256,256))
                for i in range (gen_image.shape[3]):
                        data += gen_image[0,:,:,i]
                data_array.append(data)
        new_image = image_recomposer(data_array, num_x, num_y)
        save_fits(new_image)


if __name__ == '__main__':
        main()
