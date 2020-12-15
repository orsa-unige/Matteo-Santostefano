#genera coppie di immagini semplici che rappresentino il cielo,
#una reale ed una sfocata per testare CNN.
#L'allenamento della CNN sarà poi da effettuare con immagini reali,
#queste servono solo per testare l'architettura.

#-genera 2 immagini:
#   -sfondo nero
#   -una con puntini distribuiti random
#   -una con cerchi ed ellissi centrati negli stessi punti


import numpy as np 
from PIL import Image, ImageDraw, ImageFilter
from pathlib import Path
import os


def Inizialization(dims): #crea immagine in scala di grigi a sfondo nero
    #TODO: togli sfondo completamente nero e randomizza una scala stretta di grigi scuri
    im = Image.new(mode = "L", size = dims, color='black') #L is grayscale 
    return im

def Ellipsator(dims): #crea le ellissi
    r=np.random.randint(2, (min(*dims)/25)) 
    e=np.random.rand()*2
    i=np.random.choice([True, False])
    if i:
        q=r+e
        w=r-e
    else:
        q=r-e
        w=r+e
    return r,q,w


def drawSky(dims): #naive
    img_1 = Image.new(mode = "RGB", size = dims, color='black') #stelle vere
    img_2 = Image.new(mode = "RGB", size = dims, color='black') #stelle deformate
    n = np.random.randint(0, 50) #numero di stelle
    draw_1 = ImageDraw.Draw(img_1)
    draw_2 = ImageDraw.Draw(img_2)
    for i in range(n):
        r,q,w=Ellipsator(dims)
        fill= (np.random.randint(100,255))
        fill=(fill,fill,fill)
        x, y =np.random.randint(0, dims[0]-max(r,q,w)), np.random.randint(0, dims[1]-max(r,q,w))
        draw_1.ellipse((x-r, y-r, x+r, y+r), fill) #stelle buone
        draw_2.ellipse((x-(q),y-(w),x+(q),y+(w)), fill = None , outline = fill) #ellisse base
    rad=np.random.randint(0,2)
    img_2 = img_2.filter(ImageFilter.GaussianBlur(radius = 2+ rad))
    img_1 = img_1.filter(ImageFilter.GaussianBlur(radius = 2))    
    return img_1, img_2

def PathCreator(): # crea la directory dove salvare le immagini
    path = f'{Path.cwd()}\Sky'
    path_1 = f'{path}\Test'
    path_2 = f'{path}\Train'
    path_3 = f'{path}\Val'
    Path(path).mkdir()
    Path(path_1).mkdir()
    Path(path_2).mkdir()
    Path(path_3).mkdir()
    return path_1, path_2, path_3

#def ImageSaver(im_1, im_2, path_1, path_2, i):  #salva le immagini
#    name_1 = f'SkyC{i}.jpg'
#    name_2 = f'SkyB{i}.jpg'
#    file_path_1 = os.path.join(path_1,name_1)
#    file_path_2 = os.path.join(path_2, name_2)
#    im_1.save(file_path_1, 'JPEG')
#    im_2.save(file_path_2, 'JPEG')

def accSaver(im, Path, i):
    name = f"{i+1}.jpg"
    new_path = os.path.join(Path,name)
    im.save(new_path, "JPEG")

def main():
    dims = (256,256)
    path = PathCreator()  #da spegnere se si vogliono solo visualizzare le immagini
    for j in path:
        for i in range(400):
            image_1, image_2= drawSky(dims)
            #image_1.show()  #da accendere per visualizzare le coppie di immagini
            #image_2.show()  #da accendere per visualizzare le coppie di immagini
            #ImageSaver(image_1, image_2, path_1, path_2, i) #da spegnere se si vogliono solo visualizzare le immagini
            new_im = accopiator(image_1, image_2)
            accSaver(new_im, j , i)

def accopiator(image_1, image_2): #attacca le 2 immagini
    new_im = Image.new('RGB', (512,256)) #creates a new empty image, RGB mode
    new_im.paste(image_1)
    new_im.paste(image_2, (256,0))
    return new_im


if __name__ == "__main__":
    main()


