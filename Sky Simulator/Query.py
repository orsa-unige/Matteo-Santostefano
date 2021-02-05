from astroquery.simbad import Simbad
import astropy.coordinates as coord
from astropy import wcs
from astropy import units as u
import numpy as np


def Coordinator(c, center):
    
    w = wcs.WCS(naxis=2)
    w.wcs.crpix = [1, 1]
    w.wcs.crval = [center[0], center[1]]
    w.wcs.ctype = ["RA", "DEC"]
    
    x,y = wcs.utils.skycoord_to_pixel(c, w)
    x = (-x*4389) + 512 #(deg*arc/sec*deg)*arc/sec*pixel+offset
    y = (-y*4389) + 516

    return x, y


def Data_structure():
    Coord_x = []
    Coord_y = []
    Flux_U = []
    Flux_B = []
    Flux_V = []
    Flux_R = []
    Flux_I = []
    Flux_G = []
    Flux_J = []
    Flux_H = []
    Flux_K = []
    data = [Coord_x, Coord_y, Flux_U, Flux_B, Flux_V, Flux_R, Flux_I,Flux_G,Flux_J,Flux_H,Flux_K]
    return data
    
def query(coordi = coord.SkyCoord("07 59 5.840", "+15 23 29.24", frame = 'icrs', unit=(u.hourangle, u.deg)),
          radius = (7/60)*u.deg):
    
    photo_filters = ['U', 'B', 'V', 'R', 'I','G','J','H','K']
    f =['flux({0})'.format(x) for x in photo_filters]
    [Simbad.add_votable_fields(g) for g in f]
    result_table = Simbad.query_region(coordi, radius)
    center = [0,0]
    center[0] = coordi.ra.degree
    center[1] = coordi.dec.degree

    data=Data_structure()

    Stars_number=len(result_table)

    for i in range (Stars_number):
        row = result_table[i]
        
        tot = 0
        for i in range(11,20):  #sum the fluxies
            if row[i] is np.ma.masked:
                row[i] = 0
            tot += row[i]

        if tot != 0: #if total flux is 0 the object is not passed
            ra = hms(row[1])
            de = dms(row[2])
            c = coord.SkyCoord(ra, de, frame = 'icrs')
            x, y = Coordinator(c, center)
            data[0].append(x)
            data[1].append(y)
            #the fluxies of different wavelengh 
            for f in range (11,20):
                if row[f] is np.ma.masked:
                    data[f-9].append(0)
                else:
                    data[f-9].append(row[f])
            #TODO multiply the fluxies for the efficiency spectrum
    return data

def hms(row):
    s = list(row)
    s[2] = 'h'
    for n,i in enumerate(s):
        if i == ' ':
            s[n]='m'

    s.append('s')
    if s.count('m') == 1:
        f = ''.join(s)
    else:
        s.insert(3,'m')
        s.insert(3,'0')
        s.insert(3,'0')
        f = ''.join(s)
    return f

def dms(row):
    s= list(row)
    if s[2] ==' ':
        s[2] = 'd'
    else:
        s[3] = 'd'
    for n,i in enumerate(s):
        if i == ' ':
            s[n]='m'
    s.append('s')
    if s.count('m') == 1:
        f = ''.join(s)
    else:
        s.insert(4,'m')
        s.insert(4,'0')
        s.insert(4,'0')
        f = ''.join(s)
    return f

def main():
    data=query()
    print(data)

if __name__ == '__main__':
    main()
