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
    x = (x*4389)/2 + 512 #(deg*arc/sec*deg)*arc/sec*pixel+offset
    y = (y*4389)/2 + 516

    return x, y


def Data_structure():

    Coord_x = []
    Coord_y = []
    Flux_tot = []
    data = [Coord_x, Coord_y, Flux_tot]

    return data
    
def query(coordi = "07 59 5.840  +15 23 29.24", photo_filters = ['U', 'B', 'V', 'R', 'I', 'G', 'J', 'H','K'], radius = (15.6/60)*u.deg):
    
    coordi = coord.SkyCoord(coordi, frame = 'icrs', unit=(u.hourangle, u.deg))
    f =['flux({0})'.format(x) for x in photo_filters]
    [Simbad.add_votable_fields(g) for g in f]
    result_table = Simbad.query_region(coordi, radius)
    
    center = [0,0]
    center[0] = coordi.ra.degree
    center[1] = coordi.dec.degree

    data = Data_structure()

    Stars_number = len(result_table)
    len_filter_list = len(photo_filters)
    
    for i in range (Stars_number):
        row = result_table[i]
        
        tot = 0
        for i in range(11, 11+len_filter_list):  #sum the fluxies
            if row[i] is np.ma.masked:
                row[i] = 0
            tot += row[i]
        data[2].append(tot)

        if tot != 0: #if total flux is 0 the object is not passed
            c = coord.SkyCoord(row[1], row[2], frame = 'icrs', unit=(u.hourangle, u.deg))##change here

            x, y = Coordinator(c, center)
            data[0].append(x)
            data[1].append(y)


    return data



def main():
    data = query()
    print(data)

if __name__ == '__main__':
    main()
