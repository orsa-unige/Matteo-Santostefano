import Query as Q
import SkyConstructer as SC
from numpy import random
import astropy.coordinates as coord
from astropy import units as u

def main():

    auto =['A', 'AUTO', 'AUTOMATIC']
    c = input('Automatic or Manual? (A/M) ')
    if c.upper() in auto:
        Automatics()
    else:
        SC.main()

def Automatics():

    c = input('Sky Coordinates: ')
    coordinates = coord.SkyCoord(c, frame='icrs', unit=(u.hourangle, u.deg))

    a = random.rand(5)
    ph = random.randint(360, size = 3)

    SC.datum(coordinates, a[0], a[1], a[2], a[3], a[4], ph[0], ph[1], ph[2])

if __name__ == "__main__":
    main()
