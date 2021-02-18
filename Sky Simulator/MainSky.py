import Query as Q
import SkyConstructer as SC
from numpy import random


def main():

    auto =['A', 'AUTO', 'AUTOMATIC']
    c = input('Automatic or Manual? (A/M) ')
    if c.upper() in auto:
        Automatics()
    else:
        SC.main()

def Automatics():

    coordinates = input('Sky Coordinates: ')

    SC.datum(coordinates, 6)

if __name__ == "__main__":
    main()
