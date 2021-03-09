import sys
from astropy.time import Time

def timestamp(date=False, iso=False):
    '''
    Get a time stamp for strings in iso format.
    If date=True, return only yyyy-mm-dd.
    '''

    time = Time.now()
    string = time.iso[:-4] if not iso else time.isot[:-4]
    if date:
        string = string.split()[0]

    return string

def hist(text=None):
    '''
    Add a history tag of type. Example:
    [ mainsky.fill_header.newhead ] 2020-04-29T20:06:47 > Created.
    '''
    tsp = timestamp()

    text = 'Called.' if not text else text
    routine = sys._getframe().f_back.f_code.co_name
    string = f"{tsp} [ {routine} ] > {text}"

    return string
