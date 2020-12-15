import os
import tensorflow as tf
from pathlib import Path

def pathfinder():
    URL = Path.cwd()

    #path_to_zip = tf.keras.utils.get_file('Sky.zip',
                                      #origin=URL,
                                      #extract=True)

    PATH = os.path.join(URL, 'Sky\\')
    return URL, PATH