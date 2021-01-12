import Filteretor as Fil
from PIL import Image

new_im = Fil.main()
new_im = Image.fromarray((new_im/255))
new_im.show()

#print(new_im)
