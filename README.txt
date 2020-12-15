1-Utilizzare Image_generator.py per generare il dataset da dare in pasto alla rete neurale (RN)
2-Image_elaborator è il costruttore della RN, la crea e la allena
3-Image_elaboretor lavora con modulo1 e modulo2

Per poter utilizzare Image_elaborator sarà necessario installare tensorflow (PIP install tensorflow), 
il pacchetto che gestisce le RN

# Rapido approfondimento sulla rete neurale:
#
# L'architettura di base della RN è basata su Pix2Pix, 
# una rete affinata per tradurre un immagine in un altra.
# https://arxiv.org/abs/1611.07004
#
# Questa utilizza le "conditional's GANs" ossia 2 reti neurali avversarie:
#  -la prima impara a creare nuove immagini
#  -la seconda impara a discernere quelle create da quelle originali.
# Entrambe le reti quindi devono essere generate ed allenate contemporaneamente.
#
# Ci sono vari parametri con qui si può giocare, ma la configurazione attuale è la più performante;
# proprio per questo il lavoro da fare sarà quello di downgradarla per renderla più veloce e leggera.
#
# Un'altra cosa da fare sarà salvare la RN allenata una volta affinata lei e il dataset utilizzato.
# Questo attualmente non è ancora stato implementato nel codice. 
# Una volta salvata sarà possibile richiamarla con una routine senza doverla riallenare ogni volta


NB. Image_elaborator di default tenterà di utilizzare la GPU per poter accellerare il processo
(nel mio caso con una scheda nvidia 150mx passa da ~600 a ~200 secondi per epoca,
il processo durerà 150 epoche quindi sarà necessario nella fase finale ma non in fase di test [2 epoche]).
Nel caso non si avessero i drivers giusti, il programma girerà comunque ma darà degli Alert