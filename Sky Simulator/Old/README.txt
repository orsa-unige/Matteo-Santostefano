NEWS 27/12/2020:
Caricato Pupil.py
Pupil.py è un programma autosufficiente per osservare l'effetto del defocus, 
lavora con array2d e li elabora matematicamente per ottenere le PSF al variare della distanza dal fuoco.
Dato l'approccio dovrebbe essere possibile implementare altri tipi di distorsioni ottiche.

NEWS 15/12/2020:
Il programma principare per ora è Filteretor.py:
Per ora prende solo l'immagine test.png e ci applica vari filtri salvati in Kernerls.py

i filtri sono di 3 tipi:
-da funzione come il filtro gaussiano (k2)
-disegnati a mano come la croce (k1)
-da disegno, come il defocus (k3) o l'astigmatico (k4)
sarà poi da implementarli tutti da funzione

ai filtri si possono passare diversi parametri:
k1 nessuno
k2 dimensione, deviazione standard, anisotropia
k3,k4 dimensione e anisotropia


NB. se si vuole salvare l'immagine ottenuta non si può salvare direttamente 
in quanto si possiede un array da normalizzare e convertire in immagine; 
questo viene fatto da prova.py 



