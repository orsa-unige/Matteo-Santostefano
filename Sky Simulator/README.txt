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



