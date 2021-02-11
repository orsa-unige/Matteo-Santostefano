NEWS 05/02/2021:
'TheSky.py' è ora 'SkyConstructer.py'. La creazione delle stelle viene fatta in maniera matematicamente più pulita e corretta, i dati vengono salvati su un fits con le proprietà della CCD.
'SkyConstructer.py' lavora in tandem con 'Query.py'. 'Query.py' dato un punto nel cielo (per ora in in coordinate icrs) chiama SIMBAD e trascrive le posizioni e magnitudo delle stelle in una regione dal raggio di 7 arcmin (grandezza massima della CCD), queste info sono poi passate a 'SkyConstructer.py' che le userà per costruire il cielo simulato.
Infine 'MainSky.py' è solo un codice che permette di scegliere se lavorare in maniera automatica o manuale. Il processo automatizzato richiede solo le coordinate centrali del cielo da simulare (i dati delle distorsioni sono generate randomicamente). Il processo manuale è più lungo ma permette di calibrare il cielo in maniera diretta.

NEWS 09/01/2021:
Caricato un nuovo programma 'TheSky.py' che è in grado di gestire immagini simulate in maniera completamente diversa da i precedenti.
Per questo motivo i file ormai obsoleti sono stati spostati in una nuova cartella.

'TheSky.py' può sia mostrare che salvare in fits un cielo creato semi-randomicamente da parametri impostati dall'utente:

una volta fatto partire il programma appariranno le seguenti richieste:

Build your Star!

Do you want a In-focus star? (Y/N)
n                                 #se si dice di si si visualizzerà una 'stella' in-focus e si andrà direttamente alla domanda 'Do you like your model?'
Do you want a sample star? (Y/N)
n                                 #se si dice di si si visualizzerà una 'stella'  non in-focus e si andrà direttamente alla domanda 'Do you like your model?'
How much defocus? 2               #qui vengono richiesti dei parametri numerici (i qui limiti sono ancora da definire) ed infine si visualizzerà la 'stella'
How much coma? 1                  #generata con essi con un inclinazione standar che poi sarà scelta casualmente nella creazione del cielo
How much bessel? 3
How much astigmatism? 2
The phase between coma and stigmatism? (deg) 45
Do you like your model? (Y/N)  #se si mette no, si ricomincia da capo
y
Bulding your sky, it can take few seconds #il programma gira e crea un cielo con stelle con posizioni e intensità randomiche

Do you want do visualize? (Y/N)  #se si dice di sì si visualizzerà in un grafico il cielo creato
y
Do you want to save it? (Y/N)    #se si dice di si il cielo verrà salvato come fits
y
