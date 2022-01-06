Examen machine 2022

////Mesure de temps 

Temps                      |  small_lena_gray.png  |  tiny_lena_gray.png 
Encodage image par DFT     |    121.994            |    0.4902
Sélection des p%           |    0.007053           |    0.0003159
Reconstitution de l'image  |    15.6098            |    0.06777




////OpenMP

Temps                      |  small_lena_gray.png  |  tiny_lena_gray.png 
Encodage image par DFT     |    72.9468            |    0.276648         
Sélection des p%           |    0.009053           |    0.0003159
Reconstitution de l'image  |    9.71435            |    0.03850

J'ai décidé de paralléliser l'encodage par DFT et la reconstitution de l'image et 
elle s'est faite sur les pixels.

Pour small_lena_gray on a une accélération de 1.67 pour l'encodage et 1.61 pour
la reconstitution de l'image.
Pour tiny_lena_gray on a une accélération de 1.77 pour l'encodage et 1.758 pour
la reconstitution de l'image.

En effet pour les boucles imbriquées for au niveau de l'encodage, on a une complexité 
en accès mémoire inférieur à la complexité algorithmique, ce qui donne une situation de
cpu_bound et il y'a un intéret à paralléliser et c'est pour cela qu'on a une accélération
assez interessante.
Au niveau de la reconstition de l'image, on a une logique similaire.


////MPI 1







