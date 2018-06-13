# Mode d'emploi de `reactor`

## Compilation

`gfortran -O3 irkc_m.f90 define_ivp.f90 reactorGen.f90  -o reactor`

* __irkc_m.f90__ : IRKC implicit-explicit integrator module     
  L.F. Shampine. B.P. Sommeijer and J.G. Verwer (2006)   
  _IRKC: an IMEX Solver for Stiff Diffusion-Reaction PDEs._    
  J. Comput. Appl. Math. 196:485-497.    
  (https://doi.org/10.1016/j.cam.2005.09.014)  
  
* __define_ivp.f90__ : module implementing the `FE` and `FI` 
  gradient routines for IRKC using a sparse representation of 
  the stoechiometry matrix following    
  A. Cangiani (2012)    
  _Biochemical pathways simulation._    
  (https://arxiv.org/abs/1208.4222)

* __reactorGen.f90__ : main program, described in    
  Z. Peng, N. Carrasco and P. Pernot (2014)    
  _Modeling of synchrotron-based laboratory simulations 
  of Titan's ionospheric photochemistry._    
  GeoResJ __1-2__:33-53    
  (https://doi.org/10.1016/j.grj.2014.03.002)
  

## Utilisation

### Arborescence de travail

* __ChemDBPublic__ contient les bases de données de référence pour la chimie et la photochimie
* __Projects__ contient les projets (un répertoire par projet, contenant toutes les inputs et outputs pour ce projet). La structure de chaque projet est la suivante:
    * __Scripts__ : codes pour la création de bases de données adaptées au projet et pour l'analyse des résultats
    * __Run__ : là où vit le code `reactor`
    * __MC_Input__ : contient les bases de données pour le projet
    * __MC_Output__ : contient les fichiers de résultats du code

### Créer un projet 

Dans le répertoire __Projects__

1. Copier un ancien projet et changer son nom (choisir le nom le plus explicite possible)

2. Aller dans le sous-répertoire __Scripts__ et éditer le fichier `getSamples.R` :
    + choisir le répertoire de travail
    + choisir le nombre de tirages Monte Carlo dans la base de données
    + choisir le mélange gazeux 

3. Lancer `getSamples.R`, qui va peupler __MC_Inputs__ et créer divers fichiers de configuration dans __Run__. 

4. Editer le fichier `Run/control.dat` (paramètres du réacteur, du flux de photons, composition du mélange gazeux...)

### Exécuter le code

Ouvrir un terminal

#### Run nominal
1. Se placer dans le répertoire du projet
2. Lancer la commande `./Scripts/OneRun_Loc.sh 0` (l'indice 0 réfère à la version nominale des bases de données de réactions)

#### Runs Monte Carlo
1. Se placer dans le répertoire du projet
2. Lancer la commande `./Scripts/MCRun_Loc.sh xxx`  où xxx est le nombre de runs souhaité. Le code fera le run nominal (indice 0) plus xxx runs aléatoires. Le nombre de runs demandé doit être inférieur ou égal au nombre d'échantillons dans __MC_Inputs__

## Résultats
Les résultats seront dans __MC_Outputs__:
1. `fracmol_xxxx.dat` contient les concentrations des espèces en fonction du temps (xxxx est le numéro du run: 0000 pour nominal)
2. `mc_rates_xxxx.dat` contient les constantes de vitesse dans les conditions du réacteur, pour le run xxxx.
3. `photo_rates_xxxx.dat` contient les les constantes de photolyse dans les conditions du réacteur, pour le run xxxx.

Pour traiter ces fichiers, on utilise des scripts R. Typiquement, dans Rstudio
1. choisir le script à exécuter
    1. Scripts/plotSpEvol.R : lit le run nominal et trace les courbes de fractions molaires  pour les neutres et les ions et calcule le taux de conso de CH4 et le taux d'ionisation
    2. Scripts/viewSpEvol.R : lit tous les runs et génère un fichier pdf dans lequel on a 1 figure d'évolution de concentration par espèce
    3.  Scripts/viewSample.R : lit tous les runs et génère des pseudo spectres de masse (non convolués par les fragmentation patterns) pour le temps final des simulations. Il génère aussi un fichier yStats_xxx qui est lu par le script suivant.
    4. Scripts/EI_MS.R :  lit yStats_xxx et génère un spectre de masse mar impact electronique. PAS FONCTIONNEL. A FAIRE...
2. choisir le répertoire de travail comme MC_Outputs
3. lancer le script

Selon les scripts, des fichiers de figures sont générés dans MC_Outputs, ou bien les figures apparaissent dans Rstudio. Tout ceci peut être customisé!!!
