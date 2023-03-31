Pierre ZACHARY
Rendu VTK
https://github.com/Pierre-ZACHARY/vtktp.git

Voir le pdf pour les images ou le dossier images

Explication du code :
Pour optimiser un maximum le code, j’ai repris les principes de l’out of core et du tp7 avec le
parallélisme. Dans les grandes lignes, voici comment fonctionne mon code :
- On répartit les Z entre chaque processus
- Chaque processus va charger une partie de ses Z dans le fichier raw ( il divise ses Z
  par numPasses )
- On fait une boucle pour effectuer plusieurs azimuth et/ou modifier les valeurs de
  coupe à chercher
- Chaque processus fait le rendu de ses Z en plusieurs fois : il charge les données
  pour une partie de ses Z, fait le rendu et garde le rendu finale pour ses Z via un
  zbuffer local
- Une fois que tous les processus ont fait le rendu de leurs Z, on fait un MPI Gather
  sur root pour récupérer tous les rgba et zbuffer et calculer l’image finale