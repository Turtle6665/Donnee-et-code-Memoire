# Ensemble des données et codes réalisés dans le cadre du mémoire de Diego Nève

![Python](https://img.shields.io/badge/python-3.8-blue.svg) ![License](https://img.shields.io/badge/license-MIT-green.svg)

Ce référentiel GitHub contient les données et les codes associés au [mémoire de Diego Nève](http://hdl.handle.net/2078.1/thesis:43410), intitulé "Evaluation de l’effet de la mise en place d’une rue a sens unique sur la concentration en dioxyde d’azote à Bruxelles".

## Description

La pollution de l'air est un problème de santé mondial majeur, causant des millions de décès prématurés chaque année et réduisant l'espérance de vie en bonne santé. Parmi les polluants, le dioxyde d'azote (NO2), principalement émis par les processus de combustion, reste un enjeu majeur en Région de Bruxelles-Capitale (RBC). Dans le cadre du plan de mobilité Good Move, de nombreux changements de circulation ont été mis en place à Bruxelles depuis 2020, visant notamment à améliorer la qualité de l'air.

Ce travail évalue l'effet de la mise en sens unique de l'avenue E. Cambier sur les concentrations de NO2 dans le quartier des Azalées à Bruxelles. Le modèle SIRANE [1], utilisant des estimations d'émissions, des conditions météorologiques et une géométrie simplifiée du réseau de rues, a été paramétré pour la région en 2018. Les résultats ont été comparés aux mesures du réseau télémétrique de Bruxelles Environnement, montrant une bonne adéquation aux critères d'acceptation de Chang et al. [2] et aux *model quality objectives* de Thunis et al. [3], bien que certaines améliorations soient nécessaires.

La mise en place de la rue à sens unique a entraîné une faible diminution (0,07%) des concentrations moyennes de NO2 à l'échelle du quartier. Cette diminution s'est concentrée principalement dans les zones moins polluées, tandis que des augmentations ont été observées dans les zones déjà plus exposées, notamment des rues canyons. Les implications pour la santé sont probablement très faibles.

Cependant, il est important de souligner que ce travail se concentre uniquement sur l'effet spécifique de cette mesure ponctuelle et ne peut être généralisé à l'ensemble des mesures de Good Move. Des études à plus grande échelle seraient nécessaires pour évaluer l'effet total de ces changements sur la qualité de l'air à Bruxelles.

## Contributions ?
Ce travail est fini et les contributions ne sont donc pas ouvertes. Par contre, si vous êtes intéressé de contribuer à `PySiWrap` (un wrapper Python pour SIRANE), rendez-vous dans le [dépôts prévu à cette effet](https://github.com/Turtle6665/PySiWrap). Une version plus a jour y est présente.

## Comment reproduire les résultats

Attention, une majorité des données initiales ne sont pas présentes pour des raisons de droit d'auteur.

Pour reproduire les résultats obtenus dans ce mémoire, suivez ces étapes :

0. Télécharger le git.
1. Le Preprocessing (réadaptation de certaines données, ...) est déjà réalisée et déjà introduit dans les inputs SIRANE.
2. Installer SIRANE. Si vous n'avez pas de version de SIRANE, rendez-vous sur le [site officiel](http://air.ec-lyon.fr/SIRANE/index.php) et demandez une licence.
    * Copier les executables (sirane.exe et \*.dll) dans le dossier [SIRANE](/SIRANE/)
    * Copier le fichier SETTINGS (avec les fichiers Defaut.dat) dans le dossier [SIRANE/SETTINGS](/SIRANE/SETTINGS/)
3. Installer les bonnes version des librairies python ([requirment.txt](requirment.txt)). Je recommande l'utilisation d'anaconda pour installer toutes les dépendances dans environments différent.
4. Lancer les différents scripts python depuis le répertoire principale avec par exemple~:
```python
python PythonScipts/Effet-de-bord/effet-de-bord.py
```
Les résultats de SIRANE sont directement analysé et les figures, tableaux et autres sont enregistées dans le dossier [Resultats](Resultats/)


## Auteur

Diego Nève ([@Turtle6665](https://github.com/Turtle6665))

## Licence

Ce projet est sous licence MIT. Voir le fichier [LICENSE](LICENSE) pour plus de détails.

## Références
[1] SOULHAC L. (2023). SIRANE. En ligne sur [http://air.ec-lyon.fr/SIRANE/](http://air.ec-lyon.fr/SIRANE/)

[2] Chang, J. C., & Hanna, S. R. (2004). Air quality model performance evaluation. Meteorology and Atmospheric Physics, 87(1‑3). doi:[10.1007/s00703-003-0070-7](https://doi.org/10.1007/s00703-003-0070-7)

[3] Thunis, P., Pederzoli, A., & Pernigotti, D. (2012). Performance criteria to evaluate air quality modeling applications. Atmospheric Environment, 59, 476‑482. doi:[10.1016/j.atmosenv.2012.05.043](https://doi.org/10.1016/j.atmosenv.2012.05.043)
