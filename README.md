# Homogénéisation Périodique – Résolution Numérique du Problème de Poisson

## Description du projet
Ce projet présente une étude complète et une implémentation numérique autour de l’**homogénéisation périodique appliquée au problème de Poisson** dans des matériaux à microstructure.  
L’objectif est de comparer la **solution exacte** et la **solution homogénéisée** en utilisant les méthodes des éléments finis, et d’analyser la convergence numérique selon différents paramètres de discrétisation.

Ce travail illustre non seulement la **rigueur mathématique** mais aussi la **capacité d’implémentation algorithmique** et la **validation numérique** de méthodes avancées en calcul scientifique.

---

## Objectifs pédagogiques et techniques
- Résoudre numériquement le problème de Poisson avec conditions de Dirichlet.
- Étudier l’effet de la périodicité (ϵ) de la microstructure sur la précision des solutions.
- Mettre en place une **méthode d’homogénéisation périodique** pour réduire le coût de calcul.
- Valider la solution homogénéisée par comparaison qualitative et quantitative avec la solution exacte.
- Évaluer les erreurs numériques (**normes L2 et H1**) et analyser la convergence.

---

##  Méthodologie
1. **Solution exacte par éléments finis (P1)**  
   - Implémentation d’un solveur FEM pour différentes valeurs du tenseur de diffusion A(x).  
   - Comparaison des résultats obtenus avec les solutions analytiques connues.

2. **Problème homogénéisé**  
   - Résolution des **problèmes de cellule périodiques** (méthode variationnelle + pénalisation).  
   - Calcul du **tenseur homogénéisé A_eff** à partir des solutions numériques.  
   - Intégration de A_eff dans le problème global pour résoudre une version simplifiée du problème.

3. **Comparaison et validation**  
   - Comparaison des solutions : exactes vs homogénéisées.  
   - Étude des erreurs L2 et H1 en fonction de la période ϵ.  
   - Analyse de la stabilité et des limites numériques lorsque ϵ → 0.

---

## Résultats clés
- La solution homogénéisée **approxime très fidèlement** la solution exacte pour des petites valeurs de ϵ.  
- La méthode réduit considérablement le **coût de calcul** tout en maintenant une **bonne précision**.  
- Les erreurs L2 décroissent plus rapidement que les erreurs H1, confirmant la convergence du schéma.  
- Pour des valeurs très petites de ϵ, une saturation des erreurs est observée à cause des limites de discrétisation.

---

## Compétences démontrées
- **Analyse mathématique** : équations aux dérivées partielles, homogénéisation périodique, théorème de Lax-Milgram.  
- **Programmation scientifique** : implémentation en MATLAB/Octave des solveurs FEM, discrétisation, assemblage des matrices.  
- **Validation numérique** : comparaison des solutions exactes et approchées, calcul d’erreurs, étude de convergence.  
- **Optimisation du calcul** : réduction du coût de simulation grâce à l’homogénéisation.  
- **Communication scientifique** : rédaction et présentation claire des résultats numériques et théoriques.

---

##  Organisation du projet
- `principal_dirichlet_cellule.m` → Solveur FEM pour le problème exact.  
- `periodique_cellule.m` → Résolution des problèmes de cellule périodiques.  
- `principal_dirichlet_cellule_homogene.m` → Solveur utilisant les coefficients homogénéisés.  
- **Rapport complet** → Explications détaillées, validation et résultats graphiques.

---

##  Perspectives
- Extension à des problèmes de **Neumann** ou **conditions mixtes**.  
- Application à des **matériaux composites** plus complexes.  
- Intégration avec des solveurs FEM avancés (FreeFEM++, FEniCS, Firedrake).  
- Étude de la convergence en 3D.

---

## Conclusion:
Ce projet démontre :
- Une **maîtrise solide en mathématiques appliquées et calcul scientifique**.  
- La **capacité à implémenter des algorithmes complexes** en programmation numérique.  
- Une **rigueur d’analyse et de validation** avec des résultats concrets.  
- Des compétences transférables en **modélisation numérique**, **simulation physique**, **optimisation** et **data science**.

##  Auteur
**Bréhima Samaké**  
 Octobre 2024  
Projet académique en méthodes numériques appliquées aux équations aux dérivées partielles.
