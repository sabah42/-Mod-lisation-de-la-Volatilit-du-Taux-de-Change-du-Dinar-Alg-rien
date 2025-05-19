# Modélisation de la Volatilité du Taux de Change du Dinar Algérien

## Table des matières
- [Description du projet](#description-du-projet)
- [Structure du projet](#structure-du-projet)
- [Jeu de données](#jeu-de-données)
- [Outils et bibliothèques](#outils-et-bibliothèques)
- [Nettoyage et préparation des données](#nettoyage-et-préparation-des-données)
- [Méthodologie d’analyse](#méthodologie-danalyse)
  - [Modèles SARIMA](#modèles-sarima)
  - [Modèles MAR et MAR-ARCH](#modèles-mar-et-mar-arch)
- [Principaux résultats](#principaux-résultats)
- [Recommandations](#recommandations)

##  Description du projet
Ce projet vise à modéliser et prévoir la volatilité du taux de change du **Dinar Algérien (DZD)** face à quatre devises majeures (**USD, EUR, GBP, JPY**).  
L’objectif est d’aider le **Crédit Populaire d’Algérie** à anticiper les fluctuations du taux de change pour réduire les risques financiers liés à l’import-export et à la gestion de devises.

L’approche se base sur des techniques d’analyse de séries chronologiques, incluant des modèles linéaires et non linéaires, et une comparaison de performance.

## Structure du projet
**data/**  
Contient les séries temporelles mensuelles des taux de change (source : Banque d’Algérie).

**notebooks/**  
- `Tp` : Analyse, modélisation et comparaison des modèles SARIMA, MAR et MAR-ARCH.

**results/**  
- Graphiques de séries, fonctions d’autocorrélation, diagnostics de stationnarité  
- Visualisation de la variance conditionnelle et des erreurs de prévision  
- Comparatifs RMSE, AIC, BIC

**rapport/**  
- `PFE.pdf` : Rapport complet du mémoire

## Jeu de données
Les séries temporelles utilisées couvrent plusieurs années ( 2005-2015) avec les variables suivantes :
- **Taux de change DZD/USD, DZD/EUR, DZD/GBP, DZD/JPY**
- **Données économiques macro (optionnel)**

## Outils et bibliothèques
- `MATLAB` : Modélisation MAR, MAR-ARCH, estimation par maximum de vraisemblance (EM)
- `EViews` : Analyse SARIMA et prévisions
- `LaTeX` : Rédaction du mémoire
- `Excel` : Pré-traitement et organisation des données

## Nettoyage et préparation des données
- Importation des séries temporelles
- Détection et traitement des valeurs manquantes
- Tests de stationnarité : **ADF**, **KPSS**
- Différenciation des séries pour obtenir la stationnarité
- Normalisation pour le traitement statistique

## Méthodologie d’analyse

### Modèles SARIMA
- Application de la méthode **Box & Jenkins** (identification, estimation, validation)
- Prévision de court terme
- Limites : ne capte pas la **volatilité conditionnelle**

### Modèles MAR et MAR-ARCH
- Implémentation du **modèle MAR (Wong et Li, 2000)** pour capter les changements de régime
- Extension **MAR-ARCH** pour modéliser la **variance conditionnelle** (volatilité)
- Estimation des paramètres par **maximum de vraisemblance (algorithme EM)**
- Études de simulation pour valider les performances

## Principaux résultats
- Le **modèle MAR-ARCH** offre une **réduction du RMSE de 22 %** par rapport au modèle SARIMA
- **AIC moyen inférieur de 18 %**, indiquant un meilleur ajustement
- Meilleure capture de la **volatilité stochastique**, en particulier pour le **taux DZD/USD**
- Les modèles MAR-ARCH sont plus robustes dans des contextes économiques instables

## Recommandations
- Utiliser des modèles non linéaires (MAR-ARCH) pour surveiller les devises soumises à de fortes fluctuations
- Déployer le modèle MAR-ARCH dans les outils de prévision du Crédit Populaire d’Algérie
- Adapter la stratégie de change selon les périodes de forte volatilité identifiées par le modèle
