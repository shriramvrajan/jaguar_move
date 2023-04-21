# Modeling jaguar movement ecology

## Dataset

Morato, R. G. et al. Jaguar movement database: a GPS-based movement dataset of an apex predator in the Neotropics. Ecology 99, 1691–1691 (2018). <https://doi.org/10.1002/ecy.2379>

## Model

Step-selection model based on a 'neighborhood' grid at each movement decision point. Formulated and parametrized in different ways, each fitted using MLE techniques (`optim` in R).

## Model formulations

0. Null model: step-selection has equal probability for each pixel
1. K: (Johnson, 2021) Temperature, precipitation, human footprint, forest cover, elevation, and slope influence step selection.
2. RW: Human footprint, forest cover, elevation, slope, distance to roads, and distance to water.
3. H: **Only** home range influences step selection - penalty for moving away from empirical probability distribution of steps.
4. RWH: All predictors from RW + H.

## Bibliography

- Johnson, K. Inter-individual variation in landscape resistance: jaguars as a case study. (2021). Undergraduate thesis, unpublished.
