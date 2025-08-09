# Modeling jaguar movement ecology

## Dataset

Morato, R. G. et al. Jaguar movement database: a GPS-based movement dataset of an apex predator in the Neotropics. Ecology 99, 1691–1691 (2018). <https://doi.org/10.1002/ecy.2379>

## Model

Path propagation model based on a 'neighborhood' grid at each movement decision point. Formulated and parametrized in different ways, each fitted using MLE techniques (`optim` in R). Conceived as an alternative to traditional step selection approach to animal telemetry data.