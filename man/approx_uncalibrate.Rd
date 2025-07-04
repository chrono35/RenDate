% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RenDate.R, R/RenDate_20250701.R
\name{approx_uncalibrate}
\alias{approx_uncalibrate}
\title{approx_uncalibrate}
\usage{
approx_uncalibrate(
  calibrated_result,
  grid_size = 100,
  tmin = NULL,
  tmax = NULL,
  show_progress = TRUE
)

approx_uncalibrate(
  calibrated_result,
  grid_size = 100,
  tmin = NULL,
  tmax = NULL,
  show_progress = TRUE
)
}
\arguments{
\item{calibrated_result}{Un objet de classe `RenDate`, typiquement produit par une calibration}

\item{grid_size}{Taille de la grille des mesures simulées (nombre de points)}

\item{tmin}{Optionnel : limite inférieure de la courbe de calibration à utiliser (en calBP)}

\item{tmax}{Optionnel : limite supérieure de la courbe de calibration à utiliser (en calBP)}

\item{show_progress}{Affiche une barre de progression si `TRUE` (par défaut)}
}
\value{
Une liste de type `UncalibratedMeasure` contenant, pour chaque date, la grille des mesures, la densité estimée des mesures, la courbe utilisée, et la densité calibrée d’origine

Une liste de type `UncalibratedMeasure` contenant, pour chaque date, la grille des mesures, la densité estimée des mesures, la courbe utilisée, et la densité calibrée d’origine
}
\description{
Approximation de l'opération de **"décalibration"** : reconstitue la densité de mesure à partir d'une densité calibrée (`RenDate`), en inversant le modèle d'observation par pseudo-inversion matricielle.

Approximation de l'opération de **"décalibration"** : reconstitue la densité de mesure à partir d'une densité calibrée (`RenDate`), en inversant le modèle d'observation par pseudo-inversion matricielle.
}
\details{
Cette fonction approxime l'inversion du processus de calibration, en utilisant la pseudo-inverse de la matrice de convolution entre les densités calibrées et la courbe de calibration.
Elle suppose que la calibration a été effectuée à l’aide d’une courbe sous forme de fichier `.rda` dans le dossier `data/`, contenant une matrice à trois colonnes (année calBP, mu, sigma).

L’estimation peut être bruitée et doit être interprétée avec prudence. Elle est utile à des fins exploratoires ou pédagogiques.

Cette fonction approxime l'inversion du processus de calibration, en utilisant la pseudo-inverse de la matrice de convolution entre les densités calibrées et la courbe de calibration.
Elle suppose que la calibration a été effectuée à l’aide d’une courbe sous forme de fichier `.rda` dans le dossier `data/`, contenant une matrice à trois colonnes (année calBP, mu, sigma).

L’estimation peut être bruitée et doit être interprétée avec prudence. Elle est utile à des fins exploratoires ou pédagogiques.
}
