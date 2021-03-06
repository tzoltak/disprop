% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/dysprop.R
\name{gini_disprop}
\alias{gini_disprop}
\title{Wspolczynnik Giniego (adaptacja).}
\usage{
gini_disprop(x, y)
}
\arguments{
\item{x}{wektor ważony}

\item{y}{wektor ważący}
}
\value{
Funkcja zwraca wartość liczbową.
}
\description{
Wylicza współczynnik Giniego zaadaptowany jako miara dyspoporcjonalności.
}
\details{
Uwaga, porównywane wektory nie mogą zawierać braków danych!
}
\examples{
vTG=c(10,40,50)
sTG=c(15,25,60)
gini_disprop(sTG, vTG)
}

