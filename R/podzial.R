#' @title Podzial mandatow (i ew. innych rzeczy)
#' @description
#' Funkcja wylicza podział mandatów (ew. innych niepodzielnych rzeczy) na podstawie rozkładu głosów (ew. innych wag).
#' @param x wektor z liczbą głosów (wag)
#' @param m liczba mandatów do podziału
#' @param regula ciąg znaków definiujacy metodę podziału; jeden z: "D'Hondt", "Hamilton", "Hare-Niemeyer", "Jefferson", "Saint-League", "Webster", "DH", "H", "HN", "J", "SL", "W"
#' @details
#' Równoważne sobie wartości parametru \code{regula}:
#' \itemize{
#' \item "Saint-League", "Webster", "SL", "W";
#' \item "D'Hondt", "Jefferson", "DH", "J";
#' \item "Hamilton", "Hare-Niemeyer", "H", "HN".
#' }
#' @return funkcja zwraca wektor liczb całkowitych o długości równej długości wektora \code{x}
#' @examples
#' v1 = c(47, 16, 15.8, 12, 6.1, 3.1) * 1000
#' podziel(v1, 10, "SL")
#' podziel(v1, 10, "DH")
#' podziel(v1, 10, "HN")
#' # paradoks Alabamy
#' vA = c(15, 15, 9, 5, 5, 2) * 100
#' podziel(vA, 25, "HN")
#' podziel(vA, 26, "HN")
#' @export
podziel = function(x, m, regula) {
  stopifnot(
    is.numeric(x) | is.integer(x),
    is.numeric(m) | is.integer(m),
    is.character(regula)
  )
  stopifnot(
    !any(is.na(x)),
    as.integer(m) == m,
    length(m) == 1,
    length(regula) == 1,
    regula %in% c("D'Hondt", "Hamilton", "Hare-Niemeyer", "Jefferson",
                  "Saint-League", "Webster",
                  "DH", "H", "HN", "J", "SL", "W")
  )
  stopifnot(
    !is.na(m), !is.na(regula)
  )
  if (!all(as.integer(x) == x)) warning("Argument x zawiera wartości niecałkowite.")

  if (regula %in% c("Saint-League", "Webster", "SL", "W")) {
    podzial = rep(0, length(x))
    while (sum(podzial) < m) {
      q = (x / (2 * podzial + 1))  # iloraz
      ktoremuDodac = q == max(q)
      podzial[ktoremuDodac] = podzial[ktoremuDodac] + 1
    }
  } else if (regula %in% c("D'Hondt", "Jefferson", "DH", "J")) {
    podzial = rep(0, length(x))
    while (sum(podzial) < m) {
      q = (x / (podzial + 1))  # iloraz
      ktoremuDodac = q == max(q)
      podzial[ktoremuDodac] = podzial[ktoremuDodac] + 1
    }
  } else { # "Hamilton", "Hare-Niemeyer", "H", "HN"
    q = sum(x) / m  # kwota Hare
    podzial = floor(x / q)  # część całkowita
    while (sum(podzial) < m) {  # rozdzielamy pozostałe mandaty
      r = x / q - podzial  # reszta
      ktoraNajwieksza = r == max(r)
      podzial[ktoraNajwieksza] = podzial[ktoraNajwieksza] + 1
    }
  }

  if (sum(podzial) > m) {
    r = (x / sum(x) / m)
    warning(paste0(
      "Rozdzielono więcej mandatów, niż było dostępnych!\n",
      "Elementy ex aequo najbardziej nadreprezentowane\n",
      "(miara nadrepr. - część ułamkowa z dzielenia x / kwota Hare):\n",
      "  nr el.   nadrepr.\n",
      paste0(
        "  ",
        format((1:length(x))[ktoremuDodac], width=5),
        "     ",
        round(r[ktoremuDodac], 4),
        collapse="\n"
      )
    ))
  }
  if (!is.null(names(x))) names(podzial) = names(x)
  return(podzial)
}
