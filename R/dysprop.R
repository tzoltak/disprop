#' @title Wspolczynnik D (Duncana, Pedersena, or Loosemore-Hanby'ego).
#' @description
#' Wylicza współczynnik D (Duncana, Pedersena, or Loosemore-Hanby'ego).
#' Zaimplementowano także możliwość przeprowadzenia dekompozycji na dwa komponenty
#' w postaci zaproponowanej przez Powella and Tuckera.
#' @param x wektory liczbowy
#' @param y wektory liczbowy
#' @param komponent opcjonalnie wektor definiujący podział na komponenty
#' @details
#' Uwaga, porównywane wektory nie mogą zawierać braków danych!
#' @return Funkcja zwraca wektor (jednoelementowy gdy nie jest przeprowadzana dekompozycja
#' lub (liczba_grup + 1)-elementowy, gdy jest przeprowadzana) z wartościami współczynnika (i ew. jego
#' komponentów).
#' @examples
#' vTG=c(10,40,50)
#' sTG=c(15,25,60)
#' dissim(vTG, sTG)
#' dissim(vTG, sTG, c(1, 2, 2))
#' @export
dissim = function(x, y, komponent=NULL) {
  stopifnot(
    is.numeric(x) | is.integer(x),
    is.numeric(y) | is.integer(y),
    is.null(komponent) | is.vector(komponent) | is.factor(komponent)
  )
  if (!is.null(komponent)) stopifnot(length(komponent) == length(x))
  stopifnot(
    all(x >= 0), all(y >= 0),
    any(x >  0), any(y >  0),
    all(!is.na(x)), all(!is.na(y))
  )

  if (!is.null(komponent)) {
    temp = 0.5 * abs( y / sum(y) - x / sum(x) )
    wyniki = sum(temp)
    if (!is.factor(komponent)) komponent = factor(komponent)
    for (i in levels(komponent)) {
      wyniki = c(wyniki, sum(temp[komponent == i]))
    }
    names(wyniki) = c("ogółem", levels(komponent))
    return(wyniki)
  } else {
    return(0.5 * sum( abs( y / sum(y) - x / sum(x) ) ))
  }
}
#' @title Wspolczynnik Gallaghera.
#' @description
#' Wylicza współczynnik Gallaghera.
#' Zaimplementowano także możliwość przeprowadzenia dekompozycji na komponenty
#' w postaci analogicznej do zaproponowanej przez Powella and Tuckera.
#' @param x wektory liczbowy
#' @param y wektory liczbowy
#' @param doKwadratu wartość logiczna - czy podnieść zwracaną wartość współczynnika
#' Gallaghera do kwadratu?
#' @param komponent opcjonalnie wektor definiujący podział na komponenty
#' @details
#' Uwaga, porównywane wektory nie mogą zawierać braków danych!
#'
#' Dekompozycja może być przeprowadzona tylko w sytuacji, gdy parametr \code{doKwadratu}
#' przyjmuje wartość \code{TRUE}.
#' @return Funkcja zwraca wektor (jednoelementowy gdy nie jest przeprowadzana dekompozycja
#' lub (liczba_grup + 1)-elementowy, gdy jest przeprowadzana) z wartościami współczynnika (i ew. jego
#' komponentów).
#' @examples
#' vTG=c(10,40,50)
#' sTG=c(15,25,60)
#' gallagher(vTG, sTG)
#' gallagher(vTG, sTG, komponent=c(1, 2, 2))
#' @export
gallagher = function(x, y, doKwadratu=FALSE, komponent=NULL) {
  stopifnot(
    is.numeric(x) | is.integer(x),
    is.numeric(y) | is.integer(y),
    is.logical(doKwadratu),
    is.null(komponent) | is.vector(komponent) | is.factor(komponent)
  )
  if (!is.null(komponent)) stopifnot(length(komponent) == length(x))
  stopifnot(
    all(x >= 0), all(y >= 0),
    any(x >  0), any(y >  0),
    all(!is.na(x)), all(!is.na(y)),
    length(doKwadratu) == 1
  )
  stopifnot(doKwadratu %in% c(TRUE, FALSE))

  if (!is.null(komponent) & doKwadratu != TRUE) {
    warning("Parametr 'doKwadratu' ma wartość FALSE, choć podano parametr 'komponent'! Zmieniono wartość parametru 'doKwadratu' na TRUE.",
            call. = FALSE, immediate. = TRUE)
    doKwadratu = TRUE
  }
  if (doKwadratu) {
    if (!is.null(komponent)) {
      wyniki = 0.5 * sum( ( y / sum(y) - x / sum(x) )^2 )
      if (!is.factor(komponent)) komponent = factor(komponent)
      for (i in levels(komponent)) {
        maska = komponent == i
        wyniki = c(wyniki,
                   0.5 * sum((y[maska] / sum(y) - x[maska] / sum(x) )^2) )
      }
      names(wyniki) = c("ogółem", levels(komponent))
      return(wyniki)
    } else {
      return(0.5 * sum( ( y / sum(y) - x / sum(x) )^2 ))
    }
  }
  else return(sqrt(0.5 * sum( ( y / sum(y) - x / sum(x) )^2 )))
}
#' @title Wspolczynnik 1-cos.
#' @description
#' Wylicza współczynnik 1-cos, zaproponowany przez Koppela i Diskina.
#' Zaimplementowano także możliwość przeprowadzenia dekompozycji na komponenty
#' w postaci analogicznej do zaproponowanej przez Powella and Tuckera.
#' @param x wektory liczbowy
#' @param y wektory liczbowy
#' @param pierwiastek wartość logiczna - czy wyciągnąć pierwiastek ze zwracanej wartości
#' współczynnika (uzyskując w efekcie zmodyfikowany współczynnik Gallaghera)?
#' @param komponent opcjonalnie wektor definiujący podział na komponenty
#' @details
#' Uwaga, porównywane wektory nie mogą zawierać braków danych!
#'
#' Dekompozycja może być przeprowadzona tylko w sytuacji, gdy parametr \code{pierwiastek}
#' przyjmuje wartość \code{FALSE}.
#' @return Funkcja zwraca wektor (jednoelementowy gdy nie jest przeprowadzana dekompozycja
#' lub (liczba_grup + 1)-elementowy, gdy jest przeprowadzana) z wartościami współczynnika (i ew. jego
#' komponentów).
#' @examples
#' vTG=c(10,40,50)
#' sTG=c(15,25,60)
#' cos_disprop(vTG, sTG)
#' cos_disprop(vTG, sTG, komponent=c(1, 2, 2))
#' @export
cos_disprop = function(x, y, pierwiastek=FALSE, komponent=NULL) {
  stopifnot(
    is.numeric(x) | is.integer(x),
    is.numeric(y) | is.integer(y),
    is.logical(pierwiastek),
    is.null(komponent) | is.vector(komponent) | is.factor(komponent)
  )
  if (!is.null(komponent)) stopifnot(length(komponent) == length(x))
  stopifnot(
    all(x >= 0), all(y >= 0),
    any(x >  0), any(y >  0),
    all(!is.na(x)), all(!is.na(y)),
    length(pierwiastek) == 1
  )
  stopifnot(pierwiastek %in% c(TRUE, FALSE))

  if (!is.null(komponent) & pierwiastek == TRUE) {
    warning("Parametr 'pierwiastek' ma wartość TRUE, choć podano parametr 'komponent'! Zmieniono wartość parametru 'pierwiastek' na FALSE.",
            call. = FALSE, immediate. = TRUE)
    pierwiastek = FALSE
  }
  if (pierwiastek) {
    return(sqrt( 1 - sum(x * y) / sqrt( sum(x^2) * sum(y^2) ) ))
  } else {
    if (!is.null(komponent)) {
      wyniki = 1 - sum(x * y) / sqrt( sum(x^2) * sum(y^2) )
      if (!is.factor(komponent)) komponent = factor(komponent)
      for (i in levels(komponent)) {
        maska = komponent == i
        wyniki = c(wyniki,
                   0.5 * sum( ( x[maska] / sqrt(sum(x^2)) - y[maska] / sqrt(sum(y^2)) )^2) )
      }
      names(wyniki) = c("ogółem", levels(komponent))
      return(wyniki)
    } else {
      return( 1 - sum(x * y) / sqrt( sum(x^2) * sum(y^2) ) )
    }
  }
}
#' @title Wspolczynnik Giniego (adaptacja).
#' @description
#' Wylicza współczynnik Giniego zaadaptowany jako miara dyspoporcjonalności.
#' @param x wektor ważony
#' @param y wektor ważący
#' @details
#' Uwaga, porównywane wektory nie mogą zawierać braków danych!
#' @return Funkcja zwraca wartość liczbową.
#' @examples
#' vTG=c(10,40,50)
#' sTG=c(15,25,60)
#' gini_disprop(sTG, vTG)
#' @export
gini_disprop = function(x, y) {	# x ważone y
  stopifnot(
    is.numeric(x) | is.integer(x),
    is.numeric(y) | is.integer(y)
  )
  stopifnot(
    all(x >= 0), all(y >= 0),
    any(x >  0), any(y >  0),
    all(!is.na(x)), all(!is.na(y)),
    length(x) == length(y)
  )

  porzadek = order(x / y)
  x = x[porzadek] / sum(x)
  y = y[porzadek] / sum(y)
  return( 1 - 2 * sum( y * (cumsum(x) - x / 2) ))
}
#' @title Dywergencja Kullbacka-Leiblera (entropia wzgledna).
#' @description
#' Wylicza dywergencje Kullbacka-Leiblera (entropię względną).
#' @param x wektor ważony
#' @param y wektor ważący
#' @param p podstawa logarytmu
#' @details
#' Uwaga, porównywane wektory nie mogą zawierać braków danych!
#' @return Funkcja zwraca wartość liczbową.
#' @examples
#' vTG=c(10,40,50)
#' sTG=c(15,25,60)
#' divergence_kl(sTG, vTG)
#' @export
divergence_kl = function(x, y, p=2) {	# x ważone y
  stopifnot(
    is.numeric(x) | is.integer(x),
    is.numeric(y) | is.integer(y),
    is.numeric(p) | is.integer(p)
  )
  stopifnot(
    all(x > 0), all(y > 0),
    all(!is.na(x)), all(!is.na(y)),
    length(x) == length(y),
    length(p) == 1
  )

  return(sum( y * log(y * sum(x) / (x * sum(y)), p) ) / sum(y))
}
#' @title Wariancja (wazona wariancja odsetkow).
#' @description
#' Wylicza "ważoną wariancję odsetków" - wariancję zaadaptowaną do roli miary
#' dysproporcjonalności.
#' @param x wektor ważony
#' @param y wektor ważący
#' @param unorm stała normująca lub NULL, jeśli współczynnik ma nie być normowany -
#' p. szczegóły
#' @param komponent opcjonalnie lista wektorów definiujących podział na komponenty;
#' kolejne elementy listy definiują kolejne poziomy podziału (przy czym podział na każdym
#' kolejnym poziomie będzie traktowany jako zagnieżdżony w poziomach wyższych)
#' @details
#' Uwaga, porównywane wektory nie mogą zawierać braków danych!
#'
#' W przypadku normowania wartości wskaźnika istnieją dwa typowe wybory wartości stałej
#' normującej:
#' \itemize{
#'   \item{1 - w takim wypadku zwrócony zostanie kwadrat "ważonego współczynnika
#'         zmienności odsetków", czyli wartość wskaźnika zostanie podzielona przez kwadrat
#'         ilorazu \code{sum(x) / sum(y)}. Tak przekształcony wskaźnik nie posiada
#'         ograniczenia górnego w 1 i dekomponuje się w nieco inny sposób (p. przykłady).}
#'   \item{Liczba obiektów w analizowanej zbiorowości. Jeśli wartość wskaźnika, oprócz
#'         podzielenia przez przez kwadrat ilorazu \code{sum(x) / sum(y)}, zostanie
#'         podzielona również przez liczbę obiektów w analizowanej zbiorowości, to dla
#'         liczby obiektów zbiegającej do nieskończoności będzie on posiadał ograniczenie
#'         górne w 1, wyznaczane przez sytuację, gdy cała wartość x skupiona jest w jednym
#'         obiekcie. Oczywiście unormowanie to ma sens, o ile w ramach rozpatrywanego
#'         problemu daje się sensownie wyróżnić liczbę obiektów w zbiorowości (typowo może
#'         to być suma \code{y}). Tak wyliczony wskaźnik nie poddaje się dekompozycji.}
#' }
#' Pewien problem ze wskaźnikiem polega na tym, że jego wartość bardzo szybko spada nawet
#' przy niewielkich odstępstwach od maksymalnej możliwej dysproporcjonalności, co jest
#' niezbyt intuicyjne (szczególnie wyraźnie widać to w przypadku drugiego z ww. typów
#' unormowania) - p. przykłady.
#' @return funkcja zwraca wartość liczbową
#' @examples
#' vTG=c(10,40,50)
#' sTG=c(15,25,60)
#' var_disprop(sTG, vTG)
#'
#' # przykład dekomponowalności
#' dane = data.frame(
#'   gr = rep(1:3, each=3),
#'   v  = c(150, 90, 100, 120, 80, 70, 200, 20, 30),
#'   s  = c(  3,  2,   2,   5,  1,  1,   3,  0,  0)
#' )
#' wewnGrup = c(
#'   gr1 = with(subset(dane, gr==1), var_disprop(s, v)),
#'   gr2 = with(subset(dane, gr==2), var_disprop(s, v)),
#'   gr3 = with(subset(dane, gr==3), var_disprop(s, v))
#' )
#' sumy = aggregate(dane[, -1], list(gr=dane$gr), sum)
#' wewnGrup = weighted.mean(wewnGrup, sumy$v)
#' miedzyGrup = with(sumy, var_disprop(s, v))
#' miedzyGrup + wewnGrup
#' ( ogolem = with(dane, var_disprop(s, v)) )
#' all.equal(miedzyGrup + wewnGrup, ogolem)
#'
#' # przykład własności unormowań
#' n = c(10, 100, 1000)
#' zestawienie = matrix(NA, nrow=6, ncol=length(n),
#'   dimnames=list(c("Gini n1=1", "Gini n1=2", "var  n1=1, unorm=n",
#'                   "var  n1=2, unorm=n", "var  n1=1, unorm=1", "var  n1=2, unorm=1"),
#'                 paste0("n=", n)))
#' for(i in 1:length(n)) {
#'   zestawienie[, i] = c(
#'     gini_disprop(c(1,0), c(1, n[i]-1)),
#'     gini_disprop(c(1,0), c(2, n[i]-2)),
#'     var_disprop (c(1,0), c(1, n[i]-1), n[i]),
#'     var_disprop (c(1,0), c(2, n[i]-2), n[i]),
#'     var_disprop (c(1,0), c(1, n[i]-1), 1),
#'     var_disprop (c(1,0), c(2, n[i]-2), 1)
#'   )
#' }
#' zestawienie
#' # n1=1 oznacza, że jeden obiekt posiada wszystko (a n-1 obiektów nie posiada nic)
#' # n1=2 oznacza, że dwa obiekty posiadają wszystko (a n-2 obiektów nie posiada nic)
#'
#' # dekompozycja wskaźnika wyliczonego z unorm=1
#' dane = data.frame(
#'   gr = rep(1:3, each=3),
#'   v  = c(150, 90, 100, 120, 80, 70, 200, 20, 30),
#'   s  = c(  3,  2,   2,   5,  1,  1,   3,  0,  0)
#' )
#' wewnGrup = c(
#'   gr1 = with(subset(dane, gr==1), var_disprop(s, v, unorm=1)),
#'   gr2 = with(subset(dane, gr==2), var_disprop(s, v, unorm=1)),
#'   gr3 = with(subset(dane, gr==3), var_disprop(s, v, unorm=1))
#' )
#' sumy = aggregate(dane[, -1], list(gr=dane$gr), sum)
#' wewnGrup = weighted.mean(wewnGrup * with(sumy, {(s / v)^2 / (sum(s) / sum(v))^2}),
#'                          sumy$v)
#' miedzyGrup = with(sumy, var_disprop(s, v, unorm=1))
#' miedzyGrup + wewnGrup
#' ( ogolem = with(dane, var_disprop(s, v, unorm=1)) )
#' all.equal(miedzyGrup + wewnGrup, ogolem)
#'
#' @export
var_disprop = function(x, y, unorm=NULL, komponent=NULL) {  # x ważone y
  stopifnot(
    is.numeric(x) | is.integer(x),
    is.numeric(y) | is.integer(y),
    is.numeric(unorm) | is.null(unorm),
    is.list(komponent) | is.null(komponent)
  )
  stopifnot(
    all(x >= 0), all(y > 0),
    all(!is.na(x)), all(!is.na(y)),
    length(x) == length(y)
  )
  if (!is.null(unorm)) {
    stopifnot(length(unorm) == 1)
    stopifnot(unorm > 0)
  }
  if (!is.null(komponent)) {
    stopifnot((!is.null(komponent) & unorm == 1) | is.null(unorm))
    stopifnot(all(unlist(lapply(komponent, is.vector)) |
                    unlist(lapply(komponent, is.factor))))
    stopifnot(length(unique(unlist(lapply(komponent, length)))) == 1)
    stopifnot(unique(unlist(lapply(komponent, length))) == length(x))
  }

  if (!is.null(komponent)) {
    wyniki = rep(NA, length=length(komponent))
    for (i in 1:length(komponent)) {
      sumyGr = aggregate(data.frame(x, y), komponent[1:i], sum)
      wyniki[i] = var_disprop(sumyGr$x, sumyGr$y, unorm=unorm)
      if (i > 1) wyniki[i] = wyniki[i] - sum(wyniki[1:(i - 1)])
    }
    wyniki = c(var_disprop(x, y, unorm=unorm), wyniki)
    wyniki = c(wyniki, wyniki[1] - sum(wyniki[-1]))
    names(wyniki) = c("ogółem", names(komponent), "wewnątrzgrupowa")
    return(wyniki)
  } else {
    sx = sum(x)
    sy = sum(y)
    if (!is.null(unorm)) return( ( sy * sum(x^2 / y) - sx^2      ) / sx^2 / unorm )
    else                 return( (      sum(x^2 / y) - sx^2 / sy ) / sy           )
  }
}
