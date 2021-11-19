#  Licence ----
#
#  Copyright or © or Copr. CNRS	2021
#
#
# This package is under
# The “Creative Commons Attribution-ShareAlike International License” version 4.0

#  A copy of the License is available at
#  http://www.r-project.org/Licenses/
#  _________________________________________________________________________________

# Version 2021-11-19
#

#' RenDate; inspiré du package Bchron (https://cran.r-project.org/web/packages/Bchron/index.html) développé par Andrew Parnell <Andrew.Parnell at mu.ie>, 
#' ( Haslett J, Parnell AC (2008). “A simple monotone process with application to radiocarbon-dated depth chronologies.” 
#' Journal of the Royal Statistical Society: Series C (Applied Statistics), 57(4), 399–418.
#'  http://onlinelibrary.wiley.com/doi/10.1111/j.1467-9876.2008.00623.x/full. )
#' est modifié pour pouvoir l'utiliser avec des mesures magnétiques, ou d'autres sortes de mesures nécessitants des valeurs décimales
#' @author "Philippe DUFRESNE"
#' @name RenDate
#' @docType package

#' Stocke une courbe dans le répertoire / data
#' Cette fonction est utilisée avec la fonction calibrate
#' @seealso \cite{\code{calibrate} }
#' @export
createCalCurve = function(name,
                          cal_ages,
                          uncal_ages,
                          one_sigma=rep(0,length(cal_ages)),
                          pathToCalCurves = paste0(getwd(),'/','data') )
{
  
  # This function creates a calibration curve and puts it in the appropriate place for future use with Bchron
  
  # First identify the place where the file needs to go
  file_loc = pathToCalCurves
  
  # Now interpolate so that everything is on a regular grid in calendar years
  cal_order = order(cal_ages,decreasing=TRUE)
  out = cbind(cal_ages[cal_order],uncal_ages[cal_order],one_sigma[cal_order])
  
  # Now write to an rda file
  fl = paste0(file_loc,'/',name,'.rda')
  save(out, file = fl)
  
}

# Il faut dégrader la fonction BchronCalibrate pour autoriser le passage de valeurs décimales. Il faut aussi corriger la formule

#' Fonction qui permet le calcul de la densité de probabilté d'une date
#' @param mesures mesure ou liste des mesures à calibrer
#' @param std erreur ou liste des erreurs sur les mesures à calibrer
#' @param timeScale pas de la grille de temps
#' @export
calibrate <- function (mesures, std, calCurves, ids = NULL, positions = NULL,   pathToCalCurves = paste0(getwd(),'/','data'),  timeScale = 1  ) 
{
  if (length(mesures) != length(std)) 
    stop("mesures and std should be of same length")
  
  if (length(mesures) != length(calCurves)) 
    stop("mesures and calCurves should be of same length")
  
  if (!is.null(positions)) 
    if (length(mesures) != length(positions)) 
      stop("mesures and positions should be of same length")
  
  if (is.null(ids)) 
    ids = paste("Date", 1:length(mesures), sep = "")
  
  allCalCurves = unique(calCurves)
  calCurve = calTime = calValue = calStd = timeGrid = mu = tau1 = list()
  for (i in 1:length(allCalCurves)) {
    calCurveFile = paste(pathToCalCurves, "/", allCalCurves[i], 
                         ".rda", sep = "")
    if (!file.exists(calCurveFile)) 
      stop(paste("Calibration curve file", calCurveFile, 
                 "not found"))
    x = load(calCurveFile)
    calCurve = get(x)
    calTime[[i]] = calCurve[, 1]
    calValue[[i]] = calCurve[, 2]
    calStd[[i]] = calCurve[, 3]
    
    timeGrid[[i]] = seq(min(calTime[[i]]), max(calTime[[i]]), by = timeScale)
    mu[[i]] = stats::approx(calTime[[i]], calValue[[i]], xout = timeGrid[[i]], rule = 2)$y
    tau1[[i]] = stats::approx(calTime[[i]], calStd[[i]], xout = timeGrid[[i]], rule = 2)$y
    
    if (allCalCurves[i] == "normal") {
      timeRange = range(c(calTime, mesures + 4 * std))
      timeGrid[[i]] = seq(timeRange[1], timeRange[2], by = timeScale)
      mu[[i]] = timeGrid[[i]]
      tau1[[i]] = rep(0, length(timeGrid[[i]]))
    }
  }
  matchCalCurves = match(calCurves, allCalCurves)
  out = list()
  for (i in 1:length(mesures)) {
    if (mesures[i] > max(mu[[matchCalCurves[i]]]) | mesures[i] <  min(mu[[matchCalCurves[i]]])) {
      cal_range = range(mu[[matchCalCurves[i]]])
      stop(paste("Date", ids[i], "outside of calibration range. Range of", calCurves[i], "is", cal_range[1], "to", cal_range[2]))
    }
    
    tau = std[i]^2 + tau1[[matchCalCurves[i]]]^2
    currTimeGrid = timeGrid[[matchCalCurves[i]]]
    # dens = stats::dt( (ages[i] - mu[[matchCalCurves[i]]]) /sqrt(std[i]^2 + tau1[[matchCalCurves[i]]]^2), df = dfs[i])
    dens <- NULL
    for (j in 1:length(timeGrid[[i]]) ) {
      dens[j] <- exp(-0.5 * (mesures[i] - mu[[i]][j])^2/(tau1[[i]][j]^2 + std[i]^2  ) ) /sqrt(tau1[[i]][j]^2 + std[i]^2  )
    }
    
    
    dens = dens/sum(dens)
    out[[i]] = list(mesures = mesures[i], std = std[i], 
                    calCurves = calCurves[i], timeGrid = currTimeGrid, densities = dens,
                    positions = positions[i], timeScale = timeScale)
    
    
  }
  
  names(out) = ids
  class(out) = "RenDate"
  return(out)
}

#' Fonction qui permet le calcul de la densité de probabilté de date uniforme - fonction porte
#' @param gate.min début de la porte
#' @param gate.max fin de la porte
#' @param timeGrid.min valeur minimale de la grille
#' @param timeGrid.max valeur minimale de la grille
#' @param timeScale pas de la grille de temps
#' @param ids nom ou identifiant
#' @export
date.uniform <- function( gate.min = -0, gate.max = 100, time.grid.min = -1000, time.grid.max = 2000, time.grid.scale = 1, ids = NULL, position = NULL)
{
  out = list()
  timeGrid <- seq(time.grid.min, time.grid.max, by = time.grid.scale)
  dens  <- dunif(timeGrid, min = gate.min, max = gate.max)
  
  out[[1]] = list(gate.min = gate.min, gate.max = gate.max, 
                  calCurves = NA, timeGrid = timeGrid, densities = dens,
                  positions = position, timeScale = time.grid.scale)
  
  if (is.null(ids)) 
    ids = paste("U[", gate.min, ";",  gate.max, "]", sep = "")
  names(out) = ids
  class(out) = "RenDate"
  return(out)
}

#' Fonction qui permet le calcul de la densité de probabilté de date gaussienne - fonction gauss
#' @param mean moyenne
#' @param sd standard deviations
#' @param timeGrid.min valeur minimale de la grille
#' @param timeGrid.max valeur minimale de la grille
#' @param timeScale pas de la grille de temps
#' @param ids nom ou identifiant
#' @export
date.gaussian <- function( mean = 0, sd = 10, time.grid.min = -1000, time.grid.max = 2000, time.grid.scale = 1, ids = NULL, position = NULL)
{
  out = list()
  timeGrid <- seq(time.grid.min, time.grid.max, by = time.grid.scale)
  dens  <- dnorm(timeGrid, mean = mean, sd=sd)
  
  out[[1]] = list(calCurves = NA, timeGrid = timeGrid, densities = dens,
                  positions = position, timeScale = time.grid.scale)
  
  if (is.null(ids)) 
    ids = paste("N(", mean, ";",  sd, ")", sep = "")
  names(out) = ids
  class(out) = "RenDate"
  return(out)
}
#' Calcul le hpd (hdr) sur une densité de probabilité de date
#' @param date densité produite par la fonction calibrate, générant un objet de class "RenDate" 
#' @param prob requested surface value [0, 1]
#' @export
setGeneric("hpd", package = "RenDate", valueClass = "list",
 function(date, prob = 0.95) 
{
  
  # A function to return the HPD interval for a date object which should have an timeGrid and a densities argument
  # I was previously using the hdrcde package but this has several unnecessary dependencies
  
  ag <- date$timeGrid
  de <- date$densities
  stp <- date$timeScale

  # Error checking
  if(is.null(ag)) stop('timeGrid not found in date object.')
  if(is.null(de)) stop('densities not found in date object.')
  if(findInterval(prob, c(0, 1))!=1) stop('prob value outside (0,1).')
  
  # Put the probabilities in order
  o = order(de)
  cu = cumsum(de[o])
  
  # Find which ones are above the threshold
  good_cu = which(cu>1-prob)
  good_ag = sort(ag[o][good_cu])
  
  # Pick out the extremes of each range
  breaks = diff(good_ag)>2*stp
  where_breaks = which(diff(good_ag)>2*stp)
  n_breaks = sum(breaks) + 1
  # Store output
  out = vector('list', length = n_breaks)
  low_seq = 1
  high_seq = ifelse(length(where_breaks)==0, length(breaks), where_breaks[1])
  for(i in 1:n_breaks) {
    out[[i]] = c(good_ag[low_seq], good_ag[high_seq])
    #curr_dens = round(100*sum(de[o][seq(good_cu[low_seq], good_cu[high_seq])]), 1) !! FALSE
    i_low <- match(good_ag[low_seq], ag)
    i_high <- match(good_ag[high_seq],ag)
    curr_dens = round(100*sum(de[i_low : i_high]), digits = 1)
    names(out)[[i]] = paste0(as.character(curr_dens),'%')
    low_seq = high_seq + 1
    high_seq = ifelse(i<n_breaks-1, where_breaks[i+1], length(breaks))
  }
  return(out)
  
}
)

#' Trace des courbe de densité
#' @param withHDR Calcul le hdr (hpd) et remplie la surface correspondant
#' @param dateHeigth Fixe la hauteur des densités, quand withPositions est TRUE
#' @param normalize force le maximum à la valeur 1
#' @export
plot.RenDate <- function(x, withPositions = FALSE, pause = FALSE, dateHeight = 30, normalize= FALSE, borderCol = NULL,
                         fillCols = rep('gray', length(x)), withHDR = TRUE, hdrCol = 'darkgray',
                         ...) 
{
  
  # Get extra arguments if provided
  ex = list(...)
  
  if(is.null(ex$xlab)) ex$xlab = 'Date'
  if(is.null(ex$ylab)) ex$ylab = ifelse(withPositions,'Position','Density')
  
  # First plot for individual dates
  if(length(x)==1) {
    if(normalize == TRUE) 
    {
      fac <- max(x[[1]]$densities)
    } else {
      fac <- x[[1]]$timeScale
    }
    ag = x[[1]]$timeGrid
    den = x[[1]]$densities / fac
    ex$x = ag
    ex$y = den
    ex$type = 'l'
    ex$new = new
    if(is.null(ex$main)) ex$main = names(x)
    args = utils::modifyList(ex, list(...))
    do.call("plot", args)
    #graphics::mtext(paste(x[[1]]$calCurves),side=1,line=4,adj=0,cex=0.6)
    if(withHDR) {
      my_hdr = hpd(x[[1]])
      for(j in 1:length(my_hdr)) {
        #x_seq = seq(my_hdr[[j]][1], my_hdr[[j]][2], by = 1)
        #y_lookup = match(x_seq, ag)
        #y_seq = den[y_lookup]
        y_lookup <- match(my_hdr[[j]], ag)
        x_seq <- ag[y_lookup[1]: y_lookup[2]]
        y_seq = den[y_lookup[1]: y_lookup[2]]
        graphics::polygon(c(my_hdr[[j]][1], x_seq, my_hdr[[j]][2]),
                          c(0, y_seq, 0),
                          col = hdrCol,
                          border = NA)
      }
    }
    
  }
  
  # Now for multiple dates without depths
  if(length(x)>1 & withPositions==FALSE) {
    for(i in 1:length(x)) {
      if(normalize == TRUE) 
      {
        fac <- max(x[[i]]$densities)
      } else {
        fac <- x[[i]]$timeScale
      }
      ex_curr = ex
      ag = x[[i]]$timeGrid
      den = x[[i]]$densities / fac
      ex_curr$x = ag
      ex_curr$y = den
      ex_curr$type = 'l'
      ex_curr$new = new
      if(is.null(ex_curr$main)) ex_curr$main = names(x)[i]
      args = utils::modifyList(ex_curr, list(...))
      do.call("plot", args)
      # graphics::mtext(paste(x[[i]]$calCurves),side=1,line=4,adj=0,cex=0.6)
      if(withHDR) {
        my_hdr = hpd(x[[i]])
        for(j in 1:length(my_hdr)) {
          #x_seq = seq(my_hdr[[j]][1], my_hdr[[j]][2], by = 1)
          #y_lookup = match(x_seq, ag)
          #y_seq = den[y_lookup]
          
          y_lookup <- match(my_hdr[[j]], ag)
          x_seq <- ag[y_lookup[1]: y_lookup[2]]
          y_seq = den[y_lookup[1]: y_lookup[2]]
          graphics::polygon(c(my_hdr[[j]][1], x_seq, my_hdr[[j]][2]),
                            c(0, y_seq, 0),
                            col = hdrCol,
                            border = NA)
        }
      }
      if(pause) if(i<length(x)) readline('Hit Enter for next plot...')
    }
  }
  
  # Finally for multiple dates with depths
  if(length(x)>1 & withPositions==TRUE) {
    xlimits = NULL
    ylimits = NULL
    for(i in 1:length(x)) {
      xlimits = range(c(xlimits,x[[i]]$timeGrid))
      ylimits = range(c(ylimits,x[[i]]$positions))
    }
    #dateHeight=0.2*diff(pretty(ylimits))[1]
    ylimits[1] = ylimits[1]-dateHeight
    if(is.null(ex$xlim)) ex$xlim = rev(xlimits)
    if(is.null(ex$ylim)) ex$ylim = rev(ylimits)
    ex$x = ex$y = 1
    ex$type = 'n'
    ex$new = new
    if(is.null(ex$main)) ex$main = 'Calibrated dates by position'
    args = utils::modifyList(ex, list(...))
    do.call("plot", args)
    
    for(i in 1:length(x)) {
      curr_xlimit = range(x[[i]]$timeGrid)
      curr_pos = x[[i]]$positions
      graphics::polygon(c(curr_xlimit[1], x[[i]]$timeGrid, curr_xlimit[2]),
                        c(curr_pos, x[[i]]$positions-x[[i]]$densities*dateHeight/max(x[[i]]$densities), curr_pos),
                        border=ifelse(is.null(borderCol), NA, borderCol),
                        col=fillCols[i])
    }
  }
  
}

#' Trace des courbes de densité avec leurs enveloppes d erreur
#' @param withHDR Calcul le hdr (hdp) et remplie la surface correspondant
#' @export
lines.RenDate <- function(x, withPositions=FALSE, pause=FALSE, dateHeight = 30, normalize = FALSE, borderCol = NULL,
                          fillCols = rep('gray', length(x)), withHDR = TRUE, hdrCol = 'darkgray', 
                          ...) 
{
  
  # Get extra arguments if provided
  ex = list(...)#as.list(substitute(list(...)))[-1L]
  
  if(is.null(ex$xlab)) ex$xlab = 'Date'
  if(is.null(ex$ylab)) ex$ylab = ifelse(withPositions,'Position','Density')
  
  # First plot for individual dates
  if(length(x)==1) {
    if(normalize == TRUE) 
    {
      fac <- max(x[[1]]$densities)
    } else {
      fac <- x[[1]]$timeScale
    }
    ag = x[[1]]$timeGrid
    den = x[[1]]$densities / fac
    ex$x = ag
    ex$y = den
    ex$type = 'l'
    if(is.null(ex$main)) ex$main = names(x)
    args = utils::modifyList(ex, list(...))
    do.call("lines", args)
    #graphics::mtext(paste(x[[1]]$calCurves),side=1,line=4,adj=0,cex=0.6)
    if(withHDR) {
      my_hdr = hpd(x[[1]])
      for(j in 1:length(my_hdr)) {
        #x_seq = seq(my_hdr[[j]][1], my_hdr[[j]][2], by = 1)
        #y_lookup = match(x_seq, ag)
        #y_seq = den[y_lookup]
        y_lookup <- match(my_hdr[[j]], ag)
        x_seq <- ag[y_lookup[1]: y_lookup[2]]
        y_seq = den[y_lookup[1]: y_lookup[2]]
        graphics::polygon(c(my_hdr[[j]][1], x_seq, my_hdr[[j]][2]),
                          c(0, y_seq, 0),
                          col = hdrCol,
                          border = NA)
      }
    }
    
  }
  
  # Now for multiple dates without depths
  if(length(x)>1 & withPositions==FALSE) {
    for(i in 1:length(x)) {
      ex_curr = ex
      if(normalize == TRUE) 
      {
        fac <- max(x[[i]]$densities)
      } else {
        fac <- x[[i]]$timeScale
      }
      ag = x[[i]]$timeGrid
      den = x[[i]]$densities / fac
      ex_curr$x = ag
      ex_curr$y = den
      ex_curr$type = 'l'
      if(is.null(ex_curr$main)) ex_curr$main = names(x)[i]
      args = utils::modifyList(ex_curr, list(...))
      do.call("lines", args)
      # graphics::mtext(paste(x[[i]]$calCurves),side=1,line=4,adj=0,cex=0.6)
      if(withHDR) {
        my_hdr = hpd(x[[i]])
        for(j in 1:length(my_hdr)) {
          #x_seq = seq(my_hdr[[j]][1], my_hdr[[j]][2], by = 1)
          #y_lookup = match(x_seq, ag)
          #y_seq = den[y_lookup]
          y_lookup <- match(my_hdr[[j]], ag)
          x_seq <- ag[y_lookup[1]: y_lookup[2]]
          y_seq = den[y_lookup[1]: y_lookup[2]]
          graphics::polygon(c(my_hdr[[j]][1], x_seq, my_hdr[[j]][2]),
                            c(0, y_seq, 0),
                            col = hdrCol,
                            border = NA)
        }
      }
      if(pause) if(i<length(x)) readline('Hit Enter for next plot...')
    }
  }
  
  # Finally for multiple dates with depths
  if(length(x)>1 & withPositions==TRUE) {
    xlimits = NULL
    ylimits = NULL
    for(i in 1:length(x)) {
      xlimits = range(c(xlimits,x[[i]]$timeGrid))
      ylimits = range(c(ylimits,x[[i]]$positions))
    }
    #dateHeight=0.2*diff(pretty(ylimits))[1]
    ylimits[1] = ylimits[1]-dateHeight
    if(is.null(ex$xlim)) ex$xlim = rev(xlimits)
    if(is.null(ex$ylim)) ex$ylim = rev(ylimits)
    ex$x = ex$y = 1
    ex$type = 'n'
    if(is.null(ex$main)) ex$main = 'Calibrated dates by position'
    args = utils::modifyList(ex, list(...))
    do.call("lines", args)
    
    for(i in 1:length(x)) {
      curr_xlimit = range(x[[i]]$ageGrid)
      curr_pos = x[[i]]$positions
      graphics::polygon(c(curr_xlimit[1], x[[i]]$timeGrid, curr_xlimit[2]),
                        c(curr_pos, x[[i]]$positions-x[[i]]$densities*dateHeight/max(x[[i]]$densities), curr_pos),
                        border=ifelse(is.null(borderCol),NA,borderCol),
                        col=fillCols[i])
    }
  }
  
}

#' Trace une courbe avec son enveloppe d erreur à 1 sigma et deux sigma
#' @export
courbe.enveloppe <- function(t, mean, std, col.env = "forestgreen",  xlim = NULL, ylim = NULL, new = TRUE,...)
{
  enveloppe1 <- mean + outer(std , c(1, -1))
  
  enveloppe2 <- mean + outer(std , c(2.54, -2.54)) 
  
  if (is.null(ylim))
    ylim <- range(enveloppe2)
  
  if(new == TRUE) {
    plot(
      x=t, y=mean, type="l", ylim = ylim,
      panel.first=polygon(c(t, rev(t)), c(enveloppe2[,1], rev(enveloppe2[,2])), border=NA, col=adjustcolor( col.env, alpha.f = 0.2)), xlim=xlim,  ...=... 
    )
  } else {
    lines(x = t , y=mean,
          panel.first=polygon(c(t, rev(t)), c(enveloppe2[,1], rev(enveloppe2[,2])), border=NA, col=adjustcolor( col.env, alpha.f = 0.2))
    )
  }
  
  
  lines(x = t , y=mean,
        panel.first=polygon(c(t, rev(t)), c(enveloppe1[,1], rev(enveloppe1[,2])), border=NA, col=adjustcolor( col.env, alpha.f = 0.3))
  )
  
}
#' Trace une droite représentant une mesure avec son enveloppe d erreur à 1 sigma et deux sigma
#' @export
mesure.enveloppe <- function(t, mesure, std, col.env = "gray",  col.mesure = "darkgray", ...)
{
  abline(h= mesure, col = col.mesure )
 
 # text(t[1], mesure, labels=as.character(mesure) )
  # enveloppe sur la mesure
  polygon( c(t[1], t[1], t[length(t)], t[length(t)] ), c(mesure + 2.54*std, mesure - 2.54*std,  mesure - 2.54*std,  mesure + 2.54*std), border=NA, col=adjustcolor(col.env, alpha.f = 0.3))
  polygon( c(t[1], t[1], t[length(t)], t[length(t)] ), c(mesure + std, mesure - std,  mesure - std,  mesure + std), border=NA, col=adjustcolor( col.env, alpha.f = 0.3))
  
  
}
#' Calcul la combinaison au sens produit de deux densités de class "RenDate"
#' @param timeScale permet de modifier la grille de temps
#' @export
produit.RenDate <- function(date1, date2, timeScale = 1)
{
  cal_range_i = range(date1[[1]]$timeGrid)
  cal_range_d = range(date2[[1]]$timeGrid)
  cal_range = c(max(cal_range_d[1], cal_range_i[1]), min(cal_range_d[2], cal_range_i[2]))
  
  
  currentTimeGridOut = seq(cal_range[1], cal_range[2], by = timeScale)
  Inc1 = stats::approx(date1[[1]]$timeGrid, date1[[1]]$densities, xout = currentTimeGridOut, rule = 2)$y
  Dec1 = stats::approx(date2[[1]]$timeGrid, date2[[1]]$densities, xout = currentTimeGridOut, rule = 2)$y
  if (sum(Inc1)/timeScale != 1) 
    paste("surface de date1 !=1", sum(Inc1)/timeScale , sep=" ")
  
  if (sum(Dec1)/timeScale != 1) 
    paste("surface de date2 !=1", sum(Dec1)/timeScale, sep =" " )
  
  densitiesOut = NULL
  ageGridOut = NULL
  out = list()
  
  densitiesOut <- Dec1 * Inc1
  
  densitiesOut <- densitiesOut/sum(densitiesOut)
  out[[1]] = list(  timeGrid = currentTimeGridOut, densities = densitiesOut, timeScale = timeScale)
  names(out) <- "Combinaison"
  class(out) <- "RenDate"
  return(out)
}

# en utilisant, la fonction convole de R
# les indices imin et imax sont relatifs, c.à.d qu'ils peuvent être négatifs ou nulles
wiggle.indice <- function(f, imin, imax)
{
  if (imax < imin) {
    message("imax must be greater than imin")
    return(NULL)
  }
  
  wiggle.span <- imax - imin + 1
  wiggle.span.n <- abs(wiggle.span)
  wiggle.span.gate <- rep(1/wiggle.span.n, wiggle.span.n)
  
  # condition pour convole

  if (imin >= 0) {
    ## 1 - tmin et tmax >0
    wiggle <- c(rep(0, trunc(abs(imin)) ), wiggle.span.gate)
    f0 <- f #c(f, rep(0, wiggle.span.n) )
    # Je garde le même f0= f
    wiggle.conj <- FALSE
    
  } else {
    ## 2 - tmin < 0 et tmax >0
    ## 3 - tmin et tmax < 0
    if (imax >= 0) {
      wiggle <- wiggle.span.gate
      f0 <- c(f[-c(1: trunc(abs(imin)))], rep(0, trunc(abs(imin))) )
      wiggle.conj <- FALSE
      
    } else {
      f0 <- f
      wiggle <- c(rep(0, trunc(abs(imax)) ), wiggle.span.gate)
      wiggle.conj <- TRUE # je dois inverser-conjuguer la fft
    }
    
  }
  # La longueur du produit de convolution est égale à M+N-1
  # Il faut donc réduire le tableau pour retrouver la taille originale
  conv <- convolve(f0, wiggle, conj= wiggle.conj, type = "open")[length(wiggle):(length(f0)+length(wiggle)-1)]
  # Normalisation
  conv <- conv/sum(conv)
  
  return(conv)
}

#' Calcul le wiggle (décalage) d'une densité de class "RenDate"
#' Fonction utilisée avec produit.RenDate() pour calculer le "Wiggle Matching"
#' Il s'agit d'un poroduit de convolution de la datation par une "fonction porte"
#' @param x une densité (datation) de class "RenDate"
#' @param wiggle.min la borne inférieure du décalage
#' @param wiggle.max la borne supérieure du décalage
#' @export
wiggle.uniform <- function(x, wiggle.min, wiggle.max)
{
  out <- x
  for (i in length(x)) {
    timeScale = x[[i]]$timeScale
    
    wiggle.range <- c(wiggle.min/timeScale, wiggle.max/timeScale)
    wiggle.imin <- floor(min(wiggle.range))
    wiggle.imax <- floor(max(wiggle.range))
    
    out[[i]]$densities <- wiggle.indice(x[[i]]$densities, wiggle.imin, wiggle.imax)
  }
  
  names(out) = paste("wiggle Uniforme [", wiggle.min, ";",  wiggle.max, "]", sep = "")
  return(out)
}

# Calcul de wiggle Matching en utilisant, la fonction convole() de R
# fonction utilisée dans le wiggle.gauss
wiggle.gauss.dens <- function(f, t.mean, t.sd, f.scale)
{
  f.len <- length(f)
  if (t.sd <= 0) {
    message("t.sd must be positive")
    return(NULL)
  }
  # génération de la gaussien, centrée sur l'intervalle
  xmin <- -trunc(f.len/2)
  xmax <- f.len + xmin
  x <- seq(xmin/f.scale, xmax/f.scale, length.out = f.len)
  
  wiggle.gauss.f <- dnorm(x, t.mean, t.sd)
  
  if (max(wiggle.gauss.f) == 0) {
    warning("Time scale too hight, no possible convolution")
  }
  
  # padding
  f.padding = c(f, rep(0, f.len ) )
  
  # La longueur du produit de convolution est égale à M+N-1
  # Il faut donc réduire le tableau pour retrouver la taille originale
  conv <- convolve(f.padding, wiggle.gauss.f, conj= FALSE, type = "open")[(f.len*1.5 ):( f.len*2.5 - 1)]
  # Normalisation
  conv <- conv/sum(conv)
  return(conv)
}

#' Calcul le wiggle (décalage) d'une densité de class "RenDate"
#' Fonction utilisée avec produit.RenDate() pour calculer le "Wiggle Matching"
#' Il s'agit d'un produit de convolution de la datation par une "densité Normale"
#' @param x une densité (datation) de class "RenDate"
#' @param mean la valeur moyenne
#' @param sd la valeur de l'écart type
#' @export
wiggle.gauss <- function(x, mean, sd)
{
  out <- x
  for (i in length(x)) {
    out[[i]]$densities <- wiggle.gauss.dens(x[[i]]$densities, mean, sd, x[[i]]$timeScale)
  }
  
  names(out) = paste("wiggle gauss [", mean, ";",  sd, "]", sep = "")
  return(out)
}


#' Fonction utilisée pour réduire l'intervalle support d'une densité de class "RenDate"
#' @param dens une densité (datation) de class "RenDate"
#' @param tmin la valeur moyenne
#' @param tmax la valeur de l'écart type
#' @export
period.reduction <- function(dens, tmin=0, tmax=2000)
{
  if (tmax<= tmin) {
    warning("tmax<= tmin")
  }
  
  # boucle sur liste
  res <- NULL
  for(li in length(dens)) {
    
    # recherche du imin
    imatch <- match(tmax, dens[[li]]$timeGrid, nomatch = -1);
    if ( imatch > -1 ) {
      imin <-  which(dens[[li]]$timeGrid == tmin) 
      
    } else if (tmin < dens[[li]]$timeGrid[1] ) {
      imin <- 1
      
    } else {
      for(i in seq(from=length(dens[[li]]$timeGrid), to=1, by=-1))
        if (tmin < dens[[li]]$timeGrid[i]) {
          imin <- i
        }
      
    }
    
    # recherche du imax
    imatch <- match(tmax, dens[[li]]$timeGrid, nomatch = -1);
    if ( imatch > -1 ) {
      imax <-  which(dens[[li]]$timeGrid == tmax) 
      
    } else if (tmax > max(dens[[li]]$timeGrid) ) {
      imax <- length(dens[[li]]$timeGrid) 
      
    } else {
      for(i in seq(from=1, to=length(dens[[li]]$timeGrid), by=1))
        if (tmax > dens[[li]]$timeGrid[i]) {
          imax <- i
        }
    }
    if (imax<= imin) {
      warning("imax<= imin")
    }
    res[[li]]$timeScale <- dens[[li]]$timeScale
    res[[li]]$timeGrid <- dens[[li]]$timeGrid[imin:imax]
    res[[li]]$densities <- dens[[li]]$densities[imin:imax]
    # Normalisation de la surface
    res[[li]]$densities <-res[[li]]$densities /sum(res[[li]]$densities)
    
  }
  names(res) = names(dens)
  class(res) = "RenDate"
  return(res)
}


# si on a l'erreur : la chaîne de caractères entrée 1 est incorrecte dans cet environnement linguistique
# cela correspond au mauvais encoding
#' Lecture d'un fichier de référence
#' @param encoding  Pour les fichiers du MacOS, il faut "macroman" -> difficle à connaitre, peut être "latin1" ou "utf8".
#' @return une liste de data.frame
#' @example 'curveI <- read.Ref("Calib/AM/GAL2002sph2014_I.ref")'
#' @export
read.Ref <- function(file.Ref, encoding = "macroman")
{
  # Lecture et extraction des points de la courbe
  # file.Ref="/Users/dufresne/Documents/AM/Croatie_Loron_Tar-Vabriga_2018/GAL2002sph2014_D.ref"
  lin<- NULL
  fil <- file(file.Ref, "r", encoding = "macroman") #, encoding = encoding)
  lin <- readLines(fil)
  close(fil)
  # Recherche position-ligne des points de référence
  g <- length(lin)
  # Compte le nombre de mesures
  n.measures <-0
  for (i in 1:length(lin)) {
    if (as.logical(length(grep("#", lin[i])))==FALSE  && as.logical(length(grep("/", lin[i])) )==FALSE )
      n.measures <- n.measures + 1
    
    if (length(grep("# reference points", lin[i])) || length(grep("# Reference points", lin[i])) ){
      g<- i
      break
    }
  }
  
  listtmp <- NULL
  listtmp<- read.table(file.Ref, dec='.', sep=",", header=FALSE, 
                          skip = i-n.measures-1, comment.char = "#", 
                          nrows = n.measures,
                          stringsAsFactors=FALSE) 
  colnames(listtmp)<-c("t", "value", "sigma")
  list <- NULL
  list$curve <- NULL
  list$curve <- listtmp

  listtmp <- NULL
  if (g<length(lin) ) {
    listtmp<- read.table(file.Ref, dec='.', 
                              sep=",", header=FALSE, 
                              skip = i+1, 
                              comment.char = "#",
                              stringsAsFactors=FALSE)
    list$pts.ref <- NULL
    list$pts.ref$tij1 <- listtmp[1]
    list$pts.ref$tij2 <- listtmp[2]
    list$pts.ref$tm <- listtmp[3]
    list$pts.ref$Yij <- listtmp[4]
    list$pts.ref$Err_Yijc <- listtmp[5]

  }
  return(list)
}

#' Conversion une trace en datation pour l'utiliser avec les autres fonction du package
#' @param trace  un vecteur history plot provenant de ChronoModel
#' @return une date par lissage avec un noyaux, gaussien par defaut
#' @export
trace.to.date <- function( trace,  bw = "nrd", adjust = 1, from, to, gridLength = 1024,
                           kernel = c("gaussian"),
                           ids = NULL, position = NULL)
{
  out = list()
  if (bw=="CM") {
    sigma <- sd(trace) # sans biais
    bandwidth<-1.06
    bw <- bandwidth * sigma * length(trace)^ (-1./5.) # in CM = h
    from <- min(trace) - 4. * h # in CM = a
    to <- max(trace) + 4. * h # in CM = b
  } 
  dens  <- density(x=trace, bw = bw, adjust = adjust, kernel = kernel, from = from, to = to, n = gridLength)
  dens$y <- dens$y/sum(dens$y)
  out[[1]] = list(n = dens$n, bw = dens$bw, 
                  calCurves = NA, timeGrid = dens$x, densities = dens$y,
                  positions = position, timeScale = dens$x[2] - dens$x[1])
  
  if (is.null(ids)) 
    ids = paste("trace : ", dens$data.name, sep = "")
  names(out) = ids
  class(out) = "RenDate"
  return(out)
}
