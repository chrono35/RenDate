#  Licence ----
#
#  Copyright or © or Copr. CNRS	2019
#
#
# This package is under
# The “Creative Commons Attribution-ShareAlike International License” version 4.0

#  A copy of the License is available at
#  http://www.r-project.org/Licenses/
#  _________________________________________________________________________________

# Version 2020-02-07
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
  
  # And we're done
  cat('Completed!\n')
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

#' Calcul le hpd (hdr) sur une densité de probabilité de date
#' @param date densité produite par la fonction calibrate, génértant un objet de class "RenDate" 
#' @param prob requested surface value [0, 1]
#' @export
setGeneric("hpd", package = "RenDate", valueClass = "RenDate",
 function(date, prob = 0.95) 
{
  
  # A function to return the HPD interval for a date object which should have an timeGrid and a densities argument
  # I was previously using the hdrcde package but this has several unnecessary dependencies
  
  ag = date$timeGrid
  de = date$densities

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
#' Trace une roite représentant une mesure avec son enveloppe d erreur à 1 sigma et deux sigma
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