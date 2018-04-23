"GetCurves" <- function (x, xnewdata, xtfnewdata, tgrid = NULL, ygrid = NULL, 
                         frail = NULL, CI = 0.95, PLOT = TRUE, ...) {
  if(is(x,"anovaDDP")){
    estimates = plot.anovaDDP(x, xnewdata, tgrid=tgrid, CI=CI, PLOT=PLOT, ...)
  }else if(is(x,"spCopulaDDP")){
    estimates = plot.spCopulaDDP(x, xnewdata, tgrid=tgrid, CI=CI, PLOT=PLOT, ...)
  }else if(is(x,"indeptCoxph")){
    estimates = plot.indeptCoxph(x, xnewdata, tgrid=tgrid, CI=CI, PLOT=PLOT, ...)
  }else if(is(x,"spCopulaCoxph")){
    estimates = plot.spCopulaCoxph(x, xnewdata, tgrid=tgrid, CI=CI, PLOT=PLOT, ...)
  }else if(is(x, "survregbayes")){
    estimates = plot.survregbayes(x, xnewdata, tgrid=tgrid, frail=frail, CI=CI, PLOT=PLOT, ...)
  }else if(is(x, "frailtyGAFT")){
    estimates = plot.frailtyGAFT(x, xnewdata, xtfnewdata, tgrid=tgrid, frail=frail, CI=CI, PLOT=PLOT, ...)
  }else if(is(x, "SuperSurvRegBayes")){
    estimates = plot.SuperSurvRegBayes(x, xnewdata, tgrid=tgrid, CI=CI, PLOT=PLOT, ...)
  }else if(is(x, "SpatDensReg")){
    estimates = plot.SpatDensReg(x, xnewdata, ygrid=ygrid, CI=CI, PLOT=PLOT, ...)
  }
  invisible(estimates)
}