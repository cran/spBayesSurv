
.DDPplots <- function ( xpred_, ygrid_, beta_, sigma2_, w_, probs_ ) {
  .Call("DDPplots", xpred_, ygrid_, beta_, sigma2_, w_, probs_, PACKAGE = "spBayesSurv")
}

.DistMat <- function ( si_, sj_ ) {
  .Call("DistMat", si_, sj_, PACKAGE = "spBayesSurv")
}

.CoxPHplots <- function ( xpred_, tgrid_, beta_, h_, d_, probs_ ) {
  .Call("CoxPHplots", xpred_, tgrid_, beta_, h_, d_, probs_, PACKAGE = "spBayesSurv")
}

.frailtyGAFTplots <- function(tgrid_, xcepred_, xtfpred_, betace_, betatf_, v_, sigma2_, maxL_, probs_){
  .Call("frailtyGAFTplots", tgrid_, xcepred_, xtfpred_, betace_, betatf_, v_, sigma2_, maxL_, probs_, PACKAGE = "spBayesSurv")
}

.BayesFactor <- function(betatf_, maxL_, a0_, b0_, gprior_, alpha_){
  .Call("BayesFactor", betatf_, maxL_, a0_, b0_, gprior_, alpha_, PACKAGE = "spBayesSurv")
}