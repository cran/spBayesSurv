#ifndef _spSurv_BP_tools_h
#define _spSurv_BP_tools_h

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "spSurv_common.h"
#include "spSurv_spatialtools.h"

// from conditional Ys to cumulative probs
void Ys_to_weight(const Rcpp::NumericVector& Ys, Rcpp::NumericVector& weight);

/////////////////////////////////////////////////////////////////////////
/////////////////// baseline suvival functions///////////////////////////
/////////////////////////////////////////////////////////////////////////
double S0BP(double t, double th1, double th2, Rcpp::NumericVector w, bool BP, int dist);
double F0BP(double t, double th1, double th2, Rcpp::NumericVector w, bool BP, int dist);
double logf0BP(double t, double th1, double th2, Rcpp::NumericVector w, bool BP, int dist);

/////////////////////////////////////////////////////////////////////////
//////////////////////////// AFT model //////////////////////////////////
/////////////////////////////////////////////////////////////////////////
// log densitt of t given xi
double AFT_BP_logpdf(double t, double th1, double th2, Rcpp::NumericVector w,
                     bool BP, int dist, double xibeta);
// log survival function of t given xi
double AFT_BP_logsurv(double t, double th1, double th2, Rcpp::NumericVector w, 
                      bool BP, int dist, double xibeta);
// log cdf of t given xi
double AFT_BP_logcdf(double t, double th1, double th2, Rcpp::NumericVector w, 
                     bool BP, int dist, double xibeta);
// log (S(t1|xi)-S(t2|xi))
double AFT_BP_logsurvdiff(double t1, double t2, double th1, double th2, Rcpp::NumericVector w, 
                          bool BP, int dist, double xibeta);
// log likelihood given data, frailties and parameters 
void AFT_BP_loglik(const Rcpp::NumericVector& t1, const Rcpp::NumericVector& t2, const Rcpp::NumericVector& ltr,
                   const Rcpp::IntegerVector& type, double th1, double th2, Rcpp::NumericVector w,
                   bool BP, int dist, const Rcpp::NumericVector& Xbeta, double& ll);
// log likelihood given frailties, parameters and data of block i
void AFT_BP_loglikblocki(const Rcpp::NumericVector& t1, const Rcpp::NumericVector& t2, const Rcpp::NumericVector& ltr, 
                         const Rcpp::IntegerVector& type, double th1, double th2, Rcpp::NumericVector w,
                         bool BP, int dist, const Rcpp::NumericVector& Xbeta, double& ll,
                         int ind1, int ind2, double vi);
// Calculate 1.0/likelihood for CPO
arma::vec AFT_BP_invLik(const Rcpp::NumericVector& t1, const Rcpp::NumericVector& t2, const Rcpp::NumericVector& ltr,
                        const Rcpp::IntegerVector& type, double th1, double th2, Rcpp::NumericVector w,
                        bool BP, int dist, const Rcpp::NumericVector& Xbeta);

/////////////////////////////////////////////////////////////////////////
//////////////////////////// PH model //////////////////////////////////
/////////////////////////////////////////////////////////////////////////
// log densitt of t given xi
double PH_BP_logpdf(double t, double th1, double th2, Rcpp::NumericVector w,
                    bool BP, int dist, double xibeta);
// log survival function of t given xi
double PH_BP_logsurv(double t, double th1, double th2, Rcpp::NumericVector w, 
                     bool BP, int dist, double xibeta);
// log cdf of t given xi
double PH_BP_logcdf(double t, double th1, double th2, Rcpp::NumericVector w, 
                    bool BP, int dist, double xibeta);
// log (S(t1|xi)-S(t2|xi))
double PH_BP_logsurvdiff(double t1, double t2, double th1, double th2, Rcpp::NumericVector w, 
                         bool BP, int dist, double xibeta);
// log likelihood given data, frailties and parameters 
void PH_BP_loglik(const Rcpp::NumericVector& t1, const Rcpp::NumericVector& t2, const Rcpp::NumericVector& ltr,
                  const Rcpp::IntegerVector& type, double th1, double th2, Rcpp::NumericVector w,
                  bool BP, int dist, const Rcpp::NumericVector& Xbeta, double& ll);
// log likelihood given frailties, parameters and data of block i
void PH_BP_loglikblocki(const Rcpp::NumericVector& t1, const Rcpp::NumericVector& t2, const Rcpp::NumericVector& ltr, 
                        const Rcpp::IntegerVector& type, double th1, double th2, Rcpp::NumericVector w,
                        bool BP, int dist, const Rcpp::NumericVector& Xbeta, double& ll,
                        int ind1, int ind2, double vi);
// Calculate 1.0/likelihood for CPO
arma::vec PH_BP_invLik(const Rcpp::NumericVector& t1, const Rcpp::NumericVector& t2, const Rcpp::NumericVector& ltr,
                       const Rcpp::IntegerVector& type, double th1, double th2, Rcpp::NumericVector w,
                       bool BP, int dist, const Rcpp::NumericVector& Xbeta);

/////////////////////////////////////////////////////////////////////////
//////////////////////////// PO model //////////////////////////////////
/////////////////////////////////////////////////////////////////////////
// log densitt of t given xi
double PO_BP_logpdf(double t, double th1, double th2, Rcpp::NumericVector w,
                    bool BP, int dist, double xibeta);
// log survival function of t given xi
double PO_BP_logsurv(double t, double th1, double th2, Rcpp::NumericVector w, 
                     bool BP, int dist, double xibeta);
// log cdf of t given xi
double PO_BP_logcdf(double t, double th1, double th2, Rcpp::NumericVector w, 
                    bool BP, int dist, double xibeta);
// log (S(t1|xi)-S(t2|xi))
double PO_BP_logsurvdiff(double t1, double t2, double th1, double th2, Rcpp::NumericVector w, 
                         bool BP, int dist, double xibeta);
// log likelihood given data, frailties and parameters 
void PO_BP_loglik(const Rcpp::NumericVector& t1, const Rcpp::NumericVector& t2, const Rcpp::NumericVector& ltr,
                  const Rcpp::IntegerVector& type, double th1, double th2, Rcpp::NumericVector w,
                  bool BP, int dist, const Rcpp::NumericVector& Xbeta, double& ll);
// log likelihood given frailties, parameters and data of block i
void PO_BP_loglikblocki(const Rcpp::NumericVector& t1, const Rcpp::NumericVector& t2, const Rcpp::NumericVector& ltr, 
                        const Rcpp::IntegerVector& type, double th1, double th2, Rcpp::NumericVector w,
                        bool BP, int dist, const Rcpp::NumericVector& Xbeta, double& ll,
                        int ind1, int ind2, double vi);
// Calculate 1.0/likelihood for CPO
arma::vec PO_BP_invLik(const Rcpp::NumericVector& t1, const Rcpp::NumericVector& t2, const Rcpp::NumericVector& ltr,
                       const Rcpp::IntegerVector& type, double th1, double th2, Rcpp::NumericVector w,
                       bool BP, int dist, const Rcpp::NumericVector& Xbeta);

/////////////////////////////////////////////////////////////////////////
//////////////////////////// AH model //////////////////////////////////
/////////////////////////////////////////////////////////////////////////
// log densitt of t given xi
double AH_BP_logpdf(double t, double th1, double th2, Rcpp::NumericVector w,
                    bool BP, int dist, double xibeta);
// log survival function of t given xi
double AH_BP_logsurv(double t, double th1, double th2, Rcpp::NumericVector w, 
                     bool BP, int dist, double xibeta);
// log cdf of t given xi
double AH_BP_logcdf(double t, double th1, double th2, Rcpp::NumericVector w, 
                    bool BP, int dist, double xibeta);
// log (S(t1|xi)-S(t2|xi))
double AH_BP_logsurvdiff(double t1, double t2, double th1, double th2, Rcpp::NumericVector w, 
                         bool BP, int dist, double xibeta);
// log likelihood given data, frailties and parameters 
void AH_BP_loglik(const Rcpp::NumericVector& t1, const Rcpp::NumericVector& t2, const Rcpp::NumericVector& ltr,
                  const Rcpp::IntegerVector& type, double th1, double th2, Rcpp::NumericVector w,
                  bool BP, int dist, const Rcpp::NumericVector& Xbeta, double& ll);
// log likelihood given frailties, parameters and data of block i
void AH_BP_loglikblocki(const Rcpp::NumericVector& t1, const Rcpp::NumericVector& t2, const Rcpp::NumericVector& ltr, 
                        const Rcpp::IntegerVector& type, double th1, double th2, Rcpp::NumericVector w,
                        bool BP, int dist, const Rcpp::NumericVector& Xbeta, double& ll,
                        int ind1, int ind2, double vi);
// Calculate 1.0/likelihood for CPO
arma::vec AH_BP_invLik(const Rcpp::NumericVector& t1, const Rcpp::NumericVector& t2, const Rcpp::NumericVector& ltr,
                       const Rcpp::IntegerVector& type, double th1, double th2, Rcpp::NumericVector w,
                       bool BP, int dist, const Rcpp::NumericVector& Xbeta);

/////////////////////////////////////////////////////////////////////////
/////////////////////// Super model: PH_PO_AFT //////////////////////////
/////////////////////////////////////////////////////////////////////////
// log density of t given xi
double PHPOAFT_BP_logpdf(double t, double th1, double th2, Rcpp::NumericVector w,
                         bool BP, int dist, double xibeta_h, double xibeta_o, double xibeta_q);
// log survival function of t given xi
double PHPOAFT_BP_logsurv(double t, double th1, double th2, Rcpp::NumericVector w, 
                          bool BP, int dist, double xibeta_h, double xibeta_o, double xibeta_q);
// log cdf of t given xi
double PHPOAFT_BP_logcdf(double t, double th1, double th2, Rcpp::NumericVector w, 
                         bool BP, int dist, double xibeta_h, double xibeta_o, double xibeta_q);
// log (S(t1|xi)-S(t2|xi)) 
double PHPOAFT_BP_logsurvdiff(double t1, double t2, double th1, double th2, Rcpp::NumericVector w, 
                              bool BP, int dist, double xibeta_h, double xibeta_o, double xibeta_q);
// log likelihood given data and parameters 
void PHPOAFT_BP_loglik(const Rcpp::NumericVector& t1, const Rcpp::NumericVector& t2, const Rcpp::NumericVector& ltr, 
                       const Rcpp::IntegerVector& type, double th1, double th2, Rcpp::NumericVector w,  
                       bool BP, int dist, const Rcpp::NumericVector& Xbeta_h, const Rcpp::NumericVector& Xbeta_o,
                       const Rcpp::NumericVector& Xbeta_q, double& ll);
// Calculate 1.0/likelihood for CPO
arma::vec PHPOAFT_BP_invLik(const Rcpp::NumericVector& t1, const Rcpp::NumericVector& t2, const Rcpp::NumericVector& ltr, 
                            const Rcpp::IntegerVector& type, double th1, double th2, Rcpp::NumericVector w, 
                            bool BP, int dist, const Rcpp::NumericVector& Xbeta_h, const Rcpp::NumericVector& Xbeta_o, 
                            const Rcpp::NumericVector& Xbeta_q);

#endif
