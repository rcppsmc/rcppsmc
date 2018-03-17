// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
//
// nonLinPMMH.cpp: Example 3.1 of Andrieu et al. (2010). Implementing particle marginal
// Metropolis-Hastings for a toy non-linear state space model previously described in
// Gordon et al. (1993) and Kitagawa (1996).
//
// Copyright (C) 2017         Dirk Eddelbuettel, Adam Johansen and Leah South
//
// This file is part of RcppSMC.
//
// RcppSMC is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// RcppSMC is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with RcppSMC.  If not, see <http://www.gnu.org/licenses/>.

#include "nonLinPMMH.h"

namespace nonLinPMMH {
    const double a_prior = 0.01;
    const double b_prior = 0.01;
}

using namespace std;
using namespace nonLinPMMH;


// nonLinPMMH_impl() function callable from R via Rcpp::
// [[Rcpp::export]]
Rcpp::DataFrame nonLinPMMH_impl(arma::vec data, unsigned long lNumber, unsigned long lMCMCits) {

    try {
        arma::vec sigv(lMCMCits+1), sigw(lMCMCits+1);
        long lIterates = data.n_rows;
        y = data;

        arma::vec loglike = arma::zeros(lMCMCits+1);
        arma::vec logprior = arma::zeros(lMCMCits+1);

        double loglike_prop;
        double logprior_prop;

        //Initialise and run the sampler
        smc::sampler<double,smc::nullParams> Sampler(lNumber, HistoryType::NONE);
        theta_prop.sigv = 10.0;
        theta_prop.sigw = 10.0;

        sigv(0) = theta_prop.sigv;
        sigw(0) = theta_prop.sigw;
        
        Rcpp::NumericVector sigvInnovation = Rcpp::rnorm(lMCMCits,0,0.15);
        Rcpp::NumericVector sigwInnovation = Rcpp::rnorm(lMCMCits,0,0.08);
        Rcpp::NumericVector unifRands = Rcpp::runif(lMCMCits);

        // Getting a particle filtering estimate of the log likelihood.
        smc::moveset<double,smc::nullParams> Moveset(fInitialise, fMove, NULL);
        Sampler.SetResampleParams(ResampleType::MULTINOMIAL, 0.5);
        Sampler.SetMoveSet(Moveset);
        Sampler.Initialise();
        Sampler.IterateUntil(lIterates-1);
        loglike(0) = Sampler.GetLogNCPath();

        // Inverse gamma prior
        logprior(0) = logPrior(theta_prop);

        double MH_ratio;
        for (unsigned int i = 1; i<lMCMCits+1; i++){
            // RW proposal for parameters
            theta_prop.sigv = sigv(i-1) + sigvInnovation(i-1);
            theta_prop.sigw = sigw(i-1) + sigwInnovation(i-1);

            // Getting a particle filtering estimate of the log likelihood.
            Sampler.Initialise();
            Sampler.IterateUntil(lIterates-1);
            loglike_prop = Sampler.GetLogNCPath();

            // Inverse gamma prior
            logprior_prop = logPrior(theta_prop);

            MH_ratio = exp(loglike_prop - loglike(i-1) + logprior_prop - logprior(i-1));

            if (MH_ratio>unifRands(i-1)){
                sigv(i) = theta_prop.sigv;
                sigw(i) = theta_prop.sigw;
                loglike(i) = loglike_prop;
                logprior(i) = logprior_prop;
            } else {
                sigv(i) = sigv(i-1);
                sigw(i) = sigw(i-1);
                loglike(i) = loglike(i-1);
                logprior(i) = logprior(i-1);
            }
        }

        return Rcpp::DataFrame::create(Rcpp::Named("samples_sigv") = sigv,
        Rcpp::Named("samples_sigw") = sigw,
        Rcpp::Named("loglike") = loglike,
        Rcpp::Named("logprior") = logprior);
    }
    catch(smc::exception  e) {
        Rcpp::Rcout << e;
    }
    return R_NilValue;            // to provide a return
}

namespace nonLinPMMH {

    /// A function to calculate the log prior for a proposal. The prior for this example is IG(0.01,0.01).

    /// \param proposal     The proposed values of the parameters
    double logPrior(const parameters & proposal)
    {
        return 2*a_prior*log(b_prior)-2*lgamma(a_prior)-(a_prior+1)*log(proposal.sigv)-b_prior/proposal.sigv-(a_prior+1)*log(proposal.sigw)-b_prior/proposal.sigw;
    }


    /// A function to initialise a particle.

    /// \param X            A reference to the empty particle value
    /// \param logweight    A reference to the empty particle log weight
    /// \param param        Additional algorithm parameters
    void fInitialise(double & X, double & logweight, smc::nullParams & param)
    {
        X = R::rnorm(0.0,sqrt(5.0));
        double mean = std::pow(X,2)/20.0;
        logweight = R::dnorm(y(0),mean,theta_prop.sigw,TRUE);
    }

    /// The proposal function.

    /// \param lTim     The sampler iteration.
    /// \param X            A reference to the current particle value
    /// \param logweight    A reference to the current particle log weight
    /// \param param        Additional algorithm parameters
    void fMove(long lTime, double & X, double & logweight, smc::nullParams & param)
    {
        X = X/2.0 + 25.0*X/(1+std::pow(X,2)) + 8*cos(1.2*(lTime+1)) + R::rnorm(0.0,theta_prop.sigv);
        double mean = std::pow(X,2)/20.0;
        logweight += R::dnorm(y(lTime),mean,theta_prop.sigw,TRUE);
    }

}
