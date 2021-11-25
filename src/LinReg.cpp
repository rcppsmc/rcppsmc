// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
//
// LinReg.cpp: A simple example for estimating the parameters of a
// linear regression model using data annealing SMC.
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

#include "LinReg.h"
#include <RcppArmadillo.h>

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <math.h>

#include <iostream>
#include <cmath>

namespace LinReg {
    arma::mat covRW("2500 -2.5 0.03; -2.5 130.0 0.0; 0.03 0.0 0.04");
    arma::mat cholCovRW = arma::chol(covRW);
    const double a_prior = 3.0;
    const double b_prior = std::pow(2.0*300.0*300.0,-1.0);
}

using namespace LinReg;


// LinReg() function callable from R via Rcpp::
// [[Rcpp::export]]
Rcpp::List LinReg_impl(arma::mat Data, unsigned long lNumber) {

    long lIterates;

    try {
        // Load observations -- or rather copy them in from R
        lIterates = Data.n_rows;
        data.y = Data.col(0);
        data.x = Data.col(1);
        mean_x = arma::sum(data.x)/lIterates;

        //Initialise and run the sampler
		myMove = new LinReg_move;
        smc::sampler<rad_state,smc::nullParams> Sampler(lNumber, HistoryType::RAM, myMove);

        Sampler.SetResampleParams(ResampleType::MULTINOMIAL, 0.5);
        Sampler.SetMcmcRepeats(10);
        Sampler.Initialise();
        Sampler.IterateUntil(lIterates-1);

        arma::mat theta(lNumber,3);
        arma::vec weights = Sampler.GetParticleWeight();

        for (unsigned int i = 0; i<lNumber; i++){
            theta.row(i) = Sampler.GetParticleValueN(i).theta.t();
        }

        double logNC = Sampler.GetLogNCPath();

		delete myMove;

        return Rcpp::List::create(Rcpp::Named("theta") = theta,Rcpp::Named("weights") = weights,
        Rcpp::Named("logNC") = logNC);
    }
    catch(smc::exception  e) {
        Rcpp::Rcout << e;
    }
    return R_NilValue;              // to provide a return
}

namespace LinReg {

    ///The function corresponding to the log likelihood at specified time and position (up to normalisation)

    /// \param lTime        The current time (i.e. the index of the current distribution)
    /// \param value        The state to consider
    double logWeight(long lTime, const rad_state & value){

        double mean_reg = value.theta(0) + value.theta(1)*(data.x(lTime) - mean_x);
        double sigma = std::pow(expl(value.theta(2)),0.5);
        return -log(sigma) - std::pow(data.y(lTime) - mean_reg,2.0)/(2.0*sigma*sigma) -0.5*log(2.0*M_PI);

    }

    ///The function corresponding to the (unnormalised) log posterior at specified time and position

    /// \param lTime        The current time (i.e. the index of the current distribution)
    /// \param value        The state to consider
    double logPosterior(long lTime, const rad_state & value){

        double log_prior = -log(1000.0)- std::pow(value.theta(0) - 3000.0,2.0)/(2.0*1000.0*1000.0) -log(100.0)- std::pow(value.theta(1) - 185.0,2.0)/(2.0*100.0*100.0) + value.theta(2)-1.0/b_prior/expl(value.theta(2)) -value.theta(2)*(a_prior+1.0);

        double sigma = std::pow(expl(value.theta(2)),0.5);

        double log_normpdf;

        if (lTime==0){
            double mean_reg = value.theta(0) + value.theta(1)*(data.x(0) - mean_x);
            log_normpdf = -log(sigma) - std::pow(data.y(0) - mean_reg,2.0)/(2.0*sigma*sigma) -0.5*log(2.0*M_PI);
        } else{
            arma::vec mean_reg = value.theta(0) + value.theta(1)*(data.x.rows(0,lTime) - mean_x);
            log_normpdf = arma::sum(-log(sigma) - pow(data.y.rows(0,lTime) - mean_reg,2.0)/(2.0*sigma*sigma) -0.5*log(2.0*M_PI));
        }

        return (log_normpdf + log_prior);
    }

    ///A function to initialise a particle

    /// \param value        Reference to the empty particle value
    /// \param logweight    Refernce to the empty particle log weight
    /// \param param        Additional algorithm parameters
    void LinReg_move::pfInitialise(rad_state & value, double & logweight, smc::nullParams & param)
    {
        value.theta = arma::zeros(3);
        // drawing from the prior
        value.theta(0) = R::rnorm(3000.0,1000.0);
        value.theta(1) = R::rnorm(185.0,100.0);
        value.theta(2) = log(std::pow(R::rgamma(3,std::pow(2.0*300.0*300.0,-1.0)),-1.0));

        logweight = logWeight(0, value);
    }

    ///The proposal function.

    /// \param lTime        The sampler iteration.
    /// \param value        Reference to the current particle value
    /// \param logweight    Refernce to the current particle log weight
    /// \param param        Additional algorithm parameters
    void LinReg_move::pfMove(long lTime, rad_state & value, double & logweight, smc::nullParams & param)
    {
        logweight += logWeight(lTime, value);
    }

    ///The MCMC function.

    /// \param lTime        The sampler iteration.
    /// \param value        Reference to the current particle value
    /// \param logweight    Reference to the log weight of the particle being moved
    /// \param param        Additional algorithm parameters
    bool LinReg_move::pfMCMC(long lTime, rad_state & value, double & logweight, smc::nullParams & param)
    {
        double logposterior_curr = logPosterior(lTime, value);

        rad_state value_prop;
        value_prop.theta = value.theta + cholCovRW*Rcpp::as<arma::vec>(Rcpp::rnorm(3));

        double logposterior_prop = logPosterior(lTime, value_prop);

        double MH_ratio = exp(logposterior_prop - logposterior_curr);

        if (MH_ratio>R::runif(0,1)){
            value = value_prop;
            logposterior_curr = logposterior_prop;
            return TRUE;
        }
        return FALSE;
    }
}
