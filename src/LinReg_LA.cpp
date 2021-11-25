// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
//
// LinReg_LA.cpp: A simple example for estimating the parameters of a
// linear regression model using likelihood annealing SMC.
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

#include "LinReg_LA.h"
#include <RcppArmadillo.h>

#include <cstdio> 
#include <cstdlib>
#include <cstring>
#include <math.h>

#include <iostream>
#include <cmath>

namespace LinReg_LA {
    arma::mat covRW("2500 -2.5 0.03; -2.5 130.0 0.0; 0.03 0.0 0.04");
    arma::mat cholCovRW = arma::chol(covRW);
    const double a_prior = 3.0;
    const double b_prior = pow(2.0*300.0*300.0,-1.0);
}

using namespace LinReg_LA;

// LinRegLA() function callable from R via Rcpp:: 
// [[Rcpp::export]]
Rcpp::List LinRegLA_impl(arma::mat Data, arma::vec intemps, unsigned long lNumber) {     


    try {
        temps = intemps;
        
        lIterates = Data.n_rows;
        data.y = Data.col(0);
        data.x = Data.col(1);
        mean_x = arma::sum(data.x)/lIterates;
        
        long lTemps = temps.n_rows;
        
        //Initialise and run the sampler
		myMove = new LinReg_LA_move;
        smc::sampler<rad_state,smc::nullParams> Sampler(lNumber, HistoryType::RAM, myMove);
        
        Sampler.SetResampleParams(ResampleType::SYSTEMATIC, 0.5);
        Sampler.SetMcmcRepeats(10);
        Sampler.Initialise();
        
        arma::cube theta(lNumber,3,lTemps);
        arma::mat loglike(lNumber,lTemps), logprior(lNumber,lTemps), Weights(lNumber,lTemps);
        arma::vec ESS(lTemps);
        
        for (unsigned int i=0; i<lNumber; i++){
            theta.slice(0).row(i) = Sampler.GetParticleValueN(i).theta.t();
            loglike(i,0) = Sampler.GetParticleValueN(i).loglike;
            logprior(i,0) = Sampler.GetParticleValueN(i).logprior;
        }
        
        Weights.col(0) = Sampler.GetParticleWeight();
        ESS(0) = Sampler.GetESS();
        
        for(int n=1; n < lTemps; ++n) {
            Sampler.Iterate();
            
            for (unsigned int i=0; i<lNumber; i++){
                theta.slice(n).row(i) = Sampler.GetParticleValueN(i).theta.t();
                loglike(i,n) = Sampler.GetParticleValueN(i).loglike;
                logprior(i,n) = Sampler.GetParticleValueN(i).logprior;
            }
            
            Weights.col(n) = Sampler.GetParticleWeight();
            ESS(n) = Sampler.GetESS();
        }
        
        double logNC_standard = Sampler.GetLogNCPath();
        double logNC_ps_trap2 = Sampler.IntegratePathSampling(integrand_ps,width_ps, NULL);
        double logNC_ps_rect = Sampler.IntegratePathSampling(PathSamplingType::RECTANGLE,integrand_ps,width_ps, NULL);
        double logNC_ps_trap = Sampler.IntegratePathSampling(PathSamplingType::TRAPEZOID1,integrand_ps,width_ps, NULL);
        
		delete myMove;
		
        return Rcpp::List::create(
        Rcpp::Named("theta") = theta,
        Rcpp::Named("loglike") = loglike,
        Rcpp::Named("logprior") = logprior,
        Rcpp::Named("Weights") = Weights,
        Rcpp::Named("ESS") = ESS,
        Rcpp::Named("logNC_standard") = logNC_standard,
        Rcpp::Named("logNC_ps_rect") = logNC_ps_rect,
        Rcpp::Named("logNC_ps_trap") = logNC_ps_trap,
        Rcpp::Named("logNC_ps_trap2") = logNC_ps_trap2);
    }
    catch(smc::exception  e) {
        Rcpp::Rcout << e;
    }
    return R_NilValue;              // to provide a return 
}

namespace LinReg_LA {
    
    double integrand_ps(long lTime,const rad_state & value,  void *) { return logLikelihood(value);}    

    double width_ps(long lTime, void *){
        return (temps(lTime) - temps(lTime-1));
    }   

    ///The function corresponding to the log likelihood at specified position
    /// \param value        The state to consider 
    double logLikelihood(const rad_state & value){

        double sigma = std::pow(expl(value.theta(2)),0.5);
        arma::vec mean_reg = value.theta(0) + value.theta(1)*(data.x - mean_x);
        return arma::sum(-log(sigma) - pow(data.y - mean_reg,2.0)/(2.0*sigma*sigma) -0.5*log(2.0*M_PI));  

    }
    ///The function corresponding to the (unnormalised) log prior at a specified position
    /// \param value        The state to consider 
    double logPrior(const rad_state & value){
        return -log(1000.0)- pow(value.theta(0) - 3000.0,2.0)/(2.0*1000.0*1000.0) -log(100.0)- pow(value.theta(1) - 185.0,2.0)/(2.0*100.0*100.0) + value.theta(2)-1.0/b_prior/expl(value.theta(2)) -value.theta(2)*(a_prior+1.0);
    }

    ///A function to initialise a particle

    /// \param value        Reference to the empty particle value
    /// \param logweight    Refernce to the empty particle log weight
    /// \param param        Additional algorithm parameters
    void LinReg_LA_move::pfInitialise(rad_state & value, double & logweight, smc::nullParams & param)
    {
        // drawing from the prior
        value.theta = arma::zeros(3);
        value.theta(0) = R::rnorm(3000.0,1000.0);
        value.theta(1) = R::rnorm(185.0,100.0);
        value.theta(2) = log(pow(R::rgamma(3,pow(2.0*300.0*300.0,-1.0)),-1.0));
        value.loglike = logLikelihood(value);
        value.logprior = logPrior(value);
        logweight = temps(0)*value.loglike;
    }

    ///The proposal function.

    ///\param lTime         The sampler iteration.
    /// \param value        Reference to the current particle value
    /// \param logweight    Refernce to the current particle log weight
    /// \param param        Additional algorithm parameters
    void LinReg_LA_move::pfMove(long lTime, rad_state & value, double & logweight, smc::nullParams & param)
    {
        logweight += (temps(lTime) - temps(lTime-1))*logLikelihood(value);
    }

    ///The MCMC function.

    ///\param lTime         The sampler iteration.
    ///\param value         Reference to the value of the particle being moved
    ///\param logweight     Reference to the log weight of the particle being moved
    ///\param param         Additional algorithm parameters
    bool LinReg_LA_move::pfMCMC(long lTime, rad_state & value, double & logweight, smc::nullParams & param)
    {
        rad_state value_prop;
        value_prop.theta = value.theta + cholCovRW*Rcpp::as<arma::vec>(Rcpp::rnorm(3));            
        value_prop.loglike = logLikelihood(value_prop);
        value_prop.logprior = logPrior(value_prop);
        
        double MH_ratio = exp(temps(lTime)*(value_prop.loglike - value.loglike) + value_prop.logprior - value.logprior);
        
        if (MH_ratio>R::runif(0,1)){
            value = value_prop;
            return TRUE;
        }
        return FALSE;
    }
}
