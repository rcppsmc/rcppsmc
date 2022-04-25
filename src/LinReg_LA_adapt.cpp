// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
//
// LinReg_LA_adapt.cpp: Rcpp wrapper for SMC library -- A simple example for
// estimating the parameters of a linear regression model using likelihood
// annealing SMC, with adaptation of the temperature schedule, the multivariate
// normal random walk covariance matrix and the number of MCMC repeats.
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

#include "LinReg_LA_adapt.h"

namespace LinReg_LA_adapt {
    const double a_prior = 3.0;
    const double b_prior = std::pow(2.0*300.0*300.0,-1.0);
}

using namespace LinReg_LA_adapt;

// LinRegLA_adapt_impl() function callable from R via Rcpp::
// [[Rcpp::export]]
Rcpp::List LinRegLA_adapt_impl(arma::mat Data, unsigned long lNumber, double resampTol, double tempTol) {


    try {
        rho = tempTol;

        lIterates = Data.n_rows;
        data.y = Data.col(0);
        data.x = Data.col(1);
        mean_x = arma::sum(data.x)/lIterates;

        // Initialise the sampler
        myAdapt = new rad_adapt;
        myMove = new rad_move;
        Sampler = new smc::sampler<rad_state,smc::staticModelAdapt>(lNumber, HistoryType::RAM, myMove);

        Sampler->SetResampleParams(ResampleType::SYSTEMATIC, resampTol);
        Sampler->SetAlgParam(smc::staticModelAdapt());
        Sampler->SetAdaptMethods(myAdapt);
        Sampler->Initialise();

        //Run the sampler
        while (Sampler->GetAlgParams().GetTempCurr() != 1){
            Sampler->Iterate();
        }

        //Storing the results in a sensible format
        std::vector<double> temps = Sampler->GetAlgParams().GetTemps();
        int lTemps = temps.size();

        arma::cube theta(lNumber,3,lTemps);
        arma::mat loglike(lNumber,lTemps), logprior(lNumber,lTemps), Weights(lNumber,lTemps);
        arma::vec ESS(lTemps), mcmcRepeats(lTemps);

        for(int n=0; n < lTemps; ++n) {
            for (unsigned int i=0; i<lNumber; i++){
                theta.slice(n).row(i) = Sampler->GetHistoryPopulation(n).GetValueN(i).theta.t();
                loglike(i,n) = Sampler->GetHistoryPopulation(n).GetValueN(i).loglike;
                logprior(i,n) = Sampler->GetHistoryPopulation(n).GetValueN(i).logprior;
            }
            Weights.col(n) = Sampler->GetHistoryPopulation(n).GetWeight();
            ESS(n) = Sampler->GetESS(n);
            mcmcRepeats(n) = Sampler->GetHistorymcmcRepeats(n);
        }

        double logNC_standard = Sampler->GetLogNCPath();
        double logNC_ps_trap2 = Sampler->IntegratePathSampling(integrand_ps,width_ps, NULL);
        double logNC_ps_rect = Sampler->IntegratePathSampling(PathSamplingType::RECTANGLE,integrand_ps,width_ps, NULL);
        double logNC_ps_trap = Sampler->IntegratePathSampling(PathSamplingType::TRAPEZOID1,integrand_ps,width_ps, NULL);

        delete myAdapt;
		delete myMove;

        return Rcpp::List::create(
        Rcpp::Named("theta") = theta,
        Rcpp::Named("loglike") = loglike,
        Rcpp::Named("logprior") = logprior,
        Rcpp::Named("Weights") = Weights,
        Rcpp::Named("ESS") = ESS,
        Rcpp::Named("Temps") = temps,
        Rcpp::Named("mcmcRepeats") = mcmcRepeats,
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

namespace LinReg_LA_adapt {

    double integrand_ps(long lTime,const rad_state & value,  void *) { return logLikelihood(value);}

    double width_ps(long lTime, void *){
        return (Sampler->GetAlgParams().GetTemp(lTime) - Sampler->GetAlgParams().GetTemp(lTime-1));
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
        return -log(1000.0)- std::pow(value.theta(0) - 3000.0,2.0)/(2.0*1000.0*1000.0) -log(100.0)- std::pow(value.theta(1) - 185.0,2.0)/(2.0*100.0*100.0) + value.theta(2)-1.0/b_prior/expl(value.theta(2)) -value.theta(2)*(a_prior+1.0);
    }

    ///A function to initialise a particle

    /// \param value        Reference to the empty particle value
    /// \param logweight    Refernce to the empty particle log weight
    /// \param param        Additional algorithm parameters
    void rad_move::pfInitialise(rad_state & value, double & logweight, smc::staticModelAdapt & params)
    {
        // drawing from the prior
        value.theta = arma::zeros(3);
        value.theta(0) = R::rnorm(3000.0,1000.0);
        value.theta(1) = R::rnorm(185.0,100.0);
        value.theta(2) = log(std::pow(R::rgamma(3,pow(2.0*300.0*300.0,-1.0)),-1.0));
        value.loglike = logLikelihood(value);
        value.logprior = logPrior(value);
        logweight = params.GetTemp(0)*value.loglike;
    }

    ///The proposal function.

    ///\param lTime         The sampler iteration.
    /// \param value        Reference to the current particle value
    /// \param logweight    Refernce to the current particle log weight
    /// \param param        Additional algorithm parameters
    void rad_move::pfMove(long lTime, rad_state & value, double & logweight, smc::staticModelAdapt & params)
    {
        logweight += (params.GetTemp(lTime) - params.GetTemp(lTime-1))*logLikelihood(value);
    }

    ///The MCMC function.

    ///\param lTime         The sampler iteration.
    ///\param value         Reference to the value of the particle being moved
    ///\param logweight     Reference to the log weight of the particle being moved
    ///\param param        Additional algorithm parameters
    bool rad_move::pfMCMC(long lTime, rad_state & value, double & logweight, smc::staticModelAdapt & params)
    {
        rad_state value_prop;
        value_prop.theta = value.theta + params.GetCholCov()*Rcpp::as<arma::vec>(Rcpp::rnorm(3));
        value_prop.loglike = logLikelihood(value_prop);
        value_prop.logprior = logPrior(value_prop);

        double MH_ratio = exp(params.GetTemp(lTime)*(value_prop.loglike - value.loglike) + value_prop.logprior - value.logprior);

        if (MH_ratio>R::runif(0,1)){
            value = value_prop;
            return TRUE;
        }
        return FALSE;
    }
}
