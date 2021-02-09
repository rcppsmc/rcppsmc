// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
//
// LinReg_LA_adapt.h: Rcpp wrapper for SMC library -- A simple example for
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

#include "smctc.h"

namespace LinReg_LA_adapt {
    
    class rad_state
    {
    public:
        arma::vec theta; // (alpha,beta,phi)
        double loglike;
        double logprior;
    };

    class rad_obs
    {
    public:
        arma::vec y, x;
    };

    rad_obs data;
    double mean_x;
    long lIterates;
    double rho;
    
    double integrand_ps(long,const rad_state &, void *);
    double width_ps(long, void *);
	
	double logLikelihood(const rad_state & value);
    double logPrior(const rad_state & value);
		
    //A derived class for the moves
    class rad_move:
    public smc::moveset<rad_state,smc::staticModelAdapt>
    {
    public:
	
        void pfInitialise(rad_state & value, double & logweight, smc::staticModelAdapt & param);
        void pfMove(long lTime, rad_state & value, double & logweight, smc::staticModelAdapt & param);
        bool pfMCMC(long lTime, rad_state & value, double & logweight, smc::staticModelAdapt & param);

        ~rad_move() {};

    };
	    
    //A derived class for adaptation which adapts parameters of type smc::staticModelAdapt
    class rad_adapt:
    public smc::adaptMethods<rad_state,smc::staticModelAdapt>
    {
    public:

        void updateForMove(smc::staticModelAdapt & param, const smc::population<rad_state> & pop) {
            unsigned long N = pop.GetNumber();
            arma::vec loglike(N);
            for (unsigned int i=0; i<N; i++){
                loglike(i) = pop.GetValueN(i).loglike;
            }
            param.ChooseTemp(pop.GetLogWeight(),loglike,rho*N);
        }

        void updateForMCMC(smc::staticModelAdapt & param, const smc::population<rad_state> & pop, double acceptProb, int nResampled, int & nRepeats) {
            if (nResampled == 0) {
                nRepeats = 0;
            } else {
                nRepeats = param.calcMcmcRepeats(acceptProb, 0.99, 10, 1000);
            }
            
            arma::mat thetaMat(pop.GetNumber(),3);
            for (long i=0; i<pop.GetNumber(); i++){
                thetaMat.row(i) = pop.GetValueN(i).theta.t();
            }
            param.calcCholCov(thetaMat,pop.GetLogWeight());
        }

        ~rad_adapt() {};

    };

    smc::sampler<rad_state,smc::staticModelAdapt> * Sampler;
    smc::adaptMethods<rad_state,smc::staticModelAdapt> * myAdapt;
	smc::moveset<rad_state,smc::staticModelAdapt>* myMove;
}
