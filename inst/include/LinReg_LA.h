// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
//
// LinReg_LA.h: Rcpp wrapper for SMC library -- A simple example for estimating
// the parameters of a linear regression model using likelihood annealing SMC.
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

namespace LinReg_LA {
    
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
    arma::vec temps;
	
	double logLikelihood(const rad_state & value);
    double logPrior(const rad_state & value);
	
	//A derived class for the moves
    class LinReg_LA_move:
    public smc::moveset<rad_state,smc::nullParams>
    {
    public:
    
        void pfInitialise(rad_state & value, double & logweight, smc::nullParams & param);
        void pfMove(long lTime, rad_state & value, double & logweight, smc::nullParams & param);
        bool pfMCMC(long lTime, rad_state & value, double & logweight, smc::nullParams & param);

        ~LinReg_LA_move() {};

    };
    
    double integrand_ps(long,const rad_state &, void *);
    double width_ps(long, void *);
    
	smc::moveset<rad_state,smc::nullParams>* myMove;
}
