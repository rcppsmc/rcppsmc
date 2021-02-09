// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
//
// blockpfgaussianopt.cpp: Rcpp integration of SMC library -- Block PF Gaussian
//
// Copyright (C) 2008 - 2009  Adam Johansen
// Copyright (C) 2012 - 2017  Dirk Eddelbuettel and Adam Johansen
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
#include "blockpfgaussianopt.h"

#include <cstdlib>
#include <cmath>

using namespace std;
using namespace BSPFG;

///The observations
namespace BSPFG {
    arma::vec y; 
    long lLag;
    long lIterates;
}

// [[Rcpp::export]]
Rcpp::List blockpfGaussianOpt_impl(arma::vec data, long part, long lag)
{
    long lNumber = part;
    lLag = lag;

    y = data;
    lIterates = y.size();

    //Initialise and run the sampler
	myMove = new BSPFG_move;
    smc::sampler<arma::vec,smc::nullParams> Sampler(lNumber, HistoryType::NONE, myMove);

    Sampler.SetResampleParams(ResampleType::SYSTEMATIC, 0.5);

    Sampler.Initialise();
    Sampler.IterateUntil(lIterates - 1);

    //Generate results
    arma::mat resValues(lNumber,lIterates);
    arma::vec resWeights = Sampler.GetParticleWeight();
    for(int i = 0; i < lNumber; ++i) 
    {
        resValues.row(i) = Sampler.GetParticleValueN(i).t();
    }
    
    double logNC = Sampler.GetLogNCPath();
	
	delete myMove;

    return Rcpp::List::create(Rcpp::_["weight"] = resWeights, Rcpp::_["values"] = resValues, Rcpp::_["logNC"] = logNC);
}


namespace BSPFG {

    /// The initialisation function
    ///
    /// \param value The value of the particle being moved
    /// \param logweight The log weight of the particle being moved
    /// \param param Additional algorithm parameters
    void BSPFG_move::pfInitialise(arma::vec & value, double & logweight, smc::nullParams & param)
    {
        value = arma::zeros<arma::vec>(lIterates);
        value(0) = R::rnorm(0.5 * y(0),1.0/sqrt(2.0));
        logweight = 1.0;
    }

    ///The proposal function.

    ///\param lTime The sampler iteration.
    ///\param value The value of the particle being moved
    ///\param logweight The log weight of the particle being moved
    ///\param param Additional algorithm parameters
    void BSPFG_move::pfMove(long lTime, arma::vec & value, double & logweight, smc::nullParams & param)
    {
        if(lTime == 1) {
            value(lTime) = (value(lTime-1) + y(lTime))/2.0 + R::rnorm(0.0,1.0/sqrt(2.0));
            logweight += -0.25*(y(lTime) - value(lTime-1))*(y(lTime)-value(lTime-1));
            return;
        }

        long lag = min(lTime,lLag);
    
        //These structures should really be made static 
        arma::vec mu(lag+1);
        arma::vec sigma(lag+1);
        arma::vec sigmah(lag+1);
        arma::vec mub(lag+1);

        // Forward filtering
        mu(0) = value(lTime-lag);
        sigma(0) = 0;
        for(int i = 1; i <= lag; ++i)
        {
            sigmah(i) = sigma(i-1) + 1;
            
            mu(i) = (sigmah(i) * y(lTime-lag+i) +  mu(i-1)) / (sigmah(i) + 1);
            sigma(i) = sigmah(i) / (sigmah(i) + 1);
        }
        // Backward smoothing
        mub(lag) = mu(lag);
        value(lTime) = R::rnorm(mub(lag),sqrt(sigma(lag)));
        for(int i = lag-1; i; --i)
        {
            mub(i) = (sigma(i)*value(lTime-lag+i+1) + mu(i)) / (sigma(i)+1);
            value(lTime-lag+i) = R::rnorm(mub(i),sqrt(sigma(lag)/(sigma(lag) + 1)));
        }

        // Importance weighting
        logweight += -0.5 * pow(y(lTime) - mu[lag-1],2.0) / (sigmah[lag]+1);

    }
}
