// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
//
// blockpfgaussianopt.cpp: Rcpp integration of SMC library -- Block PF Gaussian
//
// Copyright (C) 2008 - 2009  Adam Johansen
// Copyright (C) 2012 - 2017  Dirk Eddelbuettel and Adam Johansen
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
#include "rngR.h"

#include <cstdlib>
#include <cmath>
#include <vector>

using namespace std;

///The observations
namespace BSPFG {
    Rcpp::NumericVector y;
}

long lLag = 1;

extern "C" SEXP blockpfGaussianOpt_impl(SEXP dataS, SEXP partS, SEXP lagS)
{
    long lIterates;
    long lNumber = Rcpp::as<long>(partS);
    lLag = Rcpp::as<long>(lagS);

    y = Rcpp::NumericVector(dataS);
    lIterates = y.size();

    //Initialise and run the sampler
    smc::sampler<vector<double> > Sampler(lNumber, SMC_HISTORY_NONE);  
    smc::moveset<vector<double> > Moveset(fInitialiseBSPFG, fMoveBSPFG, NULL);

    Sampler.SetResampleParams(SMC_RESAMPLE_SYSTEMATIC, 0.5);
    Sampler.SetMoveSet(Moveset);

    Sampler.Initialise();
    Sampler.IterateUntil(lIterates - 1);

    //Generate results
    Rcpp::NumericMatrix resValues = Rcpp::NumericMatrix(lNumber,lIterates);
    Rcpp::NumericVector resWeights = Rcpp::NumericVector(lNumber);
    for(int i = 0; i < lNumber; ++i) 
    {
        vector<double> pValue = Sampler.GetParticleValueN(i);
        for(int j = 0; j < lIterates; ++j) {
            resValues(i,j) = pValue.at(j);
        }
        resWeights(i) = Sampler.GetParticleWeightN(i);
    }

    return Rcpp::List::create(Rcpp::_["weight"] = resWeights, Rcpp::_["values"] = resValues);
}

using namespace std;
using BSPFG::y;

/// The initialisation function
///
/// \param value The value of the particle being moved
/// \param logweight The log weight of the particle being moved
/// \param pRng A pointer to the random number generator which is to be used
void fInitialiseBSPFG(vector<double> & value, double & logweight, smc::rng *pRng)
{
    value.push_back(pRng->Normal(0.5 * y[0],1.0/sqrt(2.0)));
    logweight = 1.0;
}

///The proposal function.

///\param lTime The sampler iteration.
///\param value The value of the particle being moved
///\param logweight The log weight of the particle being moved
///\param pRng  A random number generator.
void fMoveBSPFG(long lTime, vector<double> & value, double & logweight, smc::rng *pRng)
{
    if(lTime == 1) {
        value.push_back((value.at(lTime-1) + y[int(lTime)])/2.0 + pRng->Normal(0.0,1.0/sqrt(2.0)));

        logweight += -0.25*(y[int(lTime)] - value.at(lTime-1))*(y[int(lTime)]-value.at(lTime-1));

        return;
    }

    long lag = min(lTime,lLag);

    //These structures should really be made static 
    std::vector<double> mu(lag+1);
    std::vector<double> sigma(lag+1);
    std::vector<double> sigmah(lag+1);
    std::vector<double> mub(lag+1);

    // Forward filtering
    mu[0] = value.at(lTime-lag);
    sigma[0] = 0;
    for(int i = 1; i <= lag; ++i)
    {
        sigmah[i] = sigma[i-1] + 1;
        
        mu[i] = (sigmah[i] * y[int(lTime-lag+i)] +  mu[i-1]) / (sigmah[i] + 1);
        sigma[i] = sigmah[i] / (sigmah[i] + 1);
    }
    // Backward smoothing
    mub[lag] = mu[lag];
    value.push_back(pRng->Normal(mub[lag],sqrt(sigma[lag])));
    for(int i = lag-1; i; --i)
    {
        mub[i] = (sigma[i]*value.at(lTime-lag+i+1) + mu[i]) / (sigma[i]+1);
        value.at(lTime-lag+i) = pRng->Normal(mub[i],sqrt(sigma[lag]/(sigma[lag] + 1)));
    }

    // Importance weighting
    logweight += -0.5 * pow(y[int(lTime)] - mu[lag-1],2.0) / (sigmah[lag]+1);

}
