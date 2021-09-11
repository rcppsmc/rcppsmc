// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
//
// Copyright (C) 2021 Adam Johansen, Dirk Eddelbuettel, Leah South and Ilya Zarubin
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

#include "RcppArmadillo.h"

#include "smctc.h"
#include "cSMCexamples.h"
#include <cmath>

namespace cSMCexamples {
    /// Initializing parameters
    const double varInit0 = 10;
    double VarEvol;
    double phi;
    ///The obsevations
    Measurements y;
    States X;
    Parameters params;
}

using namespace std;
using namespace cSMCexamples;

// runMonteCarloNCestimateComparisons_impl() function callable from R via Rcpp::
// [[Rcpp::export]]
Rcpp::DataFrame compareNCestimates_imp(arma::vec data,
                                       long lParticleNum,
                                       int simNum,
                                       Rcpp::List parInits)
{
    params.phi = parInits["phi"];
    params.varEvol = parInits["varEvol"];
    unsigned int tt = data.size();

    double NCestimateBPF;
    double NCestimateCSMCmultinomial;
    double NCestimateCSMCsystematic;
    double NCestimateCSMCstratified;
    double NCestimateCSMCresidual;
    double NCestimateKFgroundTruth;

    arma::mat outputNCestimates(simNum, 6);
    outputNCestimates.fill(0.5);
    try
    {
        // Create move-class object.
        MyLGSSmove = new cSMCexamples_move;
        // Initialize baseline BPF sampler.
        smc::sampler<States,smc::nullParams> SamplerBPF(lParticleNum,
                                                        HistoryType::NONE,
                                                        MyLGSSmove);
        SamplerBPF.SetResampleParams(ResampleType::RESIDUAL, 0.5);

        for(int i = 0; i < simNum; ++i) {
            SamplerBPF.Initialise();
            SamplerBPF.IterateUntil(tt - 1);
            NCestimateBPF = SamplerBPF.GetLogNCPath();

            // The following will be replaced with corresponding sampler output.
            NCestimateCSMCmultinomial = 0.0;
            NCestimateCSMCstratified  = 0.0;             NCestimateCSMCsystematic  = 0.0;
            NCestimateCSMCresidual    = 0.0;
            // The follwing will be replaced with corresponding Kalman forward Filter or backward smoothing.
            NCestimateKFgroundTruth   = 0.0;

            outputNCestimates.at(i, 0) = NCestimateBPF;
            outputNCestimates.at(i, 1) = NCestimateCSMCmultinomial;
            outputNCestimates.at(i, 2) = NCestimateCSMCsystematic;
            outputNCestimates.at(i, 3) = NCestimateCSMCstratified;
            outputNCestimates.at(i, 4) = NCestimateCSMCresidual;
            outputNCestimates.at(i, 5) = NCestimateKFgroundTruth;
        }
    }

    catch(smc::exception  e)
    {
        Rcpp::Rcout << e;
    }

    return R_NilValue; // to provide a return
}

namespace cSMCexamples
{
    ///The function corresponding to the log likelihood at specified time and position (up to normalisation)

    ///  \param lTime The current time (i.e. the index of the current distribution)
    ///  \param X     The state to consider
    double computeLogLikelihood(long lTime, const States& X)
    {
        return R::dnorm(arma::as_scalar(y.yObs[lTime]),
                        X.xState, 1.0, 1);
    }

    ///A function to initialise particles

    /// \param stateValue The value of the particle being moved
    /// \param logweight The log weight of the particle being moved
    /// \param param Additional algorithm parameters
    void cSMCexamples_move::pfInitialise(States& stateValue,
                                         double& logweight,
                                         smc::nullParams& param)
    {
        // Initialize state value
        stateValue.xState = R::rnorm(0.0,sqrt(varInit0));
        // Initilialize log-weight
        logweight = computeLogLikelihood(0, stateValue);
    }

    ///The proposal function

    ///\param lTime The sampler iteration.
    ///\param stateValue The value of the particle being moved.
    ///\param logweight The log weights of the particle being moved.
    ///\param param Additional algorithm parameters.
    void cSMCexamples_move::pfMove(long lTime,
                                   States& stateValue,
                                   double& logweight,
                                   smc::nullParams& param)
    {
        //Move/propose particles
        stateValue.xState = params.phi * stateValue.xState + R::rnorm(0.0,sqrt(params.varEvol));
        //Compute particle log-weight and increment
        logweight += computeLogLikelihood(lTime, stateValue);
    }
}