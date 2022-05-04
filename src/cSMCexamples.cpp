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

namespace cSMCexamples {
    /// Initializing parameters
    const double varInit0 = 10;
    double VarEvol;
    double phi;
    /// Class declarations for measurements, states and parameters
    Measurements y;
    States X;
    Parameters params;
}

using namespace cSMCexamples;

// compareNCestimates_impl() function callable from R via Rcpp::
// [[Rcpp::export]]
Rcpp::List compareNCestimates_imp(arma::vec data,
                                  long lParticleNum,
                                  int simNum,
                                  Rcpp::List parInits)
{
    // Initialize data class/container
    y.yObs = data;

    // Initialize parameter class
    params.phi = parInits["phi"];
    params.varEvol = parInits["varEvol"];

    // Initialize other constants needed to run the sampler
    unsigned int tt = y.yObs.size();

    // double NCestimateBPF;
    // double NCestimateCSMCmultinomial;
    // double NCestimateCSMCsystematic;
    // double NCestimateCSMCstratified;
    // double NCestimateCSMCresidual;
    // double NCestimateKFgroundTruth;
    // The follwing will be replaced with corresponding Kalman forward Filter or backward smoothing.
    // NCestimateKFgroundTruth   = 0.0;

    // Initialize output container:
    arma::mat outputNCestimatesKF(simNum, 2, arma::fill::zeros);
    arma::mat outputNCestimatesSMC(simNum, 4, arma::fill::zeros);
    arma::mat outputNCestimatesCSMC(simNum, 4, arma::fill::zeros);
    Rcpp::List outputNCestimates;
    try
    {
        // Create move-class object: same for SMC and CSMC.
        MyLGSSmove = new cSMCexamples_move;
        // Initialize baseline BPF sampler.
        smc::sampler<States,smc::nullParams> SamplerBPF(lParticleNum,
                                                        HistoryType::NONE,
                                                        MyLGSSmove);

        // smc::conditionalSampler<States,smc::nullParams> SamplerCBPF(lParticleNum, HistoryType::AL,MyLGSSmove);


        for(int i = 0; i < simNum; ++i) {
            SamplerBPF.Initialise();
            SamplerBPF.SetResampleParams(ResampleType::MULTINOMIAL, 0.5);
            SamplerBPF.IterateUntil(tt - 1);
            outputNCestimatesSMC.at(i, 0) = SamplerBPF.GetLogNCPath();

            SamplerBPF.Initialise();
            SamplerBPF.SetResampleParams(ResampleType::RESIDUAL, 0.5);
            SamplerBPF.IterateUntil(tt - 1);
            outputNCestimatesSMC.at(i, 1) = SamplerBPF.GetLogNCPath();

            SamplerBPF.Initialise();
            SamplerBPF.SetResampleParams(ResampleType::STRATIFIED, 0.5);
            SamplerBPF.IterateUntil(tt - 1);
            outputNCestimatesSMC.at(i, 2) = SamplerBPF.GetLogNCPath();

            SamplerBPF.Initialise();
            SamplerBPF.SetResampleParams(ResampleType::SYSTEMATIC, 0.5);
            SamplerBPF.IterateUntil(tt - 1);
            outputNCestimatesSMC.at(i, 3) = SamplerBPF.GetLogNCPath();

            // SamplerCBPF.Initialise();
            // SamplerCBPF.SetResampleParams(ResampleType::MULTINOMIAL, 0.5);
            // SamplerCBPF.IterateUntil(tt - 1);
            // outputNCestimatesCSMC.at(i, 0) = SamplerCBPF.GetLogNCPath();

            // SamplerCBPF.Initialise();
            // SamplerCBPF.SetResampleParams(ResampleType::RESIDUAL, 0.5);
            // SamplerCBPF.IterateUntil(tt - 1);
            // outputNCestimatesCSMC.at(i, 1) = SamplerCBPF.GetLogNCPath();

            // SamplerCBPF.Initialise();
            // SamplerCBPF.SetResampleParams(ResampleType::STRATIFIED, 0.5);
            // SamplerCBPF.IterateUntil(tt - 1);
            // outputNCestimatesCSMC.at(i, 2) = SamplerCBPF.GetLogNCPath();

            // SamplerCBPF.Initialise();
            // SamplerCBPF.SetResampleParams(ResampleType::SYSTEMATIC, 0.5);
            // SamplerCBPF.IterateUntil(tt - 1);
            // outputNCestimatesCSMC.at(i, 3) = SamplerCBPF.GetLogNCPath();

        }

        outputNCestimates["kfOut"] = outputNCestimatesKF;
        outputNCestimates["smcOut"] = outputNCestimatesSMC;
        outputNCestimates["csmcOut"] = outputNCestimatesCSMC;

        delete MyLGSSmove;
        return outputNCestimates;
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
