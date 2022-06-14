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

#include "cSMCexamples.h"

namespace cSMCexamples {
    /// Variable definitions
    double varObs;
    double stateInit;
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
                                  Rcpp::List parInits,
                                  arma::mat referenceTraj)
{
    // Initialize data class/container
    y.yObs = data;
    // Initialize parameter class
    params.statePhi = parInits["phi"];
    params.stateVarEvol = parInits["varStateEvol"];
    // Initialize other parameters: initial state and measurement variances
    stateInit = parInits["stateInit"];
    varObs = parInits["varObs"];


    // Initialize other constants needed to run the sampler
    const unsigned int tt = y.yObs.size();
    // Transform reference Trajectory of doubles to std::vector<States> type
    std::vector<States> referenceTrajectory(tt);
    copyReferenceTrajectory(referenceTraj.col(0), referenceTrajectory);

    // Initialize output container:
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
        // Initialize baseline conditional BPF sampler
        smc::conditionalSampler<States, smc::nullParams> cSamplerBPF(lParticleNum,
                                                                     HistoryType::AL,
                                                                     MyLGSSmove,
                                                                     referenceTrajectory);

        // Perform simulation: estimate (log-)likelihood/normalizing constant
        // simNum number of times under different resampling schemes for
        // conditional and standard BPFs and store the results (as returned
        // by ".GetLogNCPath()")
        for(int i = 0; i < simNum; ++i) {
            SamplerBPF.SetResampleParams(ResampleType::MULTINOMIAL, 0.5);
            SamplerBPF.Initialise();
            SamplerBPF.IterateUntil(tt - 1);
            outputNCestimatesSMC.at(i, 0) = SamplerBPF.GetLogNCPath();

            SamplerBPF.SetResampleParams(ResampleType::RESIDUAL, 0.5);
            SamplerBPF.Initialise();
            SamplerBPF.IterateUntil(tt - 1);
            outputNCestimatesSMC.at(i, 1) = SamplerBPF.GetLogNCPath();

            SamplerBPF.SetResampleParams(ResampleType::STRATIFIED, 0.5);
            SamplerBPF.Initialise();
            SamplerBPF.IterateUntil(tt - 1);
            outputNCestimatesSMC.at(i, 2) = SamplerBPF.GetLogNCPath();

            SamplerBPF.SetResampleParams(ResampleType::SYSTEMATIC, 0.5);
            SamplerBPF.Initialise();
            SamplerBPF.IterateUntil(tt - 1);
            outputNCestimatesSMC.at(i, 3) = SamplerBPF.GetLogNCPath();

            copyReferenceTrajectory(referenceTraj.col(i), referenceTrajectory);

            cSamplerBPF.SetResampleParams(ResampleType::MULTINOMIAL, 0.5);
            cSamplerBPF.SetReferenceTrajectory(referenceTrajectory);
            cSamplerBPF.Initialise();
            cSamplerBPF.IterateUntil(tt - 1);
            outputNCestimatesCSMC.at(i, 0) = cSamplerBPF.GetLogNCPath();

            cSamplerBPF.SetResampleParams(ResampleType::RESIDUAL, 0.5);
            cSamplerBPF.SetReferenceTrajectory(referenceTrajectory);
            cSamplerBPF.Initialise();
            cSamplerBPF.IterateUntil(tt - 1);
            outputNCestimatesCSMC.at(i, 1) = cSamplerBPF.GetLogNCPath();

            cSamplerBPF.SetResampleParams(ResampleType::STRATIFIED, 0.5);
            cSamplerBPF.SetReferenceTrajectory(referenceTrajectory);
            cSamplerBPF.Initialise();
            cSamplerBPF.IterateUntil(tt - 1);
            outputNCestimatesCSMC.at(i, 2) = cSamplerBPF.GetLogNCPath();

            cSamplerBPF.SetResampleParams(ResampleType::SYSTEMATIC, 0.5);
            cSamplerBPF.SetReferenceTrajectory(referenceTrajectory);
            cSamplerBPF.Initialise();
            cSamplerBPF.IterateUntil(tt - 1);
            outputNCestimatesCSMC.at(i, 3) = cSamplerBPF.GetLogNCPath();

            Rcpp::Rcout << "simulation run: " << i << ", out of total: "<< simNum << std::endl;
        }
        // add simulatio results to output, delete Move-class object and return
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
    double integrand_mean_x(const States& stateValueX, void *) {
        return stateValueX.xState;
    }
    ///The function corresponding to the log likelihood at specified time and
    ///position (up to normalisation)

    ///  \param lTime The current time (i.e. index of the current distribution)
    ///  \param X     The state to consider
    double computeLogLikelihood(long lTime, const States& stateValue)
    {
        return R::dnorm(arma::as_scalar(y.yObs(lTime)),
                        stateValue.xState, sqrt(varObs), 1);
    }

    ///A function to initialise particles

    /// \param stateValue The value of the particle being moved
    /// \param logweight The log weight of the particle being moved
    /// \param param Additional algorithm parameters
    void cSMCexamples_move::pfInitialise(States& stateValue,
                                         double& logweight,
                                         smc::nullParams& param)
    {
        // Initialize state value: X_0 set deterministically:
        stateValue.xState = stateInit;
        //Move/propose particles from t=0 to t = 1, i.e. X_0 to X_1
        stateValue.xState = params.statePhi * stateValue.xState + R::rnorm(0.0,sqrt(params.stateVarEvol));
        // Initilialize log-weight W_1
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
        stateValue.xState = params.statePhi * stateValue.xState + R::rnorm(0.0,sqrt(params.stateVarEvol));
        //Compute particle log-weight and increment
        logweight += computeLogLikelihood(lTime, stateValue);
    }
    ///Weighting function of the conditional/reference particle coordinate.
    ///
    ///Used implicitly by MyLGSSmove at the corresponding (derived) conditional
    ///SMC class.
    ///
    ///\param lTime            The sampler iteration.
    ///\param condStateValue   A reference to the current particle value
    ///\param logweight        A reference to the current particle log weight
    ///\param param            additional algorithm parameters
    void cSMCexamples_move::pfWeight(long lTime,
                                     States& condStateValue,
                                     double& logweight,
                                     smc::nullParams& param)
    {
        //Compute particle log-weight and increment
        logweight += computeLogLikelihood(lTime, condStateValue);
    }
    void copyReferenceTrajectory(const arma::vec& refArma,
                                 std::vector<States>& refStd)
    {
        int lenT = refArma.size();
        for (int i = 0; i < lenT; ++i) {
            refStd[i].xState = refArma[i];
        }
    }
}
