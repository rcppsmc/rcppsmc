// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
//
// sampler.h: Rcpp integration of SMC library -- sampler object
//
// Copyright (C) 2008 - 2009  Adam Johansen
// Copyright (C) 2017         Adam Johansen, Dirk Eddelbuettel and Leah South
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

//! \file
//! \brief Defines the overall sampler object.
//!
//! This file defines the smc::sampler class which is used to implement entire particle systems.

#ifndef __SMC_CONDITIONAL_SAMPLER_HH

#define __SMC_CONDITIONAL_SAMPLER_HH 1.0

#include "sampler.h"
#include <RcppArmadilloExtensions/sample.h>

namespace smc {
    template <class Space, class Params = nullParams> class conditionalSampler:
    public sampler<Space, Params>
    {
        using sampler<Space,Params>::MoveParticles;
        private:
            std::vector<Space> referenceTrajectory;
            arma::Col<unsigned int> referenceTrajectoryIndices;
        public:
            ///Create an particle system containing lSize uninitialised particles with the specified mode.
            conditionalSampler(long lSize, HistoryType::Enum htHistoryMode, std::vector<Space> referenceTrajectoryInit)
            : sampler<Space,Params>{lSize, htHistoryMode},
                referenceTrajectory{referenceTrajectoryInit}
            {
            }
            ///Create an particle system containing lSize uninitialised particles with the specified mode, additionally passing a moveset-object.
            conditionalSampler(long lSize, HistoryType::Enum htHistoryMode, moveset<Space,Params>* pNewMoves, std::vector<Space> referenceTrajectoryInit)
            : sampler<Space,Params>{lSize, htHistoryMode, pNewMoves},
                referenceTrajectory{referenceTrajectoryInit}
            {
            }
            ///Returns the reference trajectory
            std::vector<Space> GetReferenceTrajectory(void) const {return referenceTrajectory;}
            ///Returns a constant reference to the reference trajectory
            const std::vector<Space> & GetReferenceTrajectoryRefs(void) const {return referenceTrajectory;}
            ///Returns the n'th element of the reference trajectory
            Space GetReferenceValue(long n) const {return referenceTrajectory[n];}
            ///Returns a constant reference to the n'th element of the reference trajectory
            const Space & GetReferenceValueRefs(long n) const {return referenceTrajectory[n];}
            ///Returns index of conditional reference Trajectory
            const long GetReferenceValueIndex(long lTime) const {return(referenceTrajectoryIndices.at(lTime));
            }
            ///Sets the reference trajectory
            void SetReferenceTrajectory(const std::vector<Space>& newReferenceTrajectory) {referenceTrajectory = newReferenceTrajectory;}
            ///Sets the n'elment of the reference trajectory
            void SetReferenceValue(const Space & newReferenceValue, long n) {referenceTrajectory[n] = newReferenceValue;}
            ///Initialise the conditional sampler and its constituent particles. The reference trajectory used to set the first coordinate at T=0 is either from the last call to conditionalSampler.SetReferenceTrajectory(), or defaults to the trajectory passed to the constructor at the time the object is created.
            void Initialise(void);
            ///Initialise the sampler and its constituent particles.
            void Initialise(const std::vector<Space>& newReferenceTrajectory);
            ///Perform one conditional iteration of the simulation algorithm.
            void Iterate(void);
            ///Perform one conditional iteration of the simulation algorithm and return the resulting ess
            double IterateEss(void);
            ///Perform conditional iterations until the specified evolution time is reached
            void IterateUntil(long lTerminate);
            ///Resample the particle set using the specified resampling scheme specifically adjusted for conditional resampling
            void conditionalResample(ResampleType::Enum lMode);
    };
    /// At present this function resets the system evolution time to 0 and calls the moveset initialisor to assign each particle in the ensemble.
    ///
    /// Note that the initialisation function must be specified before calling this function.
    ///
    /// Additionally, the t=0 particle coordinate is set to the conditional reference value and re-weighted.
    template <class Space, class Params>
    void conditionalSampler<Space,Params>::Initialise(void)
    {
        sampler<Space,Params>::T = 0;
        sampler<Space,Params>::dlogNCPath = 0.0;
        sampler<Space,Params>::acceptProb = -1;

        //Set the initial values and log weights of the particles
        std::vector<Space> InitVal(sampler<Space,Params>::N);
        arma::vec InitWeights(sampler<Space,Params>::N);
        sampler<Space,Params>::pPopulation = population<Space>(InitVal,InitWeights);
        sampler<Space,Params>::pMoves->DoInit(sampler<Space,Params>::pPopulation,sampler<Space,Params>::N,sampler<Space,Params>::algParams);

        //Initialise the conditonal trajectory:
        //1. Sample uniformly initial period, T = 0, conditional index
        //   A. construct uniform weights in style of Armadillo 10.5 & earlier
        // arma::Col<double> tmpUniformWeights(sampler<Space,Params>::N, arma::fill::none);
        // tmpUniformWeights.fill(1.0/sampler<Space,Params>::N);
        Rcpp::NumericVector tmpUniformWeights(sampler<Space,Params>::N, 1.0/sampler<Space,Params>::N);
        //   B. sample uniformly from the grid
        // referenceTrajectoryIndices.at(T) = arma::sample(arma::linspace(0, sampler<Space,Params>::N - 1, sampler<Space,Params>::N), 1, false, tmpUniformWeights)(0);
        referenceTrajectoryIndices.at(sampler<Space,Params>::T) = Rcpp::sample(sampler<Space,Params>::N - 1, 1, false, tmpUniformWeights)[0];
        //2. Set first particle coordinate to conditional value at above index,
        //and re-weight using the DoConditionalMove-function (that, despite its name, works at initialization, T=0, as well as moves for subsequent T>=1 iterations)
        sampler<Space,Params>::pMoves->DoConditionalMove(sampler<Space,Params>::T,sampler<Space,Params>::pPopulation,referenceTrajectoryIndices(sampler<Space,Params>::T),sampler<Space,Params>::algParams);

        //Scaling weights by 1/N (for evidence computation)
        sampler<Space,Params>::pPopulation.SetLogWeight(sampler<Space,Params>::pPopulation.GetLogWeight() - log(static_cast<double>(sampler<Space,Params>::N)));

        //Estimate the normalising constant
        sampler<Space,Params>::dlogNCIt = sampler<Space,Params>::CalcLogNC();
        sampler<Space,Params>::dlogNCPath += sampler<Space,Params>::dlogNCIt;

        //Normalise the weights
        sampler<Space,Params>::pPopulation.SetLogWeight(sampler<Space,Params>::pPopulation.GetLogWeight() - sampler<Space,Params>::dlogNCIt);

        //Check if the ESS is below some reasonable threshold and resample if necessary.
        //A mechanism for setting this threshold is required.
        double ESS = sampler<Space,Params>::GetESS();
        if(ESS < sampler<Space,Params>::dResampleThreshold) {
            sampler<Space,Params>::nResampled = 1;
            sampler<Space,Params>::pAdapt->updateForMCMC(sampler<Space,Params>::algParams,sampler<Space,Params>::pPopulation,sampler<Space,Params>::acceptProb,sampler<Space,Params>::nResampled,sampler<Space,Params>::nRepeats);
            sampler<Space,Params>::Resample(sampler<Space,Params>::rtResampleMode);
        }
        else {
            sampler<Space,Params>::nResampled = 0;
            if(sampler<Space,Params>::htHistoryMode == HistoryType::AL) {
                sampler<Space,Params>::uRSIndices = arma::linspace<arma::Col<unsigned int>>(0, sampler<Space,Params>::N - 1, sampler<Space,Params>::N);
                //No resampling: set conditional index to previous iteration.
                referenceTrajectoryIndices.at(sampler<Space,Params>::T) = referenceTrajectoryIndices.at(sampler<Space,Params>::T - 1);
            }
            sampler<Space,Params>::pAdapt->updateForMCMC(sampler<Space,Params>::algParams,sampler<Space,Params>::pPopulation,sampler<Space,Params>::acceptProb,sampler<Space,Params>::nResampled,sampler<Space,Params>::nRepeats);
        }

        //A possible MCMC step should be included here.
        bool didMCMC =  sampler<Space,Params>::pMoves->DoMCMC(0,sampler<Space,Params>::pPopulation, sampler<Space,Params>::N, sampler<Space,Params>::nRepeats, sampler<Space,Params>::nAccepted, sampler<Space,Params>::algParams);
        if (didMCMC){
            sampler<Space,Params>::acceptProb = static_cast<double>(sampler<Space,Params>::nAccepted)/(static_cast<double>(sampler<Space,Params>::N)*static_cast<double>(sampler<Space,Params>::nRepeats));
        }
        //Normalise the weights
        sampler<Space,Params>::pPopulation.SetLogWeight(sampler<Space,Params>::pPopulation.GetLogWeight() - sampler<Space,Params>::CalcLogNC());

        //Perform any final updates to the additional algorithm parameters.
        sampler<Space,Params>::pAdapt->updateEnd(sampler<Space,Params>::algParams,sampler<Space,Params>::pPopulation);

        //Finally, the current particle set should be appended to the historical process.
        if(sampler<Space,Params>::htHistoryMode != HistoryType::NONE){
            sampler<Space,Params>::History.clear();
            historyelement<Space> histel;
            switch(sampler<Space,Params>::htHistoryMode) {
            case HistoryType::RAM:
                histel.Set(sampler<Space,Params>::N, sampler<Space,Params>::pPopulation, sampler<Space,Params>::nAccepted, sampler<Space,Params>::nRepeats, historyflags(sampler<Space,Params>::nResampled));
                break;
            case HistoryType::AL:
                histel.Set(sampler<Space,Params>::N, sampler<Space,Params>::pPopulation, sampler<Space,Params>::nAccepted, sampler<Space,Params>::nRepeats, historyflags(sampler<Space,Params>::nResampled), sampler<Space,Params>::uRSIndices);
                break;
            /// To avoid compiler warnings, HistoryType::NONE is handled
            case HistoryType::NONE:
                break;
            }
            sampler<Space,Params>::History.push_back(histel);
        }
        return;
    }
    /// At present this function resets the system evolution time to 0 and calls the moveset initialisor to assign each particle in the ensemble.
    ///
    /// Note that the initialisation function must be specified before calling this function.
    ///
    /// Additionally, a conditional particle trajectory is set and the t=0 particle coordinate is set to the conditional reference value and re-weighted.
    template <class Space, class Params>
    void conditionalSampler<Space,Params>::Initialise(const std::vector<Space> & referenceTrajectoryInit)
    {
        SetReferenceTrajectory(referenceTrajectoryInit);
        sampler<Space,Params>::T = 0;
        sampler<Space,Params>::dlogNCPath = 0.0;
        sampler<Space,Params>::acceptProb = -1;

        //Set the initial values and log weights of the particles
        std::vector<Space> InitVal(sampler<Space,Params>::N);
        arma::vec InitWeights(sampler<Space,Params>::N);
        sampler<Space,Params>::pPopulation = population<Space>(InitVal,InitWeights);
        sampler<Space,Params>::pMoves->DoInit(sampler<Space,Params>::pPopulation,sampler<Space,Params>::N,sampler<Space,Params>::algParams);

        //Initialise the conditonal trajectory:
        //1. Sample uniformly initial period, T = 0, conditional index
        //   A. construct uniform weights in style of Armadillo 10.5 & earlier
        // arma::Col<double> tmpUniformWeights(sampler<Space,Params>::N, arma::fill::none);
        // tmpUniformWeights.fill(1.0/sampler<Space,Params>::N);
        Rcpp::NumericVector tmpUniformWeights(sampler<Space,Params>::N, 1.0/sampler<Space,Params>::N);
        //   B. sample uniformly from the grid
        // referenceTrajectoryIndices.at(T) = arma::sample(arma::linspace(0, sampler<Space,Params>::N - 1, sampler<Space,Params>::N), 1, false, tmpUniformWeights)(0);
        referenceTrajectoryIndices.at(sampler<Space,Params>::T) = Rcpp::sample(sampler<Space,Params>::N - 1, 1, false, tmpUniformWeights)[0];
        //2. Set first particle coordinate to conditional value at above index,
        //and re-weight using the DoConditionalMove-function (that, despite its name, works at initialization, T=0, as well as moves for subsequent T>=1 iterations)
        sampler<Space,Params>::pMoves->DoConditionalMove(sampler<Space,Params>::T,sampler<Space,Params>::pPopulation,referenceTrajectoryIndices(sampler<Space,Params>::T),sampler<Space,Params>::algParams);

        //Scaling weights by 1/N (for evidence computation)
        sampler<Space,Params>::pPopulation.SetLogWeight(sampler<Space,Params>::pPopulation.GetLogWeight() - log(static_cast<double>(sampler<Space,Params>::N)));

        //Estimate the normalising constant
        sampler<Space,Params>::dlogNCIt = sampler<Space,Params>::CalcLogNC();
        sampler<Space,Params>::dlogNCPath += sampler<Space,Params>::dlogNCIt;

        //Normalise the weights
        sampler<Space,Params>::pPopulation.SetLogWeight(sampler<Space,Params>::pPopulation.GetLogWeight() - sampler<Space,Params>::dlogNCIt);

        //Check if the ESS is below some reasonable threshold and resample if necessary.
        //A mechanism for setting this threshold is required.
        double ESS = sampler<Space,Params>::GetESS();
        if(ESS < sampler<Space,Params>::dResampleThreshold) {
            sampler<Space,Params>::nResampled = 1;
            sampler<Space,Params>::pAdapt->updateForMCMC(sampler<Space,Params>::algParams,sampler<Space,Params>::pPopulation,sampler<Space,Params>::acceptProb,sampler<Space,Params>::nResampled,sampler<Space,Params>::nRepeats);
            sampler<Space,Params>::Resample(sampler<Space,Params>::rtResampleMode);
        }
        else {
            sampler<Space,Params>::nResampled = 0;
            if(sampler<Space,Params>::htHistoryMode == HistoryType::AL) {
                sampler<Space,Params>::uRSIndices = arma::linspace<arma::Col<unsigned int>>(0, sampler<Space,Params>::N - 1, sampler<Space,Params>::N);
                //No resampling: set conditional index to previous iteration.
                referenceTrajectoryIndices.at(sampler<Space,Params>::T) = referenceTrajectoryIndices.at(sampler<Space,Params>::T - 1);
            }
            sampler<Space,Params>::pAdapt->updateForMCMC(sampler<Space,Params>::algParams,sampler<Space,Params>::pPopulation,sampler<Space,Params>::acceptProb,sampler<Space,Params>::nResampled,sampler<Space,Params>::nRepeats);
        }

        //A possible MCMC step should be included here.
        bool didMCMC =  sampler<Space,Params>::pMoves->DoMCMC(0,sampler<Space,Params>::pPopulation, sampler<Space,Params>::N, sampler<Space,Params>::nRepeats, sampler<Space,Params>::nAccepted, sampler<Space,Params>::algParams);
        if (didMCMC){
            sampler<Space,Params>::acceptProb = static_cast<double>(sampler<Space,Params>::nAccepted)/(static_cast<double>(sampler<Space,Params>::N)*static_cast<double>(sampler<Space,Params>::nRepeats));
        }
        //Normalise the weights
        sampler<Space,Params>::pPopulation.SetLogWeight(sampler<Space,Params>::pPopulation.GetLogWeight() - sampler<Space,Params>::CalcLogNC());

        //Perform any final updates to the additional algorithm parameters.
        sampler<Space,Params>::pAdapt->updateEnd(sampler<Space,Params>::algParams,sampler<Space,Params>::pPopulation);

        //Finally, the current particle set should be appended to the historical process.
        if(sampler<Space,Params>::htHistoryMode != HistoryType::NONE){
            sampler<Space,Params>::History.clear();
            historyelement<Space> histel;
            switch(sampler<Space,Params>::htHistoryMode) {
            case HistoryType::RAM:
                histel.Set(sampler<Space,Params>::N, sampler<Space,Params>::pPopulation, sampler<Space,Params>::nAccepted, sampler<Space,Params>::nRepeats, historyflags(sampler<Space,Params>::nResampled));
                break;
            case HistoryType::AL:
                histel.Set(sampler<Space,Params>::N, sampler<Space,Params>::pPopulation, sampler<Space,Params>::nAccepted, sampler<Space,Params>::nRepeats, historyflags(sampler<Space,Params>::nResampled), sampler<Space,Params>::uRSIndices);
                break;
            /// To avoid compiler warnings, HistoryType::NONE is handled
            case HistoryType::NONE:
                break;
            }
            sampler<Space,Params>::History.push_back(histel);
        }
        return;
    }
    template <class Space, class Params>
    double conditionalSampler<Space,Params>::IterateEss()
    {

        sampler<Space,Params>::pAdapt->updateForMove(this->algParams,sampler<Space,Params>::pPopulation);

        //Move the particle set.
        MoveParticles();

        //Do add a conditional conditional move:
        //set reference particle coordinate at conditional value and re-weight.
        sampler<Space,Params>::pMoves->DoConditionalMove(sampler<Space,Params>::T,sampler<Space,Params>::pPopulation,referenceTrajectoryIndices(sampler<Space,Params>::T),sampler<Space,Params>::algParams);

        //Estimate the normalising constant.
        sampler<Space,Params>::dlogNCIt = sampler<Space,Params>::CalcLogNC();
        sampler<Space,Params>::dlogNCPath += sampler<Space,Params>::dlogNCIt;

        //Normalise the weights.
        sampler<Space,Params>::pPopulation.SetLogWeight(sampler<Space,Params>::pPopulation.GetLogWeight() - sampler<Space,Params>::dlogNCIt);

        //Check if the ESS is below some reasonable threshold and resample if necessary.
        //A mechanism for setting this threshold is required.
        double ESS = sampler<Space,Params>::GetESS();
        if(ESS < sampler<Space,Params>::dResampleThreshold) {
            sampler<Space,Params>::nResampled = 1;
            sampler<Space,Params>::pAdapt->updateForMCMC(sampler<Space,Params>::algParams,sampler<Space,Params>::pPopulation,sampler<Space,Params>::acceptProb,sampler<Space,Params>::nResampled,sampler<Space,Params>::nRepeats);
            sampler<Space,Params>::Resample(sampler<Space,Params>::rtResampleMode);
        }
        else {
            sampler<Space,Params>::nResampled = 0;
            if(sampler<Space,Params>::htHistoryMode == HistoryType::AL) {
                sampler<Space,Params>::uRSIndices = arma::linspace<arma::Col<unsigned int>>(0, sampler<Space,Params>::N - 1, sampler<Space,Params>::N);
                //No resampling: set conditional index to previous iteration.
                referenceTrajectoryIndices.at(sampler<Space,Params>::T + 1) = referenceTrajectoryIndices.at(sampler<Space,Params>::T);
            }
            sampler<Space,Params>::pAdapt->updateForMCMC(sampler<Space,Params>::algParams,sampler<Space,Params>::pPopulation,sampler<Space,Params>::acceptProb,sampler<Space,Params>::nResampled,sampler<Space,Params>::nRepeats);
        }

        //A possible MCMC step should be included here.
        bool didMCMC = sampler<Space,Params>::pMoves->DoMCMC(sampler<Space,Params>::T+1,sampler<Space,Params>::pPopulation, sampler<Space,Params>::N, sampler<Space,Params>::nRepeats, sampler<Space,Params>::nAccepted,sampler<Space,Params>::algParams);
        if (didMCMC){
            sampler<Space,Params>::acceptProb = static_cast<double>(sampler<Space,Params>::nAccepted)/(static_cast<double>(sampler<Space,Params>::N)*static_cast<double>(sampler<Space,Params>::nRepeats));
        }


        //Normalise the weights
        sampler<Space,Params>::pPopulation.SetLogWeight(sampler<Space,Params>::pPopulation.GetLogWeight() - sampler<Space,Params>::CalcLogNC());

        //Perform any final updates to the additional algorithm parameters.
        sampler<Space,Params>::pAdapt->updateEnd(sampler<Space,Params>::algParams,sampler<Space,Params>::pPopulation);

        //Finally, the current particle set should be appended to the historical process.
        if(sampler<Space,Params>::htHistoryMode != HistoryType::NONE){
            historyelement<Space> histel;
            switch(sampler<Space,Params>::htHistoryMode) {
            case HistoryType::RAM:
                histel.Set(sampler<Space,Params>::N, sampler<Space,Params>::pPopulation, sampler<Space,Params>::nAccepted, sampler<Space,Params>::nRepeats, historyflags(sampler<Space,Params>::nResampled));
                break;
            case HistoryType::AL:
                histel.Set(sampler<Space,Params>::N, sampler<Space,Params>::pPopulation, sampler<Space,Params>::nAccepted, sampler<Space,Params>::nRepeats, historyflags(sampler<Space,Params>::nResampled), sampler<Space,Params>::uRSIndices);
                break;
            /// To avoid compiler warnings, HistoryType::NONE is handled
            case HistoryType::NONE:
                break;
            }
            sampler<Space,Params>::History.push_back(histel);
        }
        // Increment the evolution time.
        sampler<Space,Params>::T++;

        return ESS;
    }
    template <class Space, class Params>
    void conditionalSampler<Space,Params>::Iterate(void)
    {
        IterateEss();
        return;
    }
    template <class Space, class Params>
    void conditionalSampler<Space,Params>::IterateUntil(long lTerminate)
    {
        while(sampler<Space,Params>::GetTime() < lTerminate)
        Iterate();
    }
    template <class Space, class Params>
    void conditionalSampler<Space, Params>::conditionalResample(ResampleType::Enum lMode)
    {
        //Conditional resampling performed following the algorithms outlined in Appendix C of the paper "...".
        int uMultinomialCount;

        switch(lMode) {
        case ResampleType::MULTINOMIAL:
        default:
            {
                //Algorithm 3
                //Step 0:
                //Sample conditional index K_t from appropriate version of the "lambda" distribution i.e. uniformly on {1,...,N} in case of Multinmial resampling.
                //    0.1. Generate uniform weights
                Rcpp::NumericVector tmpUniformWeights(sampler<Space,Params>::N, 1.0/sampler<Space,Params>::N);
                //    0.2. Sample one conditional index from uniform distribution.
                long Kt = Rcpp::sample(sampler<Space,Params>::N - 1, 1, false, tmpUniformWeights)[0]; // sample lamba(k_{t}|w_{t-1}, k_{t-1})=1/N
                referenceTrajectoryIndices.at(sampler<Space,Params>::T + 1) = Kt; // update referenceTrajectoryIndices with newly sampled K_t index.
                //Step 1:
                //Connect the "chosen" ancestor index to previous reference trajectory: A_{t - 1}^{K_t} = K_{t - 1}
                sampler<Space,Params>::uRSIndices.at(Kt) = referenceTrajectoryIndices.at(sampler<Space,Params>::T);
                //Step 2:
                //Sample remaining ancestors A_{t - 1}^{-K_t} i.i.d. from a categorical distribution.
                //    2.1. Generate weights for categorical distribution.
                sampler<Space,Params>::dRSWeights = exp(sampler<Space,Params>::pPopulation.GetLogWeight() - stableLogSumWeights(sampler<Space,Params>::pPopulation.GetLogWeight()));
                //    2.2 Sample the remaining N-1 ancestor indices from {1,...,N}.
                Rcpp::IntegerVector tmpAncestorIndices(sampler<Space,Params>::N - 1);
                tmpAncestorIndices = Rcpp::sample(sampler<Space,Params>::N - 1, sampler<Space,Params>::N - 1, true, sampler<Space,Params>::dRSWeights);
                //    2.3. Assign ancestor indices to offspring {1,...,N}\{Kt}.
                std::vector<unsigned int> tmpIterator(sampler<Space,Params>::N - 1);
                std::iota(tmpIterator.begin(), tmpIterator.end(), 0); // define appropriate tmpIterator as a sequence from 0 to N-1
                tmpIterator.erase(tmpIterator.begin() + Kt); //exclude the previoiusly sampled conditional index K_t.
                long intIncrement = 0;
                // assign sampled ancestors to children
                for (int i : tmpIterator) {
                    sampler<Space,Params>::uRSIndices.at(i) = tmpAncestorIndices(intIncrement);
                    ++intIncrement;
                }
                break;
            }
        case ResampleType::STRATIFIED:
            {
                //Algorithm 4
                //Step 0:
                //Sample conditional index K_t from appropriate version of the "lambda" distribution i.e. the distribution over the stratum with \lambda(k_t|w_{t - 1}^{1:N}, k_{t - 1})=p_{t - 1}^{k_{t - 1}(k_t)/W_{t - 1}^{k_{t - 1}} in case of stratified resampling.
                //    0.1. Calculate normalized particle weights and cumulative normalized weights.
                sampler<Space,Params>::dRSWeights = exp(sampler<Space,Params>::pPopulation.GetLogWeight() - stableLogSumWeights(sampler<Space,Params>::pPopulation.GetLogWeight()));
                arma::vec dRSWeightsCumulative = arma::cumsum(sampler<Space,Params>::dRSWeights);
                //    0.2. Calculate strata boundaries cumulative weight components for min() and max() computation parts.
                arma::Col<double> strataBoundariesAll = arma::linspace(0.0, 1.0, sampler<Space,Params>::N + 1);
                arma::Col<double> strataBoundariesUpper = strataBoundariesAll.tail(sampler<Space,Params>::N);
                arma::Col<double> strataBoundariesLower = strataBoundariesAll.head(sampler<Space,Params>::N);
                double minWeights = dRSWeightsCumulative.at(referenceTrajectoryIndices.at(sampler<Space,Params>::T));
                double maxWeights = 0;
                if(referenceTrajectoryIndices.at(sampler<Space,Params>::T) != 0){
                    maxWeights = dRSWeightsCumulative.at(referenceTrajectoryIndices.at(sampler<Space,Params>::T - 1));
                }
                //    0.3. Calculate strata weights p_{t - 1}^{k_{t - 1}}(k_t).
                arma::Col<double> strataWeights(sampler<Space,Params>::N);
                strataWeights.fill(0.0);
                for(int i = 0; i < sampler<Space,Params>::N; ++i) {
                    strataWeights.at(i) = std::min(minWeights, strataBoundariesUpper.at(i)) - std::max(maxWeights, strataBoundariesLower.at(i));
                }
                //    0.4. Calculate \lambda(k_t|.) distribution.
                Rcpp::NumericVector lambdaWeightsStratified = Rcpp::wrap(strataWeights/dRSWeightsCumulative.at(referenceTrajectoryIndices.at(sampler<Space,Params>::T)));
                //    0.5. Sample K_t from 0.4.
                long Kt = Rcpp::sample(sampler<Space,Params>::N - 1, 1, false, lambdaWeightsStratified)[0];
                //Step 1:
                //Connect the "chosen" ancestor index to previous reference trajectory: A_{t - 1}^{K_t} = K_{t - 1}
                sampler<Space,Params>::uRSIndices.at(Kt) = referenceTrajectoryIndices.at(sampler<Space,Params>::T);
                //Step 2:
                //Calculation of empirical distribution function F_{t - 1}^N(i) is done and equal to computation of cumulative normalized weights stored in dRSWeightsCumulative.
                //Step 3: Generate ancestor indices and sssign to offspring indices {1,...,N}\{Kt}.
                std::vector<unsigned int> tmpIterator(sampler<Space,Params>::N - 1);
                std::iota(tmpIterator.begin(), tmpIterator.end(), 0); // define appropriate tmpIterator as a sequence from 0 to N-1
                tmpIterator.erase(tmpIterator.begin() + Kt); // exclude the previoiusly sampled conditional index K_t
                arma::Col<unsigned int> minimalJs = arma::linspace(0, sampler<Space,Params>::N - 1, sampler<Space,Params>::N);
                // generate ancestor index and assign to offspring
                double tmpUnifRnd;
                for (int i : tmpIterator) {
                    tmpUnifRnd = R::runif(0,1);
                    tmpUnifRnd += i - 1;
                    tmpUnifRnd /= sampler<Space,Params>::N;

                    minimalJs = arma::find(dRSWeightsCumulative.elem(minimalJs) > tmpUnifRnd);
                    sampler<Space,Params>::uRSIndices.at(i) = minimalJs.at(0);
                }
                break;
            }
        case ResampleType::SYSTEMATIC:
            {
                //Algorithm 5
                //Step 0:
                //Sample conditional index K_t from appropriate version of the "lambda" distribution i.e. the distribution over the stratum with \lambda(k_t|w_{t - 1}^{1:N}, k_{t - 1})=p_{t - 1}^{k_{t - 1}(k_t)/W_{t - 1}^{k_{t - 1}} in case of stratified resampling.
                //    0.1. Calculate normalized particle weights and cumulative normalized weights.
                sampler<Space,Params>::dRSWeights = exp(sampler<Space,Params>::pPopulation.GetLogWeight() - stableLogSumWeights(sampler<Space,Params>::pPopulation.GetLogWeight()));
                arma::vec dRSWeightsCumulative = arma::cumsum(sampler<Space,Params>::dRSWeights);
                //    0.2. Calculate strata boundaries cumulative weight components for min() and max() computation parts.
                arma::Col<double> strataBoundariesAll = arma::linspace(0.0, 1.0, sampler<Space,Params>::N + 1);
                arma::Col<double> strataBoundariesUpper = strataBoundariesAll.tail(sampler<Space,Params>::N);
                arma::Col<double> strataBoundariesLower = strataBoundariesAll.head(sampler<Space,Params>::N);
                double minWeights = dRSWeightsCumulative.at(referenceTrajectoryIndices.at(sampler<Space,Params>::T));
                double maxWeights = 0;
                if(referenceTrajectoryIndices.at(sampler<Space,Params>::T) != 0){
                    maxWeights = dRSWeightsCumulative.at(referenceTrajectoryIndices.at(sampler<Space,Params>::T - 1));
                }
                //    0.3. Calculate strata weights p_{t - 1}^{k_{t - 1}}(k_t).
                arma::Col<double> strataWeights(sampler<Space,Params>::N);
                strataWeights.fill(0.0);
                for(int i = 0; i < sampler<Space,Params>::N; ++i) {
                    strataWeights.at(i) = std::min(minWeights, strataBoundariesUpper.at(i)) - std::max(maxWeights, strataBoundariesLower.at(i));
                }
                //    0.4. Calculate \lambda(k_t|.) distribution.
                Rcpp::NumericVector lambdaWeightsStratified = Rcpp::wrap(strataWeights/dRSWeightsCumulative.at(referenceTrajectoryIndices.at(sampler<Space,Params>::T)));
                //    0.5. Sample K_t from 0.4.
                long Kt = Rcpp::sample(sampler<Space,Params>::N - 1, 1, false, lambdaWeightsStratified)[0];
                //Step 1:
                //Connect the "chosen" ancestor index to previous reference trajectory: A_{t - 1}^{K_t} = K_{t - 1}
                sampler<Space,Params>::uRSIndices.at(Kt) = referenceTrajectoryIndices.at(sampler<Space,Params>::T);
                //Step 2:
                //Calculation of empirical distribution function F_{t - 1}^N(i) is done and equal to computation of cumulative normalized weights stored in dRSWeightsCumulative.
                //Step 3: Generate ancestor indices and sssign to offspring indices {1,...,N}\{Kt}.
                std::vector<unsigned int> tmpIterator(sampler<Space,Params>::N - 1);
                std::iota(tmpIterator.begin(), tmpIterator.end(), 0); // define appropriate tmpIterator as a sequence from 0 to N-1
                tmpIterator.erase(tmpIterator.begin() + Kt); // exclude the previoiusly sampled conditional index K_t
                arma::Col<unsigned int> minimalJs = arma::linspace(0, sampler<Space,Params>::N - 1, sampler<Space,Params>::N);
                //precompute necessary uniform random variable before assignment
                double tmpVUpperBound = dRSWeightsCumulative.at(Kt);
                double tmpVLowerBound;
                if(Kt == 0) {
                    tmpVLowerBound = 0;
                } else {
                    tmpVLowerBound = dRSWeightsCumulative.at(Kt - 1);
                }
                double tmpV = R::runif(tmpVLowerBound, tmpVUpperBound);
                tmpV = sampler<Space,Params>::N * tmpV - std::floor(sampler<Space,Params>::N * tmpV);
                double tmpU;
                // generate ancestor index and assign to offspring
                for (int i : tmpIterator) {
                    tmpU = tmpV + i - 1;
                    tmpU /= sampler<Space,Params>::N;

                    minimalJs = arma::find(dRSWeightsCumulative.elem(minimalJs) > tmpU);
                    sampler<Space,Params>::uRSIndices.at(i) = minimalJs.at(0);
                }
                break;
            }
        }
        //Perform the replication of the chosen.
        for(int i = 0; i < sampler<Space,Params>::N ; ++i) {
            if(sampler<Space,Params>::uRSIndices(i) != static_cast<unsigned int>(i)){
                sampler<Space,Params>::pPopulation.SetValueN(sampler<Space,Params>::pPopulation.GetValueN(static_cast<int>(sampler<Space,Params>::uRSIndices(i))), i);
            }
        }
        //After conditional resampling is implemented: a final step is to set equal normalised weights.
        sampler<Space,Params>::pPopulation.SetLogWeight(- log(static_cast<double>(sampler<Space,Params>::N))*arma::ones(sampler<Space,Params>::N));
    }
}
#endif