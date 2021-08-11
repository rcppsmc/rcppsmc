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

namespace smc {
    template <class Space, class Params = nullParams> class conditionalSampler:
    public sampler<Space, Params>
    {
        ///List of inherited members of the base sampler-class
        ///Number of particles in the system.
        using sampler<Space,Params>::N;
        ///The current evolution time of the system.
        using sampler<Space,Params>::T;
        ///The resampling mode which is to be employed.
        using sampler<Space,Params>::rtResampleMode;
        ///The effective sample size at which resampling should be used.
        using sampler<Space,Params>::dResampleThreshold;
        //Structure used internally for resampling.
        using sampler<Space,Params>::dRSWeights;
        ///Structure used internally for resampling.
        using sampler<Space,Params>::uRSIndices;

        ///The particles within the system.
        using sampler<Space,Params>::pPopulation;
        ///The set of moves available.
        using sampler<Space,Params>::pMoves;
         /// The additional algorithm parameters.
        using sampler<Space,Params>::algParams;
        /// An object for adapting additional algorithm parameters
        using sampler<Space,Params>::pAdapt;

        ///The number of MCMC moves which have been accepted during this iteration
        using sampler<Space,Params>::nAccepted;
        ///A flag which tracks whether the ensemble was resampled during this iteration
        using sampler<Space,Params>::nResampled;
        ///The number of MCMC repeats to be performed. The default is 1 if an MCMC step is supplied.
        using sampler<Space,Params>::nRepeats;
        ///The proportion of accepted MCMC proposals in the most recent MCMC step, with a default of -1 if no MCMC steps have been performed.
        using sampler<Space,Params>::acceptProb;

        ///An estimate of the log normalising constant ratio over the entire path.
        using sampler<Space,Params>::dlogNCPath;
        ///An estimate of the log normalising constant ratio over the last step.
        using sampler<Space,Params>::dlogNCIt;

        ///A mode flag which indicates whether historical information is stored
        using sampler<Space,Params>::htHistoryMode;
        ///The historical process associated with the particle system.
        using sampler<Space,Params>::History;

        ///Move the particle set by proposing and applying an appropriate move to each particle.
        using sampler<Space,Params>::MoveParticles;
        ///Returns the crude normalising constant ratio estimate implied by the weights
        using sampler<Space,Params>::CalcLogNC;
        ///Returns the Effective Sample Size of the specified particle generation.
        using sampler<Space,Params>::GetESS;
        //Returns the current evolution time of the system.
        using sampler<Space,Params>::GetTime;

        private:
            std::vector<Space> referenceTrajectory;
            long maxT = referenceTrajectory.size();
            arma::Col<unsigned int> referenceTrajectoryIndices;
        public:
            ///Create an particle system containing lSize uninitialised particles with the specified mode.
            conditionalSampler(long lSize, HistoryType::Enum htHistoryMode, std::vector<Space> referenceTrajectoryInit)
            : sampler<Space,Params>{lSize, htHistoryMode},
                referenceTrajectory{referenceTrajectoryInit},
                referenceTrajectoryIndices(maxT, arma::fill::zeros)
            {
            }
            ///Create an particle system containing lSize uninitialised particles with the specified mode, additionally passing a moveset-object.
            conditionalSampler(long lSize, HistoryType::Enum htHistoryMode, moveset<Space,Params>* pNewMoves, std::vector<Space> referenceTrajectoryInit)
            : sampler<Space,Params>{lSize, htHistoryMode, pNewMoves},
                referenceTrajectory{referenceTrajectoryInit},
                referenceTrajectoryIndices(maxT, arma::fill::zeros)
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
        T = 0;
        dlogNCPath = 0.0;
        acceptProb = -1;

        //Set the initial values and log weights of the particles
        std::vector<Space> InitVal(N);
        arma::vec InitWeights(N);
        pPopulation = population<Space>(InitVal,InitWeights);
        pMoves->DoInit(pPopulation,N,algParams);

        //Initialise the conditonal trajectory:
        //1. Sample uniformly initial period, T = 0, conditional index
        referenceTrajectoryIndices.at(T) = floor(unif_rand()*static_cast<double>(N));
        //2. Set first particle coordinate to conditional value at above index,
        //and re-weight using the DoConditionalMove-function (that, despite its name, works at initialization, T=0, as well as moves for subsequent T>=1 iterations)
        pMoves->DoConditionalMove(T,pPopulation,referenceTrajectory[T],referenceTrajectoryIndices.at(T),algParams);

        //Scaling weights by 1/N (for evidence computation)
        pPopulation.SetLogWeight(pPopulation.GetLogWeight() - log(static_cast<double>(N)));

        //Estimate the normalising constant
        dlogNCIt = CalcLogNC();
        dlogNCPath += dlogNCIt;

        //Normalise the weights
        pPopulation.SetLogWeight(pPopulation.GetLogWeight() - dlogNCIt);

        //Check if the ESS is below some reasonable threshold and resample if necessary.
        //A mechanism for setting this threshold is required.
        double ESS = GetESS();
        if(ESS < dResampleThreshold) {
            nResampled = 1;
            pAdapt->updateForMCMC(algParams,pPopulation,acceptProb,nResampled,nRepeats);
            conditionalResample(rtResampleMode);
        }
        else {
            nResampled = 0;
            if(htHistoryMode == HistoryType::AL) {
                uRSIndices = arma::linspace<arma::Col<unsigned int>>(0, N - 1, N);
                //No resampling: set conditional index to previous iteration.
                referenceTrajectoryIndices.at(T + 1) = referenceTrajectoryIndices.at(T);
            }
            pAdapt->updateForMCMC(algParams,pPopulation,acceptProb,nResampled,nRepeats);
        }

        //A possible MCMC step should be included here.
        bool didMCMC =  pMoves->DoMCMC(0,pPopulation, N, nRepeats, nAccepted, algParams);
        if (didMCMC){
            acceptProb = static_cast<double>(nAccepted)/(static_cast<double>(N)*static_cast<double>(nRepeats));
        }
        //Normalise the weights
        pPopulation.SetLogWeight(pPopulation.GetLogWeight() - CalcLogNC());

        //Perform any final updates to the additional algorithm parameters.
        pAdapt->updateEnd(algParams,pPopulation);

        //Finally, the current particle set should be appended to the historical process.
        if(htHistoryMode != HistoryType::NONE){
            History.clear();
            historyelement<Space> histel;
            switch(htHistoryMode) {
            case HistoryType::RAM:
                histel.Set(N, pPopulation, nAccepted, nRepeats, historyflags(nResampled));
                break;
            case HistoryType::AL:
                histel.Set(N, pPopulation, nAccepted, nRepeats, historyflags(nResampled), uRSIndices);
                break;
            /// To avoid compiler warnings, HistoryType::NONE is handled
            case HistoryType::NONE:
                break;
            }
            History.push_back(histel);
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
        T = 0;
        dlogNCPath = 0.0;
        acceptProb = -1;

        //Set the initial values and log weights of the particles
        std::vector<Space> InitVal(N);
        arma::vec InitWeights(N);
        pPopulation = population<Space>(InitVal,InitWeights);
        pMoves->DoInit(pPopulation,N,algParams);

        //Initialise the conditonal trajectory:
        //1. Sample uniformly initial period, T = 0, conditional index
        referenceTrajectoryIndices.at(T) = floor(unif_rand()*static_cast<double>(N));
        //2. Set first particle coordinate to conditional value at above index,
        //and re-weight using the DoConditionalMove-function (that, despite its name, works at initialization, T=0, as well as moves for subsequent T>=1 iterations)
        pMoves->DoConditionalMove(T,pPopulation,referenceTrajectory[T],referenceTrajectoryIndices.at(T),algParams);

        //Scaling weights by 1/N (for evidence computation)
        pPopulation.SetLogWeight(pPopulation.GetLogWeight() - log(static_cast<double>(N)));

        //Estimate the normalising constant
        dlogNCIt = CalcLogNC();
        dlogNCPath += dlogNCIt;

        //Normalise the weights
        pPopulation.SetLogWeight(pPopulation.GetLogWeight() - dlogNCIt);

        //Check if the ESS is below some reasonable threshold and resample if necessary.
        //A mechanism for setting this threshold is required.
        double ESS = GetESS();
        if(ESS < dResampleThreshold) {
            nResampled = 1;
            pAdapt->updateForMCMC(algParams,pPopulation,acceptProb,nResampled,nRepeats);
            conditionalResample(rtResampleMode);
        }
        else {
            nResampled = 0;
            if(htHistoryMode == HistoryType::AL) {
                uRSIndices = arma::linspace<arma::Col<unsigned int>>(0, N - 1, N);
                // No resampling: set conditional index to previous iteration.
                referenceTrajectoryIndices.at(T + 1) = referenceTrajectoryIndices.at(T);
            }
            pAdapt->updateForMCMC(algParams,pPopulation,acceptProb,nResampled,nRepeats);
        }

        //A possible MCMC step should be included here.
        bool didMCMC =  pMoves->DoMCMC(0,pPopulation, N, nRepeats, nAccepted, algParams);
        if (didMCMC){
            acceptProb = static_cast<double>(nAccepted)/(static_cast<double>(N)*static_cast<double>(nRepeats));
        }
        //Normalise the weights
        pPopulation.SetLogWeight(pPopulation.GetLogWeight() - CalcLogNC());

        //Perform any final updates to the additional algorithm parameters.
        pAdapt->updateEnd(algParams,pPopulation);

        //Finally, the current particle set should be appended to the historical process.
        if(htHistoryMode != HistoryType::NONE){
            History.clear();
            historyelement<Space> histel;
            switch(htHistoryMode) {
            case HistoryType::RAM:
                histel.Set(N, pPopulation, nAccepted, nRepeats, historyflags(nResampled));
                break;
            case HistoryType::AL:
                histel.Set(N, pPopulation, nAccepted, nRepeats, historyflags(nResampled), uRSIndices);
                break;
            /// To avoid compiler warnings, HistoryType::NONE is handled
            case HistoryType::NONE:
                break;
            }
            History.push_back(histel);
        }
        return;
    }
    template <class Space, class Params>
    double conditionalSampler<Space,Params>::IterateEss()
    {

        pAdapt->updateForMove(this->algParams,pPopulation);

        //Move the particle set.
        MoveParticles();

        //Do add a conditional conditional move:
        //set reference particle coordinate at conditional value and re-weight.
        pMoves->DoConditionalMove(T,pPopulation,referenceTrajectory[T],referenceTrajectoryIndices.at(T),algParams);

        //Estimate the normalising constant.
        dlogNCIt = CalcLogNC();
        dlogNCPath += dlogNCIt;

        //Normalise the weights.
        pPopulation.SetLogWeight(pPopulation.GetLogWeight() - dlogNCIt);

        //Check if the ESS is below some reasonable threshold and resample if necessary.
        //A mechanism for setting this threshold is required.
        double ESS = GetESS();
        if(ESS < dResampleThreshold) {
            nResampled = 1;
            pAdapt->updateForMCMC(algParams,pPopulation,acceptProb,nResampled,nRepeats);
            conditionalResample(rtResampleMode);
        }
        else {
            nResampled = 0;
            if(htHistoryMode == HistoryType::AL) {
                uRSIndices = arma::linspace<arma::Col<unsigned int>>(0, N - 1, N);
                //No resampling: set conditional index to previous iteration.
                referenceTrajectoryIndices.at(T + 1) = referenceTrajectoryIndices.at(T);
            }
            pAdapt->updateForMCMC(algParams,pPopulation,acceptProb,nResampled,nRepeats);
        }

        //A possible MCMC step should be included here.
        bool didMCMC = pMoves->DoMCMC(T+1,pPopulation, N, nRepeats, nAccepted,algParams);
        if (didMCMC){
            acceptProb = static_cast<double>(nAccepted)/(static_cast<double>(N)*static_cast<double>(nRepeats));
        }


        //Normalise the weights
        pPopulation.SetLogWeight(pPopulation.GetLogWeight() - CalcLogNC());

        //Perform any final updates to the additional algorithm parameters.
        pAdapt->updateEnd(algParams,pPopulation);

        //Finally, the current particle set should be appended to the historical process.
        if(htHistoryMode != HistoryType::NONE){
            historyelement<Space> histel;
            switch(htHistoryMode) {
            case HistoryType::RAM:
                histel.Set(N, pPopulation, nAccepted, nRepeats, historyflags(nResampled));
                break;
            case HistoryType::AL:
                histel.Set(N, pPopulation, nAccepted, nRepeats, historyflags(nResampled), uRSIndices);
                break;
            /// To avoid compiler warnings, HistoryType::NONE is handled
            case HistoryType::NONE:
                break;
            }
            History.push_back(histel);
        }
        // Increment the evolution time.
        T++;

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
        while(GetTime() < lTerminate)
        Iterate();
    }
    template <class Space, class Params>
    void conditionalSampler<Space, Params>::conditionalResample(ResampleType::Enum lMode)
    {
        //Conditional resampling performed following the algorithms outlined in Appendix C of the paper "...".

        uRSIndices = arma::zeros<arma::Col<unsigned int> >(static_cast<int>(N));

        switch(lMode) {
        case ResampleType::MULTINOMIAL:
        default:
            {
                //Algorithm 3
                //Step 0:
                //Sample conditional index K_t from appropriate version of the "lambda" distribution i.e. uniformly on {1,...,N} in case of Multinmial resampling.
                long Kt = floor(unif_rand()*static_cast<double>(N)); // sample lamba(k_{t}|w_{t-1}, k_{t-1})=1/N
                referenceTrajectoryIndices.at(T + 1) = Kt; // update referenceTrajectoryIndices with newly sampled K_t index.
                //Step 1:
                //Connect the "chosen" ancestor index to previous reference trajectory: A_{t - 1}^{K_t} = K_{t - 1}
                uRSIndices.at(Kt) = referenceTrajectoryIndices.at(T);
                //Step 2:
                //Sample remaining ancestors A_{t - 1}^{-K_t} i.i.d. from a categorical distribution.
                //    2.1. Generate weights for categorical distribution.
                dRSWeights = exp(pPopulation.GetLogWeight() - stableLogSumWeights(pPopulation.GetLogWeight()));
                //    2.2 Sample remaining N-1 ancestor indices from {1,...,N}.
                Rcpp::IntegerVector tmpAncestorIndices(N - 1);
                tmpAncestorIndices = Rcpp::sample(N, N - 1, true,  Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(dRSWeights))) - 1;
                //    2.3. Assign ancestor indices to offspring {1,...,N}\{Kt}.
                std::vector<unsigned int> tmpIterator(N);
                std::iota(tmpIterator.begin(), tmpIterator.end(), 0); // define appropriate tmpIterator as a sequence from 0 to N-1
                tmpIterator.erase(tmpIterator.begin() + Kt); //exclude the previoiusly sampled conditional index K_t.
                long intIncrement = 0;
                // assign sampled ancestors to children
                for (int i : tmpIterator) {
                    uRSIndices.at(i) = tmpAncestorIndices(intIncrement);
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
                dRSWeights = exp(pPopulation.GetLogWeight() - stableLogSumWeights(pPopulation.GetLogWeight()));
                arma::vec dRSWeightsCumulative = arma::cumsum(dRSWeights);
                //    0.2. Calculate strata boundaries cumulative weight components for min() and max() computation parts.
                arma::Col<double> strataBoundariesAll = arma::linspace(0.0, 1.0, N + 1);
                arma::Col<double> strataBoundariesUpper = strataBoundariesAll.tail(N);
                arma::Col<double> strataBoundariesLower = strataBoundariesAll.head(N);
                double minWeights = dRSWeightsCumulative.at(referenceTrajectoryIndices.at(T));
                double maxWeights = 0;
                if(T != 0){
                    maxWeights = dRSWeightsCumulative.at(referenceTrajectoryIndices.at(T - 1));
                }
                //    0.3. Calculate strata weights p_{t - 1}^{k_{t - 1}}(k_t).
                arma::Col<double> strataWeights(N);
                strataWeights.fill(0.0);
                for(int i = 0; i < N; ++i) {
                    strataWeights.at(i) = std::min(minWeights, strataBoundariesUpper.at(i)) - std::max(maxWeights, strataBoundariesLower.at(i));
                }
                //    0.4. Calculate \lambda(k_t|.) distribution.
                Rcpp::NumericVector lambdaWeightsStratified = Rcpp::wrap(strataWeights/dRSWeightsCumulative.at(referenceTrajectoryIndices.at(T)));
                //    0.5. Sample K_t from 0.4.
                long Kt = Rcpp::sample(N, 1, false, lambdaWeightsStratified)[0] - 1;
                //Step 1:
                //Connect the "chosen" ancestor index to previous reference trajectory: A_{t - 1}^{K_t} = K_{t - 1}
                uRSIndices.at(Kt) = referenceTrajectoryIndices.at(T);
                //Step 2:
                //Calculation of empirical distribution function F_{t - 1}^N(i) is done and equal to computation of cumulative normalized weights stored in dRSWeightsCumulative.
                //Step 3: Generate ancestor indices and sssign to offspring indices {1,...,N}\{Kt}.
                std::vector<unsigned int> tmpIterator(N - 1);
                std::iota(tmpIterator.begin(), tmpIterator.end(), 0); // define appropriate tmpIterator as a sequence from 0 to N-1
                tmpIterator.erase(tmpIterator.begin() + Kt); // exclude the previoiusly sampled conditional index K_t
                long minimalJ = 0;
                // generate ancestor index and assign to offspring
                double tmpUnifRnd;
                for (int i : tmpIterator) {
                    tmpUnifRnd = unif_rand();
                    tmpUnifRnd += i;//mathematical 'i' starts at 1; here at 0
                    tmpUnifRnd /= N;

                    minimalJ += arma::conv_to<long>::from(arma::find(dRSWeightsCumulative.tail(N - minimalJ) > tmpUnifRnd, 1, "first"));
                    uRSIndices.at(i) = minimalJ;
                }
            }
        case ResampleType::SYSTEMATIC:
            {
                //Algorithm 5
                //Step 0:
                //Sample conditional index K_t from appropriate version of the "lambda" distribution i.e. the distribution over the stratum with \lambda(k_t|w_{t - 1}^{1:N}, k_{t - 1})=p_{t - 1}^{k_{t - 1}(k_t)/W_{t - 1}^{k_{t - 1}} in case of stratified resampling.
                //    0.1. Calculate normalized particle weights and cumulative normalized weights.
                dRSWeights = exp(pPopulation.GetLogWeight() - stableLogSumWeights(pPopulation.GetLogWeight()));
                arma::vec dRSWeightsCumulative = arma::cumsum(dRSWeights);
                //    0.2. Calculate strata boundaries cumulative weight components for min() and max() computation parts.
                arma::Col<double> strataBoundariesAll = arma::linspace(0.0, 1.0, N + 1);
                arma::Col<double> strataBoundariesUpper = strataBoundariesAll.tail(N);
                arma::Col<double> strataBoundariesLower = strataBoundariesAll.head(N);
                double minWeights = dRSWeightsCumulative.at(referenceTrajectoryIndices.at(T));
                double maxWeights = 0;
                if(T != 0){
                    maxWeights = dRSWeightsCumulative.at(referenceTrajectoryIndices.at(T - 1));
                }
                //    0.3. Calculate strata weights p_{t - 1}^{k_{t - 1}}(k_t).
                arma::Col<double> strataWeights(N);
                strataWeights.fill(0.0);
                for(int i = 0; i < N; ++i) {
                    strataWeights.at(i) = std::min(minWeights, strataBoundariesUpper.at(i)) - std::max(maxWeights, strataBoundariesLower.at(i));
                }
                //    0.4. Calculate \lambda(k_t|.) distribution.
                Rcpp::NumericVector lambdaWeightsStratified = Rcpp::wrap(strataWeights/dRSWeightsCumulative.at(referenceTrajectoryIndices.at(T)));
                //    0.5. Sample K_t from 0.4.
                long Kt = Rcpp::sample(N, 1, false, lambdaWeightsStratified)[0];
                //Step 1:
                //Connect the "chosen" ancestor index to previous reference trajectory: A_{t - 1}^{K_t} = K_{t - 1}
                uRSIndices.at(Kt) = referenceTrajectoryIndices.at(T);
                //Step 2:
                //Calculation of empirical distribution function F_{t - 1}^N(i) is done and equal to computation of cumulative normalized weights stored in dRSWeightsCumulative.
                //Step 3: Generate ancestor indices and sssign to offspring indices {1,...,N}\{Kt}.
                std::vector<unsigned int> tmpIterator(N - 1);
                std::iota(tmpIterator.begin(), tmpIterator.end(), 0); // define appropriate tmpIterator as a sequence from 0 to N-1
                tmpIterator.erase(tmpIterator.begin() + Kt); // exclude the previoiusly sampled conditional index K_t
                //precompute necessary uniform random variable before assignment
                double tmpVUpperBound = dRSWeightsCumulative.at(Kt);
                double tmpVLowerBound;
                if(Kt == 0) {
                    tmpVLowerBound = 0;
                } else {
                    tmpVLowerBound = dRSWeightsCumulative.at(Kt - 1);
                }
                double tmpV = R::runif(tmpVLowerBound, tmpVUpperBound);
                tmpV = N * tmpV - std::floor(N * tmpV);
                double tmpU;
                long minimalJ = 0;
                // generate ancestor index and assign to offspring
                for (int i : tmpIterator) {
                    tmpU = tmpV + i;//mathematical 'i' at 1; here it starts at 0
                    tmpU /= N;

                    minimalJ += arma::conv_to<long>::from(arma::find(dRSWeightsCumulative.tail(N - minimalJ) > tmpU, 1, "first"));
                    uRSIndices.at(i) = minimalJ;
                }
            }
        // case ResampleType::RESIDUAL:
        //     {
        //         Rcpp::Rcout << "resampling: residual resampling" << std::endl;
        //         //Algorithm 6
        //         //Step 0:
        //         //Declare/initialize container for implementation; compute normalized particle weights.
        //         //    0.1. Container setup
        //         int numDeterministicOffspring = 0; //Counts the number of deterministically assigned offspring.
        //         int expectedNumberOffspring = 0;
        //         arma::Col<double> dRSWeightsResidual(N);
        //         arma::Col<unsigned int> DsetCurrent;
        //         int CardDsetCurrent = 0;
        //         arma::Col<unsigned int> DsetKtMinus1;
        //         arma::Col<unsigned int> tmpIterator = arma::linspace(0, N - 1, N);
        //         //    0.2. Calculate normalized particle weights and cumulative normalized weights.
        //         dRSWeights = exp(pPopulation.GetLogWeight() - stableLogSumWeights(pPopulation.GetLogWeight()));
        //         //Step 1:
        //         //Assign deterministic offpring indices for each i={1,...,N}.
        //         for (int i : tmpIterator) {
        //             //Compute (integer part of) expected number of offspring
        //             expectedNumberOffspring = std::floor(N * dRSWeights);
        //             //Generate D_i set:
        //             if(expectedNumberOffspring > 0) {
        //                 //Set D_i={numDeterministicOffspring, numDeterministicOffspring + 1, ..., numDeterministicOffspring + expectedNumberOffspring}
        //                 // DsetCurrent
        //                 // Set Card(D_i)=length(D_i)
        //                 // Set numDeterministicOffspring += Card(D_i);
        //             } else {
        //                 //Set D_i=EMPTYSET
        //                 //Set Card(D_i)= 0;
        //             }
        //             //Deterministic component assignment
        //             //for(int j : D_i) {
        //             //    uRSIndices.at(j) = i;
        //             //}
        //             //Convenient way of handling the conditioning path
        //             //if(i == K_{t-1}) store D_{K_{t-1}};
        //             //Compute normalized residual weights \tilde{W}_{t-1}^i
        //         }
        //         //Step 2:
        //         //Sample K_t from appropriate lambda distribution
        //         long Kt = 0;
        //         // ...........
        //         referenceTrajectoryIndices(T + 1) = Kt;
        //         //Step 3:
        //         //Compute residual ancestor indices via sampling from categorical distribution
        //         // ...........
        //         //Step 4: If not done in 1., assign conditional ancestor index:
        //         // if(K_t not in D_{K_{t  1}})
        //         uRSIndices.at(Kt) = referenceTrajectoryIndices(T);
        //         // else already assigned under 1.
        //         break;
        //     }
        }
        //Perform the replication of the chosen.
        for(int i = 0; i < N ; ++i) {
            if(uRSIndices(i) != static_cast<unsigned int>(i)){
                pPopulation.SetValueN(pPopulation.GetValueN(static_cast<int>(uRSIndices(i))), i);
            }
        }
        //After conditional resampling is implemented: a final step is to set equal normalised weights.
        pPopulation.SetLogWeight(- log(static_cast<double>(N))*arma::ones(N));
    }
}
#endif