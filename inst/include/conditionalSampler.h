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
        //Resampling is still done in place but adding a random permutation to the conditional trajectory index makes this scheme valid.
        int uMultinomialCount;

        //First obtain a count of the number of children each particle has.
        switch(lMode) {
        case ResampleType::MULTINOMIAL:
        default:{
            //Sample from a suitable multinomial vector.
            sampler<Space,Params>::dRSWeights = exp(sampler<Space,Params>::pPopulation.GetLogWeight() - stableLogSumWeights(sampler<Space,Params>::pPopulation.GetLogWeight()));
            rmultinom(static_cast<int>(sampler<Space,Params>::N), sampler<Space,Params>::dRSWeights.memptr(), static_cast<int>(sampler<Space,Params>::N), sampler<Space,Params>::uRSCount.memptr());
            break;
        }

        case ResampleType::RESIDUAL:
            {
                //Procedure for residual sampling.
                sampler<Space,Params>::dRSWeights = exp(log(static_cast<double>(sampler<Space,Params>::N)) + sampler<Space,Params>::pPopulation.GetLogWeight() - stableLogSumWeights(sampler<Space,Params>::pPopulation.GetLogWeight()));
                sampler<Space,Params>::uRSIndices = arma::zeros<arma::Col<unsigned int> >(static_cast<int>(sampler<Space,Params>::N));
                for(int i = 0; i < sampler<Space,Params>::N; ++i)
                sampler<Space,Params>::uRSIndices(i) = static_cast<unsigned int>(floor(sampler<Space,Params>::dRSWeights(i)));
                sampler<Space,Params>::dRSWeights = sampler<Space,Params>::dRSWeights - sampler<Space,Params>::uRSIndices;
                sampler<Space,Params>::dRSWeights = sampler<Space,Params>::dRSWeights/sum(sampler<Space,Params>::dRSWeights);
                uMultinomialCount = sampler<Space,Params>::N - arma::sum(sampler<Space,Params>::uRSIndices);
                rmultinom(uMultinomialCount, sampler<Space,Params>::dRSWeights.memptr(), static_cast<int>(sampler<Space,Params>::N), sampler<Space,Params>::uRSCount.memptr());
                sampler<Space,Params>::uRSCount += arma::conv_to<arma::Col<int> >::from(sampler<Space,Params>::uRSIndices);
                break;
            }

        case ResampleType::STRATIFIED:
            {
                //Procedure for stratified sampling.
                int j = 0, k = 0;
                sampler<Space,Params>::uRSCount = arma::zeros<arma::Col<int> >(static_cast<int>(sampler<Space,Params>::N));
                //Generate a vector of cumulative weights.
                arma::vec dWeightCumulative = arma::cumsum(exp(sampler<Space,Params>::pPopulation.GetLogWeight() - stableLogSumWeights(sampler<Space,Params>::pPopulation.GetLogWeight())));
                //Generate a random number between 0 and 1/N.
                double dRand = R::runif(0,1.0 / static_cast<double>(sampler<Space,Params>::N));
                while(k < sampler<Space,Params>::N) {
                    while((dWeightCumulative(k) - dRand) > static_cast<double>(j)/static_cast<double>(sampler<Space,Params>::N) && j < sampler<Space,Params>::N) {
                        sampler<Space,Params>::uRSCount(k)++;
                        j++;
                        dRand = R::runif(0,1.0 / static_cast<double>(sampler<Space,Params>::N));
                    }
                    k++;
                }
                break;
            }

        case ResampleType::SYSTEMATIC:
            {
                //Procedure for stratified sampling but with a common RV for each stratum.
                int j = 0, k = 0;
                sampler<Space,Params>::uRSCount = arma::zeros<arma::Col<int> >(static_cast<int>(sampler<Space,Params>::N));
                //Generate a vector of cumulative weights.
                arma::vec dWeightCumulative = arma::cumsum(exp(sampler<Space,Params>::pPopulation.GetLogWeight() - stableLogSumWeights(sampler<Space,Params>::pPopulation.GetLogWeight())));
                //Generate a random number between 0 and 1/N.
                double dRand = R::runif(0,1.0 / static_cast<double>(sampler<Space,Params>::N));
                while(k < sampler<Space,Params>::N) {
                    while((dWeightCumulative(k) - dRand) > static_cast<double>(j)/static_cast<double>(sampler<Space,Params>::N) && j < sampler<Space,Params>::N) {
                        sampler<Space,Params>::uRSCount(k)++;
                        j++;
                    }
                    k++;
                }
                break;
            }
        }

        sampler<Space,Params>::uRSIndices = arma::zeros<arma::Col<unsigned int> >(static_cast<int>(sampler<Space,Params>::N));
        //Map count to indices to allow in-place resampling.
        for (int i=0, j=0; i < sampler<Space,Params>::N; ++i) {
            if (sampler<Space,Params>::uRSCount(i)>0) {
                sampler<Space,Params>::uRSIndices(i) = i;
                while (sampler<Space,Params>::uRSCount(i)>1) {
                    while (sampler<Space,Params>::uRSCount(j)>0) ++j; // find next free spot
                    sampler<Space,Params>::uRSIndices(j++) = i; // assign index
                    --sampler<Space,Params>::uRSCount(i); // decrement number of remaining offsprings
                }
            }
        }

        // Performs a random permutation of indices to break above deterministic assignment and enforce the exchangeability of indices necessary for valid conditional resampling.
        sampler<Space,Params>::uRSIndices = sampler<Space,Params>::uRSIndices.elem(arma::randperm(sampler<Space,Params>::N - 1));

        //Perform the replication of the chosen particle coordinates.
        for(int i = 0; i < sampler<Space,Params>::N; ++i) {
            if(sampler<Space,Params>::uRSIndices(i) != static_cast<unsigned int>(i)){
                sampler<Space,Params>::pPopulation.SetValueN(sampler<Space,Params>::pPopulation.GetValueN(static_cast<int>(sampler<Space,Params>::uRSIndices(i))), i);
            }
        }

        //Post-hoc randomization of the conditional index i.e. sampling uniformly on {1,...,N}.
        //    1. Generate uniform weights
        Rcpp::NumericVector tmpUniformWeights(sampler<Space,Params>::N, 1.0/sampler<Space,Params>::N);
        //    2. Randomly draw conditional index
        referenceTrajectoryIndices.at(sampler<Space,Params>::T + 1) = Rcpp::sample(sampler<Space,Params>::N - 1, 1, false, tmpUniformWeights)[0];
        //    3. Adjust the ancestor/re-sampling index to the randomly drawn
        sampler<Space,Params>::uRSIndices(referenceTrajectoryIndices.at(sampler<Space,Params>::T + 1)) = referenceTrajectoryIndices.at(sampler<Space,Params>::T);

        //After conditional resampling is implemented: a final step is to set equal normalised weights.
        sampler<Space,Params>::pPopulation.SetLogWeight(- log(static_cast<double>(sampler<Space,Params>::N))*arma::ones(sampler<Space,Params>::N));
    }
}
#endif