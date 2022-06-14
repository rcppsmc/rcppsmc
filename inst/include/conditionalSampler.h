// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
//
// conditionalSampler.h: Rcpp integration of SMC library -- conditional sampler class that provides support for conditional sequential Monte Carlo algorithms.
//
// Copyright (C) 2021 Adam Johansen, Dirk Eddelbuettel, Leah South, Ilya Zarubin
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
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with RcppSMC. If not, see <http://www.gnu.org/licenses/>.

//! \file
//! \brief Defines the conditionalSampler class derived from base sampler.
//!
//! This file defines the smc::conditionalSampler class which is used to implement entire particle systems for conditional sequential Monte Carlo. The class is derived from the base class 'sampler' but currently does not implement facilities of the base class that are not striclty required for conditional SMC (such as MCMC moves or adaptation). Conditional resampling for (some) common resampling mechanisms are implemented.

#ifndef __SMC_CONDITIONAL_SAMPLER_HH

#define __SMC_CONDITIONAL_SAMPLER_HH 1.0

#include "sampler.h"

namespace smc {
    //Pre-declare the derived template class itself is necessary so that operator overloading for template (friend) functions works properlery.
    template<class Space, class Params> class conditionalSampler;
    //The compiler has to know that the overloaded operator friend function below is a template, hence the pre-declaration of operator<<.
    template<class Space, class Params> std::ostream& operator<< (std::ostream&, const conditionalSampler<Space,Params>&);

    template<class Space, class Params = nullParams> class conditionalSampler:
    public sampler<Space,Params>
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

        ///A flag which tracks whether the ensemble was resampled during this iteration
        using sampler<Space,Params>::nResampled;

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

            int digitsPrint = 6;
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
            ///Get the number of digits used in ostream<< operator printing
            int GetDigitsPrint(void) const {return digitsPrint;}
            ///Set the number of digits used in ostream<< operator printing
            void SetDigitsPrint(int newDigitsPrint) {digitsPrint = newDigitsPrint;}
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
            ///Move the reference particle by setting the coordinate to the reference value and re-weighting.
            void MoveReferenceParticle(void);
            ///Dump a specified particle to the specified output stream in a human readable form
            std::ostream & StreamParticle(std::ostream & os, long n) const;
            ///Dump the entire particle set to the specified output stream in a human readable form
            std::ostream & StreamParticles(std::ostream & os, int digits) const;
            friend std::ostream& operator << <>(std::ostream &, const conditionalSampler<Space,Params>&);

            ///Throws exception when adaptation related members of the base sampler class are used:
            void SetAdaptMethods(adaptMethods<Space,Params>* adaptMethod) {throw SMC_EXCEPTION(CSMCX_USING_ADAPTATION, "Adaptation methods not currently supported for conditional sampler class.");}
            ///Throws exception when members related to MCMC moves of the base sampler class are used:
            void SetMcmcRepeats(adaptMethods<Space,Params>* adaptMethod) {throw SMC_EXCEPTION(CSMCX_USING_MCMC, "MCMC moves not currently supported for conditional sampler class.");}
            void OstreamMCMCRecordToStream(std::ostream &os) const {throw SMC_EXCEPTION(CSMCX_USING_MCMC, "MCMC moves not currently supported for conditional sampler class.");}
            int GetAccepted(void) const {throw SMC_EXCEPTION(CSMCX_USING_MCMC, "MCMC moves not currently supported for conditional sampler class.");}
            int GetMcmcRepeats(void) const {throw SMC_EXCEPTION(CSMCX_USING_MCMC, "MCMC moves not currently supported for conditional sampler class.");}
            int GetHistorymcmcRepeats(long n) {throw SMC_EXCEPTION(CSMCX_USING_MCMC, "MCMC moves not currently supported for conditional sampler class.");}
            ///Sets the number of MCMC repeats
            void SetMcmcRepeats(int reps) {throw SMC_EXCEPTION(CSMCX_USING_MCMC, "MCMC moves not currently supported for conditional sampler class.");}
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

        //Set the initial values and log weights of the particles
        std::vector<Space> InitVal(N);
        arma::vec InitWeights(N);
        pPopulation = population<Space>(InitVal,InitWeights);
        pMoves->DoInit(pPopulation,N,algParams);

        //Initialise the conditonal trajectory:
        //1. Sample uniformly initial period, T = 0, conditional index
        // referenceTrajectoryIndices.at(T) = floor(unif_rand()*static_cast<double>(N));
        referenceTrajectoryIndices.at(T) = 0;
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
            conditionalResample(rtResampleMode);
        }
        else {
            nResampled = 0;
            if(htHistoryMode == HistoryType::AL) {
                uRSIndices = arma::linspace<arma::Col<unsigned int>>(0, N - 1, N);
                //No resampling: set conditional index to previous iteration.
                referenceTrajectoryIndices.at(T + 1) = referenceTrajectoryIndices.at(T);
            }
        }

        //Normalise the weights
        pPopulation.SetLogWeight(pPopulation.GetLogWeight() - CalcLogNC());

        //Finally, the current particle set should be appended to the historical process.
        if(htHistoryMode != HistoryType::NONE){
            History.clear();
            historyelement<Space> histel;
            switch(htHistoryMode) {
            case HistoryType::RAM:
                histel.Set(N, pPopulation, historyflags(nResampled));
                break;
            case HistoryType::AL:
                histel.Set(N, pPopulation, historyflags(nResampled), uRSIndices);
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

        //Set the initial values and log weights of the particles
        std::vector<Space> InitVal(N);
        arma::vec InitWeights(N);
        pPopulation = population<Space>(InitVal,InitWeights);
        pMoves->DoInit(pPopulation,N,algParams);

        //Initialise the conditonal trajectory:
        //1. Sample uniformly initial period, T = 0, conditional index
        // referenceTrajectoryIndices.at(T) = floor(unif_rand()*static_cast<double>(N));
        referenceTrajectoryIndices.at(T) = 0;
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
            conditionalResample(rtResampleMode);
        }
        else {
            nResampled = 0;
            if(htHistoryMode == HistoryType::AL) {
                uRSIndices = arma::linspace<arma::Col<unsigned int>>(0, N - 1, N);
                //No resampling: set conditional index to previous iteration.
                referenceTrajectoryIndices.at(T + 1) = referenceTrajectoryIndices.at(T);
            }
        }

        //Normalise the weights
        pPopulation.SetLogWeight(pPopulation.GetLogWeight() - CalcLogNC());

        //Finally, the current particle set should be appended to the historical process.
        if(htHistoryMode != HistoryType::NONE){
            History.clear();
            historyelement<Space> histel;
            switch(htHistoryMode) {
            case HistoryType::RAM:
                histel.Set(N, pPopulation, historyflags(nResampled));
                break;
            case HistoryType::AL:
                histel.Set(N, pPopulation, historyflags(nResampled), uRSIndices);
                break;
            ///To avoid compiler warnings, HistoryType::NONE is handled
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
        //Move the particle set.
        MoveParticles();

        //Do add a conditional conditional move:
        //set reference particle coordinate at conditional value and re-weight.
        MoveReferenceParticle();

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
            conditionalResample(rtResampleMode);
        }
        else {
            nResampled = 0;
            if(htHistoryMode == HistoryType::AL) {
                uRSIndices = arma::linspace<arma::Col<unsigned int>>(0, N - 1, N);
                //No resampling: set conditional index to previous iteration.
                referenceTrajectoryIndices.at(T + 1) = referenceTrajectoryIndices.at(T);
            }
        }

        //Normalise the weights
        pPopulation.SetLogWeight(pPopulation.GetLogWeight() - CalcLogNC());

        //Finally, the current particle set should be appended to the historical process.
        if(htHistoryMode != HistoryType::NONE){
            historyelement<Space> histel;
            switch(htHistoryMode) {
            case HistoryType::RAM:
                histel.Set(N, pPopulation, historyflags(nResampled));
                break;
            case HistoryType::AL:
                histel.Set(N, pPopulation, historyflags(nResampled), uRSIndices);
                break;
            //To avoid compiler warnings, HistoryType::NONE is handled
            case HistoryType::NONE:
                break;
            }
            History.push_back(histel);
        }
        //Increment the evolution time.
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
        // {
        Iterate();
        // Rcpp::Rcout << "Not failing at iteration t: " << GetTime() << std::endl;
        // }
    }
    template <class Space, class Params>
    void conditionalSampler<Space, Params>::conditionalResample(ResampleType::Enum lMode)
    {
        //Conditional resampling performed following the algorithms outlined in Appendix C of the paper "...".

        //Initialize container for ancestor indices.
        uRSIndices = arma::zeros<arma::Col<unsigned int> >(static_cast<int>(N));

        switch(lMode) {
        case ResampleType::MULTINOMIAL:
        default:
            {
                //Algorithm 3
                //Step 0:
                //Sample conditional index K_t from appropriate version of the "lambda" distribution i.e. uniformly on {1,...,N} in case of Multinmial resampling.
                // long Kt = floor(unif_rand()*static_cast<double>(N)); // sample lamba(k_{t}|w_{t-1}, k_{t-1})=1/N
                //Set conditional index K_t to zero implementing a sepcial case that works for conditional multinomial resampling (saving one draw compared to the implementation above).
                long Kt = 0;
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
                //Sample conditional index K_t from appropriate version of the "lambda" distribution i.e. the distribution over the stratum with \lambda(k_t|w_{t - 1}^{1:N}, k_{t - 1})=p_{t - 1}^{k_{t - 1}}(k_t)/W_{t - 1}^{k_{t - 1}} in case of stratified resampling.
                //0.0 Housekeeping for step 0.1 - 0.5
                //     k_{t - 1}
                int kt1 = referenceTrajectoryIndices.at(T);
                //    sum_{k=1}^{k_{t-1}}\hat{w}_{t-1}^k
                double minWeights = 0;
                //    sum_{k=1}^{k_{t-1} - 1}\hat{w}_{t-1}^k
                double maxWeights = 0;
                //    For N different strata (N being the number of particles)
                //    all strata boundaries: 0/N,...,c/N,...,N/N for c = 0,...,N
                arma::Col<double> stBoundsAll = arma::linspace(0.0, 1.0, N + 1);
                ///    upper sequence part: 1/N,2/N,...c/N,...,1
                arma::Col<double> stBoundsUpper = stBoundsAll.tail(N);
                ///   lower sequence part: 0/N,1/N,...c/N,...,(N-1)/N
                arma::Col<double> stBoundsLower = stBoundsAll.head(N);
                ///    particle k_{t-1} weight for all strata k_t=1,...,N
                arma::Col<double> strataWeights(N);
                strataWeights.fill(0.0);
                ///    weights of \lambda(k_t|.) distribution for k_t sampling
                Rcpp::NumericVector lambdaWeightsStratified(N);
                //0.1. Calculate normalized particle weights and cumulative normalized weights.
                dRSWeights = exp(pPopulation.GetLogWeight() - stableLogSumWeights(pPopulation.GetLogWeight()));
                arma::vec dRSWeightsCumulative = arma::cumsum(dRSWeights);
                //    0.2. Calculate min() and max() parts of strata weights
                minWeights = dRSWeightsCumulative.at(kt1);
                if(kt1 > 0){
                    maxWeights = dRSWeightsCumulative.at(kt1 - 1);
                } else {
                    maxWeights = 0;
                }
                //   0.3. Calculate weight of particle k_{t-1} in all strata
                //   k_t=1,...,N i.e. p_{t-1}^{k_{t-1}}(k_t)
                for(int i = 0; i < N; ++i) {
                    strataWeights.at(i) = std::max(std::min(minWeights, stBoundsUpper.at(i)) - std::max(maxWeights, stBoundsLower.at(i)), 0.0);
                }
                //    0.4. Calculate \lambda(k_t|.) distribution.
                lambdaWeightsStratified = Rcpp::wrap(strataWeights/dRSWeights.at(kt1));
                //    0.5. Sample K_t from 0.4.
                long Kt = Rcpp::sample(N, 1, false, lambdaWeightsStratified)[0] - 1;
                //Step 1:
                //Connect the "chosen" ancestor index to previous reference trajectory: A_{t - 1}^{K_t} = K_{t - 1}
                uRSIndices.at(Kt) = kt1;
                //Step 2:
                //Calculation of empirical distribution function F_{t - 1}^N(i) is done and equal to computation of cumulative normalized weights stored in dRSWeightsCumulative.
                //Step 3: Generate ancestor indices and assign to offspring indices {1,...,N}\{Kt}.
                std::vector<unsigned int> tmpIterator(N);
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
                break;
            }
        case ResampleType::SYSTEMATIC:
            {
                //Algorithm 5
                //Step 0:
                //Sample conditional index K_t from appropriate version of the "lambda" distribution i.e. the distribution over the stratum with \lambda(k_t|w_{t - 1}^{1:N}, k_{t - 1})=p_{t - 1}^{k_{t - 1}(k_t)/W_{t - 1}^{k_{t - 1}} (same as stratified resampling).
                //0.0 Housekeeping for step 0.1 - 0.5
                //      k_{t - 1}
                int kt1 = referenceTrajectoryIndices.at(T);
                //    sum_{k=1}^{k_{t-1}}\hat{w}_{t-1}^k
                double minWeights = 0;
                //    sum_{k=1}^{k_{t-1} - 1}\hat{w}_{t-1}^k
                double maxWeights = 0;
                //    For N different strata (N being the number of particles)
                //    all strata boundaries: 0/N,...,c/N,...,N/N for c = 0,...,N
                arma::Col<double> stBoundsAll = arma::linspace(0.0, 1.0, N + 1);
                ///    upper sequence part: 1/N,2/N,...c/N,...,1
                arma::Col<double> stBoundsUpper = stBoundsAll.tail(N);
                ///   lower sequence part: 0/N,1/N,...c/N,...,(N-1)/N
                arma::Col<double> stBoundsLower = stBoundsAll.head(N);
                ///    particle k_{t-1} weight for all strata k_t=1,...,N
                arma::Col<double> strataWeights(N);
                strataWeights.fill(0.0);
                ///    weights of \lambda(k_t|.) distribution for k_t sampling
                Rcpp::NumericVector lambdaWeightsStratified(N);
                //    0.1. Calculate normalized particle weights and cumulative normalized weights.
                dRSWeights = exp(pPopulation.GetLogWeight() - stableLogSumWeights(pPopulation.GetLogWeight()));
                arma::vec dRSWeightsCumulative = arma::cumsum(dRSWeights);
                //    0.2. Calculate min() and max() parts of strata weights
                minWeights = dRSWeightsCumulative.at(kt1);
                if(kt1 > 0){
                    maxWeights = dRSWeightsCumulative.at(kt1 - 1);
                } else {
                    maxWeights = 0;
                }
                //   0.3. Calculate weight of particle k_{t-1} in all strata
                //   k_t=1,...,N i.e. p_{t-1}^{k_{t-1}}(k_t)
                for(int i = 0; i < N; ++i) {
                    strataWeights.at(i) = std::max(std::min(minWeights, stBoundsUpper.at(i)) - std::max(maxWeights, stBoundsLower.at(i)), 0.0);
                }
                //    0.4. Calculate \lambda(k_t|.) distribution.
                lambdaWeightsStratified = Rcpp::wrap(strataWeights/dRSWeights.at(kt1));
                //    0.5. Sample K_t from 0.4.
                long Kt = Rcpp::sample(N, 1, false, lambdaWeightsStratified)[0] - 1;
                //Step 1:
                //Connect the "chosen" ancestor index to previous reference trajectory: A_{t - 1}^{K_t} = K_{t - 1}
                uRSIndices.at(Kt) = referenceTrajectoryIndices.at(T);
                //Step 2:
                //Calculation of empirical distribution function F_{t - 1}^N(i) is done and equal to computation of cumulative normalized weights stored in dRSWeightsCumulative.
                //Step 3: Generate ancestor indices and sssign to offspring indices {1,...,N}\{Kt}.
                std::vector<unsigned int> tmpIterator(N);
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
                break;
            }
        case ResampleType::RESIDUAL:
            {
                //Algorithm 6
                //Step 0:
                //Initialize container for implementation; pre-compute normalized/unnormalized residual weights; perform deterministic component assignment and compute other quantities
                //    0.1. Container setup
                int l = 0; //Counts the number of deterministically assigned offspring.
                double kt1 = referenceTrajectoryIndices.at(T); // k_{t-1}
                int expectedNumberOffspring = 0; //Integer part of N*W_{t - 1}.
                arma::Col<double> dRSWeightsResidual(N);
                arma::Col<double> dRSWeightsResidualNormalized(N);
                arma::Col<unsigned int> DsetCurrent; //D_i
                int CardDsetCurrent = 0;//|D_i|
                arma::Col<unsigned int> DsetKtMinus1; //D_{K_{t - 1}}
                arma::Col<unsigned int> DsetDeterministic; //The set of indices that are sampled deterministically.
                arma::Col<unsigned int> DsetStochastic; //The set of indices that are sampled stochastically.
                arma::Col<unsigned int> tmpIterator = arma::linspace<arma::Col<unsigned int> >(0, N - 1, N);
                //    0.2. Calculate normalized particle weights.
                dRSWeights = exp(pPopulation.GetLogWeight() - stableLogSumWeights(pPopulation.GetLogWeight()));

                for (int i : tmpIterator) {
                    //Compute (integer part of) expected number of offspring
                    expectedNumberOffspring = std::floor(N * dRSWeights.at(i));
                    //Generate D_i set:
                    if(expectedNumberOffspring > 0) {
                        //Set D_i={l + 1, ..., l + expectedNumberOffspring}
                        DsetCurrent = arma::linspace<arma::Col<unsigned int> >(l, l + expectedNumberOffspring - 1, expectedNumberOffspring);
                        // Set Card(D_i)=length(D_i) = floor(N*W_{t-1});
                        CardDsetCurrent = expectedNumberOffspring;
                        // Set l += Card(D_i);
                        l += CardDsetCurrent;
                    } else {
                        //Set D_i=EMPTYSET
                        //Set Card(D_i)= 0;
                        DsetCurrent.reset();
                        CardDsetCurrent = 0;
                    }
                    //Assign deterministic offpring indices for each i={1,...,N}
                    for (int j : DsetCurrent) {
                        uRSIndices.at(j) = i;
                    }
                    //Compute unnormalized residual weights \hat{W}_{t-1}^i.
                    dRSWeightsResidual.at(i) = dRSWeights.at(i) * N - static_cast<double>(CardDsetCurrent)/N;
                    //Check if convenient way of handling the conditioning path
                    //is possible: if(i == K_{t-1}) store D_{K_{t-1}}
                    //and check K_t \notin D_{K_{t-1}} later!
                    if(i == static_cast<int>(referenceTrajectoryIndices.at(T))) {
                        DsetKtMinus1 = DsetCurrent;
                    }
                }
                //Define deterministic and stochastic sets for index assignment
                DsetDeterministic = arma::linspace<arma::Col<unsigned int> >(0, l - 1, l);
                if (l < N) {
                    DsetStochastic = arma::linspace<arma::Col<unsigned int> >(l, N - 1, N - l);
                } else {
                    DsetStochastic = arma::linspace<arma::Col<unsigned int> >(N - 1, N - 1, 1);
                }
                //Compute normalized residual weights \tilde{W}_{t-1}^i.
                dRSWeightsResidualNormalized = dRSWeightsResidual/arma::sum(dRSWeightsResidual);
                //Step 1:
                //Sample K_t from appropriate lambda distribution
                //    1.1 Initializes container
                arma::Col<double> lambdaProbs(N);
                lambdaProbs.zeros();
                double prob1 = 1.0/(N * dRSWeights.at(kt1));
                double prob2 = dRSWeightsResidual.at(kt1)/(N * dRSWeights.at(kt1) * (N - l));
                //Defines the complement set of DsetKtMinus1 relative to {1,...,N} i.e. {1,...,N}\DsetKtMinus1
                //    1.2 Fills a vector with probabilities to sample from \lambda(k_t|k_{t - 1}, x_{t - 1}^{1:N})
                lambdaProbs.elem(DsetKtMinus1).fill(prob1);
                lambdaProbs.elem(DsetStochastic).fill(prob2);
                //    1.3 Sample Kt from eq. 4 on p.8 of the WP
                long Kt = Rcpp::sample(N, 1, false, Rcpp::wrap(lambdaProbs))[0] - 1;
                referenceTrajectoryIndices(T + 1) = Kt;
                //    1.4 If Kt \notin D_{K_{t-1}}
                if (DsetKtMinus1.size() == 0) {         //the set is empty
                    uRSIndices.at(Kt) = kt1;
                } else if (Kt < arma::as_scalar(DsetKtMinus1.head(1)) && Kt > arma::as_scalar(DsetKtMinus1.tail(1))) {//Kt notin (ordered) set
                    uRSIndices.at(Kt) = kt1;
                }
                //Step 2:
                //Compute residual ancestor indices via sampling from categorical distribution
                if (l < N) {
                    for(int k : DsetStochastic) {
                        if (k == Kt) continue; // skip since A^{K_t} already set
                        uRSIndices.at(k) = Rcpp::sample(N, 1, false, Rcpp::wrap(dRSWeightsResidualNormalized))[0] - 1;
                    }
                }
                break;
            }
        }
        //Copy values of pre-resampling particle set.
        std::vector<Space> populationValueCopy(pPopulation.GetValue().begin(), pPopulation.GetValue().end());
        //Perform the replication of the chosen.
        for(int i = 0; i < N ; ++i) {
            if(uRSIndices(i) != static_cast<unsigned int>(i)){
                pPopulation.SetValueN(populationValueCopy[static_cast<int>(uRSIndices(i))], i);
            }
        }
        //After conditional resampling is implemented: a final step is to set equal normalised weights.
        pPopulation.SetLogWeight(- log(static_cast<double>(N))*arma::ones(N));
    }
    ///Move the reference particle by setting the coordinate to the reference value and re-weighting.
    template<class Space, class Params>
    void conditionalSampler<Space,Params>::MoveReferenceParticle(void)
    {
        pMoves->DoConditionalMove(T + 1,pPopulation,referenceTrajectory[T + 1],referenceTrajectoryIndices.at(T + 1),algParams);
    }
    /// Produce a human-readable display of the current particle values and log weights.
    ///
    /// \param os The output stream to which the display should be made.
    template <class Space, class Params>
    std::ostream & conditionalSampler<Space,Params>::StreamParticles(std::ostream & os, int digits) const
    {
        Space val;
        double unw;
        double nw;
        int roundDigits = std::pow(10, digits);
        for(int i = 0; i < pPopulation.GetNumber() - 1; ++i){
            val = pPopulation.GetValueN(i);

            unw = pPopulation.GetLogWeightN(i);
            unw = static_cast<int>(unw * roundDigits >= 0 ? unw * roundDigits + 0.5 : unw * roundDigits - 0.5);
            unw = static_cast<double>(unw)/roundDigits;

            nw  = pPopulation.GetWeightN(i);

            std::string offsetPositiveVal = "";
            if(val>0.0) offsetPositiveVal = " ";
            std::string offsetPositiveUNW = "";
            if(unw>0.0) offsetPositiveUNW = " ";

            os << "Particle value: " << offsetPositiveVal << val << "  unnormalized weight: " << offsetPositiveUNW << unw << "  normalized weight: " << nw << std::endl;
        }
        return os;
    }
    /// Produce a human-readable display of the current n'th particle value and log weight.
    ///
    /// \param os The output stream to which the display should be made.
    /// \param n The index of the particle of interest
    template <class Space, class Params>
    std::ostream & conditionalSampler<Space,Params>::StreamParticle(std::ostream & os, long n) const
    {
        os << pPopulation.GetValueN(n) << "," << pPopulation.GetWeightN(n) << std::endl;
        return os;
    }
    /// Produce a human-readable display of the state of an smc::sampler class using the stream operator.
    /// \param os The output stream to which the display should be made.
    /// \param s  The sampler which is to be displayed.
    template <class Space, class Params>
    std::ostream & operator<<(std::ostream & os, const conditionalSampler<Space,Params> & CS)
    {
        os << "Sampler Configuration:" << std::endl;
        os << "======================" << std::endl;
        os << "Evolution Time:    " << CS.GetTime() << std::endl;
        os << "Particle Set Size: " << CS.GetNumber() << std::endl;
        os << "Effective Sample Size: " << CS.GetESS() << std::endl;
        os << std::endl;
        os << "Particle Set: " << std::endl;
        CS.StreamParticles(os, CS.GetDigitsPrint());
        os << std::endl;
        return os;
    }

}
#endif
