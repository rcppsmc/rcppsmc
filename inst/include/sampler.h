// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
//
// sampler.h: Rcpp integration of SMC library -- sampler object
//
// Copyright (C) 2008 - 2009  Adam Johansen
// Copyright (C) 2017 - 2020  Adam Johansen, Dirk Eddelbuettel and Leah South
// Copyright (C) 2021         Adam Johansen, Dirk Eddelbuettel, Leah South and Ilya Zarubin
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
//! \brief Defines the overall sampler object.
//!
//! This file defines the smc::sampler class which is used to implement entire particle systems.

#ifndef __SMC_SAMPLER_HH

#define __SMC_SAMPLER_HH 1.0

#include <algorithm>
#include <cstdlib>
#include <iostream>

#include "population.h"
#include "history.h"
#include "moveset.h"
#include "adaptMethods.h"
#include "smc-exception.h"
#include "staticModelAdapt.h"

///Specifiers for various resampling algorithms:
namespace ResampleType
{
    enum Enum { MULTINOMIAL = 0,
        RESIDUAL,
        STRATIFIED,
        SYSTEMATIC };
}

///Specifiers for various path sampling methods:
namespace PathSamplingType
{
    enum Enum { TRAPEZOID2 = 0,
        TRAPEZOID1,
        RECTANGLE};
}


///Storage types for the history of the particle system.
namespace HistoryType
{
    enum Enum { NONE = 0,
        RAM,
        AL};
}

namespace smc {

    /// An empty class for use when additional algorithm parameters are not required.
    class nullParams{};

    /// A template class for an interacting particle system suitable for SMC sampling
    template <class Space, class Params = nullParams>
    class sampler
    {
    protected:
        ///Number of particles in the system.
        long N;
        ///The current evolution time of the system.
        long T;

        ///The resampling mode which is to be employed.
        ResampleType::Enum rtResampleMode;
        ///The effective sample size at which resampling should be used.
        double dResampleThreshold;

        ///Structure used internally for resampling.
        arma::vec dRSWeights;
        ///Structure used internally for resampling.
        arma::Col<int> uRSCount;
        ///Structure used internally for resampling.
        arma::Col<unsigned int> uRSIndices;

        ///The particles within the system.
        population<Space> pPopulation;
        ///The set of moves available.
        moveset<Space,Params>* pMoves;
        ///A flag to track whether the moveset object needs to be included in this destructor.
        bool movesetBelong;
        /// The additional algorithm parameters.
        Params algParams;
        /// An object for adapting additional algorithm parameters
        adaptMethods<Space,Params>* pAdapt;
        ///A flag to track whether the adaptation object needs to be included in this destructor.
        bool adaptBelong;

        ///The number of MCMC moves which have been accepted during this iteration
        int nAccepted;
        ///A flag which tracks whether the ensemble was resampled during this iteration
        int nResampled;
        ///The number of MCMC repeats to be performed. The default is 1 if an MCMC step is supplied.
        int nRepeats;
        ///The proportion of accepted MCMC proposals in the most recent MCMC step, with
        /// a default of -1 if no MCMC steps have been performed.
        double acceptProb;

        ///An estimate of the log normalising constant ratio over the entire path.
        double dlogNCPath;
        ///An estimate of the log normalising constant ratio over the last step.
        double dlogNCIt;

        ///A mode flag which indicates whether historical information is stored
        HistoryType::Enum htHistoryMode;
        ///The historical process associated with the particle system.
        std::vector<historyelement<Space> > History;

        ///Helper function for copy constructor and assignment overloading
        void _copy(const sampler<Space,Params> & sFrom);

    public:
        ///Create an particle system containing lSize uninitialised particles with the specified mode.
        sampler(long lSize, HistoryType::Enum htHistoryMode);
        ///Create an particle system containing lSize uninitialised particles with the specified mode.
        sampler(long lSize, HistoryType::Enum htHistoryMode, moveset<Space,Params>*);
        ///Dispose of a sampler.
        ~sampler();
        ///Copy constructor
        sampler(const sampler<Space,Params> & sFrom);
        ///Assignment overloading
        sampler<Space,Params> & operator=(const sampler<Space,Params> & sFrom);
        ///Calculates and Returns the Effective Sample Size.
        double GetESS(void) const;
        ///Returns the number of accepted proposals from the most recent MCMC iteration
        int GetAccepted(void) const {return nAccepted;}
        ///Returns a flag for whether the ensemble was resampled during the most recent iteration
        int GetResampled(void) const {return nResampled;}
        ///Returns the number of MCMC repeats used in the most recent iteration
        int GetMcmcRepeats(void) const {return nRepeats;}
        ///Returns the History of the particle system
        const std::vector<historyelement<Space> > & GetHistory(void) const { return History; }
        ///Returns the number of particles within the system.
        long GetNumber(void) const {return N;}
        ///Returns the number of evolution times stored in the history.
        long GetHistoryLength(void) const {return History.size();}
        ///Returns the current particle set stored in the history.
        population<Space> GetHistoryPopulation(long n) const {return History[n].GetValues();}
        ///Returns a reference to the particle set stored in the history.
        const population<Space> & GetHistoryPopulationRefs(long n) const {return History[n].GetValues();}
        ///Returns the history flags
        historyflags GetHistoryFlags(long n) const {return History[n].GetFlags();}
        /// Returns the Effective Sample Size of the specified particle generation.
        double GetESS(long n) const {return History[n].GetESS();}
        /// Returns a vector of population objects re-ordered to match ancestral lines of the particle system
        std::vector<population<Space>> GetAPop(void) const;
        /// Returns indices for ancestral line of particle 'n'
        arma::Col<unsigned int> GetALineInd(long n) const;
        /// Returns ancestral line of particle 'n'
        std::vector<Space> GetALineSpace(long n) const;
        /// Returns current iteration resampled indices
        arma::Col<unsigned int> GetuRSIndices(void) const {return uRSIndices;}
        /// Returns current iteration resampled index number 'n'
        int GetuRSIndex(long n) const {return uRSIndices(n);}
        ///Returns the history number of MCMC iterations performed during this iteration.
        int GetHistorymcmcRepeats(long n) {return History[n].mcmcRepeats();}
        ///Returns the additional algorithm parameters.
        const Params & GetAlgParams(void) const {return algParams;}
        ///Return the value of particle n
        const Space &  GetParticleValueN(long n) const { return pPopulation.GetValueN(n); }
        ///Return the logarithmic unnormalized weight of particle n
        double GetParticleLogWeightN(long n) const { return pPopulation.GetLogWeightN(n); }
        ///Return the unnormalized weight of particle n
        double GetParticleWeightN(long n) const { return pPopulation.GetWeightN(n); }
        ///Return the unnormalized particle weights
        arma::vec GetParticleWeight(void) const { return pPopulation.GetWeight(); }
        ///Returns the current evolution time of the system.
        long GetTime(void) const {return T;}
        ///Returns the current estimate of the log normalising constant ratio over the entire path
        double GetLogNCPath(void) const { return dlogNCPath; }
        ///Returns the current estimate of the log normalising constant ratio over the last step
        double GetLogNCStep(void) const { return dlogNCIt; }
        ///Returns the current estimate of the normalising constant ratio over the entire path
        double GetNCPath(void) const { return exp(dlogNCPath); }
        ///Returns the current estimate of the normalising constant ratio over the last step
        double GetNCStep(void) const { return exp(dlogNCIt); }
        ///Initialise the sampler and its constituent particles.
        void Initialise(void);
        ///Integrate the supplied function with respect to the current particle set.
        double Integrate(double(*pIntegrand)(const Space &,void*), void* pAuxiliary);
        ///Integrate the supplied function over the path using the supplied width function and integration method.
        double IntegratePathSampling(PathSamplingType::Enum, double (*pIntegrand)(long,const Space &, void*), double (*pWidth)(long,void*), void* pAuxiliary);
        ///Integrate the supplied function over the path using the supplied width function and the default integration method (the corrected trapezoid rule).
        double IntegratePathSampling(double (*pIntegrand)(long,const Space &,void*), double (*pWidth)(long,void*), void* pAuxiliary) {return IntegratePathSampling(PathSamplingType::TRAPEZOID2, pIntegrand, pWidth, pAuxiliary);}
        ///Perform one iteration of the simulation algorithm.
        void Iterate(void);
        ///Cancel one iteration of the simulation algorithm.
        void IterateBack(void);
        ///Perform one iteration of the simulation algorithm and return the resulting ess
        double IterateEss(void);
        ///Perform iterations until the specified evolution time is reached
        void IterateUntil(long lTerminate);
        ///Move the particle set by proposing and applying an appropriate move to each particle.
        void MoveParticles(void);
        ///Resample the particle set using the specified resampling scheme.
        void Resample(ResampleType::Enum lMode);
        ///Sets the entire moveset to the one which is supplied
        void SetMoveSet(moveset<Space,Params>* pNewMoves) {delete pMoves; pMoves = pNewMoves; movesetBelong = 0;}
        ///Sets the entire moveset to the one which is supplied - backwards compatibility
        void SetMoveSet(moveset<Space,Params> & NewMoves) {delete pMoves; pMoves = &NewMoves; movesetBelong = 0;}
        ///Set Resampling Parameters
        void SetResampleParams(ResampleType::Enum rtMode, double dThreshold);
        ///Set additional algorithm parameters
        void SetAlgParam(Params parameters) {algParams = parameters;}
        ///Set the methods to adapt the additional algorithm parameters
        void SetAdaptMethods(adaptMethods<Space,Params>* adaptMethod) {delete pAdapt; pAdapt = adaptMethod; adaptBelong = 0;}
        ///Sets the number of MCMC repeats
        void SetMcmcRepeats(int reps) {nRepeats = reps;}
        ///Dump a specified particle to the specified output stream in a human readable form
        std::ostream & StreamParticle(std::ostream & os, long n) const;
        ///Dump the entire particle set to the specified output stream in a human readable form
        std::ostream & StreamParticles(std::ostream & os) const;
        ///Output a vector indicating the number of accepted MCMC moves at each time instance
        void OstreamMCMCRecordToStream(std::ostream &os) const;
        ///Output a 0-1 value vector indicating the times at which resampling occured to an output stream
        void OstreamResamplingRecordToStream(std::ostream &os) const;

    protected:
        ///Returns the crude normalising constant ratio estimate implied by the weights.
        double CalcLogNC(void) const {return stableLogSumWeights(pPopulation.GetLogWeight());}
    };

    /// The constructor prepares a sampler for use but does not assign any moves to the moveset, initialise the particles
    /// or otherwise perform any sampling related tasks. Its main function is to allocate a region of memory in which to
    /// store the particle set.
    ///
    /// \param lSize The number of particles present in the ensemble (at time 0 if this is a variable quantity)
    /// \param htHM The history mode to use: set this to HistoryType::RAM to store the whole history of the system and HistoryType::NONE to avoid doing so.
    /// \param pNewMoves Pointer to a moveset.
    /// \tparam Space The class used to represent a point in the sample space.
    /// \tparam Params (optional) The class used for any additional parameters.
    template <class Space, class Params>
    sampler<Space,Params>::sampler(long lSize, HistoryType::Enum htHM, moveset<Space,Params>* pNewMoves)
    {
        N = lSize;
        uRSCount = arma::zeros<arma::Col<int> >(static_cast<int>(N));

        //Some workable defaults.
        htHistoryMode = htHM;
        rtResampleMode = ResampleType::STRATIFIED;;
        dResampleThreshold = 0.5 * N;

        //Create an empty adaptation object by default
        pAdapt = new adaptMethods<Space,Params>;
        adaptBelong = 1;
        nRepeats = 1;

		// Adding the moveset
		pMoves = pNewMoves;
		movesetBelong = 0;
    }

    /// The constructor prepares a sampler for use but does not assign any moves to the moveset, initialise the particles
    /// or otherwise perform any sampling related tasks. Its main function is to allocate a region of memory in which to
    /// store the particle set.
    ///
    /// \param lSize The number of particles present in the ensemble (at time 0 if this is a variable quantity)
    /// \param htHM The history mode to use: set this to HistoryType::RAM to store the whole history of the system and HistoryType::NONE to avoid doing so.
    /// \tparam Space The class used to represent a point in the sample space.
    /// \tparam Params (optional) The class used for any additional parameters.
    template <class Space, class Params>
    sampler<Space,Params>::sampler(long lSize, HistoryType::Enum htHM)
    {
        N = lSize;
        uRSCount = arma::zeros<arma::Col<int> >(static_cast<int>(N));

        //Some workable defaults.
        htHistoryMode = htHM;
        rtResampleMode = ResampleType::STRATIFIED;;
        dResampleThreshold = 0.5 * N;

        //Create an empty adaptation object by default
        pAdapt = new adaptMethods<Space,Params>;
        adaptBelong = 1;
        nRepeats = 1;

		// Creating a default moveset
        pMoves = new moveset<Space,Params>;
        movesetBelong = 1;
    }

    template <class Space, class Params>
    sampler<Space,Params>::~sampler()
    {
        if(adaptBelong)
        delete pAdapt;
        if(movesetBelong)
        delete pMoves;
    }

    // deep-copy, to be used both for copy constructor and assignment overload.
    template <class Space, class Params>
    void sampler<Space, Params>::_copy(const sampler<Space,Params> & sFrom)
    {
        ///Number of particles in the system.
        N = sFrom.N;
        ///The current evolution time of the system.
        T = sFrom.T;

        ///The resampling mode which is to be employed.
        rtResampleMode = sFrom.rtResampleMode;
        ///The effective sample size at which resampling should be used.
        dResampleThreshold = sFrom.dResampleThreshold;

        ///Structure used internally for resampling.
        dRSWeights = sFrom.dRSWeights;
        ///Structure used internally for resampling.
        uRSCount = sFrom.uRSCount;
        ///Structure used internally for resampling.
        uRSIndices = sFrom.uRSIndices;

        ///The particles within the system.
        pPopulation = sFrom.pPopulation;
        ///The set of moves.
        if(sFrom.movesetBelong) {
            // this can only happen if the default moveset was used,
            // i.e., no call to SetMoveSet
            pMoves = new moveset<Space,Params>;
            movesetBelong = 1;
        } else {
            // this can only happen if SetMoveSet was called,
            // i.e., pMoves points to an external object which should not be deleted with sampler
            pMoves = sFrom.pMoves;
            movesetBelong = 0;
        }
        // /// Moveset object
        // pMoves = sFrom.pMoves;
        // ///A flag to track whether the moveset object needs to be included in this destructor.
        // movesetBelong = sFrom.movesetBelong sFrom.movesetBelong;


        /// The additional algorithm parameters.
        algParams = sFrom.algParams;
        if(sFrom.adaptBelong) {
            // this can only happen if the default adaptMethods was used,
            // i.e., no call to SetAdaptMethods
            pAdapt = new adaptMethods<Space,Params>;
            adaptBelong = 1;
        } else {
            // this can only happen if SetAdaptMethods was called,
            // i.e., pAdapt points to an external object which should not be deleted with sampler
            pAdapt = sFrom.pAdapt;
            adaptBelong = 0;
        }
        // /// An object for adapting additional algorithm parameters
        // pAdapt = sFrom.pAdapt;
        // ///A flag to track whether the adaptation object needs to be included in this destructor.
        // adaptBelong = sFrom.adaptBelong sFrom.adaptBelong;

        ///The number of MCMC moves which have been accepted during this iteration
        nAccepted = sFrom.nAccepted;
        ///A flag which tracks whether the ensemble was resampled during this iteration
        nResampled = sFrom.nResampled;
        ///The number of MCMC repeats to be performed. The default is 1 if an MCMC step is supplied.
        nRepeats = sFrom.nRepeats;
        ///The proportion of accepted MCMC proposals in the most recent MCMC step, with
        /// a default of -1 if no MCMC steps have been performed.
        acceptProb = sFrom.acceptProb;

        ///An estimate of the log normalising constant ratio over the entire path.
        dlogNCPath = sFrom.dlogNCPath;
        ///An estimate of the log normalising constant ratio over the last step.
        dlogNCIt = sFrom.dlogNCIt;

        ///A mode flag which indicates whether historical information is stored
        htHistoryMode = sFrom.htHistoryMode;
        ///The historical process associated with the particle system.
        History = sFrom.History;
    }

    template <class Space, class Params>
    sampler<Space,Params>::sampler(const sampler<Space,Params> & sFrom)
    {
        _copy(sFrom);
    }

    template <class Space, class Params>
    sampler<Space,Params> & sampler<Space,Params>::operator=(const sampler<Space,Params> & sFrom)
    {
        if (this != &sFrom) {
          if (adaptBelong) {
            delete pAdapt;
          }
          if (movesetBelong) {
            delete pMoves;
          }
          _copy(sFrom);
        }
        return *this;
    }

    template <class Space, class Params>
    double sampler<Space,Params>::GetESS(void) const
    {
        return expl(2*stableLogSumWeights(pPopulation.GetLogWeight())-stableLogSumWeights(2.0*pPopulation.GetLogWeight()));
    }
    template <class Space, class Params>
    std::vector<population<Space>> sampler<Space,Params>::GetAPop(void) const
    {
        long N = GetNumber();
        long historyLength = GetHistoryLength();
        /// Output container for particle population ancestor lines
        std::vector<population<Space> > populationAL(historyLength);
        /// Get final iteration indices
        arma::Col<unsigned int> tmpAI = History[historyLength-1].GetAIndices();
        /// Get final iteration particle values
        populationAL[historyLength-1] = GetHistoryPopulation(historyLength-1);
        /// Iterating from last iteration to initiliaization period whereby:
        /// 1. Re-ordering particle population elements for each t=0,...,T-1 accoridng to ancestral lines
        /// 2. Storing these particle population classes in output: populationAL
        if (historyLength > 1){
            for (int t = historyLength - 2; t >= 0; --t){
                for(int i = 0; i < N; ++i){
                    /// Setting each particle population element at indices i=0,...,N-1 to its ancesor index in the output vector populationAL
                    populationAL[t].SetValueN(GetHistoryPopulation(t).GetValueN(static_cast<int>(tmpAI(i))), i);
                }
                /// Re-mapping ancestor indices to retrieve ancestor indices a period before
                tmpAI = (History[t].GetAIndices()).elem(tmpAI);
            }
        }
        /// The output container stores population classes re-ordered according to their ancestor lines for t=0,...,T-1. So, for example, the first each element of populationAL[t] is the direct ancestor population class element of populationAL[t+1].
        return populationAL;
    }
    template <class Space, class Params>
    arma::Col<unsigned int> sampler<Space, Params>::GetALineInd(long n) const
    {
        long historyLength = GetHistoryLength();

        /// Get final iteration indices and final iteration particle values
        unsigned int tmpAI = (History[historyLength-1].GetAIndices()).at(static_cast<int>(n));
        arma::Col<unsigned int> ALineInd(historyLength);
        ALineInd[historyLength - 1] = static_cast<int>(n);

        for(int t = historyLength - 2; t >= 0; --t){
            tmpAI = (History[t].GetAIndices()).at(tmpAI);
            ALineInd[t] = tmpAI;
        }
        return ALineInd;
    }
    template <class Space, class Params>
    std::vector<Space> sampler<Space,Params>::GetALineSpace(long n) const
    {
        long historyLength = GetHistoryLength();

        std::vector<Space> ALineSpace(historyLength);

        /// Get final iteration indices and final iteration particle values
        arma::Col<unsigned int> ALineInd(historyLength);
        ALineInd = this->GetALineInd(static_cast<int>(n));

        for(int t = 0; t < historyLength; ++t){
            ALineSpace[t] = GetHistoryPopulation(t).GetValueN(static_cast<int>(ALineInd(t)));
        }
        return ALineSpace;
    }
    /// At present this function resets the system evolution time to 0 and calls the moveset initialisor to assign each
    /// particle in the ensemble.
    ///
    /// Note that the initialisation function must be specified before calling this function.
    template <class Space, class Params>
    void sampler<Space,Params>::Initialise(void)
    {
        T = 0;
        dlogNCPath = 0.0;
        acceptProb = -1;

        //Set the initial values and log weights of the particles
        std::vector<Space> InitVal(N);
        arma::vec InitWeights(N);
        pPopulation = population<Space>(InitVal,InitWeights);
        pMoves->DoInit(pPopulation,N,algParams);

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
            Resample(rtResampleMode);
        }
        else {
            nResampled = 0;
            if(htHistoryMode == HistoryType::AL) {
                uRSIndices = arma::linspace<arma::Col<unsigned int>>(0, N - 1, N);
            }
            pAdapt->updateForMCMC(algParams,pPopulation,acceptProb,nResampled,nRepeats);
        }

        //A possible MCMC step should be included here.
        bool didMCMC =  pMoves->DoMCMC(0,pPopulation, N, nRepeats, nAccepted,algParams);
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

    /// This function returns the result of integrating the supplied function under the empirical measure associated with the
    /// particle set at the present time. The final argument of the integrand function is a pointer which will be supplied
    /// with pAuxiliary to allow for arbitrary additional information to be passed to the function being integrated.
    ///
    /// \param pIntegrand The function to integrate with respect to the particle set
    /// \param pAuxiliary A pointer to any auxiliary data which should be passed to the function

    template <class Space, class Params>
    double sampler<Space,Params>::Integrate(double(*pIntegrand)(const Space&,void*), void * pAuxiliary)
    {
        long double rValue = 0;
        long double wSum = expl(stableLogSumWeights(pPopulation.GetLogWeight()));
        for(int i =0; i < N; i++)
        {
            rValue += expl(pPopulation.GetLogWeightN(i)) * pIntegrand(pPopulation.GetValueN(i), pAuxiliary);
        }

        rValue /= wSum;
        return static_cast<double>(rValue);
    }


    /// This function is intended to be used to estimate integrals of the sort which must be evaluated to determine the
    /// normalising constant of a distribution obtained using a sequence of potential functions proportional to densities with respect
    /// to the initial distribution to define a sequence of distributions leading up to the terminal, interesting distribution.
    ///
    /// In this context, the particle set at each time is used to make an estimate of the path sampling integrand, and
    /// numerical integration is then performed to obtain an estimate of the path sampling integral which is the natural logarithm
    /// of the ratio of normalising densities.
    ///
    /// The integrand is integrated at every time point in the population history. The results of this integration are
    /// taken to be point-evaluations of the path sampling integrand which are spaced on a grid of intervals given by the
    /// width function. The path sampling integral is then calculated by performing a suitable numerical integration and
    /// the results of this integration is returned.
    ///
    /// pAuxiliary is passed to both of the user specified functions to allow the user to pass additional data to either or
    /// both of these functions in a convenient manner. It is safe to use NULL if no such data is used by either function.
    ///
    /// \param PStype  The numerical integration method to use
    /// \param pIntegrand  The function to integrated. The first argument is evolution time, the second the particle value at which the function is to be evaluated and the final argument is always pAuxiliary.
    /// \param pWidth      The function which returns the width of the path sampling grid at the specified evolution time. The final argument is always pAuxiliary
    /// \param pAuxiliary  A pointer to auxiliary data to pass to both of the above functions
    /// \tparam Space The class used to represent a point in the sample space.
    /// \tparam Params (optional) The class used for any additional parameters.
    ///
    /// The PStype parameter should be set to one of the following:
    /// -# PathSamplingType::RECTANGLE to use the rectangle rule for integration.
    /// -# PathSamplingType::TRAPEZOID1 to use the trapezoidal rule for integration.
    /// -# PathSamplingType::TRAPEZOID2 to use the trapezoidal rule for integration with a second order correction.

    template <class Space, class Params>
    double sampler<Space,Params>::IntegratePathSampling(PathSamplingType::Enum PStype, double (*pIntegrand)(long,const Space &, void*), double (*pWidth)(long,void*), void* pAuxiliary)
    {
        if(htHistoryMode == HistoryType::NONE)
        throw SMC_EXCEPTION(SMCX_MISSING_HISTORY, "The path sampling integral cannot be computed as the history of the system was not stored.");

        // historyelement<Space> histel;
        // histel.Set(N, pPopulation, nAccepted, nRepeats, historyflags(nResampled));
        // History.push_back(histel);


        long lTime = 1;
        long double rValue = 0.0;
        typename std::vector<historyelement<Space> >::const_iterator it;

        switch(PStype) {
        case PathSamplingType::RECTANGLE:
            {
                for(it = ++History.begin(); it!=History.end(); it++){
                    rValue += it->Integrate(lTime, pIntegrand, pAuxiliary) * static_cast<long double>(pWidth(lTime,pAuxiliary));
                    lTime++;
                }
                break;
            }


        case PathSamplingType::TRAPEZOID1:
            {

                long double previous_expt = History.begin()->Integrate(0,pIntegrand,pAuxiliary);
                long double current_expt;
                for(it = ++History.begin(); it!=History.end(); it++){
                    current_expt = it->Integrate(lTime, pIntegrand, pAuxiliary);
                    rValue += static_cast<long double>(pWidth(lTime,pAuxiliary))/2.0 * (previous_expt + current_expt) ;
                    lTime++;
                    previous_expt = current_expt;
                }

                break;
            }

        case PathSamplingType::TRAPEZOID2:
        default:
            {
                long double previous_expt = History.begin()->Integrate(0,pIntegrand,pAuxiliary);
                long double previous_var = History.begin()->Integrate_Var(0,pIntegrand,previous_expt,pAuxiliary);
                long double current_expt;
                long double current_var;
                long double width = 0.0;
                for(it = ++History.begin(); it!=History.end(); it++){
                    current_expt = it->Integrate(lTime, pIntegrand, pAuxiliary);
                    current_var = it->Integrate_Var(lTime, pIntegrand, current_expt, pAuxiliary);
                    width = static_cast<long double>(pWidth(lTime,pAuxiliary));
                    rValue += width/2.0 * (previous_expt + current_expt) - std::pow(width,2.0)/12.0*(current_var - previous_var);
                    lTime++;
                    previous_expt = current_expt;
                    previous_var = current_var;
                }

                break;
            }

        }

        // History.pop_back();

        return static_cast<double>(rValue);
    }

    /// The iterate function:
    ///         -# moves the current particle set
    ///         -# checks the effective sample size and resamples if necessary
    ///         -# performs a mcmc step if required
    ///         -# appends the current particle set to the history if desired
    ///         -# increments the current evolution time
    template <class Space, class Params>
    void sampler<Space,Params>::Iterate(void)
    {
        IterateEss();
        return;
    }

    template <class Space, class Params>
    void sampler<Space,Params>::IterateBack(void)
    {
        if(htHistoryMode == HistoryType::NONE)
        throw SMC_EXCEPTION(SMCX_MISSING_HISTORY, "An attempt to undo an iteration was made; unfortunately, the system history has not been stored.");

        History.pop_back();
        historyelement<Space> recent = History.back();
        pPopulation = recent.GetRefs();
        N =recent.GetNumber();
        nAccepted = recent.AcceptCount();
        nResampled = recent.WasResampled();
        nRepeats = recent.mcmcRepeats();
        T--;
        return;
    }

    template <class Space, class Params>
    double sampler<Space,Params>::IterateEss(void)
    {
        pAdapt->updateForMove(algParams,pPopulation);

        //Move the particle set.
        MoveParticles();

        //Estimate the normalising constant
        dlogNCIt = CalcLogNC();
        dlogNCPath += dlogNCIt;

        //Normalise the weights
        pPopulation.SetLogWeight(pPopulation.GetLogWeight()  - dlogNCIt);

        //Check if the ESS is below some reasonable threshold and resample if necessary.
        //A mechanism for setting this threshold is required.
        double ESS = GetESS();
        if(ESS < dResampleThreshold) {
            nResampled = 1;
            pAdapt->updateForMCMC(algParams,pPopulation,acceptProb,nResampled,nRepeats);
            Resample(rtResampleMode);
        }
        else {
            nResampled = 0;
            if(htHistoryMode == HistoryType::AL) {
                uRSIndices = arma::linspace<arma::Col<unsigned int>>(0, N - 1, N);
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
    void sampler<Space,Params>::IterateUntil(long lTerminate)
    {
        while(T < lTerminate)
        Iterate();
    }

    template <class Space, class Params>
    void sampler<Space,Params>::MoveParticles(void)
    {
        pMoves->DoMove(T+1,pPopulation, N,algParams);
    }

    template <class Space, class Params>
    void sampler<Space,Params>::Resample(ResampleType::Enum lMode)
    {
        //Resampling is done in place.
        int uMultinomialCount;

        //First obtain a count of the number of children each particle has.
        switch(lMode) {
        case ResampleType::MULTINOMIAL:
            //Sample from a suitable multinomial vector
            dRSWeights = exp(pPopulation.GetLogWeight() - stableLogSumWeights(pPopulation.GetLogWeight()));
            rmultinom(static_cast<int>(N), dRSWeights.memptr(), static_cast<int>(N), uRSCount.memptr());
            break;

        case ResampleType::RESIDUAL:
            dRSWeights = exp(log(static_cast<double>(N)) + pPopulation.GetLogWeight() - stableLogSumWeights(pPopulation.GetLogWeight()));
            uRSIndices = arma::zeros<arma::Col<unsigned int> >(static_cast<int>(N));
            for(int i = 0; i < N; ++i)
            uRSIndices(i) = static_cast<unsigned int>(floor(dRSWeights(i)));
            dRSWeights = dRSWeights - uRSIndices;
            dRSWeights = dRSWeights/sum(dRSWeights);
            uMultinomialCount = N - arma::sum(uRSIndices);
            rmultinom(uMultinomialCount, dRSWeights.memptr(), static_cast<int>(N), uRSCount.memptr());
            uRSCount += arma::conv_to<arma::Col<int> >::from(uRSIndices);
            break;

        case ResampleType::STRATIFIED:
        default:
            {
                // Procedure for stratified sampling
                //Generate a random number between 0 and 1/N
                double dRand = R::runif(0,1.0 / static_cast<double>(N));
                arma::vec dWeightCumulative = arma::cumsum(exp(pPopulation.GetLogWeight() - stableLogSumWeights(pPopulation.GetLogWeight())));
                int j = 0, k = 0;
                uRSCount = arma::zeros<arma::Col<int> >(static_cast<int>(N));
                //while(j < N) {
                while(k < N) {
                    while((dWeightCumulative(k) - dRand) > static_cast<double>(j)/static_cast<double>(N) && j < N) {
                        uRSCount(k)++;
                        j++;
                        dRand = R::runif(0,1.0 / static_cast<double>(N));
                    }
                    k++;
                }
                break;
            }

        case ResampleType::SYSTEMATIC:
            {
                // Procedure for stratified sampling but with a common RV for each stratum
                //Generate a random number between 0 and 1/N
                double dRand = R::runif(0,1.0 / static_cast<double>(N));
                int j = 0, k = 0;
                uRSCount = arma::zeros<arma::Col<int> >(static_cast<int>(N));
                arma::vec dWeightCumulative = arma::cumsum(exp(pPopulation.GetLogWeight() - stableLogSumWeights(pPopulation.GetLogWeight())));
                //while(j < N) {
                while(k < N) {
                    while((dWeightCumulative(k) - dRand) > static_cast<double>(j)/static_cast<double>(N) && j < N) {
                        uRSCount(k)++;
                        j++;
                    }
                    k++;
                }
                break;
            }
        }

        uRSIndices = arma::zeros<arma::Col<unsigned int> >(static_cast<int>(N));
        //Map count to indices to allow in-place resampling
        for (int i=0, j=0; i<N; ++i) {
            if (uRSCount(i)>0) {
                uRSIndices(i) = i;
                while (uRSCount(i)>1) {
                    while (uRSCount(j)>0) ++j; // find next free spot
                    uRSIndices(j++) = i; // assign index
                    --uRSCount(i); // decrement number of remaining offsprings
                }
            }
        }

        //Perform the replication of the chosen.
        for(int i = 0; i < N ; ++i) {
            if(uRSIndices(i) != static_cast<unsigned int>(i)){
                pPopulation.SetValueN( pPopulation.GetValueN(static_cast<int>(uRSIndices(i))) ,i);
            }
        }

        //Set equal normalised weights
        pPopulation.SetLogWeight(- log(static_cast<double>(N))*arma::ones(N));
    }

    /// This function configures the resampling parameters, allowing the specification of both the resampling
    /// mode and the threshold at which resampling is used.
    ///
    /// \param rtMode The resampling mode to be used.
    /// \param dThreshold The threshold at which resampling is deemed necesary.
    ///
    /// The rtMode parameter should be set to one of the following:
    /// -# ResampleType::MULTINOMIAL to use multinomial resampling
    /// -# ResampleType::RESIDUAL to use residual resampling
    /// -# ResampleType::STRATIFIED to use stratified resampling
    /// -# ResampleType::SYSTEMATIC to use systematic resampling
    ///
    /// The dThreshold parameter can be set to a value in the range [0,1) corresponding to a fraction of the size of
    /// the particle set or it may be set to an integer corresponding to an actual effective sample size.

    template <class Space, class Params>
    void sampler<Space,Params>::SetResampleParams(ResampleType::Enum rtMode, double dThreshold)
    {
        rtResampleMode = rtMode;
        if(dThreshold < 1)
        dResampleThreshold = dThreshold * N;
        else
        dResampleThreshold = dThreshold;
    }

    /// Produce a human-readable display of the current nth particle value and log weight.
    ///
    /// \param os The output stream to which the display should be made.
    /// \param n The index of the particle of interest
    template <class Space, class Params>
    std::ostream & sampler<Space,Params>::StreamParticle(std::ostream & os, long n) const
    {
        os << pPopulation.GetValueN(n) << "," << pPopulation.GetWeightN(n) << std::endl;
        return os;
    }

    /// Produce a human-readable display of the current particle values and log weights.
    ///
    /// \param os The output stream to which the display should be made.
    template <class Space, class Params>
    std::ostream & sampler<Space,Params>::StreamParticles(std::ostream & os) const
    {
        os << pPopulation << std::endl;
        return os;
    }

    /// This function records the MCMC acceptance history to the specified output stream as a list of
    /// the number of moves accepted at each time instant.
    ///
    /// \param os The output stream to send the data to.
    template <class Space, class Params>
    void sampler<Space,Params>:: OstreamMCMCRecordToStream(std::ostream &os) const
    {
        os << "Accepted MCMC proposals history:" << std::endl;
        os << "======================" << std::endl;
        for(typename std::vector<historyelement<Space> >::const_iterator it = History.begin(); it!=History.end(); it++){
            os << it->AcceptCount() << std::endl;
        }
    }
    /// This function records the resampling history to the specified output stream as a 0-1 valued list which takes
    /// the value 1 for those time instances when resampling occured and 0 otherwise.
    ///
    /// \param os The output stream to send the data to.
    template <class Space, class Params>
    void sampler<Space,Params>:: OstreamResamplingRecordToStream(std::ostream &os) const
    {
        os << "Resampling history:" << std::endl;
        os << "======================" << std::endl;
        os << "Flag\t" << "ESS\t" << std::endl;
        for(typename std::vector<historyelement<Space> >::const_iterator it = History.begin(); it!=History.end(); it++){
            if(it->WasResampled())
            os << "1\t";
            else
            os << "0\t";

            os << it->GetESS() << std::endl;
        }
    }

}

namespace std {
    /// Produce a human-readable display of the state of an smc::sampler class using the stream operator.

    /// \param os The output stream to which the display should be made.
    /// \param s  The sampler which is to be displayed.
    template <class Space, class Params>
    std::ostream & operator<< (std::ostream & os, smc::sampler<Space,Params> & s)
    {
        os << "Sampler Configuration:" << std::endl;
        os << "======================" << std::endl;
        os << "Evolution Time:   " << s.GetTime() << std::endl;
        os << "Particle Set Size:" << s.GetNumber() << std::endl;
        os << "Effective Sample Size:" << s.GetESS() << std::endl;
        os << std::endl;
        os << "Particle Set:" << std::endl;
        s.StreamParticles(os);
        os << std::endl;
        return os;
    }
}
#endif
