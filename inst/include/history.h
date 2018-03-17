// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
//
// history.h: Rcpp integration of SMC library -- sampler history
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
//

//! \file
//! \brief Classes and function related to the history of the sampler.
//!
//! This file contains template definitions for the classes used to store the history of an SMCTC sampler.
//! It defines smc::historyflags and smc::historyelement.

#ifndef __SMC_HISTORY_HH
#define __SMC_HISTORY_HH 1.0

#include "population.h"

namespace smc {
    /// The historyflags class holds a set of flags which describe various properties of the particle system at a given time.
    class historyflags
    {
    private:
        /// true if the particle system was resampled during the described iteration.
        unsigned int Resampled : 1;
    public:
      // ///Create a new set of history flags corresponding to the specified properties
      // historyflags(int wasResampled);

      /// This constructor produces an initialised historyflags instance.
      ///
      /// \param wasResampled An indicator which should be nonzero if the particle
      /// system was resampled during the iteration being described
      historyflags(int wasResampled)
      {
	if(wasResampled)
	  Resampled = 1;
	else
	  Resampled = 0;
      }

        ///This function returns true if the flag set indicates that the ensemble was resampled during the described iteration.
        int WasResampled(void) {return Resampled;}
    };

    /// A template class for the elements of a linked list to be used for the history of the sampler.
    template <class Space>class historyelement
    {
    private:
        long number; //!< The number of particles (presently redundant as this is not a function of iteration)
        int nAccepted; //!< Number of MCMC moves accepted during this iteration.
        int nRepeat; //!< Number of MCMC iterations performed at this iteration (per particle)
        population<Space> pop; //!< The particles themselves (values and weights)
        historyflags flags; //!< Flags associated with this iteration.

    public:
        /// The null constructor creates an empty history element.
        historyelement();
        /// A constructor with four arguments initialises the particle set.
        historyelement(long lNumber, population<Space> pNew, int nAccepts, int nRepeats, historyflags hf);

        /// The destructor tidies up.
        ~historyelement();

        /// Returns the effective sample size of this particle generation.
        double GetESS(void) const;
        /// Returns the flags
        historyflags GetFlags(void) const{return flags;}
        /// Returns the number of particles present.
        long GetNumber(void) const {return number;}
        /// Returns the current particle set.
        population<Space> GetValues(void) const { return pop; }
        /// Returns a reference to the current particle set.
        population<Space> & GetRefs(void) { return pop; }
        /// Monte Carlo estimate of the expectation of the supplied function with respect to the empirical measure of the particle ensemble.
        long double Integrate(long lTime, double (*pIntegrand)(long,const Space&,void*), void* pAuxiliary) const;
        /// Monte Carlo estimate of the variance of the supplied function with respect to the empirical measure of the particle ensemble (to be used in second order trapezoidal correction).
        long double Integrate_Var(long lTime, double (*pIntegrand)(long,const Space&,void*), double Expectation, void* pAuxiliary) const;
        /// Sets the particle set to the specified values.
        void Set(long lNumber, const population<Space> &New, int inAccepted, int nRepeats, const historyflags &histflags){number = lNumber; pop = New; nAccepted = inAccepted; nRepeat = nRepeats; flags = histflags;};

        /// Returns the number of MCMC moves accepted during this iteration.
        int AcceptCount(void) {return nAccepted; }
        /// Returns the number of MCMC iterations performed during this iteration.
        int mcmcRepeats(void) {return nRepeat; }
        /// Returns true if the particle set
        int WasResampled(void) {return flags.WasResampled(); }




    };

    template <class Space>
    historyelement<Space>::historyelement(): flags(0)
    {
        number = 0;
        nAccepted = 0;
        nRepeat = 0;
    }


    /// \param lNumber The number of particles present in the particle generation
    /// \param New    The array of particles which are present in the particle generation
    /// \param nAccepts The number of MCMC moves that were accepted during this particle generation
    /// \param nRepeats The number of MCMC iterations that were performed during this particle generation
    /// \param hf      The historyflags associated with the particle generation

    template <class Space>
    historyelement<Space>::historyelement(long lNumber, population<Space> New, int nAccepts, int nRepeats, historyflags hf) :
    number(lNumber), nAccepted(nAccepts), nRepeat(nRepeats), pop(New), flags(hf)
    {
    }

    template <class Space>
    historyelement<Space>::~historyelement(void)
    {
    }

    template <class Space>
    double historyelement<Space>::GetESS(void) const
    {
        return expl(2*stableLogSumWeights(pop.GetLogWeight())-stableLogSumWeights(2.0*pop.GetLogWeight()));
    }

    /// \param lTime The timestep at which the integration is to be carried out
    /// \param pIntegrand The function which is to be integrated
    /// \param pAuxiliary A pointer to additional information which is passed to the integrand function

    template <class Space>
    long double historyelement<Space>::Integrate(long lTime, double (*pIntegrand)(long,const Space&,void*), void* pAuxiliary) const
    {
        long double rValue = 0;
        long double wSum = expl(stableLogSumWeights(pop.GetLogWeight()));
        for(long i =0; i < number; i++)
        {
            rValue += expl(pop.GetLogWeightN(i)) * static_cast<long double>(pIntegrand(lTime, pop.GetValueN(i), pAuxiliary));
        }

        rValue /= wSum;
        return rValue;
    }

    /// \param lTime The timestep at which the integration is to be carried out
    /// \param pIntegrand The function which is to be integrated
    /// \param pAuxiliary A pointer to additional information which is passed to the integrand function

    template <class Space>
    long double historyelement<Space>::Integrate_Var(long lTime, double (*pIntegrand)(long,const Space &,void*), double Expectation, void* pAuxiliary) const
    {
        long double rValue = 0;
        long double wSum = expl(stableLogSumWeights(pop.GetLogWeight()));
        for(long i =0; i < number; i++)
        {
            rValue += expl(pop.GetLogWeightN(i)) * std::pow(static_cast<long double>(pIntegrand(lTime, pop.GetValueN(i), pAuxiliary) - Expectation),2.0);
        }

        rValue /= wSum;
        return rValue;
    }

}

#endif
