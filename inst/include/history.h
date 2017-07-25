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
//! It defines smc::history, smc::historyelement and smc::history.

#ifndef __SMC_HISTORY_HH
#define __SMC_HISTORY_HH 1.0

#include <RcppArmadillo.h>

namespace smc {
    /// The historyflags class holds a set of flags which describe various properties of the particle system at a given time.
    class historyflags
    {
    private:
        /// true if the particle system was resampled during the described iteration.
        unsigned int Resampled : 1;
    public:
        ///Create a new set of history flags corresponding to the specified properties
        historyflags(int wasResampled);

        ///This function returns true if the flag set indicates that the ensemble was resampled during the described iteration.
        int WasResampled(void) {return Resampled;}
    };

    /// A template class for the elements of a linked list to be used for the history of the sampler.
    template <class Population>class historyelement
    {
    private:
        long number; //!< The number of particles (presently redundant as this is not a function of iteration)
        int nAccepted; //!< Number of MCMC moves accepted during this iteration.
        Population value; //!< The particles themselves (values and weights)
        historyflags flags; //!< Flags associated with this iteration.
        
    public:
        /// The null constructor creates an empty history element.
        historyelement();
        /// A constructor with four arguments initialises the particle set.
        historyelement(long lNumber, Population pNew, int nAccepts, historyflags hf);

        /// The destructor tidies up.
        ~historyelement();

        /// Returns the effective sample size of this particle generation.
        double GetESS(void) const;
        /// Returns the flags
        historyflags GetFlags(void) const{return flags;}
        /// Returns the number of particles present.
        long GetNumber(void) const {return number;} 
        /// Returns the current particle set.
        Population GetValues(void) const { return value; }
		/// Returns a reference to the current particle set.
		Population & GetRefs(void) { return value; }
        /// Integrate the supplied function according to the empirical measure of the particle ensemble.
        long double Integrate(long lTime, double (*pIntegrand)(long,const Population&,long,void*), void* pAuxiliary) const;
		/// Sets the particle set to the specified values.  
		void Set(long lNumber, const Population &New, int inAccepted, const historyflags &histflags){number = lNumber; value = New; nAccepted = inAccepted; flags = histflags;};

        /// Returns the number of MCMC moves accepted during this iteration.
        int AcceptCount(void) {return nAccepted; }
        /// Returns true if the particle set 
        int WasResampled(void) {return flags.WasResampled(); }
        
        

    };

	template <class Population>
	historyelement<Population>::historyelement(): flags(0)
	{
		number = 0;
		nAccepted = 0;
	}


    /// \param lNumber The number of particles present in the particle generation
    /// \param New    The array of particles which are present in the particle generation
    /// \param nAccepts The number of MCMC moves that were accepted during this particle generation
    /// \param hf      The historyflags associated with the particle generation

    template <class Population>
	historyelement<Population>::historyelement(long lNumber, Population New, int nAccepts, historyflags hf) :
	number(lNumber), nAccepted(nAccepts), value(New), flags(hf)
	{
	}

    template <class Population>
    historyelement<Population>::~historyelement(void)
    {
    }

    template <class Population>
    double historyelement<Population>::GetESS(void) const
    {
        double sum = arma::sum(exp(value.GetLogWeight()));
		double sumsq = arma::sum(exp(2.0*value.GetLogWeight()));
		return expl(-log(sumsq) + 2*log(sum));
    }

    /// \param lTime The timestep at which the integration is to be carried out
    /// \param pIntegrand The function which is to be integrated
    /// \param pAuxiliary A pointer to additional information which is passed to the integrand function

    template <class Population>
	long double historyelement<Population>::Integrate(long lTime, double (*pIntegrand)(long,const Population&,long,void*), void* pAuxiliary) const
	{
		long double rValue = 0;
		long double wSum = 0;
		for(long i =0; i < number; i++)
		{
			rValue += expl(value.GetLogWeightN(i)) * (long double)pIntegrand(lTime, value,i, pAuxiliary); //may want to change input type for this pIntegrand
			wSum  += expl(value.GetLogWeightN(i));
		}

		rValue /= wSum;
		return rValue;
	}

}

#endif
