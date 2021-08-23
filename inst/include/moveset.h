// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
//
// moveset.h: Rcpp integration of SMC library -- sampler proposal moves
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
//! \brief Classes and functions which deal with collections of sampler proposal "moves".
//!
//! This file contains definitions of smc::moveset.
//! It deals with the collections of proposal moves (including initialisation and MCMC moves) which must be dealt with by the sampler.

#ifndef __SMC_MOVESET_HH
#define __SMC_MOVESET_HH 1.0

#include "population.h"

namespace smc {

	/// A template class for a set of moves for use in an SMC samplers framework.
    template <class Space, class Params> class moveset {

    private:

        ///Default functions (only needed so that they can be overriden for backwards compatibility)
        ///The function which initialises a single particle.
        void (*defaultInitialise)(Space &, double &, Params &);
        ///The functions which perform actual moves on a single particle.
        void (*defaultMove)(long, Space &, double &, Params &);
        ///One iteration of a Markov Chain Monte Carlo move for a single particle.
        bool (*defaultMCMC)(long, Space &, double &, Params &);
        ///The function which weights a single (reference) particle coordinate.
        void (*defaultWeight)(long, Space &, double &, Params &);

    public:

	    ///Create a completely unspecified moveset
        moveset();

        ///An alternative constructor for backwards compatibility
        moveset(void (*pfInit)(Space &, double &, Params &),
        void (*pfNewMove)(long, Space &, double &, Params &),
        bool (*pfNewMCMC)(long, Space &, double &, Params &));

        ///An alternative constructor used for conditional SMC specifically.
        moveset(void (*pfInit)(Space &, double &, Params &),
        void (*pfNewMove)(long, Space &,double &, Params &),
        bool (*pfNewMCMC)(long, Space &,double &, Params &),
        void (*pfNewWeight)(long, Space &, double &, Params &));

        /// Free the workspace allocated for the algorithm parameters.
        virtual ~moveset() {
        }

        /// Holder function for updates to be done before the move step.
        virtual void pfInitialise(Space & value, double & weight, Params & myParams) {(*defaultInitialise)(value, weight, myParams);}

        /// Holder function for updates to be done before the MCMC step.
        virtual void pfMove(long time, Space & value, double & weight, Params & myParams) {(*defaultMove)(time, value, weight, myParams);}

        /// Holder function for updates to be done at the end of each iteration.
        virtual bool pfMCMC(long time, Space & value,double & weight, Params & myParams) {
            if(defaultMCMC){
                return (*defaultMCMC)(time, value, weight, myParams);
            } else{
                return 0;
            }
        }
        /// Holder function for weighting of conditioning reference trajectory.
        virtual void pfWeight(long time, Space & referenceValue, double & referenceWeight, Params & myParams) {(*defaultWeight)(time, referenceValue, referenceWeight, myParams);}

        ///Initialise the population of particles
        virtual void DoInit(population<Space> & pFrom, long N, Params &);
        ///Perform an MCMC move on the particles
        virtual bool DoMCMC(long lTime, population<Space> & pFrom, long N, int nRepeats, int & nAccepted, Params &);
        ///Select an appropriate move at time lTime and apply it to pFrom
        virtual void DoMove(long lTime, population<Space> & pFrom,long N, Params &);
        ///Performs a conditional move: sets conditional reference value and re-weights corresponding particle coordinate
        virtual void DoConditionalMove(long lTime, population<Space> & pFrom, const Space & referenceValue, long lReferenceIndex, Params & params);
    };


    /// The argument free smc::moveset constructor simply sets the number of available moves to zero and sets
    /// all of the associated function pointers to NULL.
    template <class Space, class Params>
    moveset<Space,Params>::moveset()
    {
        defaultInitialise = NULL;
        defaultMove = NULL;
        defaultMCMC = NULL;
        defaultWeight = NULL;
    }

    template <class Space, class Params>
    moveset<Space,Params>::moveset(void (*pfInit)(Space &, double &, Params &),
    void (*pfNewMove)(long, Space &,double &, Params &),
    bool (*pfNewMCMC)(long,Space &,double &, Params &))
    {
        defaultInitialise = pfInit;
        defaultMove = pfNewMove;
        defaultMCMC = pfNewMCMC;
    }

    template <class Space, class Params>
    moveset<Space,Params>::moveset(void (*pfInit)(Space &, double &, Params &),
    void (*pfNewMove)(long, Space &, double &, Params &),
    bool (*pfNewMCMC)(long, Space &, double &, Params &),
    void (*pfNewWeight)(long, Space &, double &, Params &))
    {
        defaultInitialise = pfInit;
        defaultMove = pfNewMove;
        defaultMCMC = pfNewMCMC;
        defaultWeight = pfNewWeight;
    }

    template <class Space, class Params>
    void moveset<Space,Params>::DoInit(population<Space> & pFrom, long N, Params & params) {
        for (long i=0; i<N; i++){
            pfInitialise(pFrom.GetValueRefN(i),pFrom.GetLogWeightRefN(i),params);
        }
    }

    template <class Space, class Params>
    bool moveset<Space,Params>::DoMCMC(long lTime, population<Space> & pFrom, long N, int nRepeats, int & nAccepted, Params & params)
    {
		// NOTE: previously this checked for the existence of pfMCMC but now it will always exist.
		// Need to check behaviour of this and add a warning that the interpretation
		// won't make sense if there is no MCMC function in the derived moveset class.
		//if(pfMCMC && nRepeats>0) {
		if(nRepeats>0) {
            nAccepted = 0;
            for (int j=0; j<nRepeats; j++){
                for (long i=0; i<N; i++){
                    nAccepted += pfMCMC(lTime,pFrom.GetValueRefN(i),pFrom.GetLogWeightRefN(i),params);
                }
            }
            return TRUE;
        }
        else {
            nAccepted = 0;
            return FALSE; // an MCMC step was not performed because there was no user defined MCMC step or because the number of MCMC repeats was zero.
        }
    }

    template <class Space, class Params>
    void moveset<Space,Params>::DoMove(long lTime, population<Space> & pFrom, long N, Params & params)
    {
        for (long i=0; i<N; i++){
            pfMove(lTime,pFrom.GetValueRefN(i),pFrom.GetLogWeightRefN(i),params);
        }
    }

    template <class Space, class Params>
    void moveset<Space,Params>::DoConditionalMove(long lTime, population<Space> & pFrom, const Space & referenceValue, long lReferenceIndex, Params & params) {
        /// Sets conditional reference value of lTime-coordinate.
        pFrom.SetValueN(referenceValue,
                        lReferenceIndex);
        /// Re-weights this coordinate.
        pfWeight(lTime,
                 pFrom.GetValueRefN(lReferenceIndex),
                 pFrom.GetLogWeightRefN(lReferenceIndex),
                 params);
    }
}
#endif
