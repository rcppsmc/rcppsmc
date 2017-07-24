// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
//
// moveset.h: Rcpp integration of SMC library -- sampler proposal moves 
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
//! \brief Classes and functions which deal with collections of sampler proposal "moves".
//!
//! This file contains definitions of smc::moveset.
//! It deals with the collections of proposal moves (including initialisation and MCMC moves) which must be dealt with by the sampler.

#ifndef __SMC_MOVESET_HH
#define __SMC_MOVESET_HH 1.0

#include "population.h"

namespace smc {

    /// A template class for a set of moves for use in an SMC samplers framework.
    template <class Space> class moveset
    {
    private:
        ///The number of moves which are present in the set
        long number;
        ///The function which initialises a single particle.
        void (*pfInitialise)(Space &, double &, rng*);
        ///The function which selects a move for the population of particles at a given time.
        long (*pfMoveSelect)(long , const Space &, const double &, rng*);
        ///The functions which perform actual moves on a single particle.
        void (**pfMoves)(long, Space &, double &, rng*);
        ///A Markov Chain Monte Carlo move for a single particle.
        int (*pfMCMC)(long, Space &,double &, rng*);

    public:
        ///Create a completely unspecified moveset
        moveset();
        ///Create a reduced moveset with a single move
        moveset(void (*pfInit)(Space &, double &, rng*),
        void (*pfNewMoves)(long, Space &,double &, rng*),
        int (*pfNewMCMC)(long,Space &,double &, rng*));
        ///Create a fully specified moveset
        moveset(void (*pfInit)(Space &, double &, rng*),long (*pfMoveSelector)(long , const population<Space> &, rng*), 
        long nMoves, void (**pfNewMoves)(long, Space &,double &, rng*),
        int (*pfNewMCMC)(long,Space &, double &, rng*));
        
        ///Initialise the population of particles
        void DoInit(rng * pRng, population<Space> & pFrom, long N);
        ///Perform an MCMC move on the particles
        int DoMCMC(long lTime, population<Space> & pFrom, rng* pRng, long N);
        ///Select an appropriate move at time lTime and apply it to pFrom
        void DoMove(long lTime, population<Space> & pFrom, rng* pRng,long N);
        
        ///Free the memory used for the array of move pointers when deleting
        ~moveset();

        /// \brief Set the initialisation function.
        /// \param pfInit is a function which returns a particle generated according to the initial distribution 
        void SetInitialisor( void (*pfInit)(Space &, double &, rng*))
        {pfInitialise = pfInit;}

        /// \brief Set the MCMC function
        /// \param pfNewMCMC  The function which performs an MCMC move
        void SetMCMCFunction(int (*pfNewMCMC)(long,Space &,double &,rng*)) {pfMCMC = pfNewMCMC;}

        /// \brief Set the move selection function
        /// \param pfMoveSelectNew returns the index of move to perform at the specified time given specified particles
        void SetMoveSelectionFunction(long (*pfMoveSelectNew)(long , const Space &, const double &, rng*))
        {pfMoveSelect = pfMoveSelectNew;};

        ///Set the individual move functions to the supplied array of such functions
        void SetMoveFunctions(long nMoves, void (**pfNewMoves)(long, Space &, double &, rng*));
        
        ///Moveset assignment should allocate buffers and deep copy all members.
        moveset<Space> & operator= (moveset<Space> & pFrom);
    };


    /// The argument free smc::moveset constructor simply sets the number of available moves to zero and sets
    /// all of the associated function pointers to NULL.
    template <class Space>
    moveset<Space>::moveset()
    {
        number = 0;

        pfInitialise = NULL;
        pfMoveSelect = NULL;
        pfMoves = NULL;
        pfMCMC = NULL;
    }

    /// The three argument moveset constructor creates a new set of moves and sets all of the relevant function
    /// pointers to the supplied values. Only a single move should exist if this constructor is used.
    /// \param pfInit The function which should be used to initialise a particle when the system is initialised
    /// \param pfNewMoves An functions which moves particles at a specified time to a new location
    /// \param pfNewMCMC The function which should be called to apply an MCMC move (if any)
    template <class Space>
    moveset<Space>::moveset(void (*pfInit)(Space &, double &, rng*),
    void (*pfNewMoves)(long, Space &,double &, rng*),
    int (*pfNewMCMC)(long,Space &,double &, rng*))
    {
        SetInitialisor(pfInit);
        SetMoveSelectionFunction(NULL);
        pfMoves = NULL; //This presents an erroneous deletion by the following call
        SetMoveFunctions(1, &pfNewMoves);
        SetMCMCFunction(pfNewMCMC);
    }

    /// The five argument moveset constructor creates a new set of moves and sets all of the relevant function
    /// pointers to the supplied values.
    /// \param pfInit The function which should be used to initialise particles when the system is initialised
    /// \param pfMoveSelector The function which selects a move to apply, at a specified time, to specified particles
    /// \param nMoves The number of moves which are defined in general
    /// \param pfNewMoves An array of functions which moves particles at a specified time to a new location
    /// \param pfNewMCMC The function which should be called to apply an MCMC move (if any)
    template <class Space>
    moveset<Space>::moveset(void (*pfInit)(Space &, double &, rng*),long (*pfMoveSelector)(long ,const population<Space> &, rng*), 
    long nMoves, void (**pfNewMoves)(long, Space &,double &, rng*),
    int (*pfNewMCMC)(long,Space &,double &, rng*))
    {
        SetInitialisor(pfInit);
        SetMoveSelectionFunction(pfMoveSelector);
        pfMoves = NULL; //This presents an erroneous deletion by the following call
        SetMoveFunctions(nMoves, pfNewMoves);
        SetMCMCFunction(pfNewMCMC);
    }

    template <class Space>
    moveset<Space>::~moveset()
    {    if(pfMoves)
        delete [] pfMoves;
    }


    template <class Space>
    void moveset<Space>::DoInit(rng * pRng, population<Space> & pFrom, long N) {
        for (long i=0; i<N; i++){
            (*pfInitialise)(pFrom.GetValueRefN(i),pFrom.GetLogWeightRefN(i),pRng);
        }
    }

    template <class Space>
    int moveset<Space>::DoMCMC(long lTime, population<Space> & pFrom, rng *pRng, long N)
    {    if(pfMCMC) {
            int count = 0;
            for (long i=0; i<N; i++){
                count += pfMCMC(lTime,pFrom.GetValueRefN(i),pFrom.GetLogWeightRefN(i),pRng);
            }
            return count;
        }
        else {
            return 0;
        }
    }

    template <class Space>
    void moveset<Space>::DoMove(long lTime, population<Space> & pFrom, rng *pRng, long N)
    {
        if(number > 1)
        for (long i=0; i<N; i++){
            pfMoves[pfMoveSelect(lTime,pFrom.GetValueRefN(i),pFrom.GetLogWeightRefN(i),pRng)](lTime,pFrom.GetValueRefN(i),pFrom.GetLogWeightRefN(i),pRng);
        }
        else
        
        for (long i=0; i<N; i++){
            pfMoves[0](lTime,pFrom.GetValueRefN(i),pFrom.GetLogWeightRefN(i),pRng);
        }
    }

    /// \param nMoves The number of moves which are defined in general.
    /// \param pfNewMoves An array of functions which moves particles at a specified time to a new location
    ///
    /// The move functions accept two arguments, the first of which corresponds to the system evolution time and the
    /// second to an initial particle position and the second to a weighted starting position. It returns a new 
    /// weighted position corresponding to the moved particle.
    template <class Space>
    void moveset<Space>::SetMoveFunctions(long nMoves,void (**pfNewMoves)(long,Space &,double &, rng*))
    {
        number = nMoves;
        if(pfMoves)
        delete [] pfMoves;
        pfMoves = (void (**)(long int, Space &, double &, rng*)) new void* [nMoves];
        for(int i = 0; i < nMoves; i++)
        pfMoves[i] = pfNewMoves[i];
        return;
    }

    template <class Space>
    moveset<Space> & moveset<Space>::operator= (moveset<Space> & pFrom)
    {
        SetInitialisor(pFrom.pfInitialise);
        SetMCMCFunction(pFrom.pfMCMC);
        SetMoveSelectionFunction(pFrom.pfMoveSelect);
        SetMoveFunctions(pFrom.number, pFrom.pfMoves);       
        return *this;
    }
}
#endif
