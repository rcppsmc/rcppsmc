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
    template <class Space, class Params> class moveset
    {
    private:
        ///The number of moves which are present in the set
        long number;
        ///The function which initialises a single particle.
        void (*pfInitialise)(Space &, double &, Params &);
        ///The function which selects a move for a given particle at a given time.
        long (*pfMoveSelect)(long , const Space &, const double &, Params &);
        ///The functions which perform actual moves on a single particle.
        void (**pfMoves)(long, Space &, double &, Params &);
        ///One iteration of a Markov Chain Monte Carlo move for a single particle.
        bool (*pfMCMC)(long, Space &,double &, Params &);

    public:
        ///Create a completely unspecified moveset
        moveset();
        ///Create a reduced moveset with a single move
        moveset(void (*pfInit)(Space &, double &, Params &),
        void (*pfNewMoves)(long, Space &,double &, Params &),
        bool (*pfNewMCMC)(long,Space &,double &, Params &));
        ///Create a fully specified moveset
        moveset(void (*pfInit)(Space &, double &, Params &),long (*pfMoveSelector)(long , const Space &, const double &, Params &), 
        long nMoves, void (**pfNewMoves)(long, Space &,double &, Params &),
        bool (*pfNewMCMC)(long,Space &, double &, Params &));
        
        ///Initialise the population of particles
        void DoInit(population<Space> & pFrom, long N, Params &);
        ///Perform an MCMC move on the particles
        bool DoMCMC(long lTime, population<Space> & pFrom, long N, int nRepeats, int & nAccepted, Params &);
        ///Select an appropriate move at time lTime and apply it to pFrom
        void DoMove(long lTime, population<Space> & pFrom,long N, Params &);
        
        ///Free the memory used for the array of move pointers when deleting
        ~moveset();

        /// \brief Set the initialisation function.
        /// \param pfInit is a function which returns a particle generated according to the initial distribution 
        void SetInitialisor( void (*pfInit)(Space &, double &, Params &))
        {pfInitialise = pfInit;}

        /// \brief Set the MCMC function
        /// \param pfNewMCMC  The function which performs an MCMC move
        void SetMCMCFunction(bool (*pfNewMCMC)(long,Space &,double &, Params &)) {pfMCMC = pfNewMCMC;}

        /// \brief Set the move selection function
        /// \param pfMoveSelectNew returns the index of move to perform at the specified time given the specified particle
        void SetMoveSelectionFunction(long (*pfMoveSelectNew)(long , const Space &, const double &, Params &))
        {pfMoveSelect = pfMoveSelectNew;};

        ///Set the individual move functions to the supplied array of such functions
        void SetMoveFunctions(long nMoves, void (**pfNewMoves)(long, Space &, double &, Params &));
        
        ///Moveset assignment should allocate buffers and deep copy all members.
        moveset<Space,Params> & operator= (const moveset<Space,Params> & pFrom);
    };


    /// The argument free smc::moveset constructor simply sets the number of available moves to zero and sets
    /// all of the associated function pointers to NULL.
    template <class Space, class Params>
    moveset<Space,Params>::moveset()
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
    template <class Space, class Params>
    moveset<Space,Params>::moveset(void (*pfInit)(Space &, double &, Params &),
    void (*pfNewMoves)(long, Space &,double &, Params &),
    bool (*pfNewMCMC)(long,Space &,double &, Params &))
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
    /// \param pfMoveSelector The function which selects a move to apply, at a specified time, to a specified particle
    /// \param nMoves The number of moves which are defined in general
    /// \param pfNewMoves An array of functions which moves particles at a specified time to a new location
    /// \param pfNewMCMC The function which should be called to apply an MCMC move (if any)
    template <class Space, class Params>
    moveset<Space,Params>::moveset(void (*pfInit)(Space &, double &, Params &),long (*pfMoveSelector)(long ,const Space &, const double &, Params &), 
    long nMoves, void (**pfNewMoves)(long, Space &,double &, Params &),
    bool (*pfNewMCMC)(long,Space &,double &, Params &))
    {
        SetInitialisor(pfInit);
        SetMoveSelectionFunction(pfMoveSelector);
        pfMoves = NULL; //This presents an erroneous deletion by the following call
        SetMoveFunctions(nMoves, pfNewMoves);
        SetMCMCFunction(pfNewMCMC);
    }

    template <class Space, class Params>
    moveset<Space,Params>::~moveset()
    {    if(pfMoves)
        delete [] pfMoves;
    }


    template <class Space, class Params>
    void moveset<Space,Params>::DoInit(population<Space> & pFrom, long N, Params & params) {
        for (long i=0; i<N; i++){
            (*pfInitialise)(pFrom.GetValueRefN(i),pFrom.GetLogWeightRefN(i),params);
        }
    }

    template <class Space, class Params>
    bool moveset<Space,Params>::DoMCMC(long lTime, population<Space> & pFrom, long N, int nRepeats, int & nAccepted, Params & params)
    {
        if(pfMCMC && nRepeats>0) {
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
        if(number > 1)
        for (long i=0; i<N; i++){
            pfMoves[pfMoveSelect(lTime,pFrom.GetValueRefN(i),pFrom.GetLogWeightRefN(i),params)](lTime,pFrom.GetValueRefN(i),pFrom.GetLogWeightRefN(i),params);
        }
        else
        
        for (long i=0; i<N; i++){
            pfMoves[0](lTime,pFrom.GetValueRefN(i),pFrom.GetLogWeightRefN(i),params);
        }
    }

    /// \param nMoves The number of moves which are defined in general.
    /// \param pfNewMoves An array of functions which moves particles at a specified time to a new location
    ///
    /// The move functions accept two arguments, the first of which corresponds to the system evolution time and the
    /// second to an initial particle position and the second to a weighted starting position. It returns a new 
    /// weighted position corresponding to the moved particle.
    template <class Space, class Params>
    void moveset<Space,Params>::SetMoveFunctions(long nMoves,void (**pfNewMoves)(long,Space &,double &, Params &))
    {
        number = nMoves;
        if(pfMoves)
        delete [] pfMoves;
        pfMoves = (void (**)(long int, Space &, double &, Params &)) new void* [nMoves];
        for(int i = 0; i < nMoves; i++)
        pfMoves[i] = pfNewMoves[i];
        return;
    }

    template <class Space, class Params>
    moveset<Space,Params> & moveset<Space,Params>::operator= (const moveset<Space,Params> & pFrom)
    {
        SetInitialisor(pFrom.pfInitialise);
        SetMCMCFunction(pFrom.pfMCMC);
        SetMoveSelectionFunction(pFrom.pfMoveSelect);
        SetMoveFunctions(pFrom.number, pFrom.pfMoves);       
        return *this;
    }
}
#endif
