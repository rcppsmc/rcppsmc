// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
//
// rareEvents.cpp: Rcpp wrapper for SMC library -- second example of Johansen (2009)
//
// Copyright (C) 2012         Dirk Eddelbuettel
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
// along with RInside.  If not, see <http://www.gnu.org/licenses/>.
// from examples/pf/pfexample.cc; pffuncs.cc and pffuncs.hh also used

// RcppSMC builds on top of, and wrap, SMCTC which is
//
//   Copyright Adam Johansen, 2008.
//
// and released under GPL-3, see the copyright headers in inst/include/ 

#include <Rcpp.h>

#include <iostream>
#include <cmath>
#include "simfunctions.hh"

///Length of Markov Chain
long lChainLength = 15;
///Number of distributions used
long lIterates;
///Rare event threshold
double dThreshold = 5.0;
///Annealing schedule constant
double dSchedule = 30.0;

// rareEvents() function callable from R via Rcpp -- which is essentially 
// the same as main() from SMCTC's examples/rare-events/main.cc 
extern "C" SEXP rareEvents(SEXP numberS, SEXP iteratesS, SEXP thresholdS, SEXP scheduleS) { 	

    long lNumber = Rcpp::as<long>(numberS);			// Number of particles
    lIterates  = Rcpp::as<long>(iteratesS);			// Number of iterations
    dThreshold = Rcpp::as<double>(thresholdS);			// Rare event threshold
    dSchedule  = Rcpp::as<double>(scheduleS);			// Annealing schedule constant

    try{
        ///An array of move function pointers
        void (*pfMoves[])(long, smc::particle<mChain<double> > &,smc::rng*) = {fMove1, fMove2};
        smc::moveset<mChain<double> > Moveset(fInitialiseMC, fSelect, sizeof(pfMoves) / sizeof(pfMoves[0]), pfMoves, fMCMC);
        smc::sampler<mChain<double> > Sampler(lNumber, SMC_HISTORY_RAM);
        
        Sampler.SetResampleParams(SMC_RESAMPLE_STRATIFIED,0.5);
        Sampler.SetMoveSet(Moveset);

        Sampler.Initialise();
        Sampler.IterateUntil(lIterates);
      
        ///Estimate the normalising constant of the terminal distribution
        double zEstimate = Sampler.IntegratePathSampling(pIntegrandPS, pWidthPS, NULL) - log(2.0);
        ///Estimate the weighting factor for the terminal distribution
        double wEstimate = Sampler.Integrate(pIntegrandFS, NULL);
      
        // cout << zEstimate << " " << log(wEstimate) << " " << zEstimate + log(wEstimate) << endl;
        Rcpp::NumericVector res = Rcpp::NumericVector::create(Rcpp::Named("zEstimate") = zEstimate,
                                                              Rcpp::Named("log(wEstimate)") = log(wEstimate),
                                                              Rcpp::Named("sum") = zEstimate + log(wEstimate));
        return res;
    }
    catch(smc::exception  e) {
        Rcpp::Rcout << e;       	// not cerr, as R doesn't like to mix with i/o 
        //exit(e.lCode);		// we're just called from R so we should not exit
    }
    return R_NilValue;          	// to provide a return 

}

