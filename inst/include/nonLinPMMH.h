// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
//
// nonLinPMMH.h: Example 3.1 of Andrieu et al. (2010). Implementing particle marginal
// Metropolis-Hastings for a toy non-linear state space model previously described in
// Gordon et al. (1993) and Kitagawa (1996).
//
// Copyright (C) 2017         Dirk Eddelbuettel, Adam Johansen and Leah South
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

#include "smctc.h"

namespace nonLinPMMH {

    class parameters 
    {
    public:
        double sigv, sigw;
    };
    arma::vec y; //data
    
    double logPrior(const parameters & proposal);
	
    //A derived class for the moves
    class nonLinPMMH_move:
    public smc::moveset<double,smc::nullParams>
    {
    public:
	    
        void pfInitialise(double & value, double & logweight, smc::nullParams & param);
        void pfMove(long lTime, double & value, double & logweight, smc::nullParams & param);

        ~nonLinPMMH_move() {};

    };
    
    parameters theta_prop;
	smc::moveset<double,smc::nullParams>* myMove;
}

