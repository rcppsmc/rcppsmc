// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
//
// pflineart.h: Rcpp wrapper for SMC library -- first example of Johansen (2009)
//
// Copyright (C) 2008 - 2009  Adam Johansen
// Copyright (C) 2012         Dirk Eddelbuettel and Adam Johansen
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

#include "smctc.h"

namespace pflineart {
    class cv_state
    {
    public:
        double x_pos, y_pos;
        double x_vel, y_vel;
    };

    class cv_obs
    {
    public:
        arma::vec x_pos, y_pos;
    };
	
	double logLikelihood(long lTime, const cv_state & X);
	
	//A derived class for the moves
    class pflineart_move:
    public smc::moveset<cv_state,smc::nullParams>
    {
    public:
        
		void pfInitialise(cv_state & value, double & logweight, smc::nullParams & param);
        void pfMove(long lTime, cv_state & value, double & logweight, smc::nullParams & param);

        ~pflineart_move() {};

    };

    double integrand_mean_x(const cv_state&, void*);
    double integrand_mean_y(const cv_state&, void*);
    double integrand_var_x(const cv_state&, void*);
    double integrand_var_y(const cv_state&, void*);
	
	smc::moveset<cv_state,smc::nullParams>* myMove;

}
