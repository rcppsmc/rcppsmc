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

class cv_state 
{
public:
    double x_pos, y_pos;
    double x_vel, y_vel;
};

class cv_obs
{
public:
    double x_pos, y_pos;
};

double logLikelihood(long lTime, const cv_state & X);

void fInitialise(cv_state & value, double & logweight, smc::rng *pRng);
void fMove(long lTime, cv_state & value, double & logweight, smc::rng *pRng);

extern double nu_x;
extern double nu_y;
extern double Delta;

extern std::vector<cv_obs> y;
