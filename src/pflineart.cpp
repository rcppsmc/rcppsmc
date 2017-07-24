// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
//
// pflineart.cpp: Rcpp wrapper for SMC library -- first example of Johansen (2009)
//
// Copyright (C) 2008 - 2009  Adam Johansen
// Copyright (C) 2012 - 2017  Dirk Eddelbuettel and Adam Johansen
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

#include <RcppArmadillo.h>

#include "smctc.h"
#include "pflineart.h"
#include "rngR.h"

#include <cstdio> 
#include <cstdlib>
#include <cstring>

using namespace std;

///The observations
std::vector<cv_obs> y;

double integrand_mean_x(const cv_state&, void*);
double integrand_mean_y(const cv_state&, void*);
double integrand_var_x(const cv_state&, void*);
double integrand_var_y(const cv_state&, void*);

// pf() function callable from R via Rcpp:: essentially the same as main() from pf.cc 
// minor interface change to pass data down as matrix, rather than a filename
extern "C" SEXP pfLineartBS_impl(SEXP dataS, SEXP partS, SEXP usefS, SEXP funS) { 	

    long lIterates;

    try {
        unsigned long lNumber = Rcpp::as<unsigned long>(partS);
        bool useF = Rcpp::as<bool>(usefS);
        Rcpp::Function f(funS);

        // Load observations -- or rather copy them in from R
        Rcpp::NumericMatrix dat = Rcpp::NumericMatrix(dataS); // so we expect a matrix
        lIterates = dat.nrow();
        y.reserve(lIterates);
        for (long i = 0; i < lIterates; ++i) {
            y[i].x_pos = dat(i,0);
            y[i].y_pos = dat(i,1);
        }

        //Initialise and run the sampler
        smc::sampler<cv_state> Sampler(lNumber, SMC_HISTORY_NONE);  
        smc::moveset<cv_state> Moveset(fInitialise, fMove, NULL);

        Sampler.SetResampleParams(SMC_RESAMPLE_RESIDUAL, 0.5);
        Sampler.SetMoveSet(Moveset);
        Sampler.Initialise();

        Rcpp::NumericVector Xm(lIterates), Xv(lIterates), Ym(lIterates), Yv(lIterates);

        Xm(0) = Sampler.Integrate(integrand_mean_x, NULL);
        Xv(0) = Sampler.Integrate(integrand_var_x, (void*)&Xm(0));
        Ym(0) = Sampler.Integrate(integrand_mean_y, NULL);
        Yv(0) = Sampler.Integrate(integrand_var_y, (void*)&Ym(0));

        for(int n=1; n < lIterates; ++n) {
            Sampler.Iterate();
            
            Xm(n) = Sampler.Integrate(integrand_mean_x, NULL);
            Xv(n) = Sampler.Integrate(integrand_var_x, (void*)&Xm(n));
            Ym(n) = Sampler.Integrate(integrand_mean_y, NULL);
            Yv(n) = Sampler.Integrate(integrand_var_y, (void*)&Ym(n));

            if (useF) f(Xm, Ym);
        }

        return Rcpp::DataFrame::create(Rcpp::Named("Xm") = Xm,
        Rcpp::Named("Xv") = Xv,
        Rcpp::Named("Ym") = Ym,
        Rcpp::Named("Yv") = Yv);
    }
    catch(smc::exception  e) {
        Rcpp::Rcout << e;
    }
    return R_NilValue;          	// to provide a return 
}


double integrand_mean_x(const cv_state& s, void *)
{
    return s.x_pos;
}

double integrand_var_x(const cv_state& s, void* vmx)
{
    double* dmx = (double*)vmx;
    double d = (s.x_pos - (*dmx));
    return d*d;
}

double integrand_mean_y(const cv_state& s, void *)
{
    return s.y_pos;
}

double integrand_var_y(const cv_state& s, void* vmy)
{
    double* dmy = (double*)vmy;
    double d = (s.y_pos - (*dmy));
    return d*d;
}
#include <iostream>
#include <cmath>


using namespace std;

double var_s0 = 4;
double var_u0 = 1;
double var_s  = 0.02;
double var_u  = 0.001;

double scale_y = 0.1;
double nu_y = 10.0;
double Delta = 0.1;

///The function corresponding to the log likelihood at specified time and position (up to normalisation)

///  \param lTime The current time (i.e. the index of the current distribution)
///  \param X     The state to consider 
double logLikelihood(long lTime, const cv_state & X)
{
    return - 0.5 * (nu_y + 1.0) * (log(1 + pow((X.x_pos - y[lTime].x_pos)/scale_y,2) / nu_y) + log(1 + pow((X.y_pos - y[lTime].y_pos)/scale_y,2) / nu_y));
}

///A function to initialise particles

/// \param value The value of the particle being moved
/// \param logweight The log weight of the particle being moved
/// \param pRng A pointer to the random number generator which is to be used
void fInitialise(cv_state & value, double & logweight, smc::rng *pRng)
{
    value.x_pos = pRng->Normal(0,sqrt(var_s0));
    value.y_pos = pRng->Normal(0,sqrt(var_s0));
    value.x_vel = pRng->Normal(0,sqrt(var_u0));
    value.y_vel = pRng->Normal(0,sqrt(var_u0));

    logweight = logLikelihood(0,value);
}

///The proposal function.

///\param lTime The sampler iteration.
///\param value The value of the particle being moved
///\param logweight The log weight of the particle being moved
///\param pRng  A random number generator.
void fMove(long lTime, cv_state & value, double & logweight, smc::rng *pRng)
{
    value.x_pos += value.x_vel * Delta + pRng->Normal(0,sqrt(var_s));
    value.x_vel += pRng->Normal(0,sqrt(var_u));
    value.y_pos += value.y_vel * Delta + pRng->Normal(0,sqrt(var_s));
    value.y_vel += pRng->Normal(0,sqrt(var_u));

    logweight += logLikelihood(lTime, value);
}
