// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
//
// pf.cpp: Rcpp wrapper for SMC library -- first example of Johansen (2009)
//
// Copyright (C) 2008         Adam Johansen
// Copyright (C) 2012         Dirk Eddelbuettel and Adam Johansen
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
// from examples/pf/pfexample.cc; pffuncs.cc and pffuncs.hh also used

#include <Rcpp.h>

#include "smctc.h"
#include "pffuncs.h"
#include "rngR.h"

#include <cstdio> 
#include <cstdlib>
#include <cstring>

using namespace std;

///The observations
//cv_obs * y;
//Rcpp::NumericMatrix y;
std::vector<cv_obs> y;
//long load_data(char const * szName, cv_obs** y);

double integrand_mean_x(const cv_state&, void*);
double integrand_mean_y(const cv_state&, void*);
double integrand_var_x(const cv_state&, void*);
double integrand_var_y(const cv_state&, void*);

// pf() function callable from R via Rcpp:: essentially the same as main() from pf.cc 
// minor interface change to pass data down as matrix, rather than a filename
extern "C" SEXP pf(SEXP dataS, SEXP partS, SEXP usefS, SEXP funS) { 	

    long lIterates;

    try {

        //std::string filename = Rcpp::as<std::string>(fileS);
        unsigned long lNumber = Rcpp::as<unsigned long>(partS);
        bool useF = Rcpp::as<bool>(usefS);
        Rcpp::Function f(funS);

        // Load observations -- or rather copy them in from R
        //lIterates = load_data(filename.c_str(), &y);
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

        for(int n=0; n < lIterates; ++n) {
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
        Rcpp::Rcout << e;       	// not cerr, as R doesn't like to mix with i/o 
        //exit(e.lCode);		// we're just called from R so we should not exit
    }
    return R_NilValue;          	// to provide a return 
}

#if 0
long load_data(char const * szName, cv_obs** yp)
{
  FILE * fObs = fopen(szName,"rt");
  if (!fObs)
    throw SMC_EXCEPTION(SMCX_FILE_NOT_FOUND, "Error: pf assumes that the current directory contains an appropriate data file called data.csv\nThe first line should contain a constant indicating the number of data lines it contains.\nThe remaining lines should contain comma-separated pairs of x,y observations.");
  char szBuffer[1024];
  char* rc = fgets(szBuffer, 1024, fObs);
  if (rc==NULL)
    throw SMC_EXCEPTION(SMCX_FILE_NOT_FOUND, "Error: no data found.");
  long lIterates = strtol(szBuffer, NULL, 10);

  *yp = new cv_obs[lIterates];
  
  for(long i = 0; i < lIterates; ++i)
    {
      rc = fgets(szBuffer, 1024, fObs);
      (*yp)[i].x_pos = strtod(strtok(szBuffer, ",\r\n "), NULL);
      (*yp)[i].y_pos = strtod(strtok(szBuffer, ",\r\n "), NULL);
    }
  fclose(fObs);

  return lIterates;
}
#endif

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
//#include <gsl/gsl_randist.h>


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
    //return - 0.5 * (nu_y + 1.0) * (log(1 + pow((X.x_pos - y[lTime].x_pos)/scale_y,2) / nu_y) + log(1 + pow((X.y_pos - y[lTime].y_pos)/scale_y,2) / nu_y));

    //double y_xpos = y(lTime,0); 
    //double y_ypos = y(lTime,1); 
    //return - 0.5 * (nu_y + 1.0) * (log(1 + pow((X.x_pos - y_xpos)/scale_y,2) / nu_y) + log(1 + pow((X.y_pos - y_ypos)/scale_y,2) / nu_y));
    
    return - 0.5 * (nu_y + 1.0) * (log(1 + pow((X.x_pos - y[lTime].x_pos)/scale_y,2) / nu_y) + log(1 + pow((X.y_pos - y[lTime].y_pos)/scale_y,2) / nu_y));
}

///A function to initialise particles

/// \param pRng A pointer to the random number generator which is to be used
smc::particle<cv_state> fInitialise(smc::rng *pRng)
{
  cv_state value;
  
  value.x_pos = pRng->Normal(0,sqrt(var_s0));
  value.y_pos = pRng->Normal(0,sqrt(var_s0));
  value.x_vel = pRng->Normal(0,sqrt(var_u0));
  value.y_vel = pRng->Normal(0,sqrt(var_u0));

  return smc::particle<cv_state>(value,logLikelihood(0,value));
}

///The proposal function.

///\param lTime The sampler iteration.
///\param pFrom The particle to move.
///\param pRng  A random number generator.
void fMove(long lTime, smc::particle<cv_state > & pFrom, smc::rng *pRng)
{
  cv_state * cv_to = pFrom.GetValuePointer();

  cv_to->x_pos += cv_to->x_vel * Delta + pRng->Normal(0,sqrt(var_s));
  cv_to->x_vel += pRng->Normal(0,sqrt(var_u));
  cv_to->y_pos += cv_to->y_vel * Delta + pRng->Normal(0,sqrt(var_s));
  cv_to->y_vel += pRng->Normal(0,sqrt(var_u));
  pFrom.AddToLogWeight(logLikelihood(lTime, *cv_to));
}
