// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
//
// pf.cpp: Rcpp wrapper for SMC library -- first example
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

// RcppSMC builds on and wrap SMCTC which is
//
//   Copyright Adam Johansen, 2008.
//
// and released under GPL-3, see the copyright headers in inst/include/ 

#include <Rcpp.h>

#include "smctc.hh"
#include "pffuncs.hh"
#include <cstdio> 
#include <cstdlib>
#include <cstring>

using namespace std;

///The observations
cv_obs * y;
long load_data(char const * szName, cv_obs** y);

double integrand_mean_x(const cv_state&, void*);
double integrand_mean_y(const cv_state&, void*);
double integrand_var_x(const cv_state&, void*);
double integrand_var_y(const cv_state&, void*);

// pf() function callable from R via Rcpp:: essentially the same as main() from pf.cc 
extern "C" SEXP pf(SEXP fileS, SEXP partS) { 	

    long lIterates;

    try {

        std::string filename = Rcpp::as<std::string>(fileS);
        unsigned long lNumber = Rcpp::as<unsigned long>(partS);

        //Load observations
        lIterates = load_data(filename.c_str(), &y);

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

long load_data(char const * szName, cv_obs** yp)
{
  FILE * fObs = fopen(szName,"rt");
  if (!fObs)
    throw SMC_EXCEPTION(SMCX_FILE_NOT_FOUND, "Error: pf assumes that the current directory contains an appropriate data file called data.csv\nThe first line should contain a constant indicating the number of data lines it contains.\nThe remaining lines should contain comma-separated pairs of x,y observations.");
  char* szBuffer = new char[1024];
  char* rc = fgets(szBuffer, 1024, fObs);
  if (rc==NULL)
    throw SMC_EXCEPTION(SMCX_FILE_NOT_FOUND, "Error: no data found.");
  long lIterates = strtol(szBuffer, NULL, 10);

  *yp = new cv_obs[lIterates];
  
  for(long i = 0; i < lIterates; ++i)
    {
      rc = fgets(szBuffer, 1024, fObs);
      (*yp)[i].x_pos = strtod(strtok(szBuffer, ",\r\n "), NULL);
      (*yp)[i].y_pos = strtod(strtok(NULL, ",\r\n "), NULL);
    }
  fclose(fObs);

  delete [] szBuffer;

  return lIterates;
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
