// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
//
// rngR.h: Rcpp wrapper for SMC library -- R-only RNG class
//
// Copyright (C) 2008 - 2009  Adam Johansen
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
// along with RInside.  If not, see <http://www.gnu.org/licenses/>.


//! \file 
//! \brief Random number generation.
//!
//! This file contains the definitions for the smc::rng and smc::rnginfo class.
//! It wraps the random number generation facilities provided by the GSL and provides a convenient interfaces to access several of its more commonly-used features.

#ifndef __SMC_RNGR_H
#define __SMC_RNGR_H 1.0

#include <Rcpp.h>

namespace smc {
    /// A GSL-RNG information handling class (not templated)
    /// Now essentially empty as using R RNGs requires less state
    class gslrnginfo {
    private:
        ///This is a null terminated array of the available random number generators.
        //const gsl_rng_type** typePtArray;
        ///The number of available random number generators.
        //int nNumber;

    protected:
        gslrnginfo() { 
            // Nothing to do for us
        }

    public:
        ///Returns a reference to the sole static instance of this class.
        static gslrnginfo* GetInstance();

        ///Returns the number of available random number generators.
        //int GetNumber();
        ///Returns the name of random number generator number nIndex.
        //const char* GetNameByIndex(int nIndex);
        ///Returns a pointer to random number generator nIndex.
        //const gsl_rng_type* GetPointerByIndex(int nIndex);
        ///Returns a pointer to the random number generator with name szName (or null if it doesn't exist).
        //const gsl_rng_type* GetPointerByName(const char* szName);
    };

    ///The global application instance of the gslrnginfo class:
    extern gslrnginfo rngset;


    ///A random number generator class.

    ///    At present this serves as a wrapper for the gsl random number generation code.
    ///    Rewritten for R package using R's RNGs
    class rng {
    private:
        ///This is the type of random number generator underlying the class.
        //const gsl_rng_type* type;
        ///This is a pointer to the internal workspace of the rng including its current state.
        //gsl_rng* pWorkspace;

        Rcpp::RNGScope *scopeptr; 	// RNGScope saves and later restores R's RNG state 

    public:
        ///Initialise the random number generator using default settings
        rng(void) {
            scopeptr = new Rcpp::RNGScope; 	// needed for state saving/restoring of R RNGs
        }

        ///Initialise the random number generator using the default seed for the type
        //rng(const gsl_rng_type* Type);
        ///Initialise the random number generator using specified type and seed
        //rng(const gsl_rng_type* Type,unsigned long int lSeed);
        
        ///Free the workspace allocated for random number generation
        ~rng() {
            delete scopeptr;
        }

        ///Provide access to the raw random number generator
        //gsl_rng* GetRaw(void);

        ///Generate a multinomial random vector with parameters (n,w[1:k]) and store it in X
        void Multinomial(unsigned n, unsigned k, double* w, unsigned* X) {
            Rcpp::IntegerVector v(k);
            double sum = 0.0;
            unsigned int i;
            for (i=0; i<k; i++) sum += w[i]; 
            for (i=0; i<k; i++) w[i] = w[i] / sum;

            // R sources:  rmultinom(int n, double* prob, int K, int* rN);
            rmultinom(static_cast<int>(n), const_cast<double*>(w), static_cast<int>(k), &(v[0]));

            for (i=0; i<k; i++) {
                X[i] = static_cast<unsigned int>(v[i]);
            }
        }

        ///Returns a random integer generated uniformly between the minimum and maximum values specified
        long UniformDiscrete(long lMin, long lMax) {
            return ::Rf_runif(static_cast<double>(lMin), static_cast<double>(lMax));
        }

        ///Returns a random number generated from a Beta distribution with the specified parameters.
        double Beta(double da, double db) {
            return ::Rf_rbeta(da, db);
        }

        ///Returns a random number generated from a Cauchy distribution with the specified scale parameter.
        double Cauchy(double dScale) {
            return ::Rf_rcauchy(0.0, dScale);
        }

        ///Returns a random number generated from an exponential distribution with the specified mean.
        double Exponential(double dMean) {
            return ::Rf_rexp(dMean);
        }

        ///Return a random number generated from a gamma distribution with shape alpha and scale beta.
        double Gamma(double dAlpha, double dBeta) {
            return ::Rf_rgamma(dAlpha, dBeta);
        }

        ///Returns a random number generated from a Laplacian distribution with the specified scale.
        // "ported" from the gsl_ran_laplace functon
        double Laplacian(double dScale) {
            double u;
            do {
                u = 2 * ::Rf_runif(0.0,1.0) - 1.0;
            } while (u == 0.0);

            if (u < 0) {
                return dScale * log (-u);
            } else {
                return -dScale * log (u);
            }
        }

        ///Returns a random number generated from a Lognormal distribution of location mu and scale sigma
        double Lognormal(double dMu, double dSigma) {
            return ::Rf_rlnorm(dMu, dSigma);
        }

        ///Return a random number generated from a normal distribution with a specified mean and standard deviation
        double Normal(double dMean, double dStd) {
            return ::Rf_rnorm(dMean, dStd);
        }

        ///Return a random number generated from a standard normal distribution
        double NormalS(void) {
            return ::Rf_rnorm(0.0, 1.0);
        }

        ///Returns a random number from a normal distribution, conditional upon it exceeding the specified threshold.
        // no Rcpp yet -- double NormalTruncated(double dMean, double dStd, double dThreshold);
#if 0
        ///This function simply calls gsl_ran_gaussian_tail with the specified parameters and performs appropriate shifting.
        ///     \param dMean The mean of the distribution.
        ///     \param dStd  The standard deviation of the distribution
        ///     \param dThreshold The lower truncation threshold.
        double rng::NormalTruncated(double dMean, double dStd, double dThreshold)
        {
            return dMean + gsl_ran_gaussian_tail(pWorkspace, dThreshold - dMean, dStd);
        }
#endif

        ///Return a student-t random number generated with a specified number of degrees of freedom
        double StudentT(double dDF) {
            return ::Rf_rt(dDF);
        }

        ///Return a random number generated uniformly between dMin and dMax
        double Uniform(double dMin, double dMax) {
            return ::Rf_runif(dMin, dMax);
        }

        ///Returns a random number generated from the standard uniform[0,1) distribution
        double UniformS(void) {
            return ::Rf_runif(0.0, 1.0);
        }

    };
}
 
 
#endif
