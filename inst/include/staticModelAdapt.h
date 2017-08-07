// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
//
// staticModelAdapt.h: A class containing parameters and functions to update
// these parameters in the context of static Bayesian models. The methods to
// estimate the empirical covariance and its Cholesky decomposition are applicable
// for applications where the particle values are a vector of doubles. The methods
// to compute the next 'temperature' are relevant to likelihood annealing SMC where the
// power posteriors are defined by P_t(theta|y) \propto P(y|theta)^\gamma_t P(theta).
// Here y is observed data, theta denotes the parameters of the model and gamma
// denotes the 'temperatures' which start at 0 (the prior) and finish at 1
// (the posterior).
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

#ifndef __SMC_STATICMODELADAPT_HH
#define __SMC_STATICMODELADAPT_HH 1.0

#include <population.h>

namespace smc {

    /// A class containing parameters and functions to update these parameters in the context of static Bayesian models.
    class staticModelAdapt{
    private:
        /// The temperature schedule;
        std::vector<double> temp;
        /// The empirical covariance matrix estimated from the population of particles.
        arma::mat empCov;
        /// The Cholesky decomposition of the empirical covariance matrix.
        arma::mat cholCov;

        /// Computes the difference between the conditional ESS given the specified temperature difference and the desired conditional ESS.
        double CESSdiff(const arma::vec & logweight, const arma::vec & loglike, double tempDiff, double desiredCESS);
        ///Performs the bisection method to find the temperature within (temp_curr,1) which gives the desired conditional ESS.
        double bisection(double curr, const arma::vec & logweight, const arma::vec & loglike, double desiredCESS, double epsilon);

    public:
        ~staticModelAdapt() {}
        ///The class constructor which sets the current and previous temperatures to zero.
        staticModelAdapt() {temp.push_back(0.0);}

        /// Chooses the next temperature such that a desired conditional ESS is maintained.
        void ChooseTemp(const arma::vec & logweight, const arma::vec & loglike, double desiredCESS, double epsilon = 0.01);
        /// Calculates the empirical covariance matrix based on the current weighted particle set.
        void calcEmpCov(const arma::mat & theta, const arma::vec & logweight);
        /// Calculates the Cholesky decomposition of the empirical covariance matrix based on the current weighted particle set.
        void calcCholCov(const arma::mat & theta, const arma::vec logweight);
        /// Calculates the number of MCMC repeats based on the results of the most recent set of MCMC moves.
        int calcMcmcRepeats(double acceptProb, double desiredAcceptProb, int initialN, int maxReps);

        /// Returns the t-th temperature.
        double GetTemp(int t) const {return temp[t];}
        /// Returns the current temperature.
        double GetTempCurr(void) const {return temp.back();}
        /// Returns the previous temperature.
        double GetTempPrevious(void) const {return temp.rbegin()[1];}
        /// Returns the vector of temperatures
        std::vector<double> GetTemps(void) const {return temp;}
        /// Returns the Cholesky decomposition of the empirical covariance matrix based on the current weighted particle set.
        arma::mat GetCholCov(void) const {return cholCov;}
        /// Returns the empirical covariance matrix based on the current weighted particle set.
        arma::mat GetEmpCov(void) const {return empCov;}
    };
}
#endif
