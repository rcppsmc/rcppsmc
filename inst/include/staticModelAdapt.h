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

      /// Computes the difference between the conditional ESS given the specified temperature difference and the desired conditional ESS.
    inline double staticModelAdapt::CESSdiff(const arma::vec & logweight, const arma::vec & loglike, double tempDiff, double desiredCESS){
        double logsum1 = stableLogSumWeights(logweight + tempDiff*loglike);
        double logsum2 = stableLogSumWeights(logweight + 2*tempDiff*loglike);

        return exp(log(static_cast<double>(logweight.n_rows)) + 2*logsum1 - logsum2) - desiredCESS;
    }
    
    ///Performs the bisection method to find the temperature within (temp_curr,1) which gives the desired conditional ESS.
    inline double staticModelAdapt::bisection(double curr, const arma::vec & logweight, const arma::vec & loglike, double desiredCESS, double epsilon){
        double a = curr;
        double b = 1.0;
        double f_a = CESSdiff(logweight,loglike,a-curr,desiredCESS);
        double f_b = CESSdiff(logweight,loglike,b-curr,desiredCESS);
        if (f_a*f_b>0){
            Rcpp::stop("Bisection method to choose the next temperature failed");
        } else{
            double m, f_m, err;
            m = (a+b)/2.0;
            f_m = CESSdiff(logweight,loglike,m-curr,desiredCESS);
            err = 10.0;
            while (err > epsilon){
                if (f_m<0.0){
                    b = m;
                    f_b = f_m;
                } else{
                    a = m;
                    f_a = f_m;
                }
                m = (a+b)/2.0;
                f_m = CESSdiff(logweight,loglike,m-curr,desiredCESS);
                err = std::abs(f_m);
            }
            return m;
        }
    }

    /// Chooses the next temperature such that a desired conditional ESS is maintained.
    ///
    /// \param logweight An armadillo vector containing the logarithm of the current particle weights.
    /// \param loglike An armadillo vector containing the log likelihood of the current particle values.
    /// \param desiredCESS The target conditional ESS for the next temperature (generally fixed).
    /// \param epsilon The tolerance for the bisection method (maximum difference between desired and actual conditional ESS).
    inline void staticModelAdapt::ChooseTemp(const arma::vec & logweight, const arma::vec & loglike, double desiredCESS, double epsilon) {
        double temp_curr = temp.back();
        if (CESSdiff(logweight,loglike,1.0-temp_curr,desiredCESS)>=-epsilon){
            temp.push_back(1.0);
        } else {
            temp.push_back(bisection(temp_curr, logweight, loglike, desiredCESS, epsilon));
        }
    }

    /// Calculates the empirical covariance matrix based on the current weighted particle set.
    ///
    /// \param theta An [Nxd] armadillo matrix of doubles for the current particle values, where N is
    /// the number of particles and d is the dimension of the parameter.
    /// \param logweight An armadillo vector containing the logarithm of the current particle weights.
    inline void staticModelAdapt::calcEmpCov(const arma::mat & theta, const arma::vec & logweight){
        int N = logweight.n_rows;
        arma::vec normWeights = exp(logweight - stableLogSumWeights(logweight));

        arma::mat diff = theta - arma::ones(N,1)*arma::mean(theta,0);
        empCov = diff.t()*diagmat(normWeights)*diff;
    }

    /// Calculates the Cholesky decomposition of the empirical covariance matrix based on the current weighted particle set.
    ///
    /// \param theta An [Nxd] armadillo matrix of doubles for the current particle values, where N is
    /// the number of particles and d is the dimension of the parameter
    /// \param logweight An armadillo vector of the logarithm of the current particle weights
    inline void staticModelAdapt::calcCholCov(const arma::mat & theta, const arma::vec logweight){
        calcEmpCov(theta,logweight);
        cholCov = arma::chol(empCov);
    }

    /// Calculates the number of MCMC repeats based on the results of the most recent set of MCMC moves.
    ///
    /// \param acceptProb The proportion of accepted MCMC steps in the most recent iteration.
    /// \param desiredAcceptProb The desired probability of a successful move for each particle.
    /// \param initialN The initial number of MCMC repeats.
    /// \param maxReps The maximum number of MCMC repeats.
    inline int staticModelAdapt::calcMcmcRepeats(double acceptProb, double desiredAcceptProb, int initialN, int maxReps){
        if (acceptProb + 1.0 <= 1e-9){
            return initialN;
        } else if (acceptProb - 1.0 >= -1e-9){
            return 1;
        } else if (acceptProb <= 1e-9){
            return maxReps;
        } else {
            return std::min(maxReps,static_cast<int>(std::ceil(log(1.0-desiredAcceptProb)/log(1.0-acceptProb))));
        }
    }
}
#endif
