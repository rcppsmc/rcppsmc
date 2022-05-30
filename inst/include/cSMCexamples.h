// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
//
// Copyright (C) 2021         Adam Johansen, Dirk Eddelbuettel, Leah South and Ilya Zarubin
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

#include "RcppSMC.h"

namespace cSMCexamples {

    // I. Class declarations

    // Creating a class for a single double may be unneccessary but this example
    // could be extended later on to a D-dimensional state space similar to the
    // velocity example. Idea: to test printing facilities I'd like to
    // implement and could not get to during GSoC.
    class States
    {
    public:
        double xState;
    };
    // Same as above: to be extened to a multivariate LGSSM, hence the class.
    class Measurements
    {
    public:
        arma::vec yObs;
    };
    // Same as above: to be extened for parameters in a multivariate LGSS.
    class Parameters
    {
    public:
        double statePhi;
        double stateVarEvol;
    };

    class cSMCexamples_move:
    public smc::moveset<States,smc::nullParams>
    {
    public:
        void pfInitialise(States& stateValue,
                          double& logweight,
                          smc::nullParams& param);
        void pfMove(long lTime,
                    States& stateValue,
                    double& logweight,
                    smc::nullParams& param);
        void pfWeight(long lTime,
                      States& condStateValue,
                      double& logweight,
                      smc::nullParams& param);
        ~cSMCexamples_move() {};
    };

    // II. Helper function declaration:
    // 1. log-likelihood computation and copying
    // of reference trajectory from arma::vec to std::vect<States> container.
    double computeLogLikelihood(long lTime, const States& stateValue);
    // 2. Transform reference Trajectory of doubles to std::vector<States> type
    void copyReferenceTrajectory(const arma::vec& refArma,
                                 std::vector<States>& refStd);

    // III. Pointer declaration for moveset-class
    smc::moveset<States, smc::nullParams>* MyLGSSmove;
}
