// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
//
// helper.cpp: Rcpp integration of SMC library -- helper functions of general interest
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

//! \file
//! \brief This file contains the untemplated functions of general interest.

namespace smc {
    /// This function performs a stable calculation of the log sum of the weights, which is useful for
    /// normalising weights, calculating the effective sample size and estimating the normalising constant.
    ///
    /// \param logw The log weights of interest.
    double stableLogSumWeights(const arma::vec & logw){
        double dMaxWeight = arma::max(logw);
        double sum = arma::sum(exp(logw - dMaxWeight));
        return (dMaxWeight + log(sum));
    }
}
