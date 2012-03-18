// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
//
// blockpfgaussianopt.h: Rcpp integration of SMC library -- Block PF Gaussian
//
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

#include "smctc.h"
#include <vector>

using namespace std;

smc::particle<vector<double> > fInitialiseBSPFG(smc::rng *pRng);
//long fSelectBSPFG(long lTime, const smc::particle<vector<double> > & p, smc::rng *pRng);
void fMoveBSPFG(long lTime, smc::particle<vector<double> > & pFrom, smc::rng *pRng);

namespace BSPFG {
    extern Rcpp::NumericVector y; 
}
using BSPFG::y;

extern long lLag;
