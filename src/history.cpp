// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
//
// history.cpp: Rcpp integration of SMC library -- sampler history
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
// along with RcppSMC.  If not, see <http://www.gnu.org/licenses/>.

#include "smctc.h"

//! \file
//! \brief This file contains the untemplated functions used for storing the history of the system.

namespace smc {
  /// This constructor produces an initialised historyflags instance.
  ///
  /// \param wasResampled An indicator which should be nonzero if the particle 
  /// system was resampled during the iteration being described

  historyflags::historyflags(int wasResampled)
  {
    if(wasResampled)
      Resampled = 1;
    else
      Resampled = 0;
  }
}
