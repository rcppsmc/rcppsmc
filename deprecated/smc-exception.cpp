// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
//
// smc-exceptions.cpp: Rcpp integration of SMC library -- exception handling
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

#include "smc-exception.h"

//! \file
//! \brief The untemplated smc::exception class is implemented here.

namespace smc
{
  //! Generate an SMCTC Exception class with the specified initialisation.

  //! This constructor fills the four elements of the class with their specified values.
  //! It is used to allow a single-line command to create and throw an exception.
  //!
  //! \param szN The name of the source file generating the exception.
  //! \param lL The line in that file responsible for the exception.
  //! \param lC The numerical code identifying the exception.
  //! \param szM An textual explanation of the problem.
  exception::exception(char const * szN, long lL, long lC, char const * szM)
  {
    szFile = szN; 
    lLine = lL; 
    lCode = lC; 
    szMessage = szM;
  }
}

namespace std {
  ///Display a human-readable version of an SMC exception.

  /// \param os The stream to write to.
  /// \param e The exception class to display.
  /// \return os
  std::ostream & operator<< (std::ostream & os, smc::exception & e)
  {
    os << "SMC Exception: #" << e.lCode << endl;
    os << "Error occured in file " << e.szFile << " at line " << e.lLine << "." << endl;
    os << "Details: " << endl << '\t' << e.szMessage << endl;
    return os;
  }
}

