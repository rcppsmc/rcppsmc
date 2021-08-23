// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
//
// exception.h: Rcpp integration of SMC library -- handling exceptions
//
// Copyright (C) 2008 - 2020  Adam Johansen
// Copyright (C) 2021         Adam Johansen and Ilya Zarubin
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
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with RcppSMC. If not, see <http://www.gnu.org/licenses/>.

//! \file
//! \brief This file defines exception-handling facilities.
//!
//! The smc::exception class, which is used for exception handling by SMCTC, is defined.

#ifndef __SMC_EXCEPT_HH
#define __SMC_EXCEPT_HH 1.0

#include <iostream>

///A macro which autocompletes the housekeeping components of an smc::exception
#define SMC_EXCEPTION(code,error) smc::exception(__FILE__, __LINE__, code, error)

///Exception thrown if a file cannot be accessed.
#define SMCX_FILE_NOT_FOUND 0x0020
///Exception thrown if the sampler attempts to access history data which wasn't stored.
#define SMCX_MISSING_HISTORY 0x0010
///Exception thrown if an attempt is made to instantiate a class of which a single instance is permitted more than once.
#define SMCX_MULTIPLE_INSTANTIATION 0x1000
///Exception thrown if an attempt is made to use facilities of the base sampler class related to MCMC moves.
#define CSMCX_USING_MCMC 0x0100
///Exception thrown if an attempt is made to use facilities of the base sampler class related to adaptation.
#define CSMCX_USING_ADAPTATION 0x0001

namespace smc {
  ///SMC Exception class

  /// This class holds details of unrecoverable errors which occur within the SMCTC library.
  /// An instance of it is thrown whenever such an error occurs.
  class exception {
  public:
    char const * szFile; //!< The source file from which the code generating the exception was generated.
    long lLine;   //!< The line of that source file which generates the exception.
    long lCode;   //!< A numerical code indicating the nature of the exception generated.
    char const * szMessage; //!< A human-readable explanation of the cause of the exception.

    //! Generate an SMCTC Exception class with the specified initialisation.

    //! This constructor fills the four elements of the class with their specified values.
    //! It is used to allow a single-line command to create and throw an exception.
    //!
    //! \param szN The name of the source file generating the exception.
    //! \param lL The line in that file responsible for the exception.
    //! \param lC The numerical code identifying the exception.
    //! \param szM An textual explanation of the problem.
    exception(char const * szN, long lL, long lC, char const * szM)
    {
      szFile = szN;
      lLine = lL;
      lCode = lC;
      szMessage = szM;
    }
  };
}



namespace std {
  /// Produce a human-readable display of the state of an smc::exception class using the stream operator.

  ///Display a human-readable version of an SMC exception.

  /// \param os The stream to write to.
  /// \param e The exception class to display.
  /// \return os
  inline std::ostream & operator<< (std::ostream & os, smc::exception & e)
  {
    os << "SMC Exception: #" << e.lCode << endl;
    os << "Error occured in file " << e.szFile << " at line " << e.lLine << "." << endl;
    os << "Details: " << endl << '\t' << e.szMessage << endl;
    return os;
  }
}

#endif
