#include "smctc.hh"
#include "blockpfgaussianopt.hh"
#include <cstdio> 
#include <cstdlib>
#include <cstring>

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

using namespace std;

///The observations
//double *y;
namespace BSPFG {
Rcpp::NumericVector y;
};

long lLag = 1;
long load_data(char const * szName, double** y);

extern "C" SEXP blockpfGaussianOpt(SEXP dataS, SEXP partS, SEXP lagS)
//int main(int argc, char** argv)
{
    //Load observations
    long lIterates;
    long lNumber = Rcpp::as<long>(partS);
    lLag = Rcpp::as<long>(lagS);

    y = Rcpp::NumericVector(dataS);
    lIterates = y.size();

    //Initialise and run the sampler
    smc::sampler<vector<double> > Sampler(lNumber, SMC_HISTORY_NONE);  
    smc::moveset<vector<double> > Moveset(fInitialiseBSPFG, fMoveBSPFG, NULL);
	
    Sampler.SetResampleParams(SMC_RESAMPLE_SYSTEMATIC, 0.5);
    Sampler.SetMoveSet(Moveset);
    
    Sampler.Initialise();
    Sampler.IterateUntil(lIterates - 1);

    //Generate results
    Rcpp::NumericMatrix resValues = Rcpp::NumericMatrix(lNumber,lIterates);
    Rcpp::NumericVector resWeights = Rcpp::NumericVector(lNumber);
    for(int i = 0; i < lNumber; ++i) 
    {
	vector<double> pValue = Sampler.GetParticleValue(i);
	for(int j = 0; j < lIterates; ++j) {
	    resValues(i,j) = pValue.at(j);
	}
	resWeights(i) = Sampler.GetParticleWeight(i);
    }

    return Rcpp::List::create(Rcpp::_["weight"] = resWeights, Rcpp::_["values"] = resValues);
}

