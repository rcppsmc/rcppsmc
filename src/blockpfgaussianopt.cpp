#include "smctc.h"
#include "blockpfgaussianopt.h"
#include "rngR.h"

#include <cstdlib>
#include <cmath>
#include <vector>

using namespace std;

///The observations
namespace BSPFG {
Rcpp::NumericVector y;
}

long lLag = 1;

extern "C" SEXP blockpfGaussianOpt(SEXP dataS, SEXP partS, SEXP lagS)
{
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

using namespace std;
using BSPFG::y;

/// \param pRng A pointer to the random number generator which is to be used
smc::particle<vector<double> > fInitialiseBSPFG(smc::rng *pRng)
{
  vector<double> value;
  
  value.push_back(pRng->Normal(0.5 * y[0],1.0/sqrt(2.0)));

  return smc::particle<vector<double> >(value,1.0);
}

///The proposal function.

///\param lTime The sampler iteration.
///\param pFrom The particle to move.
///\param pRng  A random number generator.
void fMoveBSPFG(long lTime, smc::particle<vector<double> > & pFrom, smc::rng *pRng)
{
    std::vector<double> * cv_to = pFrom.GetValuePointer();

    if(lTime == 1) {
	cv_to->push_back((cv_to->at(lTime-1) + y[int(lTime)])/2.0 + pRng->Normal(0.0,1.0/sqrt(2.0)));

	pFrom.AddToLogWeight(-0.25*(y[int(lTime)] - cv_to->at(lTime-1))*(y[int(lTime)]-cv_to->at(lTime-1)));

	return;
    }

    long lag = min(lTime,lLag);

    //These structures should really be made static 
    std::vector<double> mu(lag+1);
    std::vector<double> sigma(lag+1);
    std::vector<double> sigmah(lag+1);
    std::vector<double> mub(lag+1);

    // Forward filtering
    mu[0] = cv_to->at(lTime-lag);
    sigma[0] = 0;
    for(int i = 1; i <= lag; ++i)
    {
	sigmah[i] = sigma[i-1] + 1;
	
	mu[i] = (sigmah[i] * y[int(lTime-lag+i)] +  mu[i-1]) / (sigmah[i] + 1);
	sigma[i] = sigmah[i] / (sigmah[i] + 1);
    }
    // Backward smoothing
    mub[lag] = mu[lag];
    cv_to->push_back(pRng->Normal(mub[lag],sqrt(sigma[lag])));
    for(int i = lag-1; i; --i)
    {
	mub[i] = (sigma[i]*cv_to->at(lTime-lag+i+1) + mu[i]) / (sigma[i]+1);
	cv_to->at(lTime-lag+i) = pRng->Normal(mub[i],sqrt(sigma[lag]/(sigma[lag] + 1)));
    }
    
    // Importance weighting
    pFrom.AddToLogWeight(-0.5 * pow(y[int(lTime)] - mu[lag-1],2.0) / (sigmah[lag]+1) );

}
