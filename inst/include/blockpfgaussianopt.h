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
