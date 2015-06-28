//
// File: BetaDiscreteDistribution.cpp
// Created by: Laurent Guéguen
// Created on: lundi 31 mai 2010, à 11h 15
//

/*
Copyright or © or Copr. CNRS, (November 17, 2004)

This software is a computer program whose purpose is to provide classes
for numerical calculus.

This software is governed by the CeCILL  license under French law and
abiding by the rules of distribution of free software.  You can  use, 
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info". 

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability. 

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or 
data to be ensured and,  more generally, to use and operate it in the 
same conditions as regards security. 

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.
*/

#include "BetaDiscreteDistribution.h"
#include "../Random/RandomTools.h"
#include "../NumConstants.h"
#include "../../Utils/MapTools.h"

using namespace bpp;

// From the STL:
#include <cmath>

using namespace std;

/** Constructor: **************************************************************/

BetaDiscreteDistribution::BetaDiscreteDistribution(size_t n, double alpha, double beta) :
  AbstractParameterAliasable("Beta."),
  AbstractDiscreteDistribution(n,NumConstants::VERY_TINY(),"Beta."), alpha_(alpha), beta_(beta), diffln_(0)
{
  addParameter_(new Parameter("Beta.alpha", alpha, new IntervalConstraint(1, 0.0001, true), true));

  // For precision issues, beta cannot be too low
  addParameter_(new Parameter("Beta.beta", beta, new IntervalConstraint(1, 0.1, true), true));
  intMinMax_.setLowerBound(0, true);
  intMinMax_.setUpperBound(1, true);

  diffln_ = exp(RandomTools::lnBeta(alpha_ + 1, beta_) - RandomTools::lnBeta(alpha_, beta_));
  discretize();
}

BetaDiscreteDistribution::BetaDiscreteDistribution(const BetaDiscreteDistribution& bdd) :
  AbstractParameterAliasable(bdd),
  AbstractDiscreteDistribution(bdd), alpha_(bdd.alpha_), beta_(bdd.beta_), diffln_(bdd.diffln_)
{
}

BetaDiscreteDistribution& BetaDiscreteDistribution::operator=(const BetaDiscreteDistribution& bdd)
{
  AbstractParameterAliasable::operator=(bdd);
  AbstractDiscreteDistribution::operator=(bdd);

  alpha_=bdd.alpha_;
  beta_=bdd.beta_;
  diffln_=bdd.diffln_;
  
  return *this;
}
  
/******************************************************************************/

void BetaDiscreteDistribution::fireParameterChanged(const ParameterList& parameters)
{
  AbstractDiscreteDistribution::fireParameterChanged(parameters);
  alpha_=getParameterValue("alpha");
  beta_=getParameterValue("beta");

  if (alpha_<=1 && intMinMax_.getLowerBound()==0)
    intMinMax_.setLowerBound(precision(),false);

  if (beta_<=1 && intMinMax_.getUpperBound()==1)
    intMinMax_.setUpperBound(1-precision(),false);

  diffln_=exp(RandomTools::lnBeta(alpha_+1,beta_)-RandomTools::lnBeta(alpha_,beta_));
  discretize();
}

/******************************************************************************/

double BetaDiscreteDistribution::qProb(double x) const
{
  return RandomTools::qBeta(x, alpha_, beta_);
}

double BetaDiscreteDistribution::pProb(double x) const 
{
  return RandomTools::pBeta(x, alpha_, beta_);
}

double BetaDiscreteDistribution::Expectation(double a) const
{
  return RandomTools::pBeta(a,alpha_+1,beta_)*diffln_;
}

    
