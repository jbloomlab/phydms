//
// File: GaussianDiscreteDistribution.cpp
// Created by: Laurent Guéguen
// Created on: April 2010
//

/*
Copyright or © or Copr. Bio++ Development Team, (November 17, 2004)

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

#include "GaussianDiscreteDistribution.h"
#include "../Random/RandomTools.h"
#include "../NumConstants.h"
#include "../../Utils/MapTools.h"

using namespace bpp;

// From the STL:
#include <cmath>

using namespace std;

/** Constructor: **************************************************************/

GaussianDiscreteDistribution::GaussianDiscreteDistribution(size_t n, double mu, double sigma) :
  AbstractParameterAliasable("Gaussian."),
  AbstractDiscreteDistribution(n,"Gaussian."), mu_(mu), sigma_(sigma)
{
  addParameter_(new Parameter("Gaussian.mu", mu));
  addParameter_(new Parameter("Gaussian.sigma", sigma, &Parameter::R_PLUS_STAR));
  discretize();
}

GaussianDiscreteDistribution::GaussianDiscreteDistribution(const GaussianDiscreteDistribution& gdd) :
  AbstractParameterAliasable(gdd),
  AbstractDiscreteDistribution(gdd),
  mu_(gdd.mu_),
  sigma_(gdd.sigma_)
{
}

GaussianDiscreteDistribution& GaussianDiscreteDistribution::operator=(const GaussianDiscreteDistribution& gdd) 
{
  AbstractParameterAliasable::operator=(gdd);
  AbstractDiscreteDistribution::operator=(gdd);
  mu_=gdd.mu_;
  sigma_=gdd.sigma_;

  return *this;
}

GaussianDiscreteDistribution::~GaussianDiscreteDistribution() {}

/******************************************************************************/

void GaussianDiscreteDistribution::fireParameterChanged(const ParameterList& parameters)
{
  AbstractDiscreteDistribution::fireParameterChanged(parameters);
  mu_ = getParameterValue("mu");
  sigma_ = getParameterValue("sigma");
  discretize();  
}

/******************************************************************************/

double GaussianDiscreteDistribution::qProb(double x) const
{
  return RandomTools::qNorm(x, mu_, sigma_);
}

double GaussianDiscreteDistribution::pProb(double x) const
{
  return RandomTools::pNorm(x, mu_, sigma_);
}

double GaussianDiscreteDistribution::Expectation(double a) const
{
  return -sigma_/sqrt(2*M_PI)*exp(-pow((a-mu_)/sigma_,2)/2)
    +mu_*RandomTools::pNorm(a,mu_,sigma_);
}
