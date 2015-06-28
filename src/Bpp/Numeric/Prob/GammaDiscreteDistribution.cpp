//
// File: GammaDiscreteDistribution.cpp
// Created by: Julien Dutheil
// Created on: Sun Oct 26 20:36:12 2003
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

#include "GammaDiscreteDistribution.h"
#include "../Random/RandomTools.h"
#include "../NumConstants.h"
#include "../../Utils/MapTools.h"

using namespace bpp;

// From the STL:
#include <cmath>

using namespace std;

/** Constructor: **************************************************************/

GammaDiscreteDistribution::GammaDiscreteDistribution(size_t n, double alpha, double beta, double minimumAlpha, double minimumBeta, bool paramOffset, double offset) :
  AbstractParameterAliasable("Gamma."),
  AbstractDiscreteDistribution(n, "Gamma."),
  alpha_(alpha),
  beta_(beta),
  offset_(offset),
  ga1_(1)
{
  // We use a lower bound of 0.0001 for alpha and beta to prohibe errors due to computer
  // floating precision: if alpha is quite low (gamma -> constant), some classes
  // may have the same category value, leading to a classe number lower than expected.
  // NB: if this is the case, then a warning is shown. This may happen in optimization
  // algorithms.
  addParameter_(new Parameter("Gamma.alpha", alpha, new IntervalConstraint(1, minimumAlpha, true), true));
  addParameter_(new Parameter("Gamma.beta", beta, new IntervalConstraint(1, minimumBeta, true), true));
  if (paramOffset)
    addParameter_(new Parameter("Gamma.offset", offset));
  
  ga1_ = exp(RandomTools::lnGamma(alpha_ + 1) - RandomTools::lnGamma(alpha_));

  intMinMax_.setLowerBound(offset_, true);
  discretize();
}

GammaDiscreteDistribution::GammaDiscreteDistribution(const GammaDiscreteDistribution& gdd) :
  AbstractParameterAliasable(gdd),
  AbstractDiscreteDistribution(gdd),
  alpha_(gdd.alpha_),
  beta_(gdd.beta_),
  offset_(gdd.offset_),
  ga1_(gdd.ga1_)
{
}

GammaDiscreteDistribution& GammaDiscreteDistribution::operator=(const GammaDiscreteDistribution& gdd)
{
  AbstractParameterAliasable::operator=(gdd);
  AbstractDiscreteDistribution::operator=(gdd);
  alpha_=gdd.alpha_;
  beta_=gdd.beta_;
  offset_=gdd.offset_;
  ga1_=gdd.ga1_;

  return *this;
}

GammaDiscreteDistribution::~GammaDiscreteDistribution() {}

/******************************************************************************/

void GammaDiscreteDistribution::fireParameterChanged(const ParameterList& parameters)
{
  AbstractDiscreteDistribution::fireParameterChanged(parameters);
  alpha_ = getParameterValue("alpha");
  beta_ = getParameterValue("beta");
  if (hasParameter("offset"))
      offset_ = getParameterValue("offset");
  ga1_ = exp(RandomTools::lnGamma(alpha_ + 1) - RandomTools::lnGamma(alpha_));

  discretize();
}

/******************************************************************************/

// Adapted from function DiscreteGamma of Yang

double GammaDiscreteDistribution::qProb(double x) const
{
  return offset_ + RandomTools::qGamma(x, alpha_, beta_);
}


double GammaDiscreteDistribution::pProb(double x) const
{
  return RandomTools::pGamma(x-offset_, alpha_, beta_);
}

double GammaDiscreteDistribution::Expectation(double a) const
{
  return RandomTools::pGamma(a-offset_, alpha_ + 1, beta_) / beta_ * ga1_ + (offset_ ? offset_ * RandomTools::pGamma(a - offset_, alpha_, beta_):0);
}

