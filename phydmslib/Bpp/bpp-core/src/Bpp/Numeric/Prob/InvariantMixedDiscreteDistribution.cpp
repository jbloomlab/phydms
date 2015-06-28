//
// File: InvariantMixedDiscreteDistribution.cpp
// Created by: Julien Dutheil
// Created on: Mon Dec 24 12:02 2007
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

#include "InvariantMixedDiscreteDistribution.h"
#include "../../Utils/MapTools.h"


using namespace bpp;
using namespace std;

/******************************************************************************/

InvariantMixedDiscreteDistribution::InvariantMixedDiscreteDistribution(
  DiscreteDistribution* dist, double p, double invariant) :
  AbstractParameterAliasable("Invariant."),
  AbstractDiscreteDistribution(1, "Invariant."),
  dist_(dist),
  invariant_(invariant),
  p_(p),
  nestedPrefix_(dist->getNamespace())
{
  // We first change the namespace of the nested distribution:
  dist_->setNamespace("Invariant." + nestedPrefix_);
  addParameters_(dist_->getIndependentParameters());
  addParameter_(new Parameter("Invariant.p", p, &Parameter::PROP_CONSTRAINT_IN));

  updateDistribution();
}

/******************************************************************************/

void InvariantMixedDiscreteDistribution::fireParameterChanged(const ParameterList& parameters)
{
  AbstractDiscreteDistribution::fireParameterChanged(parameters);
  p_ = getParameterValue("p");
  dist_->matchParametersValues(parameters);

  updateDistribution();
}

/******************************************************************************/

void InvariantMixedDiscreteDistribution::updateDistribution()
{
  distribution_.clear();
  bounds_.clear();

  size_t distNCat = dist_->getNumberOfCategories();
  vector<double> probs = dist_->getProbabilities();
  vector<double> cats  = dist_->getCategories();

  distribution_[invariant_] = p_;
  for (size_t i = 0; i < distNCat; i++)
  {
    if (cats[i] == invariant_)
      distribution_[invariant_] += (1. - p_) * probs[i];
    else
      distribution_[cats[i]] = (1. - p_) * probs[i];
  }

  intMinMax_.setLowerBound(dist_->getLowerBound(), !dist_->strictLowerBound());
  intMinMax_.setUpperBound(dist_->getUpperBound(), !dist_->strictUpperBound());

  if (invariant_ <= intMinMax_.getLowerBound())
    intMinMax_.setLowerBound(invariant_, true);
  if (invariant_ >= intMinMax_.getUpperBound())
    intMinMax_.setUpperBound(invariant_, true);

  numberOfCategories_ = distribution_.size();

  // bounds_

  // if invariant_ is between 2 values of dist_, bounds_ are set in the
  // middle of the 3 values

  bool nv = true;

  double a = dist_->getCategory(0), b;
  if (nv && (invariant_ < a))
  {
    bounds_.push_back((a + invariant_) / 2);
    nv = false;
  }

  for (size_t i = 1; i < distNCat; i++)
  {
    b = dist_->getCategory(i);
    if (nv && (invariant_ < b))
    {
      bounds_.push_back((a + invariant_) / 2);
      bounds_.push_back((invariant_ + b) / 2);
      nv = false;
    }
    else
      bounds_.push_back(dist_->getBound(i - 1));
    a = b;
  }

  if (nv)
    bounds_.push_back((a + invariant_) / 2);
}

/******************************************************************************/

void InvariantMixedDiscreteDistribution::setNamespace(const string& prefix)
{
  AbstractDiscreteDistribution::setNamespace(prefix);
  // We also need to update the namespace of the nested distribution:
  dist_->setNamespace(prefix + nestedPrefix_);
}

/******************************************************************************/

void InvariantMixedDiscreteDistribution::restrictToConstraint(const Constraint& c)
{
  if (!c.isCorrect(invariant_))
    throw ConstraintException("Impossible to restrict to Constraint", &getParameter_("p"), invariant_);

  dist_->restrictToConstraint(c);
  updateDistribution();
}
