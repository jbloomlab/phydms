//
// File: MixtureOfDiscreteDistributions
// Created by: Laurent Guéguen
// Created on: mercredi 9 juin 2010, à 23h 09
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

#include "MixtureOfDiscreteDistributions.h"
#include "../NumConstants.h"
#include "../../Utils/MapTools.h"
#include "../../Text/TextTools.h"

using namespace bpp;
using namespace std;

MixtureOfDiscreteDistributions::MixtureOfDiscreteDistributions(const vector<DiscreteDistribution*>& distributions,
                                                               const vector<double>& probas ) :
  AbstractParameterAliasable("Mixture."),
  AbstractDiscreteDistribution(1, "Mixture."),
  vdd_(),
  probas_(),
  vNestedPrefix_()
{
  if (distributions.size() != probas.size())
  {
    throw Exception("MixtureOfDiscreteDistributions. Distributions and probabilities vectors must have the same size (" + TextTools::toString(distributions.size()) + " != " + TextTools::toString(probas.size()) + ").");
  }

  size_t size = distributions.size();
  for (size_t i = 0; i < size; i++)
  {
    if (distributions[i] == 0)
      throw Exception("MixtureOfDiscreteDistributions. Null DiscreteDistribution* in argument vector at index " + TextTools::toString(i));
  }

  for (size_t i = 0; i < size; i++)
  {
    probas_.push_back(probas[i]);
  }

  double sum = VectorTools::sum(probas);
  if (fabs(1. - sum) > precision())
    throw Exception("MixtureOfDiscreteDistributions. Probabilities must equal 1 (sum =" + TextTools::toString(sum) + ").");

  double y = 1;
  for (size_t i = 0; i < size - 1; i++)
  {
    addParameter_(new Parameter("Mixture.theta" + TextTools::toString(i + 1), probas[i] / y, &Parameter::PROP_CONSTRAINT_IN));
    y -= probas[i];
  }


  for (size_t i = 0; i < size; i++)
  {
    vdd_.push_back(distributions[i]->clone());
  }

  //  Parameters

  for (size_t i = 0; i < size; i++)
  {
    vNestedPrefix_.push_back(TextTools::toString(i + 1) + "_" + distributions[i]->getNamespace());
  }

  for (size_t i = 0; i < size; i++)
  {
    vdd_[i]->setNamespace("Mixture." + vNestedPrefix_[i]);
  }

  for (size_t i = 0; i < size; i++)
  {
    addParameters_(vdd_[i]->getParameters());
  }

  updateDistribution();
}

MixtureOfDiscreteDistributions::MixtureOfDiscreteDistributions(const MixtureOfDiscreteDistributions& mdd) :
  AbstractParameterAliasable(mdd),
  AbstractDiscreteDistribution(mdd),
  vdd_(),
  probas_(),
  vNestedPrefix_()
{
  for (size_t i = 0; i < mdd.vdd_.size(); i++)
  {
    probas_.push_back(mdd.probas_[i]);
    vdd_.push_back(mdd.vdd_[i]->clone());
    vNestedPrefix_.push_back(mdd.vNestedPrefix_[i]);
  }
}

MixtureOfDiscreteDistributions& MixtureOfDiscreteDistributions::operator=(const MixtureOfDiscreteDistributions& mdd)
{
  AbstractParameterAliasable::operator=(mdd);
  AbstractDiscreteDistribution::operator=(mdd);
  vdd_.clear();
  probas_.clear();
  vNestedPrefix_.clear();

  for (size_t i = 0; i < mdd.vdd_.size(); i++)
  {
    probas_.push_back(mdd.probas_[i]);
    vdd_.push_back(mdd.vdd_[i]->clone());
    vNestedPrefix_.push_back(mdd.vNestedPrefix_[i]);
  }

  return *this;
}

MixtureOfDiscreteDistributions::~MixtureOfDiscreteDistributions()
{
  for (size_t i = 0; i < vdd_.size(); i++)
  {
    delete vdd_[i];
  }

  vdd_.clear();
}

void MixtureOfDiscreteDistributions::setNumberOfCategories(size_t nbClasses)
{
  for (size_t i = 0; i < vdd_.size(); i++)
  {
    vdd_[i]->setNumberOfCategories(nbClasses);
  }

  updateDistribution();
}


void MixtureOfDiscreteDistributions::fireParameterChanged(const ParameterList& parameters)
{
  AbstractDiscreteDistribution::fireParameterChanged(parameters);
  size_t size = vdd_.size();
  double x = 1.0;
  for (size_t i = 0; i < size - 1; i++)
  {
    probas_[i] = getParameterValue("theta" + TextTools::toString(i + 1)) * x;
    x *= 1 - getParameterValue("theta" + TextTools::toString(i + 1));
  }

  probas_[size - 1] = x;


  for (size_t i = 0; i < size; i++)
  {
    vdd_[i]->matchParametersValues(parameters);
  }

  updateDistribution();
}

void MixtureOfDiscreteDistributions::updateDistribution()
{
  size_t size = vdd_.size();
  distribution_.clear();
  // calculation of distribution

  for (size_t i = 0; i < size; i++)
  {
    vector<double> values = vdd_[i]->getCategories();
    for (size_t j = 0; j < values.size(); j++)
    {
      distribution_[values[j]] = 0;
    }
  }

  for (size_t i = 0; i < size; i++)
  {
    vector<double> values = vdd_[i]->getCategories();
    vector<double> probas2 = vdd_[i]->getProbabilities();
    for (size_t j = 0; j < values.size(); j++)
    {
      distribution_[values[j]] += probas2[j] * probas_[i];
    }
  }

  numberOfCategories_ = distribution_.size();

  // intMinMax_

  double uB, lB;
  uB = -NumConstants::VERY_BIG();
  lB = NumConstants::VERY_BIG();

  bool suB = true, slB = true;

  for (size_t i = 0; i < size; i++)
  {
    if (vdd_[i]->getLowerBound() <= lB)
    {
      lB = vdd_[i]->getLowerBound();
      slB = vdd_[i]->strictLowerBound();
    }
    if (vdd_[i]->getUpperBound() >= uB)
    {
      uB = vdd_[i]->getUpperBound();
      suB = vdd_[i]->strictUpperBound();
    }
  }

  intMinMax_.setLowerBound(lB, slB);
  intMinMax_.setUpperBound(uB, suB);

  // Compute midpoint bounds_:
  vector<double> values = MapTools::getKeys<double, double, AbstractDiscreteDistribution::Order>(distribution_);

  bounds_.resize(numberOfCategories_ - 1);

  // Fill from 0 to numberOfCategories_-1 with midpoints:
  for (size_t i = 0; i < numberOfCategories_ - 1; i++)
  {
    bounds_[i] = (values[i] + values[i + 1]) / 2.;
  }
}

void MixtureOfDiscreteDistributions::setMedian(bool median)
{
  if (median_ != median)
  {
    median_ = median;
    for (size_t i = 0; i < vdd_.size(); i++)
    {
      vdd_[i]->setMedian(median);
    }
    updateDistribution();
  }
}
void MixtureOfDiscreteDistributions::discretize()
{
  for (size_t i = 0; i < vdd_.size(); i++)
  {
    vdd_[i]->discretize();
  }

  updateDistribution();
}

double MixtureOfDiscreteDistributions::pProb(double x) const
{
  double s = 0;
  for (size_t i = 0; i < vdd_.size(); i++)
  {
    s += probas_[i] * vdd_[i]->pProb(x);
  }
  return s;
}

double MixtureOfDiscreteDistributions::qProb(double x) const
{
  throw Exception("MixtureOfDiscreteDistributions::qProb to difficult to compute: not implemented");
  return 0;
}

double MixtureOfDiscreteDistributions::Expectation(double a) const
{
  double s = 0;
  for (size_t i = 0; i < vdd_.size(); i++)
  {
    s += probas_[i] * vdd_[i]->Expectation(a);
  }
  return s;
}

void MixtureOfDiscreteDistributions::restrictToConstraint(const Constraint& c)
{
  for (size_t i = 0; i < vdd_.size(); i++)
  {
    vdd_[i]->restrictToConstraint(c);
  }

  updateDistribution();
}

/******************************************************************************/

void MixtureOfDiscreteDistributions::setNamespace(const string& prefix)
{
  AbstractDiscreteDistribution::setNamespace(prefix);
  // We also need to update the namespace of the nested distributions:
  for (size_t i = 0; i < vdd_.size(); i++)
  {
    vdd_[i]->setNamespace(prefix + vNestedPrefix_[i]);
  }
}

