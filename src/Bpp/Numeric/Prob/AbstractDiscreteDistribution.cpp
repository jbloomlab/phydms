//
// File: AbstractDiscreteDistribution.cpp
// Created by: Julien Dutheil
// Created on: ?
//

/*
   Copyright or © or Copr. Bio++ Development Team, (November 19, 2004)

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

#include "AbstractDiscreteDistribution.h"

#include "../Random/RandomTools.h"
#include "../VectorTools.h"

using namespace bpp;
using namespace std;


AbstractDiscreteDistribution::AbstractDiscreteDistribution(size_t nbClasses, const std::string& prefix) :
  AbstractParameterAliasable(prefix),
  numberOfCategories_(nbClasses),
  distribution_(),
  bounds_(nbClasses-1),
  intMinMax_(-NumConstants::VERY_BIG(), NumConstants::VERY_BIG(), true, true),
  median_(false)
{}

AbstractDiscreteDistribution::AbstractDiscreteDistribution(size_t nbClasses, double delta, const std::string& prefix) :
  AbstractParameterAliasable(prefix),
  numberOfCategories_(nbClasses),
  distribution_(Order(delta)),
  bounds_(nbClasses-1),
  intMinMax_(-NumConstants::VERY_BIG(), NumConstants::VERY_BIG(),true, true),
  median_(false)
{}

AbstractDiscreteDistribution::AbstractDiscreteDistribution(const AbstractDiscreteDistribution& adde) :
  AbstractParameterAliasable(adde),
  numberOfCategories_(adde.numberOfCategories_),
  distribution_(adde.distribution_),
  bounds_(adde.bounds_),
  intMinMax_(adde.intMinMax_),
  median_(adde.median_)
{
}

AbstractDiscreteDistribution& AbstractDiscreteDistribution::operator=(const AbstractDiscreteDistribution& adde) 
{
  AbstractParameterAliasable::operator=(adde);
  numberOfCategories_=adde.numberOfCategories_;
  distribution_=adde.distribution_;
  bounds_=adde.bounds_;
  intMinMax_=adde.intMinMax_;
  median_=adde.median_;

  return *this;
}

/******************************************************************************/

size_t AbstractDiscreteDistribution::getNumberOfCategories() const
{
  return numberOfCategories_;
}

void AbstractDiscreteDistribution::setNumberOfCategories(size_t nbClasses)
{
  if (nbClasses <= 0)
    cerr << "DEBUG: ERROR!!! Number of categories is <= 0 in AbstractDiscreteDistribution::setNumberOfCategories()." << endl;

  if (numberOfCategories_ != nbClasses)
  {
    numberOfCategories_ = nbClasses;
    discretize();
  }
}


/******************************************************************************/

double AbstractDiscreteDistribution::getCategory(size_t categoryIndex) const
{
  map<double, double>::const_iterator it = distribution_.begin();
  for (unsigned int i = 0; i < categoryIndex; i++)
  {
    it++;
  }
  return it->first;
}

/******************************************************************************/

double AbstractDiscreteDistribution::getProbability(size_t categoryIndex) const
{
  map<double, double>::const_iterator it = distribution_.begin();
  for (unsigned int i = 0; i < categoryIndex; i++)
  {
    it++;
  }
  return it->second;
}

/******************************************************************************/

double AbstractDiscreteDistribution::getProbability(double category) const
{
  return distribution_.find(category)->second;
}

/******************************************************************************/

Vdouble AbstractDiscreteDistribution::getCategories() const
{
  Vdouble result(distribution_.size());
  unsigned int i = 0;
  for (map<double, double>::const_iterator it = distribution_.begin();
       it != distribution_.end();
       it++)
  {
    result[i] = it->first;
    i++;
  }
  return result;
}

/******************************************************************************/

Vdouble AbstractDiscreteDistribution::getProbabilities() const
{
  Vdouble result(distribution_.size());
  size_t i = 0;
  for (map<double, double>::const_iterator it = distribution_.begin();
       it != distribution_.end();
       it++)
  {
    result[i] = it->second;
    i++;
  }
  return result;
}

/******************************************************************************/

void AbstractDiscreteDistribution::set(double category, double probability)
{
  distribution_[category] = probability;
}

/******************************************************************************/

void AbstractDiscreteDistribution::add(double category, double probability)
{
  if (distribution_.find(category) == distribution_.end())
  {
    // new category
    distribution_[category] = probability;
  }
  else
  {
    // existing category
    distribution_[category] += probability;
  }
}

/******************************************************************************/

double AbstractDiscreteDistribution::rand() const
{
  double r = RandomTools::giveRandomNumberBetweenZeroAndEntry(1);
  double cumprob = 0;
  for (map<double, double>::const_iterator i = distribution_.begin();
       i != distribution_.end();
       i++)
  {
    cumprob += i->second;
    if (r <= cumprob)
      return i->first;
  }
  // This line can't be reached:
  return -1.;
}

/******************************************************************************/

double AbstractDiscreteDistribution::getInfCumulativeProbability(double category) const
{
  double prob = 0;
  map<double, double>::const_iterator it = distribution_.find(category);
  for (map<double, double>::const_iterator i = distribution_.begin();
       i != it;
       i++)
  {
    prob += i->second;
  }
  return prob;
}

/******************************************************************************/

double AbstractDiscreteDistribution::getIInfCumulativeProbability(double category) const
{
  double prob = 0;
  map<double, double>::const_iterator it = distribution_.find(category);
  if (it == distribution_.end())
    return 0;
  for (map<double, double>::const_iterator i = ++it;
       i != distribution_.end();
       i++)
  {
    prob += i->second;
  }
  return 1. - prob;
}

/******************************************************************************/

double AbstractDiscreteDistribution::getSupCumulativeProbability(double category) const
{
  double prob = 0;
  map<double, double>::const_iterator it = distribution_.find(category);
  if (it == distribution_.end())
    return 0;
  for (map<double, double>::const_iterator i = ++it;
       i != distribution_.end();
       i++)
  {
    prob += i->second;
  }
  return prob;
}

/******************************************************************************/

double AbstractDiscreteDistribution::getSSupCumulativeProbability(double category) const
{
  double prob = 0;
  map<double, double>::const_iterator it = distribution_.find(category);
  for (map<double, double>::const_iterator i = distribution_.begin();
       i != it;
       i++)
  {
    prob += i->second;
  }
  return 1. - prob;
}

/******************************************************************************/

void AbstractDiscreteDistribution::print(OutputStream& out) const
{
  for (map<double, double>::const_iterator i = distribution_.begin(); i != distribution_.end(); i++)
  {
    (out << "Pr(" << (i->first) << ") = " << (i->second)).endLine();
  }
}

/******************************************************************************/

double AbstractDiscreteDistribution::getValueCategory(double value) const
{
  if (!(intMinMax_.isCorrect(value)))
    throw Exception("AbstractDiscreteDistribution::getValueCategory out of bounds:" + TextTools::toString(value));

  map<double, double>::const_iterator it = distribution_.begin();
  for (unsigned int i=1;i<bounds_.size();i++)
    if (value<bounds_[i])
      break;
    else
      it++;

  return it->first;
}

/***********************************************************************/


void AbstractDiscreteDistribution::discretize()
{
  /* discretization of distribution with equal proportions in each
     category
   */

  distribution_.clear();
  bounds_.resize(numberOfCategories_ - 1);

  double minX = pProb(intMinMax_.getLowerBound());
  double maxX = pProb(intMinMax_.getUpperBound());

  double ec;
  size_t i;
  vector<double> values(numberOfCategories_);

  if (maxX != minX)
  {
    // divide the domain into equiprobable intervals
    ec = (maxX - minX) / static_cast<double>(numberOfCategories_);

    for (i = 1; i < numberOfCategories_; i++)
    {
      bounds_[i-1] = qProb(minX + static_cast<double>(i) * ec);
    }

    // for each category, sets the value v as the median, adjusted
    //      such that the sum of the values = 1
    if (median_)
    {
      double t=0;
      for (i = 0; i < numberOfCategories_; i++)
        values[i] = qProb(minX + (static_cast<double>(i) + 0.5) * ec);

      for (i = 0, t = 0; i < numberOfCategories_; i++)
        t += values[i];

      double mean = Expectation(intMinMax_.getUpperBound()) - Expectation(intMinMax_.getLowerBound());

      for (i = 0; i < numberOfCategories_; i++)
        values[i] *= mean / t / ec;
    }
    else
      // for each category, sets the value v such that
      //      v * length_of_the_interval = the surface of the category
      {
      double a = Expectation(intMinMax_.getLowerBound()), b;
      for (i = 0; i < numberOfCategories_-1; i++)
        {
          b = Expectation(bounds_[i]);
          values[i] = (b - a) / ec;
          a = b;
        }
      values[numberOfCategories_-1] = (Expectation(intMinMax_.getUpperBound())-a) / ec;
    }
  }
  else 
    // if maxX==minX, uniform discretization of the range
  {
    ec = (intMinMax_.getUpperBound() - intMinMax_.getLowerBound()) / static_cast<double>(numberOfCategories_);
    for (i = 1; i < numberOfCategories_; i++)
      bounds_[i-1] = intMinMax_.getLowerBound() + static_cast<double>(i) * ec;

    values[0] = (intMinMax_.getLowerBound() + bounds_[0]) / 2;
    
    for (i = 1; i < numberOfCategories_-1; i++)
      values[i] = (bounds_[i-1] + bounds_[i]) / 2;

    values[numberOfCategories_-1] = (intMinMax_.getUpperBound()+bounds_[numberOfCategories_ - 1]) / 2;
  }

  // adjustments near the boundaries of the domain, according to the precision chosen
  if (intMinMax_.strictLowerBound())
  {
    for (i = 0; i < numberOfCategories_; i++)
    {
      if (values[i] < intMinMax_.getLowerBound() + precision())
        values[i] = intMinMax_.getLowerBound() + precision();
      else
        break;
    }
  }
  else
  {
    for (i = 0; i < numberOfCategories_; i++)
    {
      if (values[i] < intMinMax_.getLowerBound())
        values[i] = intMinMax_.getLowerBound() + precision();
      else
        break;
    }
  }

  if (intMinMax_.strictUpperBound())
  {
    for (i = numberOfCategories_; i > 0; i--)
    {
      if (values[i-1] > intMinMax_.getUpperBound() - precision())
        values[i-1] = intMinMax_.getUpperBound() - precision();
      else
        break;
    }
  }
  else
  {
    for (i = numberOfCategories_; i > 0; i--)
    {
      if (values[i-1] > intMinMax_.getUpperBound())
        values[i-1] = intMinMax_.getUpperBound() - precision();
      else
        break;
    }
  }

  // now the distribution_ map, taking care that all values are different
  
  double p = 1. / static_cast<double>(numberOfCategories_);
  for (i = 0; i < numberOfCategories_; i++)
  {
    if (distribution_.find(values[i]) != distribution_.end())
    {
      int j = 1;
      int f = ((values[i] + NumConstants::TINY()) >= intMinMax_.getUpperBound()) ? -1 : 1;
      while (distribution_.find(values[i] + f * j * precision()) != distribution_.end())
      {
        j++;
        f = ((values[i] + f * j * precision()) >= intMinMax_.getUpperBound()) ? -1 : 1;
      }
      distribution_[values[i] + f * j * precision()] = p;
    }
    else
      distribution_[values[i]] = p;
  }

  return;
}

Vdouble AbstractDiscreteDistribution::getBounds() const
{
  Vdouble vb(numberOfCategories_ + 1);
  vb[0] = getLowerBound();
  for (unsigned int i = 0; i < numberOfCategories_ - 1; i++)
    vb[i + 1] = getBound(i);
  vb[numberOfCategories_] = getUpperBound();
  return vb;
}

void AbstractDiscreteDistribution::restrictToConstraint(const Constraint& c)
{
  const IntervalConstraint* pi = dynamic_cast<const IntervalConstraint*>(&c);

  if (!pi)
    throw Exception("AbstractDiscreteDistribution::restrictToConstraint: the constraint is not an interval");

  if (!(intMinMax_ <= (*pi)))
  {
    intMinMax_ &= c;
    discretize();
  }
}
