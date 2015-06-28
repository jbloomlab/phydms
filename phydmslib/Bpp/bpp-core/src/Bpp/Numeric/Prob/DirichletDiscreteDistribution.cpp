//
// File: BetaDiscreteDistribution.cpp
// Created by: Laurent Guéguen
// Created on: jeudi 2 septembre 2010, à 17h 03
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

#include "DirichletDiscreteDistribution.h"
#include "BetaDiscreteDistribution.h"
#include "../Random/RandomTools.h"
#include "../NumConstants.h"
#include "../../Text/TextTools.h"

using namespace bpp;

// From the STL:
#include <cmath>

using namespace std;

/** Constructor: **************************************************************/

DirichletDiscreteDistribution::DirichletDiscreteDistribution(vector<size_t> vn, Vdouble valpha) :
  AbstractParameterAliasable("Dirichlet."),
  vpBDD_()
{
  if (vn.size() <= 0 || vn.size() != valpha.size() - 1)
    throw Exception("Wrong number of categories for Dirichlet distribution: " + TextTools::toString(vn.size()));

  for (size_t j = 0; j < valpha.size(); j++)
  {
    addParameter_(new Parameter("Dirichlet.alpha_" + TextTools::toString(j + 1), valpha[j], new IntervalConstraint(1, 0.0001, true), true));
  }

  for (size_t j = 0; j < vn.size(); j++)
  {
    if (vn[j] <= 0)
      throw Exception("Wrong number of categories in Dirichlet distribution constructor: " + TextTools::toString(vn[j]));
  }

  for (size_t j = 0; j < vn.size(); j++)
  {
    vpBDD_.push_back(new BetaDiscreteDistribution(vn[j], 1, 1));
  }

  discretize(valpha);
}

DirichletDiscreteDistribution::~DirichletDiscreteDistribution()
{
  for (unsigned int i = 0; i < vpBDD_.size(); i++)
  {
    delete vpBDD_[i];
  }
}

/******************************************************************************/

void DirichletDiscreteDistribution::fireParameterChanged(const ParameterList& parameters)
{
  AbstractParameterAliasable::fireParameterChanged(parameters);
  applyParameters();
}

/******************************************************************************/

void DirichletDiscreteDistribution::applyParameters()
{
  Vdouble valpha;

  for (unsigned int j = 0; j < vpBDD_.size() + 1; j++)
  {
    valpha.push_back(getParameterValue("alpha_" + TextTools::toString(j + 1)));
  }

  discretize(valpha);
}

/******************************************************************************/

void DirichletDiscreteDistribution::discretize(Vdouble& valpha)
{
  /* discretization of Dirichlet distribution with equal proportions in
     each category */

  double x;
  ParameterList pL;
  pL.addParameter(Parameter("Beta.alpha", 1));
  pL.addParameter(Parameter("Beta.beta", 1));
  for (unsigned int j = 0; j < vpBDD_.size(); j++)
  {
    x = 0;
    for (unsigned int i = j + 1; i < valpha.size(); i++)
    {
      x += valpha[i];
    }

    pL.setParameterValue("Beta.alpha", valpha[j]);
    pL.setParameterValue("Beta.beta", x);

    vpBDD_[j]->matchParametersValues(pL);
  }
}


/******************************************************************************/

size_t DirichletDiscreteDistribution::getNumberOfCategories() const
{
  size_t n = 1;

  for (size_t j = 0; j < vpBDD_.size(); j++)
  {
    n *= vpBDD_[j]->getNumberOfCategories();
  }

  return n;
}

/******************************************************************************/

Vdouble DirichletDiscreteDistribution::getValueCategory(Vdouble& value) const
{
  if (value.size() != vpBDD_.size() + 1)
    throw Exception("Bad Vdouble parameter in DirichletDiscreteDistribution::getValueCategory");

  Vdouble vd;
  double y, sumc = 0;

  for (size_t j = 0; j < vpBDD_.size(); j++)
  {
    if (1 - sumc < NumConstants::VERY_TINY())
      y = vpBDD_[j]->getValueCategory(NumConstants::VERY_TINY());
    else
      y = vpBDD_[j]->getValueCategory(value[j] / (1 - sumc)) * (1 - sumc);
    sumc += y;
    vd.push_back(y);
  }
  vd.push_back(1 - sumc);
  return vd;
}

/******************************************************************************/

double DirichletDiscreteDistribution::getProbability(Vdouble& category) const
{
  if (category.size() != vpBDD_.size() + 1)
    throw Exception("Bad Vdouble parameter in DirichletDiscreteDistribution::getProbability");

  double sumc = 0;
  double p = 1;

  for (unsigned int j = 0; j < vpBDD_.size(); j++)
  {
    p *= vpBDD_[j]->getProbability(category[j] / (1 - sumc));
    sumc += category[j];
  }
  return p;
}

/******************************************************************************/

VVdouble DirichletDiscreteDistribution::getCategories() const
{
  VVdouble vvd1, vvd2;
  Vdouble vdj, vd;
  double sumc = 0;

  vdj = vpBDD_[0]->getCategories();
  for (unsigned int k = 0; k < vdj.size(); k++)
  {
    vd.push_back(vdj[k]);
    vvd1.push_back(vd);
    vd.pop_back();
  }

  for (unsigned int j = 1; j < vpBDD_.size(); j++)
  {
    vdj = vpBDD_[j]->getCategories();
    vvd2.clear();
    for (unsigned int i = 0; i < vvd1.size(); i++)
    {
      vd = vvd1[i];
      sumc = 0;
      for (unsigned int k = 0; k < vd.size(); k++)
      {
        sumc += vd[k];
      }
      for (unsigned int k = 0; k < vdj.size(); k++)
      {
        vd.push_back(vdj[k] * (1 - sumc));
        vvd2.push_back(vd);
        vd.pop_back();
      }
    }
    vvd1 = vvd2;
  }

  vvd2.clear();
  for (unsigned int i = 0; i < vvd1.size(); i++)
  {
    vd = vvd1[i];
    sumc = 0;
    for (unsigned int k = 0; k < vd.size(); k++)
    {
      sumc += vd[k];
    }
    vd.push_back(1 - sumc);
    vvd2.push_back(vd);
  }

  return vvd2;
}

/******************************************************************************/

Vdouble DirichletDiscreteDistribution::rand() const
{
  Vdouble vd;
  double x, sumc = 0;
  for (unsigned int j = 0; j < vpBDD_.size(); j++)
  {
    x = vpBDD_[j]->rand() * (1 - sumc);
    sumc += x;
    vd.push_back(x);
  }

  vd.push_back(1 - sumc);
  return vd;
}

/******************************************************************************/

Vdouble DirichletDiscreteDistribution::randC() const
{
  Vdouble vd;
  double x, sumc = 0;
  for (unsigned int j = 0; j < vpBDD_.size(); j++)
  {
    x = vpBDD_[j]->randC() * (1 - sumc);
    sumc += x;
    vd.push_back(x);
  }

  vd.push_back(1 - sumc);
  return vd;
}

/******************************************************************************/

