//
// File: Simplex.cpp
// Created by: Laurent Guéguen
// Created on: mardi 31 mai 2011, à 13h 16
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

#include "Simplex.h"
#include "../NumConstants.h"

#include "../VectorTools.h"

using namespace bpp;
using namespace std;

Simplex::Simplex(const std::vector<double>& probas, unsigned short method, bool allowNull, const std::string& name) : AbstractParameterAliasable(name),
  dim_(probas.size()),
  method_(method),
  vProb_(),
  valpha_()
{
  if  (dim_==0)
    return;

  double sum = VectorTools::sum(probas);
  if (fabs(1. - sum) > NumConstants::SMALL())
    throw Exception("Simplex. Probabilities must equal 1 (sum =" + TextTools::toString(sum) + ").");

  const Constraint* pc = (allowNull ? &Parameter::PROP_CONSTRAINT_IN : &Parameter::PROP_CONSTRAINT_EX);

  for (unsigned int i = 0; i < dim_; i++)
  {
    vProb_.push_back(probas[i]);
  }

  double y = 1;
  switch (method_)
  {
  case 1:
    for (unsigned int i = 0; i < dim_ - 1; i++)
    {
      addParameter_(new Parameter(name + "theta" + TextTools::toString(i + 1), vProb_[i] / y, pc));
      y -= vProb_[i];
    }
    break;
  case 2:
    for (unsigned int i = 0; i < dim_ - 1; i++)
    {
      addParameter_(new Parameter(name + "theta" + TextTools::toString(i + 1), vProb_[i] / (vProb_[i] + vProb_[i + 1]), pc));
    }
    for (unsigned int i = 0; i < dim_ - 1; i++)
    {
      valpha_.push_back(vProb_[i + 1] / vProb_[i]);
    }
    break;
  case 3:
    for (size_t i = 1; i < dim_; i++)
    {
      size_t o = i;
      size_t li2 = 0; // rank of the strongest bit
      while (o)
      {
        li2++;
        o = o >> 1;
      }

      double i1 = 0, i0 = 0;
      size_t j = 0;
      size_t pi = i &  ~(1 << (li2 - 1));
      while (j < dim_)
      {
        size_t t = (j << li2) + pi;
        if (t >= dim_)
          break;
        else
          i0 += vProb_[t];
        t += (1 << (li2 - 1));
        if (t < dim_)
          i1 += vProb_[t];
        j++;
      }
      addParameter_(new Parameter(name + "theta" + TextTools::toString(i), i1 / (i0 + i1), pc));
    }
    break;
  }
}

Simplex::Simplex(size_t dim, unsigned short method, bool allowNull, const std::string& name) :
  AbstractParameterAliasable(name),
  dim_(dim),
  method_(method),
  vProb_(),
  valpha_()
{
  if  (dim_==0)
    return;
  
  for (size_t i = 0; i < dim_; i++)
  {
    vProb_.push_back(1. / static_cast<double>(dim_));
  }

  const Constraint* pc = (allowNull ? &Parameter::PROP_CONSTRAINT_IN : &Parameter::PROP_CONSTRAINT_EX);

  double y = 1;
  switch (method_)
  {
  case 1:
    for (unsigned int i = 0; i < dim_ - 1; i++)
    {
      addParameter_(new Parameter(name + "theta" + TextTools::toString(i + 1), vProb_[i] / y, pc));
      y -= vProb_[i];
    }
    break;
  case 2:
    for (unsigned int i = 0; i < dim_ - 1; i++)
    {
      addParameter_(new Parameter(name + "theta" + TextTools::toString(i + 1), 0.5, pc));
    }
    for (unsigned int i = 0; i < dim_ - 1; i++)
    {
      valpha_.push_back(1.);
    }
    break;
  case 3:
    for (unsigned int i = 0; i < dim_ - 1; i++)
    {
      addParameter_(new Parameter(name + "theta" + TextTools::toString(i + 1), 0.5, pc));
    }
    setFrequencies(vProb_);
    
    break;
  }
}

void Simplex::fireParameterChanged(const ParameterList& parameters)
{
  if  (dim_==0)
    return;

  AbstractParameterAliasable::fireParameterChanged(parameters);

  double x = 1.0;
  switch (method_)
  {
  case 1:
    double th;
    for (unsigned int i = 0; i < dim_ - 1; i++)
    {
      th = getParameterValue("theta" + TextTools::toString(i + 1));
      vProb_[i] = th * x;
      x *= 1 - th;
    }
    vProb_[dim_ - 1] = x;
    break;
  case 2:
    for (unsigned int i = 0; i < dim_ - 1; i++)
    {
      th = getParameterValue("theta" + TextTools::toString(i + 1));
      valpha_[i] = (1 - th) / th;
    }
    th = 1;
    vProb_[0] = 1;
    x = 1.0;
    for (unsigned int i = 0; i < dim_ - 1; i++)
    {
      th *= valpha_[i];
      vProb_[i + 1] = th;
      x += vProb_[i + 1];
    }
    for (unsigned int i = 0; i < dim_; i++)
    {
      vProb_[i] /= x;
    }

    break;
  case 3:
    size_t o = dim_;
    size_t ld2 = 0; // rank of the strongest bit
    while (o)
    {
      ld2++;
      o = o >> 1;
    }
    for (size_t i = 0; i < dim_; i++)
    {
      x = 1;
      size_t ld = ld2;
      size_t k = i;
      while (ld)
      {
        if (k >> (ld - 1))
          x *= getParameterValue("theta" + TextTools::toString(k));
        else
        {
          if ((k + (1 << (ld - 1))) < dim_)
            x *= 1 - getParameterValue("theta" + TextTools::toString(k + (1 << (ld - 1))));
        }
        k &= ~(1 << (--ld));
      }
      vProb_[i] = x;
    }
    break;
  }
}


void Simplex::setFrequencies(const std::vector<double>& probas)
{
  if  (dim_==0)
    return;

  double sum = VectorTools::sum(probas);
  if (fabs(1. - sum) > NumConstants::SMALL())
    throw Exception("Simplex::setFrequencies. Probabilities must equal 1 (sum =" + TextTools::toString(sum) + ").");

  double y = 1;

  ParameterList pl;
  switch (method_)
  {
  case 1:
    for (unsigned int i = 0; i < dim_ - 1; i++)
    {
      pl.addParameter(Parameter(getNamespace() + "theta" + TextTools::toString(i + 1), probas[i] / y));
      y -= probas[i];
    }
    break;
  case 2:
    for (unsigned int i = 0; i < dim_ - 1; i++)
    {
      pl.addParameter(Parameter(getNamespace() + "theta" + TextTools::toString(i + 1), probas[i] / (probas[i] + probas[i + 1])));
      valpha_[i] = probas[i + 1] / probas[i];
    }
    break;
  case 3:
    for (size_t i = 1; i < dim_; i++)
    {
      size_t o = i;
      size_t li2 = 0; // rank of the strongest bit
      while (o)
      {
        li2++;
        o = o >> 1;
      }

      double i1 = 0, i0 = 0;
      size_t j = 0;
      size_t pi = i &  ~(1 << (li2 - 1));
      while (j < dim_)
      {
        size_t t = (j << li2) + pi;
        if (t >= dim_)
          break;
        else
          i0 += probas[t];
        t += (1 << (li2 - 1));
        if (t < dim_)
          i1 += probas[t];
        j++;
      }
      pl.addParameter(Parameter(getNamespace() + "theta" + TextTools::toString(i), i1 / (i0 + i1)));
    }
    break;
  }

  matchParametersValues(pl);
}

