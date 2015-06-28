//
// File: GranthamAAChemicalDistance.cpp
// Created by: Julien Dutheil
// Created on: Mon Feb 21 17:42 2005
//

/*
   Copyright or © or Copr. Bio++ Development Team, (November 17, 2004)

   This software is a computer program whose purpose is to provide classes
   for sequences analysis.

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

// from the STL:
#include <string>

#include "GranthamAAChemicalDistance.h"
#include "../Alphabet/AlphabetTools.h"
#include <Bpp/Numeric/NumTools.h>

using namespace bpp;

short int GranthamAAChemicalDistance::SIGN_NONE = 0;
short int GranthamAAChemicalDistance::SIGN_ARBITRARY = 1;
short int GranthamAAChemicalDistance::SIGN_PC1 = 2;

GranthamAAChemicalDistance::GranthamAAChemicalDistance() :
  distanceMatrix_(20, 20),
  signMatrix_(20, 20),
  alpha_(&AlphabetTools::PROTEIN_ALPHABET),
  sign_(0)
{
  // Load the matrix:
  #include "__GranthamMatrixCode"
}
GranthamAAChemicalDistance::~GranthamAAChemicalDistance() {}

double GranthamAAChemicalDistance::getIndex(int state1, int state2) const
throw (BadIntException)
{
  size_t stateIndex1 = alpha_->getStateIndex(state1);
  size_t stateIndex2 = alpha_->getStateIndex(state2);
  double d = distanceMatrix_(stateIndex1, stateIndex2);
  if (sign_ == SIGN_NONE)
    return NumTools::abs<double>(d);
  if (sign_ == SIGN_PC1)
    return signMatrix_(stateIndex1, stateIndex2) * NumTools::abs<double>(d);
  return d;
}

double GranthamAAChemicalDistance::getIndex(const std::string& state1, const std::string& state2) const
throw (BadCharException)
{
  return getIndex(alpha_->charToInt(state1), alpha_->charToInt(state2));
}

Matrix<double>* GranthamAAChemicalDistance::getIndexMatrix() const
{
  RowMatrix<double>* m = new RowMatrix<double>(distanceMatrix_);
  if (sign_ == SIGN_NONE)
  {
    for (size_t i = 0; i < 20; ++i)
    {
      for (size_t j = 0; j < 20; ++j)
      {
        (*m)(i, j) = NumTools::abs<double>((*m)(i, j));
      }
    }
  }
  else if (sign_ == SIGN_PC1)
  {
    for (size_t i = 0; i < 20; ++i)
    {
      for (size_t j = 0; j < 20; ++j)
      {
        (*m)(i, j) = signMatrix_(i, j) * NumTools::abs<double>((*m)(i, j));
      }
    }
  }
  return m;
}

