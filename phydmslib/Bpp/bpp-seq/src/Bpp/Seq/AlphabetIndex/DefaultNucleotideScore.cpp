//
// File: DefaultNucleotideScore.cpp
// Created by: Julien Dutheil
// Created on: Fri Jan 19 10:30 2007
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

#include "DefaultNucleotideScore.h"

// from the STL:
#include <string>

using namespace std;
using namespace bpp;

DefaultNucleotideScore::DefaultNucleotideScore(const NucleicAlphabet* alphabet) :
  distanceMatrix_(4, 4),
  alpha_(alphabet)
{
  // Load the matrix:
  distanceMatrix_(0, 0) = 10;
  distanceMatrix_(0, 1) = -3;
  distanceMatrix_(0, 2) = -1;
  distanceMatrix_(0, 3) = -4;

  distanceMatrix_(1, 0) = -3;
  distanceMatrix_(1, 1) = 9;
  distanceMatrix_(1, 2) = -5;
  distanceMatrix_(1, 3) = 0;

  distanceMatrix_(2, 0) = -1;
  distanceMatrix_(2, 1) = -5;
  distanceMatrix_(2, 2) = 7;
  distanceMatrix_(2, 3) = -3;

  distanceMatrix_(3, 0) = -4;
  distanceMatrix_(3, 1) = 0;
  distanceMatrix_(3, 2) = -3;
  distanceMatrix_(3, 3) = 8;
}

double DefaultNucleotideScore::getIndex(int state1, int state2) const
throw (BadIntException)
{
  if (alpha_->isGap(state1) || !alpha_->isIntInAlphabet(state1))
    throw BadIntException(state1, "DefaultNucleotideScore::getIndex(). Invalid state1.", alpha_);
  if (alpha_->isGap(state2) || !alpha_->isIntInAlphabet(state2))
    throw BadIntException(state2, "DefaultNucleotideScore::getIndex(). Invalid state1.", alpha_);
  if (!alpha_->isUnresolved(state1) && !alpha_->isUnresolved(state2))
    return distanceMatrix_(
        static_cast<size_t>(state1),
        static_cast<size_t>(state2));
  vector<int> states1 = alpha_->getAlias(state1);
  vector<int> states2 = alpha_->getAlias(state2);
  double score = -5;
  double tmp_score;
  for (size_t i = 0; i < states1.size(); i++)
  {
    for (size_t j = 0; j < states2.size(); j++)
    {
      tmp_score = getIndex(states1[i], states2[j]);
      if (tmp_score > score)
        score = tmp_score;
    }
  }
  return score / static_cast<double>(states1.size() + states2.size() - 1);
}

double DefaultNucleotideScore::getIndex(const std::string& state1, const std::string& state2) const
throw (BadCharException)
{
  return distanceMatrix_(
      static_cast<size_t>(alpha_->charToInt(state1)),
      static_cast<size_t>(alpha_->charToInt(state2)));
}

LinearMatrix<double>* DefaultNucleotideScore::getIndexMatrix() const
{
  return new LinearMatrix<double>(distanceMatrix_);
}

