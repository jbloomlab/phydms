//
// File: NumericAlphabet.cpp
// Created by: Laurent Gueguen
// Created on: March 2010
//

/*
   Copyright or © or Copr. Bio++ Development Team, (November 17, 2004)

   This software is a computer program whose purpose is to provide
   classes for sequences analysis.

   This software is governed by the CeCILL license under French law and
   abiding by the rules of distribution of free software. You can use,
   modify and/ or redistribute the software under the terms of the CeCILL
   license as circulated by CEA, CNRS and INRIA at the following URL
   "http://www.cecill.info".

   As a counterpart to the access to the source code and rights to copy,
   modify and redistribute granted by the license, users are provided
   only with a limited warranty and the software's author, the holder of
   the economic rights, and the successive licensors have only limited
   liability.

   In this respect, the user's attention is drawn to the risks associated
   with loading, using, modifying and/or developing or reproducing the
   software by the user in light of its specific status of free software,
   that may mean that it is complicated to manipulate, and that also
   therefore means that it is reserved for developers and experienced
   professionals having in-depth computer knowledge. Users are therefore
   encouraged to load and test the software's suitability as regards
   their requirements in conditions enabling the security of their
   systems and/or data to be ensured and, more generally, to use and
   operate it in the same conditions as regards security.

   The fact that you are presently reading this means that you have had
   knowledge of the CeCILL license and that you accept its terms.
 */

#include "NumericAlphabet.h" // class's header file

#include <Bpp/Text/TextTools.h>

#include <iostream>
#include <typeinfo>

using namespace std;
using namespace bpp;

/****************************************************************************************/

NumericAlphabet::NumericAlphabet(const UniformDiscreteDistribution& pd) :
  AbstractAlphabet(), pdd_(pd.clone()), values_()
{
  
  // Alphabet size definition
  size_t size = pdd_->getNumberOfCategories();

  Vdouble vd = pdd_->getCategories();
  
  for (size_t i = 0; i < size; ++i){
    registerState(new AlphabetNumericState(static_cast<int>(i), vd[i], TextTools::toString(vd[i]), TextTools::toString(vd[i])));
  }
}

NumericAlphabet::NumericAlphabet(const NumericAlphabet& na) :
  AbstractAlphabet(na),
  pdd_(na.pdd_->clone()),
  values_(na.values_)
{}

NumericAlphabet& NumericAlphabet::operator=(const NumericAlphabet& na) 
{
  AbstractAlphabet::operator=(na);
  pdd_ = na.pdd_->clone();
  return *this;
}

/****************************************************************************************/

void NumericAlphabet::setState(size_t stateIndex, AlphabetState* st) throw (Exception)
{
  try {
    AbstractAlphabet::setState(stateIndex, st);
    double x = dynamic_cast<AlphabetNumericState*>(st)->getValue();
    if (values_.find(x) == values_.end())
    values_[x] = stateIndex;  
  } catch(std::bad_cast&) {
    throw Exception("NumericAlphabet::setState. Incorrect alphabet type.");
  }
}

void NumericAlphabet::registerState(AlphabetState* st) throw (Exception)
{
  try {
    AbstractAlphabet::registerState(st);
    double x = dynamic_cast<AlphabetNumericState*>(st)->getValue();
    if (values_.find(x) == values_.end())
      values_[x] = getSize();
  } catch(std::bad_cast&) {
    throw Exception("NumericAlphabet::registerState. Incorrect alphabet type.");
  }
}


vector<int> NumericAlphabet::getAlias(int state) const throw (BadIntException)
{
  vector<int> v;
  v.push_back(state);
  return v;
}

vector<string> NumericAlphabet::getAlias(const string& state) const throw (BadCharException)
{
  vector<string> v;
  v.push_back(state);
  return v;
}

/****************************************************************************************/

bool NumericAlphabet::isGap(int state) const
{
  return false;
}

bool NumericAlphabet::containsGap(const string& state) const throw (BadCharException)
{
  return false;
}

bool NumericAlphabet::isUnresolved(const string& state) const
{
  return false;
}

bool NumericAlphabet::isUnresolved(int state) const
{
  return false;
}

unsigned int NumericAlphabet::getSize() const
{
  return static_cast<unsigned int>(values_.size());
}

unsigned int NumericAlphabet::getNumberOfTypes() const
{
  return static_cast<unsigned int>(values_.size());
}  

void NumericAlphabet::remap()
{
  AbstractAlphabet::remap();
  values_.clear();
  for (size_t pos = 0 ; pos < getSize() ; pos++) {
    double x = getStateAt(pos).getValue();
    if (values_.find(x) == values_.end())
      values_[x] = pos;
  }
}

/************************************************************/

double NumericAlphabet::getDelta() const
{
  return (pdd_->getUpperBound() - pdd_->getLowerBound()) / static_cast<double>(pdd_->getNumberOfCategories());
}

double NumericAlphabet::intToValue(int state) const throw (BadIntException)
{
  return static_cast<const AlphabetNumericState& >(getState(state)).getValue();
}

size_t NumericAlphabet::getValueIndex(double value) const
{
  map<double, size_t>::const_iterator it = values_.find(pdd_->getValueCategory(value));
  return it->second;
}

AlphabetNumericState& NumericAlphabet::getStateAt(size_t stateIndex) throw (IndexOutOfBoundsException)
{
  return static_cast<AlphabetNumericState&>(AbstractAlphabet::getStateAt(stateIndex));
}

const AlphabetNumericState& NumericAlphabet::getStateAt(size_t stateIndex) const throw (IndexOutOfBoundsException)
{
  return static_cast<const AlphabetNumericState&>(AbstractAlphabet::getStateAt(stateIndex));
}
