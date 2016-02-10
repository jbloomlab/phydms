//
// File: AbstractAlphabet.cpp
// Authors: Guillaume Deuchst
//          Julien Dutheil
//          Sylvain Gaillard
// Created on: Tue Jul 22 2003
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

#include "AbstractAlphabet.h"
#include <Bpp/Text/TextTools.h>
#include <Bpp/Utils/MapTools.h>

using namespace bpp;

// From the STL:
#include <ctype.h>
#include <map>
#include <iostream>

using namespace std;

/******************************************************************************/

void AbstractAlphabet::updateMaps_(size_t pos, const AlphabetState& st) {
  if (letters_.find(st.getLetter()) == letters_.end())
    letters_[st.getLetter()] = pos;
  else
    throw Exception("AbstractAlphabet::updateMaps_. A state with the same character code already exists! " + st.getLetter() + ".");
  if (nums_.find(st.getNum()) == nums_.end())
    nums_[st.getNum()] = pos;
  else
    nums_[st.getNum()] = min(pos, nums_[st.getNum()]);
}

/******************************************************************************/

void AbstractAlphabet::registerState(AlphabetState* st) throw (Exception) {
  // Add the state to the vector
  alphabet_.push_back(st);
  // Update the maps
  updateMaps_(alphabet_.size() - 1, *st);
}

/******************************************************************************/

void AbstractAlphabet::setState(size_t pos, AlphabetState* st)
  throw (Exception, IndexOutOfBoundsException) {
    if (pos > alphabet_.size())
      throw IndexOutOfBoundsException("AbstractAlphabet::setState: incorect position", pos, 0, alphabet_.size());
    // Delete the state if not empty
    if (alphabet_[pos] != 0)
      delete alphabet_[pos];
    // Put the state in the vector
    alphabet_[pos] = st;
    // Update the maps
    updateMaps_(pos, *st);
  }

/******************************************************************************/

const AlphabetState& AbstractAlphabet::getState(const std::string& letter) const throw (BadCharException) {
  map<string, size_t>::const_iterator it = letters_.find(letter);
  if (it == letters_.end())
    throw BadCharException(letter, "AbstractAlphabet::getState(string): Specified base unknown", this);
  return * (alphabet_[it->second]);
}

/******************************************************************************/

size_t AbstractAlphabet::getStateIndex(const std::string& letter) const throw (BadCharException) {
  map<string, size_t>::const_iterator it = letters_.find(letter);
  if (it == letters_.end())
    throw BadCharException(letter, "AbstractAlphabet::getStateIndex(string): Specified base unknown", this);
  return it->second;
}

/******************************************************************************/

const AlphabetState& AbstractAlphabet::getState(int num) const throw (BadIntException) {
  map<int, size_t>::const_iterator it = nums_.find(num);
  if (it == nums_.end())
    throw BadIntException(num, "AbstractAlphabet::getState(int): Specified base unknown", this);
  return *(alphabet_[it->second]);
}

/******************************************************************************/

size_t AbstractAlphabet::getStateIndex(int num) const throw (BadIntException) {
  map<int, size_t>::const_iterator it = nums_.find(num);
  if (it == nums_.end())
    throw BadIntException(num, "AbstractAlphabet::getStateIndex(int): Specified base unknown", this);
  return it->second;
}

/******************************************************************************/

AlphabetState& AbstractAlphabet::getState(const std::string& letter) throw (BadCharException) {
  map<string, size_t>::iterator it = letters_.find(letter);
  if (it == letters_.end())
    throw BadCharException(letter, "AbstractAlphabet::getState(string): Specified base unknown", this);
  return * (alphabet_[it->second]);
}

/******************************************************************************/

AlphabetState& AbstractAlphabet::getState(int num) throw (BadIntException) {
  map<int, size_t>::iterator it = nums_.find(num);
  if (it == nums_.end())
    throw BadIntException(num, "AbstractAlphabet::getState(int): Specified base unknown", this);
  return * (alphabet_[it->second]);
}

/******************************************************************************/

AlphabetState& AbstractAlphabet::getStateAt(size_t pos) throw (IndexOutOfBoundsException) {
  if (pos > alphabet_.size())
    throw IndexOutOfBoundsException("AbstractAlphabet::getStateAt: incorect position", pos, 0, alphabet_.size());
  return * (alphabet_[pos]);
}

/******************************************************************************/

const AlphabetState& AbstractAlphabet::getStateAt(size_t pos) const throw (IndexOutOfBoundsException) {
  if (pos > alphabet_.size())
    throw IndexOutOfBoundsException("AbstractAlphabet::getStateAt: incorect position", pos, 0, alphabet_.size());
  return * (alphabet_[pos]);
}

/******************************************************************************/

std::string AbstractAlphabet::getName(const std::string& state) const throw (BadCharException)
{
  return (getState(state)).getName();
}

/******************************************************************************/

std::string AbstractAlphabet::getName(int state) const throw (BadIntException)
{
  return (getState(state)).getName();
}

/******************************************************************************/

int AbstractAlphabet::charToInt(const std::string& state) const throw (BadCharException)
{
  return getState(state).getNum();
}

/******************************************************************************/

std::string AbstractAlphabet::intToChar(int state) const throw (BadIntException)
{
  return (getState(state)).getLetter();
}

/******************************************************************************/

bool AbstractAlphabet::isIntInAlphabet(int state) const
{
  map<int, size_t>::const_iterator it = nums_.find(state);
  if (it != nums_.end())
    return true;
  return false;
}

/******************************************************************************/

bool AbstractAlphabet::isCharInAlphabet(const std::string& state) const
{
  map<string, size_t>::const_iterator it = letters_.find(state);
  if (it != letters_.end())
    return true;
  return false;
}	

/******************************************************************************/

std::vector<int> AbstractAlphabet::getAlias(int state) const throw (BadIntException) 
{
  if (!isIntInAlphabet(state)) throw BadIntException(state, "AbstractAlphabet::getAlias(int): Specified base unknown.");
  vector<int> v(1);
  v[0] = state;
  return v;
}

/******************************************************************************/

std::vector<std::string> AbstractAlphabet::getAlias(const std::string& state) const throw (BadCharException) 
{
  if (!isCharInAlphabet(state)) throw BadCharException(state, "AbstractAlphabet::getAlias(char): Specified base unknown.");
  vector<string> v(1);
  v[0] = state;
  return v;
}

/******************************************************************************/

int AbstractAlphabet::getGeneric(const std::vector<int>& states) const throw (BadIntException) {
  map<int, int> m;
  for (unsigned int i = 0 ; i < states.size() ; ++i) {
    vector<int> tmp_s = this->getAlias(states[i]); // get the states for generic characters
    for (unsigned int j = 0 ; j < tmp_s.size() ; ++j) {
      m[tmp_s[j]] ++; // add each state to the list
    }
  }
  vector<int> ve = MapTools::getKeys(m);

  string key;
  for (unsigned int i = 0 ; i < ve.size() ; ++i) {
    if (!isIntInAlphabet(ve[i])) throw BadIntException(ve[i], "AbstractAlphabet::getGeneric(const vector<int>): Specified base unknown.");
    key += "_" + TextTools::toString(ve[i]);
  }
  int v;
  if (ve.size() == 1) {
    v = ve[0];
  } else {
    v = this->getUnknownCharacterCode();
  }
  return v;
}

/******************************************************************************/

std::string AbstractAlphabet::getGeneric(const std::vector<std::string>& states) const throw (AlphabetException) {
  map <string, int> m;
  for (unsigned int i = 0 ; i < states.size() ; ++i) {
    vector<string> tmp_s = this->getAlias(states[i]); // get the states for generic characters
    for (unsigned int j = 0 ; j < tmp_s.size() ; ++j) {
      m[tmp_s[j]] ++; // add each state to the list
    }
  }
  vector<string> ve = MapTools::getKeys(m);

  string key;
  for (unsigned int i = 0 ; i < ve.size() ; ++i) {
    if (!isCharInAlphabet(ve[i])) throw BadCharException(ve[i], "AbstractAlphabet::getAlias(const vector<string>): Specified base unknown.");
    key += TextTools::toString(ve[i]);
  }
  string v;
  if (ve.size() == 1) {
    v = ve[0];
  } else {
    throw CharStateNotSupportedException("AbstractAlphabet::getAlias(const vector<string>): No generic char state.");
  }
  return v;
}

/******************************************************************************/

const std::vector<int>& AbstractAlphabet::getSupportedInts() const
{
  if(intList_.size() != alphabet_.size())
  {
    intList_.resize(alphabet_.size());
    charList_.resize(alphabet_.size());
    for (size_t i = 0; i < alphabet_.size(); ++i)
    {
      intList_[i]  = alphabet_[i]->getNum();
      charList_[i] = alphabet_[i]->getLetter();
    }
  }
  return intList_;
}

/******************************************************************************/

const std::vector<std::string>& AbstractAlphabet::getSupportedChars() const
{
  if(charList_.size() != alphabet_.size())
  {
    intList_.resize(alphabet_.size());
    charList_.resize(alphabet_.size());
    for (size_t i = 0; i < alphabet_.size(); ++i)
    {
      intList_[i]  = alphabet_[i]->getNum();
      charList_[i] = alphabet_[i]->getLetter();
    }
  }
  return charList_;
}

/******************************************************************************/

const std::vector<std::string> & AbstractAlphabet::getResolvedChars() const
{
  charList_.clear();
  for(size_t i = 0; i < alphabet_.size(); ++i)
    // well, non-gap chars also
    if(!isGap(alphabet_[i]->getLetter()) and !isUnresolved(alphabet_[i]->getLetter()))
      charList_.push_back(alphabet_[i]->getLetter());

  return charList_;
}
