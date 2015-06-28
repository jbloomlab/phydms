//
// File: WordAlphabet.h
// Authors: Laurent Gueguen
//          Sylvain Gaillard
// Created on: Sun Dec 28 2008
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

#include "WordAlphabet.h"
#include <Bpp/Text/TextTools.h>

using namespace bpp;

// From the STL:
#include <iostream>

using namespace std;

WordAlphabet::WordAlphabet(const vector<const Alphabet*>& vAlpha) :
  AbstractAlphabet(),
  vAbsAlph_(vAlpha)
{
  build_();
}

WordAlphabet::WordAlphabet(const Alphabet* pAlpha, unsigned int num) :
  AbstractAlphabet(),
  vAbsAlph_(0)
{
  for (unsigned int i = 0; i < num; i++)
  {
    vAbsAlph_.push_back(pAlpha);
  }

  build_();
}

void WordAlphabet::build_()
{
  size_t size = 1;

  for (size_t i = 0; i < vAbsAlph_.size(); ++i)
  {
    size *= vAbsAlph_[i]->getSize();
  }

  vector<AlphabetState*> states(size + 2);

  string s = "";
  for (size_t i = 0; i < vAbsAlph_.size(); ++i)
  {
    s += "-";
  }

  states[0] = new AlphabetState(-1, s, "gap");

  for (size_t i = 0; i < size; ++i)
  {
    states[i + 1] = new AlphabetState(static_cast<int>(i), "", "");
  }

  size_t lr = size;
  char c;
  for (size_t na = 0; na < vAbsAlph_.size(); ++na)
  {
    lr /= vAbsAlph_[na]->getSize();
    size_t j = 1;
    int i = 0;
    while (j <= size)
    {
      c = vAbsAlph_[na]->intToChar(i)[0];

      for (size_t k = 0; k < lr; k++)
      {
        states[j]->setLetter(states[j]->getLetter() + c);
        j++;
        // alphabet[j++].letter += c;
      }

      if (++i == static_cast<int>(vAbsAlph_[na]->getSize()))
        i = 0;
    }
  }

  s = "";
  for (size_t i = 0; i < vAbsAlph_.size(); ++i)
  {
    s += "N";
  }

  states[size + 1] = new AlphabetState(static_cast<int>(size), s, "Unresolved");

  //Now register all states once for all:
  for (size_t i = 0; i < states.size(); ++i) {
    registerState(states[i]);
  }
  //jdutheil on 24/07/14: this should not be necessary anymore.
  //remap();
}

/******************************************************************************/

std::string WordAlphabet::getAlphabetType() const
{
  string s = "Word alphabet:";
  for (unsigned int i = 0; i < vAbsAlph_.size(); i++)
  {
    s += " " +  vAbsAlph_[i]->getAlphabetType();
  }

  return s;
}

bool WordAlphabet::hasUniqueAlphabet() const
{
  string s = vAbsAlph_[0]->getAlphabetType();
  for (unsigned int i = 1; i < vAbsAlph_.size(); i++)
  {
    if (vAbsAlph_[i]->getAlphabetType() != s)
      return false;
  }
  return true;
}

bool WordAlphabet::containsUnresolved(const std::string& state) const throw (BadCharException)
{
  size_t s = vAbsAlph_.size();
  if (state.length() != s)
    throw BadCharException(state, "WordAlphabet::containsUnresolved", this);

  for (size_t i = 0; i < vAbsAlph_.size(); i++)
  {
    if (vAbsAlph_[i]->isUnresolved(state.substr(i, 1)))
    {
      return true;
    }
  }
  return false;
}

/******************************************************************************/

bool WordAlphabet::containsGap(const std::string& state) const throw (BadCharException)
{
  size_t s = vAbsAlph_.size();
  if (state.length() != s)
    throw BadCharException(state, "WordAlphabet::containsGap", this);

  for (size_t i = 0; i < vAbsAlph_.size(); i++)
  {
    if (vAbsAlph_[i]->isGap(state.substr(i, 1)))
      return true;
  }

  return false;
}

/******************************************************************************/

std::string WordAlphabet::getName(const std::string& state) const throw (BadCharException)
{
  if (state.size() != vAbsAlph_.size())
    throw BadCharException(state, "WordAlphabet::getName", this);
  if (containsUnresolved(state))
    return getStateAt(getSize() + 1).getName();
  if (containsGap(state))
    return getStateAt(0).getName();
  else
    return AbstractAlphabet::getName(state);
}

/******************************************************************************/

std::vector<int> WordAlphabet::getAlias(int state) const throw (BadIntException)
{
  if (!isIntInAlphabet(state))
    throw BadIntException(state, "WordAlphabet::getAlias(int): Specified base unknown.");
  vector<int> v;
  size_t s = getSize();

  if (static_cast<size_t>(state) == s)
  {
    v.resize(s);
    for (size_t i = 0; i < s; ++i)
    {
      v[i] = static_cast<int>(i);
    }
  }
  else
  {
    v.resize(1); v[0] = state;
  }
  return v;
}

/******************************************************************************/

std::vector<std::string> WordAlphabet::getAlias(const std::string& state) const throw (BadCharException)
{
  string locstate = TextTools::toUpper(state);
  if (!isCharInAlphabet(locstate))
    throw BadCharException(locstate, "WordAlphabet::getAlias(string): Specified base unknown.");
  vector<string> v;

  size_t s = getSize();

  string st = "";
  for (size_t i = 0; i < vAbsAlph_.size(); ++i)
  {
    st += "N";
  }

  if (locstate == st)
  {
    v.resize(s);
    for (size_t i = 0; i < s; ++i)
    {
      v[i] = intToChar(static_cast<int>(i));
    }
  }
  else
  {
    v.resize(1); v[0] = state;
  }
  return v;
}

/******************************************************************************/

int WordAlphabet::getGeneric(const std::vector<int>& states) const throw (BadIntException)
{
  return states[0];
}

/******************************************************************************/

std::string WordAlphabet::getGeneric(const std::vector<std::string>& states) const throw (BadCharException)
{
  return states[0];
}

/******************************************************************************/

int WordAlphabet::getWord(const Sequence& seq, size_t pos) const throw (IndexOutOfBoundsException)
{
  if (seq.size() < pos + vAbsAlph_.size())
    throw IndexOutOfBoundsException("WordAlphabet::getWord", pos, 0, seq.size() - vAbsAlph_.size());

  vector<string> vs;
  for (size_t i = 0; i < vAbsAlph_.size(); i++)
  {
    vs.push_back(vAbsAlph_[i]->intToChar(seq[i + pos]));
  }

  return charToInt(getWord(vs)); // This can't throw a BadCharException!
}


/******************************************************************************/

int WordAlphabet::getWord(const std::vector<int>& vint, size_t pos) const throw (IndexOutOfBoundsException)
{
  if (vint.size() < pos + vAbsAlph_.size())
    throw IndexOutOfBoundsException("WordAlphabet::getWord", pos, 0, vint.size() - vAbsAlph_.size());

  vector<string> vs;
  for (size_t i = 0; i < vAbsAlph_.size(); i++)
  {
    vs.push_back(vAbsAlph_[i]->intToChar(vint[i + pos]));
  }

  return charToInt(getWord(vs)); // This can't throw a BadCharException!
}

/****************************************************************************************/

std::string WordAlphabet::getWord(const std::vector<string>& vpos, size_t pos) const throw (IndexOutOfBoundsException, BadCharException)
{
  if (vpos.size() < pos + vAbsAlph_.size())
    throw IndexOutOfBoundsException("WordAlphabet::getWord", pos, 0, vpos.size() - vAbsAlph_.size());

  string s = "";
  for (size_t i = 0; i < vAbsAlph_.size(); i++)
  {
    s += vpos[pos + i];
  }
  // test
  charToInt(s);
  return s;
}

/****************************************************************************************/

Sequence* WordAlphabet::translate(const Sequence& sequence, size_t pos) const throw (AlphabetMismatchException, Exception)
{
  if ((!hasUniqueAlphabet()) or
      (sequence.getAlphabet()->getAlphabetType() != vAbsAlph_[0]->getAlphabetType()))
    throw AlphabetMismatchException("No matching alphabets", sequence.getAlphabet(), vAbsAlph_[0]);

  vector<int> content;

  size_t s = sequence.size();
  unsigned int l = getLength();
  size_t i = pos;

  while (i + l <= s)
  {
    content.push_back(getWord(sequence, i));
    i += l;
  }

  return new BasicSequence(sequence.getName(), content, this);
}

/****************************************************************************************/

Sequence* WordAlphabet::reverse(const Sequence& sequence) const throw (AlphabetMismatchException, Exception)
{
  if ((!hasUniqueAlphabet()) or
      (sequence.getAlphabet()->getAlphabetType() != getAlphabetType()))
    throw AlphabetMismatchException("No matching alphabets");

  Sequence* pseq = new BasicSequence(sequence.getName(), "", getNAlphabet(0));

  size_t s = sequence.size();
  for (size_t i = 0; i < s; i++)
  {
    pseq->append(getPositions(sequence[i]));
  }

  return pseq;
}

/****************************************************************************************/

