//
// File: Sequence.cpp
// Created by: Guillaume Deuchst
//             Julien Dutheil
// Created on: Tue Aug 21 2003
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

#include "Sequence.h" // class's header file

#include "Alphabet/AlphabetTools.h"
#include "StringSequenceTools.h"
#include <Bpp/Text/TextTools.h>

using namespace bpp;

 // From the STL:
#include <iostream>

using namespace std;

/* Constructors: **************************************************************/

BasicSequence::BasicSequence(const Alphabet* alpha):
  AbstractCoreSequence(),
  BasicSymbolList(alpha)
{}

BasicSequence::BasicSequence(
  const std::string& name,
  const std::string& sequence,
  const Alphabet* alpha
) throw (BadCharException) :
  AbstractCoreSequence(name),
  BasicSymbolList(alpha)
{
  if (sequence!="")
    setContent(sequence);
}

BasicSequence::BasicSequence(
  const std::string& name,
  const std::string& sequence,
  const Comments& comments,
  const Alphabet* alpha
) throw (BadCharException) :
  AbstractCoreSequence(name, comments),
  BasicSymbolList(alpha)
{
  if (sequence != "")
    setContent(sequence);
}

BasicSequence::BasicSequence(
  const std::string& name,
  const std::vector<std::string>& sequence,
  const Alphabet* alpha
) throw (BadCharException) :
  AbstractCoreSequence(name),
  BasicSymbolList(sequence, alpha)
{}

BasicSequence::BasicSequence(
  const std::string& name,
  const std::vector<std::string>& sequence,
  const Comments& comments,
  const Alphabet* alpha
) throw (BadCharException) :
  AbstractCoreSequence(name, comments),
  BasicSymbolList(sequence, alpha)
{}

BasicSequence::BasicSequence(
  const std::string& name,
  const std::vector<int>& sequence,
  const Alphabet* alpha
) throw (BadIntException) :
  AbstractCoreSequence(name),
  BasicSymbolList(sequence, alpha)
{}

BasicSequence::BasicSequence(
  const std::string& name,
  const std::vector<int>& sequence,
  const Comments& comments,
  const Alphabet* alpha
) throw (BadIntException) :
  AbstractCoreSequence(name, comments),
  BasicSymbolList(sequence, alpha)
{}

/* Copy constructors: *********************************************************/

BasicSequence::BasicSequence(const Sequence& s) :
  AbstractCoreSequence(s),
  BasicSymbolList(s)
{
}

BasicSequence::BasicSequence(const BasicSequence& s) :
  AbstractCoreSequence(s),
  BasicSymbolList(s)
{
}

/* Assignation operator: ******************************************************/

BasicSequence& BasicSequence::operator=(const Sequence& s)
{
  AbstractCoreSequence::operator=(s);
  BasicSymbolList::operator=(s);
  return *this;
}

BasicSequence& BasicSequence::operator=(const BasicSequence& s)
{
  AbstractCoreSequence::operator=(s);
  BasicSymbolList::operator=(s);
  return *this;
}

/******************************************************************************/

void BasicSequence::setContent(const std::string& sequence) throw (BadCharException)
{
  // Remove blanks in sequence
  content_ = StringSequenceTools::codeSequence(TextTools::removeWhiteSpaces(sequence), getAlphabet());
  //Warning, an exception may be thrown here!
}

/******************************************************************************/

void BasicSequence::setToSizeR(size_t newSize)
{
  // Size verification
  size_t seqSize = content_.size();
  if (newSize == seqSize) return;

  if (newSize < seqSize)
  {
    content_.resize(newSize);
    return;
  }

  // Add gaps up to specified size
  int gap = getAlphabet()->getGapCharacterCode();
  while (content_.size() < newSize) content_.push_back(gap);
}

/******************************************************************************/

void BasicSequence::setToSizeL(size_t newSize)
{
  // Size verification
  size_t seqSize = content_.size();
  if (newSize == seqSize) return;

  if (newSize < seqSize)
  {
    //We must truncate sequence from the left.
    //This is a very unefficient method!
    content_.erase(content_.begin(), content_.begin() + static_cast<ptrdiff_t>(seqSize - newSize));
    return;
  }

  // Add gaps up to specified size
  int gap = getAlphabet()->getGapCharacterCode();
  content_.insert(content_.begin(), newSize - seqSize, gap);
}

/******************************************************************************/

void BasicSequence::append(const Sequence& seq) throw (AlphabetMismatchException)
{
  if (seq.getAlphabet()->getAlphabetType() != getAlphabet()->getAlphabetType())
    throw AlphabetMismatchException("BasicSequence::append");
  // Check list for incorrect characters
  for (size_t i = 0; i < seq.size(); i++)
    content_.push_back(seq[i]);
}

void BasicSequence::append(const std::vector<int>& content) throw (BadIntException)
{
  // Check list for incorrect characters
  for (size_t i = 0; i < content.size(); i++)
    if(!getAlphabet()->isIntInAlphabet(content[i]))
      throw BadIntException(content[i], "BasicSequence::append", getAlphabet());
  //BasicSequence is valid:
  for (size_t i = 0; i < content.size(); i++)
    content_.push_back(content[i]);
}

void BasicSequence::append(const std::vector<std::string>& content) throw (BadCharException)
{
  // Check list for incorrect characters
  for (size_t i = 0; i < content.size(); i++)
    if(!getAlphabet()->isCharInAlphabet(content[i]))
      throw BadCharException(content[i], "BasicSequence::append", getAlphabet());
  
  //BasicSequence is valid:
  for (size_t i = 0; i < content.size(); i++)
    content_.push_back(getAlphabet()->charToInt(content[i]));
}

void BasicSequence::append(const std::string& content) throw (BadCharException)
{
  append(StringSequenceTools::codeSequence(content, getAlphabet()));
}

/******************************************************************************/

