//
// File: LexicalAlphabet.h
// Created by: Laurent Gueguen
// Created on: mercredi 13 juillet 2016, à 11h 26
//

/*
  Copyright or © or Copr. Bio++ Development Team, (November 17, 2004)

  This software is a computer program whose purpose is to provide classes
  for sequences analysis.

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

#ifndef _LEXICAL_ALPHABET_H_
#define _LEXICAL_ALPHABET_H_

#include "AbstractAlphabet.h"

// From the STL:
#include <string>
#include <vector>

namespace bpp
{
/**
 * @brief Alphabet made of given words.
 *
 */
  class LexicalAlphabet :
    public AbstractAlphabet
  {
  public:
    // Constructor and destructor.
    /**
     * @brief Builds a new word alphabet from a vector of given words.
     *
     * Words unicity is checked. 
     *
     * @param vocab The vector of words to be used.
     */
  
    LexicalAlphabet(const std::vector<std::string>& vocab);
  
    LexicalAlphabet(const LexicalAlphabet& bia) : AbstractAlphabet(bia) {}
  
    LexicalAlphabet& operator=(const LexicalAlphabet& bia)
    {
      AbstractAlphabet::operator=(bia);
      return *this;
    }

    LexicalAlphabet* clone() const
    {
      return new LexicalAlphabet(*this);
    }

    virtual ~LexicalAlphabet() {}

  public:
    /**
     * @name Methods redefined from Alphabet
     *
     * @{
     */
    /**
     * @brief Get the complete name of a state given its string description.
     *
     * In case of undefined characters (i.e. N and X for nucleic alphabets),
     * this method will return the name of the undefined word.
     *
     * @param state The string description of the given state.
     * @return The name of the state.
     * @throw BadCharException When state is not a valid char description.
     */
    
    unsigned int getSize() const
    {
      return getNumberOfChars() - 2;
    }

    /** @} */


    /**
     * @brief Returns the number of resolved states + one for unresolved
     *
     */
    unsigned int getNumberOfTypes() const
    {
      return getNumberOfChars() - 1;
    }

    std::string getAlphabetType() const;

    int getUnknownCharacterCode() const
    {
      return static_cast<int>(getSize());
    }

    bool isUnresolved(int state) const { return state == getUnknownCharacterCode(); }
    bool isUnresolved(const std::string& state) const { return charToInt(state) == getUnknownCharacterCode(); }

  };
} // end of namespace bpp.

#endif  // _LEXICAL_ALPHABET_H_

