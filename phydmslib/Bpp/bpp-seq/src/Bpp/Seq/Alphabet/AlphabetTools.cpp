//
// File: AlphabetTools.cpp
// Created by: Julien Dutheil
// Created on: Fri Oct 10 17:27:39 2003
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

#include "AlphabetTools.h"
#include <Bpp/Text/TextTools.h>

using namespace bpp;

// From the STL:
#include <ctype.h>
#include <iostream>

using namespace std;

/**********************************************************************************************/

const DNA AlphabetTools::DNA_ALPHABET;
const RNA AlphabetTools::RNA_ALPHABET;
const ProteicAlphabet AlphabetTools::PROTEIN_ALPHABET;
const DefaultAlphabet AlphabetTools::DEFAULT_ALPHABET;

/**********************************************************************************************/

int AlphabetTools::getType(char state)
{
  if (state == '-')
    return -1;
  state = static_cast<char>(toupper(static_cast<int>(state))); // toupper works on int
  bool d = DNA_ALPHABET.isCharInAlphabet(TextTools::toString(state));
  bool r = RNA_ALPHABET.isCharInAlphabet(TextTools::toString(state));
  bool p = PROTEIN_ALPHABET.isCharInAlphabet(TextTools::toString(state));

  if (!d && !r && !p)
    return 0;  // Unknown character
  else if (d && !r && !p)
    return 1;  // DNA specific
  else if (!d && r && !p)
    return 2;  // RNA specific
  else if (!d && !r && p)
    return 3;  // Protein specific
  else if (d && r && !p)
    return 4;  // Nucleotide specific
  else if (d && !r && p)
    return 5;  // DNA or Protein specific
  else if (!d && r && p)
    return 6;  // RNA or Protein specific
  else
    return 7;  // Non-specific character
}

/**********************************************************************************************/

bool AlphabetTools::checkAlphabetCodingSize(const Alphabet& alphabet) throw (AlphabetException)
{
  if (alphabet.getNumberOfChars() == 0)
    return true;  // Will this really happen?
  size_t size = alphabet.intToChar(0).size();
  for (int i = 1; i < static_cast<int>(alphabet.getNumberOfTypes()); ++i)
  {
    if (alphabet.intToChar(i).size() != size)
      return false;
  }
  return true;
}

/**********************************************************************************************/

bool AlphabetTools::checkAlphabetCodingSize(const Alphabet* alphabet) throw (AlphabetException)
{
  return checkAlphabetCodingSize(*alphabet);
}

/**********************************************************************************************/

unsigned int AlphabetTools::getAlphabetCodingSize(const Alphabet& alphabet) throw (AlphabetException)
{
  if (!checkAlphabetCodingSize(alphabet))
    throw AlphabetException("Bad alphabet in function Alphabet::getAlphabetCodingSize().");
  return static_cast<unsigned int>(alphabet.intToChar(0).size());
}

/**********************************************************************************************/

unsigned int AlphabetTools::getAlphabetCodingSize(const Alphabet* alphabet) throw (AlphabetException)
{
  return getAlphabetCodingSize(*alphabet);
}
