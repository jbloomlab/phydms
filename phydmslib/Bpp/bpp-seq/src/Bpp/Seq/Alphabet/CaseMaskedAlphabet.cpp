//
// File: CaseMaskedAlphabet.cpp
// Created by: Julien Dutheil
// Created on: Sun Sep 05 2010
//

/*
Copyright or © or Copr. CNRS, (November 17, 2004)

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

#include "CaseMaskedAlphabet.h"
#include <Bpp/Text/TextTools.h>

using namespace bpp;

//From the STL:
#include <vector>
#include <iostream>

using namespace std;

CaseMaskedAlphabet::CaseMaskedAlphabet(const LetterAlphabet* nocaseAlphabet) :
  LetterAlphabet(true),
  nocaseAlphabet_(nocaseAlphabet)
{
  vector<string> chars = nocaseAlphabet_->getSupportedChars();
  for (size_t i = 0; i < chars.size(); ++i) {
    AlphabetState* state = nocaseAlphabet_->getState(chars[i]).clone();
    registerState(state);
    char c = *chars[i].c_str();
    if (isalpha(c)) {
      if (isupper(c)) {
        registerState(new AlphabetState(state->getNum() + 100, TextTools::toLower(state->getLetter()), string("Masked ") + state->getName()));
      }
    }
  }
}

int CaseMaskedAlphabet::getMaskedEquivalentState(int state) const
  throw (BadIntException)
{
  if (!isIntInAlphabet(state))
    throw BadIntException(state, "CaseMaskedAlphabet::getMaskedEquivalentState. Unsupported state code.");
  if (state >= 100) return state;
  else {
    state += 100;
    if (!isIntInAlphabet(state))
      throw BadIntException(state, "CaseMaskedAlphabet::getMaskedEquivalentState. State has masked equivalent.");
    return state;
  }
}

const string CaseMaskedAlphabet::getMaskedEquivalentState(const string& state) const
  throw (BadCharException, BadIntException)
{
  if (!isCharInAlphabet(state))
    throw BadCharException(state, "CaseMaskedAlphabet::getMaskedEquivalentState. Unsupported state code.");
  int code = charToInt(state);
  if (code >= 100) return state;
  else {
    code += 100;
    if (!isIntInAlphabet(code))
      throw BadIntException(code, "CaseMaskedAlphabet::getMaskedEquivalentState. State has masked equivalent.");
    return intToChar(code);
  }
}


