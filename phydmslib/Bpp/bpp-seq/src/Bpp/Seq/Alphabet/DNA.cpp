//
// File: DNA.cpp
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

#include "DNA.h"
#include "AlphabetState.h"
#include <Bpp/Text/TextTools.h>
#include <Bpp/Utils/MapTools.h>

using namespace bpp;

// From STL:
#include <map>

using namespace std;

/******************************************************************************/

DNA::DNA(bool exclamationMarkCountsAsGap)
{
	// Alphabet content definition
	// all unresolved bases use n°14
  registerState(new NucleicAlphabetState(-1, "-",  0, "Gap"));
  registerState(new NucleicAlphabetState( 0, "A",  1, "Adenine"));
  registerState(new NucleicAlphabetState( 1, "C",  2, "Cytosine"));
  registerState(new NucleicAlphabetState( 2, "G",  4, "Guanine"));
  registerState(new NucleicAlphabetState( 3, "T",  8, "Thymine"));
  registerState(new NucleicAlphabetState( 4, "M",  3, "Adenine or Cytosine"));
  registerState(new NucleicAlphabetState( 5, "R",  5, "Purine (Adenine or Guanine)"));
  registerState(new NucleicAlphabetState( 6, "W",  9, "Adenine or Thymine"));
  registerState(new NucleicAlphabetState( 7, "S",  6, "Cytosine or Guanine"));
  registerState(new NucleicAlphabetState( 8, "Y", 10, "Pyrimidine (Cytosine or Thymine)"));
  registerState(new NucleicAlphabetState( 9, "K", 12, "Guanine or Thymine"));
  registerState(new NucleicAlphabetState(10, "V",  7, "Adenine or Cytosine or Guanine"));
  registerState(new NucleicAlphabetState(11, "H", 11, "Adenine or Cytosine or Thymine"));
  registerState(new NucleicAlphabetState(12, "D", 13, "Adenine or Guanine or Thymine"));
  registerState(new NucleicAlphabetState(13, "B", 14, "Cytosine or Guanine or Thymine"));
  registerState(new NucleicAlphabetState(14, "N", 15, "Unresolved base"));
  registerState(new NucleicAlphabetState(14, "X", 15, "Unresolved base"));
  registerState(new NucleicAlphabetState(14, "O", 15, "Unresolved base"));
  registerState(new NucleicAlphabetState(14, "0", 15, "Unresolved base"));
  registerState(new NucleicAlphabetState(14, "?", 15, "Unresolved base"));
  if (exclamationMarkCountsAsGap)
    registerState(new NucleicAlphabetState(-1, "!", 0, "Frameshift"));
  else
    registerState(new NucleicAlphabetState(14, "!", 15, "Unresolved base"));
}

/******************************************************************************/

std::vector<int> DNA::getAlias(int state) const throw (BadIntException) 
{
	if (!isIntInAlphabet(state))
    throw BadIntException(state, "DNA::getAlias(int): Specified base unknown.");
	vector<int> v;
  const NucleicAlphabetState& st = getState(state);
  if (state == -1)
    v.push_back(-1);
  if (st.getBinaryCode() & 1)
    v.push_back(0);
  if (st.getBinaryCode() & 2)
    v.push_back(1);
  if (st.getBinaryCode() & 4)
    v.push_back(2);
  if (st.getBinaryCode() & 8)
    v.push_back(3);
	return v;
}


/******************************************************************************/

std::vector<std::string> DNA::getAlias(const std::string& state) const throw (BadCharException) 
{
  string locstate = TextTools::toUpper(state);
	if(!isCharInAlphabet(locstate)) throw BadCharException(locstate, "DNA::getAlias(int): Specified base unknown.");
  vector<int> vi = this->getAlias(this->charToInt(state));
	vector<string> v;
  for (unsigned int i = 0 ; i < vi.size() ; i++)
    v.push_back(this->intToChar(vi[i]));
	return v;
}

/******************************************************************************/

int DNA::getGeneric(const std::vector<int>& states) const throw (BadIntException)
{
  int v = 0;
  for (size_t i = 0 ; i < states.size() ; ++i) {
    if (!isIntInAlphabet(states[i])) throw BadIntException(states[i], "DNA::getGeneric(const vector<int>& states): Specified base unknown.");
    v |= getState(states[i]).getBinaryCode();
  }
  return getStateByBinCode(v).getNum();
}

/******************************************************************************/

std::string DNA::getGeneric(const std::vector<std::string>& states) const throw (BadCharException)
{
  vector<int> vi;
  for (unsigned int i = 0 ; i < states.size() ; ++i) {
    if (!isCharInAlphabet(states[i])) throw BadCharException(states[i], "DNA::getGeneric(const vector<string>& states): Specified base unknown.");
    vi.push_back(this->charToInt(states[i]));
  }
  return intToChar(getGeneric(vi));
}

/******************************************************************************/

