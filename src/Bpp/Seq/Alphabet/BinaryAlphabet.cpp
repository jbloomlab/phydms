//
// File: BinaryAlphabet.cpp
// Authors: Laurent Gueguen
// Created on: 2009
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


#include "BinaryAlphabet.h"
#include "AlphabetState.h"

// From Utils:
#include <Bpp/Text/TextTools.h>

using namespace bpp;

BinaryAlphabet::BinaryAlphabet()
{
  // Alphabet content definition
  registerState(new AlphabetState(-1, "-", "Gap"));
  registerState(new AlphabetState(2, "?", "Unresolved state"));

  for (int i = 0; i < 2; i++)
  {
    registerState(new AlphabetState(i, TextTools::toString(i), ""));
  }
}

/******************************************************************************/

std::vector<int> BinaryAlphabet::getAlias(int state) const throw (BadIntException) 
{
  if (!isIntInAlphabet(state)) throw BadIntException(state, "BinaryAlphabet::getAlias(int): Specified base unknown.");
  std::vector<int> v;
  switch(state)
  {
  case -1:
    v.push_back(-1);
    break;
  case 0:
    v.push_back(0);
    break;
  case 1:
    v.push_back(1);
    break;
  case 2:
    v.push_back(0);
    v.push_back(1);
    break;
  }
  
  return v;
}

/******************************************************************************/

std::vector<std::string> BinaryAlphabet::getAlias(const std::string& state) const throw (BadCharException) 
{
  if (!isCharInAlphabet(state)) throw BadCharException(state, "BinaryAlphabet::getAlias(char): Specified base unknown.");
  
  std::vector<std::string> v(1);
  if (state=="?")
  {
    v.push_back("0");
    v.push_back("1");
  }
  else
    v.push_back(state);
  
  return v;
}

