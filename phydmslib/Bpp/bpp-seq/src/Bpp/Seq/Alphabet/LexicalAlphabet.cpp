//
// File: LexicalAlphabet.h
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

#include "LexicalAlphabet.h"
#include <Bpp/Text/TextTools.h>

using namespace bpp;

// From the STL:
#include <iostream>

using namespace std;

LexicalAlphabet::LexicalAlphabet(const vector<std::string>& vocab) :
  AbstractAlphabet()
{
  if (vocab.size()==0)
    throw Exception("LexicalAlphabet::LexicalAlphabet: not constructible from empty vocabulary.");
  
  size_t t=vocab[0].size();
  
  string s="";  
  for (size_t i=0; i<t; i++)
    s+="-";
  
  registerState(new AlphabetState(-1, s, "gap"));


  for (size_t i = 0; i < vocab.size(); ++i)
  {
    if (t!=vocab[i].size())
      throw Exception("LexicalAlphabet: several lengths in vocabulary.");

    try
    {
      string s2=getName(vocab[i]);
      throw Exception("LexicalAlphabet : " + vocab[i] + " defined twice.");
    }
    catch (BadCharException& e)
    {
    }
    
    registerState(new AlphabetState(static_cast<int>(i), vocab[i], vocab[i]));
  }
  
  s="";  
  for (size_t i=0; i<t; i++)
    s+="?";
  
  registerState(new AlphabetState(static_cast<int>(vocab.size()), s, "Unresolved word"));
}


/******************************************************************************/

std::string LexicalAlphabet::getAlphabetType() const
{
  string s = "Lexicon(";

  for (size_t i=1; i<getNumberOfStates()-1; i++)
  {
    if (i!=1)
      s+=",";
    
    s += getStateAt(i).getLetter();
  }

  s+=")";
  
  return s;
}


