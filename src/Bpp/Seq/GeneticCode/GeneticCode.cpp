//
// File: GeneticCode.cpp
// Created by: Julien Dutheil
// Created on: Mon Oct 13 15:37:25 2003
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

#include "GeneticCode.h"
#include "../SequenceTools.h"
#include "../Alphabet/AlphabetTools.h"

using namespace bpp;
using namespace std;

/**********************************************************************************************/

StopCodonException::StopCodonException(const std::string& text, const std::string& codon) :
  Exception("StopCodonException: " + text + "(" + codon + ")"),
  codon_(codon) {}

/**********************************************************************************************/

int GeneticCode::translate(int state) const throw (BadIntException, Exception)
{
  if (isStop(state))
    throw StopCodonException("GeneticCode::translate().", codonAlphabet_.intToChar(state)); 
    
  map<int, int>::const_iterator it = tlnTable_.find(state);
  if (it == tlnTable_.end())
    throw BadIntException(state, "GeneticCode::translate().");
  
  return it->second;
}		

/**********************************************************************************************/

std::string GeneticCode::translate(const std::string& state) const throw (BadCharException, Exception)
{
  int x = codonAlphabet_.charToInt(state);
  return proteicAlphabet_.intToChar(translate(x));
}

/**********************************************************************************************/

vector<int> GeneticCode::getSynonymous(int aminoacid) const throw (BadIntException)
{
  // test:
  proteicAlphabet_.intToChar(aminoacid);

  vector<int> synonymes;
  for (int i = 0; i < static_cast<int>(codonAlphabet_.getSize()); ++i)
  {
    try
    {
      if (translate(i) == aminoacid)
        synonymes.push_back(i);
    }
    catch (StopCodonException)
    { }
  }
  return synonymes;
}

/**********************************************************************************************/

std::vector<std::string> GeneticCode::getSynonymous(const std::string& aminoacid) const throw (BadCharException)
{
  // test:
  int aa = proteicAlphabet_.charToInt(aminoacid);

  vector<string> synonymes;
  for (int i = 0; i < static_cast<int>(codonAlphabet_.getSize()); ++i)
  {
    try
    {
      if (translate(i) == aa)
        synonymes.push_back(codonAlphabet_.intToChar(i));
    }
    catch (StopCodonException)
    { }
  }
  return synonymes;
}

/**********************************************************************************************/

bool GeneticCode::isFourFoldDegenerated(int val) const
{
  if (isStop(val))
    return false;

  vector<int> codon = codonAlphabet_.getPositions(val);
  int acid = translate(val);

  // test all the substitution on third codon position
  for (int an = 0; an < 4; an++)
  {
    if (an == codon[2])
      continue;
    vector<int> mutcodon = codon;
    mutcodon[2] = an;
    int intcodon = codonAlphabet_.getCodon(mutcodon[0], mutcodon[1], mutcodon[2]);
    if (isStop(intcodon))
      return false;
    int altacid = translate(intcodon);
    if (altacid != acid)   // if non-synonymous
    {
      return false;
    }
  }

  return true;
}

/**********************************************************************************************/

Sequence* GeneticCode::getCodingSequence(const Sequence& sequence, bool lookForInitCodon, bool includeInitCodon) const throw (Exception)
{
  size_t initPos = 0;
  size_t stopPos = sequence.size();
  if (AlphabetTools::isCodonAlphabet(sequence.getAlphabet()))
  {
    // Look for AUG(or ATG) codon:
    if (lookForInitCodon)
    {
      for (size_t i = 0; i < sequence.size(); i++)
      {
        vector<int> pos = codonAlphabet_.getPositions(sequence[i]);
        if (pos[0] == 0 && pos[1] == 3 && pos[2] == 2)
        {
          initPos = includeInitCodon ? i : i + 1;
          break;
        }
      }
    }
    // Look for stop codon:
    for (size_t i = initPos; i < sequence.size(); i++)
    {
      if (isStop(sequence[i]))
      {
        stopPos = i;
        break;
      }
    }
  }
  else if (AlphabetTools::isNucleicAlphabet(sequence.getAlphabet()))
  {
    // Look for AUG(or ATG) codon:
    if (lookForInitCodon)
    {
      for (size_t i = 0; i < sequence.size() - 2; i++)
      {
        if (sequence[i] == 0 && sequence[i + 1] == 3 && sequence[i + 2] == 2)
        {
          initPos = includeInitCodon ? i : i + 3;
          break;
        }
      }
    }
    // Look for stop codon:
    const NucleicAlphabet* nucAlpha = codonAlphabet_.getNucleicAlphabet();
    for (size_t i = initPos; i < sequence.size() - 2; i += 3)
    {
      string codon = nucAlpha->intToChar(sequence[i])
                     + nucAlpha->intToChar(sequence[i + 1])
                     + nucAlpha->intToChar(sequence[i + 2]);
      if (isStop(codon))
      {
        stopPos = i;
        break;
      }
    }
  }
  else
    throw AlphabetMismatchException("Sequence must have alphabet of type nucleic or codon in GeneticCode::getCodingSequence.", 0, sequence.getAlphabet());

  return SequenceTools::subseq(sequence, initPos, stopPos - 1);
}

/**********************************************************************************************/

