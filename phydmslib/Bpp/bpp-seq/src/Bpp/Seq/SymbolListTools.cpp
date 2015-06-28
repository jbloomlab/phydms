//
// File: SymbolListTools.cpp
// Created by: Julien Dutheil
// Created on: Wed Apr 9 2004
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

#include "SymbolListTools.h"
#include "Alphabet/AlphabetTools.h"
#include <Bpp/Numeric/Random/RandomTools.h>

//From the STL:
#include <algorithm>

using namespace std;

using namespace bpp;

void SymbolListTools::getCounts(const SymbolList& list, map<int, double>& counts, bool resolveUnknowns)
{
  if (!resolveUnknowns)
  {
    for (size_t i = 0; i < list.size(); ++i)
      counts[list[i]]++;
  }
  else
  {
    for (size_t i = 0; i < list.size(); ++i)
    {
      vector<int> alias = list.getAlphabet()->getAlias(list[i]);
      double n = static_cast<double>(alias.size());
      for (size_t j = 0; j < alias.size(); j++) counts[alias[j]] += 1./n ;
    }
  }
}

void SymbolListTools::getCounts(const SymbolList& list1, const SymbolList& list2,  map< int, map<int, double> >& counts, bool resolveUnknowns) throw (DimensionException)
{
  if (list1.size() != list2.size()) throw DimensionException("SymbolListTools::getCounts: the two sites must have the same size.", list1.size(), list2.size());
  if (!resolveUnknowns)
  {
    for (size_t i = 0; i < list1.size(); i++)
      counts[list1[i]][list2[i]]++;
  }
  else
  {
    for (size_t i = 0; i < list1.size(); i++)
    {
      vector<int> alias1 = list1.getAlphabet()->getAlias(list1[i]);
      vector<int> alias2 = list2.getAlphabet()->getAlias(list2[i]);
      double n1 = (double)alias1.size();
      double n2 = (double)alias2.size();
      for (size_t j = 0; j < alias1.size(); j++)
        for (size_t k = 0; k < alias2.size(); k++)
          counts[alias1[j]][alias2[k]] += 1./(n1*n2) ;
    }
  }
}

void SymbolListTools::getFrequencies(const SymbolList& list, map<int, double>& frequencies, bool resolveUnknowns)
{
	double n = (double)list.size();
  map<int, double> counts;
  getCounts(list, counts, resolveUnknowns);
  for (map<int, double>::iterator i = counts.begin(); i != counts.end(); i++)
  {
    frequencies[i->first] = i->second / n;
  }
}

void SymbolListTools::getFrequencies(const SymbolList& list1, const SymbolList& list2, map<int, map<int, double> >& frequencies, bool resolveUnknowns) throw (DimensionException)
{
	double n2 = (double)list1.size() * (double)list1.size();
  map<int, map<int, double> > counts;
  getCounts(list1, list2, counts, resolveUnknowns);
  for (map<int, map<int, double> >::iterator i = counts.begin(); i != counts.end(); i++)
    for (map<int, double>::iterator j = i->second.begin(); j != i->second.end(); j++)
  {
    frequencies[i->first][j->first] = j->second / n2;
  }
}

double SymbolListTools::getGCContent(const SymbolList& list, bool ignoreUnresolved, bool ignoreGap) throw (AlphabetException)
{
  const Alphabet * alphabet = list.getAlphabet();
  if (!AlphabetTools::isNucleicAlphabet(alphabet))
    throw AlphabetException("SymbolListTools::getGCContent. Method only works on nucleotides.", alphabet);
  double gc = 0;
  double total = 0;
  for (size_t i = 0; i < list.size(); i++) {
    int state = list.getValue(i);
    if (state > -1) { // not a gap
      if (state == 1 || state == 2) { // G or C
        gc++;
        total++;
      } else if (state == 0 || state == 3) { // A, T or U
        total++;
      } else { // Unresolved character
        if (!ignoreUnresolved) {
          total++;
          switch(state) {
            case(7): gc++; break;// G or C
            case(4): gc+=0.5; break;// A or C
            case(5): gc+=0.5; break;// A or G
            case(6): gc+=0.5; break;// C or T
            case(9): gc+=0.5; break;// G or T
            case(10): gc+=2./3.; break;// A or C or G
            case(11): gc+=1./3.; break;// A or C or T
            case(12): gc+=1./3.; break;// A or G or T
            case(13): gc+=2./3.; break;// C or G or T
            case(14): gc+=0.5; break;// A or C or G or T
          }
        }
      }
    } else {
      if (!ignoreGap) total++;
    }
  }
  return total != 0 ? gc/total : 0;
}

size_t SymbolListTools::getNumberOfDistinctPositions(const SymbolList& l1, const SymbolList& l2) throw (AlphabetMismatchException)
{
	if (l1.getAlphabet()->getAlphabetType() != l2.getAlphabet()->getAlphabetType()) throw AlphabetMismatchException("SymbolListTools::getNumberOfDistinctPositions.", l1.getAlphabet(), l2.getAlphabet());
	size_t n = min(l1.size(), l2.size());
	size_t count = 0;
	for (size_t i = 0; i < n; i++) {
		if (l1[i] != l2[i]) count++;
	}
	return count;
}

size_t SymbolListTools::getNumberOfPositionsWithoutGap(const SymbolList& l1, const SymbolList& l2) throw (AlphabetMismatchException)
{
	if (l1.getAlphabet() -> getAlphabetType() != l2.getAlphabet() -> getAlphabetType()) throw AlphabetMismatchException("SymbolListTools::getNumberOfDistinctPositions.", l1.getAlphabet(), l2.getAlphabet());
	size_t n = min(l1.size(), l2.size());
	size_t count = 0;
	for (size_t i = 0; i < n; i++) {
		if (l1[i] != -1 && l2[i] != -1) count++;
	}
	return count;
}

void SymbolListTools::changeGapsToUnknownCharacters(SymbolList& l)
{
  int unknownCode = l.getAlphabet()->getUnknownCharacterCode();
  for (size_t i = 0; i < l.size(); i++)
  {
    if (l.getAlphabet()->isGap(l[i])) l[i] = unknownCode;
  }
}

void SymbolListTools::changeUnresolvedCharactersToGaps(SymbolList& l)
{
  int gapCode = l.getAlphabet()->getGapCharacterCode();
  for (size_t i = 0; i < l.size(); i++)
  {
    if (l.getAlphabet()->isUnresolved(l[i])) l[i] = gapCode;
  }
}

