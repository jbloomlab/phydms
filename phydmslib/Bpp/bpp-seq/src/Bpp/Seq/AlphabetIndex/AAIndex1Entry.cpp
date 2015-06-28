//
// File: AAIndex1Entry.cpp
// Created by: Julien Dutheil
// Created on: Fri Jan 19 17:07 2007
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

#include "AAIndex1Entry.h"
#include "../Alphabet/AlphabetTools.h"
#include <Bpp/Io/FileTools.h>
#include <Bpp/Text/TextTools.h>
#include <Bpp/Text/StringTokenizer.h>

using namespace bpp;
using namespace std;

AAIndex1Entry::AAIndex1Entry(std::istream& input) throw (IOException) :
  property_(20),
  alpha_(&AlphabetTools::PROTEIN_ALPHABET)
{
  // Parse entry:
  string line;
  bool ok = false;
  do
  {
    line = FileTools::getNextLine(input);
    if (line[0] == 'I')
    {
      string line1 = FileTools::getNextLine(input);
      string line2 = FileTools::getNextLine(input);
      StringTokenizer st1(line1, " ");
      StringTokenizer st2(line2, " ");
      if (st1.numberOfRemainingTokens() != 10 || st1.numberOfRemainingTokens() != 10)
        break;
      // Amino acids are in the same order in the AAIndex1 database than in the ProteicAlphabet class:
      for (unsigned int i = 0; i < 10; i++)
      {
        property_[i] = TextTools::toDouble(st1.nextToken());
      }
      for (unsigned int i = 10; i < 20; i++)
      {
        property_[i] = TextTools::toDouble(st2.nextToken());
      }
      // Jump to next entry...
      FileTools::getNextLine(input);
      ok = true;
    }
  }
  while (!ok);
  if (!ok)
    throw IOException("AAIndex1Entry: invalid AAIndex1 entry.");
}

