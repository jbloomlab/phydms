//
// File: DCSE.cpp
// Created by: Julien Dutheil
// Created on: Wed Mar 3 2004
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

#include "Dcse.h"
#include "AbstractIAlignment.h"
#include "../Sequence.h"
#include "../Container/SequenceContainer.h"
#include "../Container/VectorSequenceContainer.h"
#include <Bpp/Text/TextTools.h>
#include <Bpp/Text/StringTokenizer.h>
#include <Bpp/Io/FileTools.h>

using namespace bpp;
using namespace std;

void DCSE::appendAlignmentFromStream(istream& input, SiteContainer& sc) const throw (Exception)
{
  // Checking the existence of specified file
  if (!input) { throw IOException ("DCSE::read : fail to open file"); }

  // Initialization
  const Alphabet * alpha = sc.getAlphabet();
  string line, name, sequence = "";

  line = FileTools::getNextLine(input); // Copy current line in temporary string
  //StringTokenizer st(line);
  //st.nextToken();
  //First line ignored for now!
  //int n1 = TextTools::toInt(st.nextToken());
  //int n2 = TextTools::toInt(st.nextToken());
  //int nbSites = n2 - n1
  //cout << nbSpecies << " species and " << nbSites << " sites." << endl;

  // Main loop : for all file lines
  while (!input.eof())
  {
    line = FileTools::getNextLine(input); // Copy current line in temporary string
    if(line == "") break;
    string::size_type endOfSeq = line.find("     ");
    if(endOfSeq == line.npos) break;
    sequence = string(line.begin(), line.begin() + static_cast<ptrdiff_t>(endOfSeq));
    sequence = TextTools::removeWhiteSpaces(sequence);
    sequence = TextTools::removeChar(sequence, '{');
    sequence = TextTools::removeChar(sequence, '}');
    sequence = TextTools::removeChar(sequence, '[');
    sequence = TextTools::removeChar(sequence, ']');
    sequence = TextTools::removeChar(sequence, '(');
    sequence = TextTools::removeChar(sequence, ')');
    sequence = TextTools::removeChar(sequence, '^');
    name     = string(line.begin() + static_cast<ptrdiff_t>(endOfSeq + 1), line.end()),
    name     = TextTools::removeFirstWhiteSpaces(name);
    if(name.find("Helix numbering") == name.npos
    && name.find("mask") == name.npos)
      sc.addSequence(BasicSequence(name, sequence, alpha), true);
  }
}

const string DCSE::getFormatName() const { return "DCSE"; }

const string DCSE::getFormatDescription() const { return "RNA structure format"; }

