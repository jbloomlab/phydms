//
// File: GenBank.cpp
// Created by: Julien Dutheil
// Created on: Tue Oct 2 2007
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

#include "GenBank.h"

#include <Bpp/Text/TextTools.h>
#include <Bpp/Text/StringTokenizer.h>

using namespace bpp;
using namespace std;

/****************************************************************************************/

void GenBank::appendSequencesFromStream(std::istream& input, SequenceContainer& vsc) const throw (Exception)
{
  if (!input) { throw IOException ("GenBank::read: fail to open file"); }

  string temp, name, sequence = "";  // Initialization

  // Main loop : for all file lines
  while (!input.eof())
  {
    getline(input, temp, '\n');  // Copy current line in temporary string

    if(temp.size() >= 9 && temp.substr(0,9) == "ACCESSION")
    {
      name = TextTools::removeSurroundingWhiteSpaces(temp.substr(10));
      StringTokenizer st(name, " ");
      name = st.nextToken();
      //cout << name << endl;
    }
    if (temp.size() >=6 && temp.substr(0,6) == "ORIGIN")
    {
      sequence = "";
      getline(input, temp, '\n');  // Copy current line in temporary string
      while (!input.eof() && temp.size() > 2 && temp.substr(0,2) != "//")
      {
        sequence += TextTools::removeWhiteSpaces(temp.substr(10));
        getline(input, temp, '\n');  // Copy current line in temporary string
      }
      if(name == "") throw Exception("GenBank::read(). Sequence with no ACCESSION number!");
      Sequence* seq = new BasicSequence(name, sequence, vsc.getAlphabet());
      vsc.addSequence(*seq, true);
      name = "";
    }
  }  
}

/****************************************************************************************/

