//
// File: Stockholm.cpp
// Authors: Julien Dutheil
// Created: Thu Apr 15 2010
//

/*
Copyright or © or Copr. Bio++ Development Team (2010)

Julien.Dutheil@univ-montp2.fr

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

#include "Stockholm.h"

#include "../StringSequenceTools.h"
#include <Bpp/Text/TextTools.h>
#include <Bpp/Text/StringTokenizer.h>
#include <Bpp/Io/FileTools.h>

using namespace bpp;
using namespace std;

/******************************************************************************/

void Stockholm::writeAlignment(ostream& output, const SiteContainer& sc) const throw (Exception)
{
	if (!output)
    throw IOException("Stockholm::writeAlignment: can't write to ostream output");

  output << "# STOCKHOLM 1.0" << endl; 
  // Loop for all general comments
  for (size_t i = 0; i < sc.getGeneralComments().size(); ++i)
  {
    output << "#=GF CC " << sc.getGeneralComments()[i] << endl;
  }

	// Main loop : for all sequences in vector container
	vector<string> names = sc.getSequencesNames();
  size_t maxSize = 0; 
  for(unsigned int i = 0; i < names.size(); ++i)
  {
    names[i] = TextTools::removeWhiteSpaces(names[i]);
    if (names[i].size() > maxSize) maxSize = names[i].size();
  }
  if (maxSize > 255) maxSize = 255;
  for (size_t i = 0; i < sc.getNumberOfSequences(); ++i)
  {
    output << TextTools::resizeRight(names[i], maxSize) << " " << sc.getSequence(i).toString() << endl;
	}
  output << "//" << endl;
}

/******************************************************************************/

