//
// File: PhredPoly.cpp
// Created by: Sylvain Gaillard
// Created on: Fri Oct 31 2008
//

/*
Copyright or © or Copr. CNRS, (October 31, 2008)

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

#include "PhredPoly.h"
#include <Bpp/Text/TextTools.h>
#include <Bpp/Text/StringTokenizer.h>
#include <Bpp/Numeric/NumTools.h>

using namespace bpp;
using namespace std;

/******************************************************************************/

PhredPoly::PhredPoly(double ratio) : ratio_(ratio) {}

/******************************************************************************/

bool PhredPoly::nextSequence(istream& input, Sequence& seq) const throw (Exception) {
  if (!input) { throw IOException ("PhredPoly::read: fail to open stream"); }

  string temp, name, sequence = "";  // Initialization
  bool flag = false;

  // Read first line
  if (!input.eof()) {
    getline(input, temp, '\n');  // Copy current line in temporary string
    StringTokenizer st(temp, " ");
    name = st.getToken(0);
  }

  const Alphabet* alpha = seq.getAlphabet();

  // Main loop : for all other lines
  while (!input.eof()) {
    getline(input, temp, '\n');  // Copy current line in temporary string
    StringTokenizer st(temp, " ");
    if (st.numberOfRemainingTokens() == 12) {
      double a = TextTools::toDouble(st.getToken(3));
      double b = TextTools::toDouble(st.getToken(7));
      if (a < b) {
        NumTools::swap(a, b);
      }
      vector<string> v;
      v.push_back(st.getToken(0)); // Get the called base
      if (b / a > this->ratio_) {
        v.push_back(st.getToken(4)); // Get the uncalled base if relative picks areas are similar
      }
      sequence += alpha->getGeneric(v);
    }
  }
  if(name == "") {
    throw Exception("PhredPoly::read: sequence without name!");
  } else {
    seq.setName(name);
    seq.setContent(sequence);
    flag = true;
  }
  return flag;
}

/******************************************************************************/
