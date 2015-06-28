//
// File: StatTools.cpp
// Created by: Julien Dutheil
// Created on: Sun Jan 30 19:10 2011
//

/*
Copyright or © or Copr. Bio++ Development Team, (November 17, 2004)

This software is a computer program whose purpose is to provide classes
for numerical calculus.

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

#include "StatTools.h"

//From the STL:
#include <algorithm>

using namespace bpp;
using namespace std;

vector<double> StatTools::computeFdr(const vector<double>& pvalues) {
  size_t n = pvalues.size();
  vector<PValue_> sortedPValues;
  for (size_t i = 0; i < n; ++i) {
    sortedPValues.push_back(PValue_(pvalues[i], i));  
  }
  sort(sortedPValues.begin(), sortedPValues.end());
  vector<double> fdr(pvalues.size());
  for (size_t i = 0; i < sortedPValues.size(); ++i) {
    fdr[sortedPValues[i].index_] = sortedPValues[i].pvalue_ * static_cast<double>(n) / ( static_cast<double>(sortedPValues[i].index_ + 1));
  }
  return fdr;
}
