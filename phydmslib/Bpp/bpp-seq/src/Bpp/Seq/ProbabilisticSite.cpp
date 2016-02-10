//
// File ProbabilisticSite.cpp
// Author: Murray Patterson
// Created on: Tue Oct 13 2015
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

#include "ProbabilisticSite.h"

using namespace bpp;

/****************************************************************************************/

BasicProbabilisticSite::BasicProbabilisticSite(const Alphabet* alpha) :
  AbstractCoreSite(), BasicProbabilisticSymbolList(alpha) {}

BasicProbabilisticSite::BasicProbabilisticSite(const Alphabet* alpha, int position) :
  AbstractCoreSite(position), BasicProbabilisticSymbolList(alpha) {}

BasicProbabilisticSite::BasicProbabilisticSite(const DataTable & site, const Alphabet* alpha) throw (Exception) :
  AbstractCoreSite(), BasicProbabilisticSymbolList(site, alpha) {}

BasicProbabilisticSite::BasicProbabilisticSite(const DataTable & site, const Alphabet* alpha, int position) throw (Exception) :
  AbstractCoreSite(position), BasicProbabilisticSymbolList(site, alpha) {}

/****************************************************************************************/

BasicProbabilisticSite::BasicProbabilisticSite(const ProbabilisticSite & site) :
  AbstractCoreSite(site), BasicProbabilisticSymbolList(site) {}

BasicProbabilisticSite::BasicProbabilisticSite(const BasicProbabilisticSite & site) :
  AbstractCoreSite(site), BasicProbabilisticSymbolList(site) {}

BasicProbabilisticSite & BasicProbabilisticSite::operator=(const ProbabilisticSite & site)
{
  AbstractCoreSite::operator=(site);
  BasicProbabilisticSymbolList::operator=(site);
  return *this;
}

BasicProbabilisticSite & BasicProbabilisticSite::operator=(const BasicProbabilisticSite & site)
{
  AbstractCoreSite::operator=(site);
  BasicProbabilisticSymbolList::operator=(site);
  return *this;
}
