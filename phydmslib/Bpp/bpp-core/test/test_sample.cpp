//
// File: test_sample.cpp
// Created by: Julien Dutheil
// Created on: Wed Sep 28 14:36 2011
//

/*
Copyright or © or Copr. Bio++ Development Team, (November 17, 2004)

This software is a computer program whose purpose is to provide classes
for numerical calculus. This file is part of the Bio++ project.

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

#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Numeric/Random/RandomTools.h>
#include <vector>
#include <string>
#include <iostream>
#include <cmath>

using namespace bpp;
using namespace std;

int main() {
  //Create vector:
  vector<string> pop;
  pop.push_back("A");
  pop.push_back("B");
  pop.push_back("C");
  pop.push_back("D");
  pop.push_back("E");
  unsigned int n = 10000;

  cout << "-*- Check without replacement -*-" << endl;
  for (unsigned int k = 1; k < 5; ++k) {
    map<string, unsigned int> counts;
    for (unsigned int i = 0; i < n; ++i) {
      vector<string> sample(k);
      RandomTools::getSample(pop, sample, false);
      for (size_t j = 0; j < sample.size(); ++j) {
        counts[sample[j]]++;
      }
    }
    for (map<string, unsigned int>::iterator it = counts.begin(); it != counts.end(); ++it) {
      double fobs = static_cast<double>(it->second) / static_cast<double>(n);
      double fexp = static_cast<double>(k) / 5.;
      cout << it->first << "\t" << it->second << "\t" << fobs << "\t" << fexp << endl;
      if (abs(fobs - fexp) > 0.1)
        return 1;
    }
    cout << "---------------------------------------" << endl;
  }

  cout << "-*- Check with replacement -*-" << endl;
  for (unsigned int k = 1; k < 5; ++k) {
    map<string, unsigned int> counts;
    for (unsigned int i = 0; i < n; ++i) {
      vector<string> sample(k);
      RandomTools::getSample(pop, sample, true);
      for (size_t j = 0; j < sample.size(); ++j) {
        counts[sample[j]]++;
      }
    }
    for (map<string, unsigned int>::iterator it = counts.begin(); it != counts.end(); ++it) {
      double fobs = static_cast<double>(it->second) / static_cast<double>(n);
      double fexp = static_cast<double>(k) / 5.;
      cout << it->first << "\t" << it->second << "\t" << fobs << "\t" << fexp << endl;
      if (abs(fobs - fexp) > 0.1)
        return 1;
    }
    cout << "---------------------------------------" << endl;
  }

  cout << "-*- Check with replacement and weights -*-" << endl;
  vector<double> weights;
  weights.push_back(2);
  weights.push_back(3);
  weights.push_back(8);
  weights.push_back(2);
  weights.push_back(1);
  double sumw = VectorTools::sum(weights);
  vector<double> fexp = weights / sumw;
  for (unsigned int k = 1; k < 5; ++k) {
    map<string, unsigned int> counts;
    for (unsigned int i = 0; i < n; ++i) {
      vector<string> sample(k);
      RandomTools::getSample(pop, weights, sample, true);
      for (size_t j = 0; j < sample.size(); ++j) {
        counts[sample[j]]++;
      }
    }
    for (size_t i = 0; i < pop.size(); ++i) {
      double fobs = static_cast<double>(counts[pop[i]]) / static_cast<double>(n*k);
      cout << pop[i] << "\t" << counts[pop[i]] << "\t" << fobs << "\t" << fexp[i] << endl;
      if (abs(fobs - fexp[i]) > 0.1)
        return 1;
    }
    cout << "---------------------------------------" << endl;
  }



  return 0;
}

