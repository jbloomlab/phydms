//
// File: test_distirbutions.cpp
// Created by: Julien Dutheil
// Created on: Fri Jun 8 20:54 2012
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

#include <Bpp/Numeric/Prob/ExponentialDiscreteDistribution.h>
#include <Bpp/Numeric/Prob/TruncatedExponentialDiscreteDistribution.h>

using namespace bpp;
using namespace std;

void testSumProbs(const DiscreteDistribution& dist) throw (Exception) {
  cout << "Test sum of probabilities: ";
  cout.flush();
  double s = VectorTools::sum(dist.getProbabilities());
  if (abs(s - 1.0) > 0.000001)
    throw Exception("Probabilities sum to " + TextTools::toString(s));
  cout << "ok" << endl;
}

void testEqualProbs(const DiscreteDistribution& dist) throw (Exception) {
  cout << "Test equality of probabilities: ";
  cout.flush();
  for (size_t i = 1; i < dist.getNumberOfCategories(); ++i)
    if (abs(dist.getProbability(i) - dist.getProbability(i - 1)) > 0.000001)
      throw Exception("Probability " + TextTools::toString(i) + " and " + TextTools::toString(i - 1) + " differ.");
  cout << "ok" << endl;
}

int main() {
  try {

    cout << "Testing Exponential distribution with lambda = 0.2:" << endl;
    for (unsigned int i = 2; i <= 100; i += 5) {
      cout << "Testing " << i << " classes." << endl;
      ExponentialDiscreteDistribution expDist(i, 0.2);
      testSumProbs(expDist);
      testEqualProbs(expDist);
    }

    cout << "Testing Exponential distribution with lambda = 5:" << endl;
    for (unsigned int i = 2; i <= 100;  i += 5) {
      cout << "Testing " << i << " classes." << endl;
      ExponentialDiscreteDistribution expDist(i, 5);
      testSumProbs(expDist);
      testEqualProbs(expDist);
    }

    cout << "Testing Truncated Exponential distribution with lambda = 0.2 and t = 0.1:" << endl;
    for (unsigned int i = 2; i <= 100; i += 5) {
      cout << "Testing " << i << " classes." << endl;
      TruncatedExponentialDiscreteDistribution trExpDist(i, 0.2, 0.1);
      testSumProbs(trExpDist);
      testEqualProbs(trExpDist);
    }

    cout << "Testing Truncated Exponential distribution with lambda = 5 and t = 0.1:" << endl;
    for (unsigned int i = 2; i <= 100; i += 5) {
      cout << "Testing " << i << " classes." << endl;
      TruncatedExponentialDiscreteDistribution trExpDist(i, 5, 0.1);
      testSumProbs(trExpDist);
      testEqualProbs(trExpDist);
    }

    cout << "Testing Truncated Exponential distribution with lambda = 0.2 and t = 4:" << endl;
    for (unsigned int i = 2; i <= 100; i += 5) {
      cout << "Testing " << i << " classes." << endl;
      TruncatedExponentialDiscreteDistribution trExpDist(i, 0.2, 4);
      testSumProbs(trExpDist);
      testEqualProbs(trExpDist);
    }

    cout << "Testing Truncated Exponential distribution with lambda = 5 and t = 4:" << endl;
    for (unsigned int i = 2; i <= 100; i += 5) {
      cout << "Testing " << i << " classes." << endl;
      TruncatedExponentialDiscreteDistribution trExpDist(i, 5, 4);
      testSumProbs(trExpDist);
      testEqualProbs(trExpDist);
    }

    //Now makes some more precise calculations:

    cout << "Check values for exponential distribution:" << endl;
    // x = -log(1-p)/lambda
    ExponentialDiscreteDistribution expDist(4, 0.2);
    cout << expDist.getLowerBound() << "\t" << "0" << endl;;
    cout << expDist.getBound(0)     << "\t" << "1.438410" << endl;;
    cout << expDist.getBound(1)     << "\t" << "3.465736" << endl;;
    cout << expDist.getBound(2)     << "\t" << "6.931472" << endl;;
    cout << expDist.getUpperBound() << "\t" << "+inf" << endl;;

    if (abs(expDist.getLowerBound() - 0) > 0.0001) throw Exception("Unvalid bound.");
    if (abs(expDist.getBound(0) - 1.438410) > 0.0001) throw Exception("Unvalid bound.");
    if (abs(expDist.getBound(1) - 3.465736) > 0.0001) throw Exception("Unvalid bound.");
    if (abs(expDist.getBound(2) - 6.931472) > 0.0001) throw Exception("Unvalid bound.");

    cout << "Check values for truncated exponential distribution:" << endl;
    // x = -log(1 - (1 - lambda*t) * p) / lambda
    TruncatedExponentialDiscreteDistribution trExpDist(4, 0.2, 1);
    cout << trExpDist.getLowerBound() << "\t" << "0" << endl;;
    cout << trExpDist.getBound(0)     << "\t" << "0.2318813" << endl;;
    cout << trExpDist.getBound(1)     << "\t" << "0.4750416" << endl;;
    cout << trExpDist.getBound(2)     << "\t" << "0.7306344" << endl;;
    cout << trExpDist.getUpperBound() << "\t" << "1" << endl;;

    if (abs(trExpDist.getLowerBound() - 0) > 0.0001) throw Exception("Unvalid bound.");
    if (abs(trExpDist.getBound(0) - 0.2318813) > 0.0001) throw Exception("Unvalid bound.");
    if (abs(trExpDist.getBound(1) - 0.4750416) > 0.0001) throw Exception("Unvalid bound.");
    if (abs(trExpDist.getBound(2) - 0.7306344) > 0.0001) throw Exception("Unvalid bound.");
    if (abs(trExpDist.getUpperBound() - 1) > 0.0001) throw Exception("Unvalid bound.");




    return 0;
  } catch(Exception& ex) {
    cout << "failed :(" << endl;
    return 1;
  }
}

