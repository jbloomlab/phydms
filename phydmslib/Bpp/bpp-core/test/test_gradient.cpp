//
// File: test_gradient.cpp
// Created by: Julien Dutheil
// Created on: Thu Oct 28 10:28 2010
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

#include <Bpp/Numeric/Function/ConjugateGradientMultiDimensions.h>
#include <vector>
#include <iostream>
#include "PolynomialFunction.h"

using namespace bpp;
using namespace std;

int main() {
  PolynomialFunction1Der1 f;
  cout << f.getValue() << endl;
  ConjugateGradientMultiDimensions optimizer(&f);
  optimizer.init(f.getParameters());
  optimizer.optimize();
  double minf = optimizer.getFunctionValue();
  double x = f.getParameterValue("x");
  double y = f.getParameterValue("y");
  double z = f.getParameterValue("z");
  cout << "x=" << x << endl;
  cout << "y=" << y << endl;
  cout << "z=" << z << endl;
  cout << "f=" << minf << endl;
  cout << setprecision(20) << (abs(minf) + abs(x - 5) + abs(y + 2) + abs(z - 3)) << endl;
  bool test = abs(minf) + abs(x - 5) + abs(y + 2) + abs(z - 3) < optimizer.getStopCondition()->getTolerance();
  return (test ? 0 : 1);
}
