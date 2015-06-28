//
// File: test_matrices.cpp
// Created by: Julien Dutheil
// Created on: Nov 2011
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

#include <Bpp/Numeric/Matrix/Matrix.h>
#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <vector>
#include <iostream>

using namespace bpp;
using namespace std;

int main() {
  RowMatrix<double> m(2,2);
  m(0,0) = 2.3;
  m(0,1) = 1.4;
  m(1,0) = 5.0;
  m(1,1) = -0.9;
  MatrixTools::print(m);
  RowMatrix<double> n(2, 2);
  MatrixTools::transpose(m, n),
  MatrixTools::print(n);
  RowMatrix<double> m2(2, 2);
  MatrixTools::transpose(n, m2);
  RowMatrix<double> o(2, 2);
  MatrixTools::mult(m, n, o);
  MatrixTools::print(o);
 
  bool test = m.equals(m2, 0.000001);
  ApplicationTools::displayBooleanResult("Test passed", test);
  return (test ? 0 : 1);
}
