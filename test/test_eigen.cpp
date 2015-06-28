//
// File: test_eigen.cpp
// Created by: Julien Dutheil
// Created on: Thu Feb 5 07:50 2009
//

/*
Copyright or © or Copr. Bio++Development Team, (November 17, 2004)

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
  EigenValue<double> eigen(m);
  RowMatrix<double> D  = eigen.getD();
  const vector<double> L  = eigen.getRealEigenValues();
  RowMatrix<double> V1 = eigen.getV();
  RowMatrix<double> V2;
  MatrixTools::inv(V1, V2);
  cout << "M=" << endl;
  MatrixTools::print(m);
  cout << "D=" << endl;
  MatrixTools::print(D);
  cout << "V1=" << endl;
  MatrixTools::print(V1);
  cout << "V2=" << endl;
  MatrixTools::print(V2);
  RowMatrix<double> test;
  MatrixTools::mult(V1, L, V2, test);
  cout << "V1 . D . V2=" << endl;
  MatrixTools::print(test);
  return (test.equals(m) ? 0 : 1);
}
