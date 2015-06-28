//
// File: test_stats.cpp
// Created by: Julien Dutheil
// Created on: Thu Dec 9 15:38 2010
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

#include <Bpp/Numeric/Stat/ContingencyTableTest.h>
#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Numeric/Random/ContingencyTableGenerator.h>
#include <Bpp/Numeric/Random/RandomTools.h>
#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <vector>
#include <iostream>
#include <cmath>

using namespace bpp;
using namespace std;

//tbl<-rbind(c(6,12,16,20),c(9,34,28,12))
//chisq.test(tbl);
int main() {
  vector< vector<size_t> > table;
  vector<size_t> row1;
  row1.push_back(6);
  row1.push_back(12);
  row1.push_back(16);
  row1.push_back(20);
  table.push_back(row1);
  vector<size_t> row2;
  row2.push_back(9);
  row2.push_back(34);
  row2.push_back(28);
  row2.push_back(12);
  table.push_back(row2);
  ContingencyTableTest test(table);
  VectorTools::print(test.getMarginRows());
  VectorTools::print(test.getMarginColumns());
  ContingencyTableGenerator ctRand(test.getMarginRows(), test.getMarginColumns());
  RowMatrix<size_t> rtable = ctRand.rcont2();
  MatrixTools::print(rtable);

  cout << test.getStatistic() << " \t" << test.getPValue() << endl;
  if (abs(test.getPValue() - 0.01324) > 0.0001)
    return 1;

  //Now test permutations:
  ContingencyTableTest test2(table, 20000);
  cout << test2.getStatistic() << " \t" << test2.getPValue() << endl;
  if (abs(test2.getPValue() - 0.01324) > 0.01)
    return 1;

  return 0;
}

