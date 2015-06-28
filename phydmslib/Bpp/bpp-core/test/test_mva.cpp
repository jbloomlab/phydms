//
// File: test_mva.cpp
// Created by: Matheu Groussin
// Created on: Tue Mar 22 14:54 2011
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

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <ctime>
#include <set>
#include <map>
#include <cmath>

using namespace std;

#include <Bpp/Numeric/Matrix/EigenValue.h>
#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/Numeric/Matrix/Matrix.h>
#include <Bpp/Numeric/Stat/Mva/DualityDiagram.h>
#include <Bpp/Numeric/Stat/Mva/PrincipalComponentAnalysis.h>
#include <Bpp/Numeric/Stat/Mva/CorrespondenceAnalysis.h>
#include <Bpp/Numeric/VectorTools.h>

using namespace bpp;

/*This program performs some tests for the Pca, Coa and Dudi classes*/

int main(int args, char ** argv)
{
	unsigned int nRows = 3;
	unsigned int nCols = 3;
	RowMatrix<double> matrix(nRows, nCols);
	RowMatrix<double> matrix2(nRows, nCols+1);
	RowMatrix<double> matrix3(nRows, nCols);
	
	matrix(0,0) = 10;
	matrix(1,0) = 20;
	matrix(2,0) = 30;
	matrix(0,1) = 20;
	matrix(1,1) = 10;
	matrix(2,1) = 40;
	matrix(0,2) = 30;
	matrix(1,2) = 40;
	matrix(2,2) = 10;
	
	cout << endl;
	cout << "First test for the Pca class, with a square matrix : " << endl;
	cout << endl;
	cout << "Here's the input matrix : " << endl;
	MatrixTools::print(matrix, cout);
	
	vector<double> rowW(nRows);
	vector<double> colW(nCols);
	VectorTools::fill(rowW, 1./static_cast<double>(nRows));
	VectorTools::fill(colW, 1.);
	
	//The constructor with row and column weights is called
	PrincipalComponentAnalysis* pca1 = new PrincipalComponentAnalysis(matrix, 3, rowW, colW, true, true);
	
	cout << "The matrix of Row Coordinates : " << endl;
	MatrixTools::print(pca1->getRowCoordinates(),cout);
	cout << endl;
	cout << endl;
	
			
	matrix2(0,0) = 10;
	matrix2(1,0) = 20;
	matrix2(2,0) = 30;
	matrix2(0,1) = 20;
	matrix2(1,1) = 10;
	matrix2(2,1) = 40;
	matrix2(0,2) = 30;
	matrix2(1,2) = 40;
	matrix2(2,2) = 10;
	matrix2(0,3) = 50;
	matrix2(1,3) = 10;
	matrix2(2,3) = 10;
	
	cout << endl;
	cout << "Second test for the Pca class, with a matrix containing more columns than rows: " << endl;
	cout << endl;
	cout << "Here's the input matrix : " << endl;
	MatrixTools::print(matrix2,cout);
	
	//The constructor without row and column weigths is called. Default weights will be created.
	PrincipalComponentAnalysis* pca2 = new PrincipalComponentAnalysis(matrix2, 3, true, true);
	
	cout << "The matrix of Principal Axes : " << endl;
	MatrixTools::print(pca2->getPrincipalAxes(),cout);
	cout << endl;
	cout << endl;
	
	
	matrix3(0,0) = 0.10;
	matrix3(1,0) = 0.20;
	matrix3(2,0) = 0.30;
	matrix3(0,1) = 0.40;
	matrix3(1,1) = 0.50;
	matrix3(2,1) = 0.60;
	matrix3(0,2) = 0.50;
	matrix3(1,2) = 0.30;
	matrix3(2,2) = 0.10;	
	
	cout << endl;
	cout << "Test for the Coa class, with a square matrix : " << endl;
	cout << endl;
	cout << "Here's the input matrix : " << endl;
	MatrixTools::print(matrix3,cout);
	
	//The Coa constructor is called.
	CorrespondenceAnalysis* coa = new CorrespondenceAnalysis(matrix3, 3);

	cout << "The matrix of Principal Components : " << endl;
	MatrixTools::print(coa->getPrincipalComponents(),cout);
	cout << endl;
	cout << endl;
	
	
	delete pca1;
	delete pca2;
	delete coa;
  return 0;	
}
















