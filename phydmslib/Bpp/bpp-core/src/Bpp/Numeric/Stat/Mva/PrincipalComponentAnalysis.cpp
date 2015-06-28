//
// File: PrincipalComponentAnalysis.cpp
// Created by: Mathieu Groussin
//

/*
   Copyright or © or Copr. Bio++ Development Team, (November 16, 2004)

   This software is a computer program whose purpose is to provide classes
   for phylogenetic data analysis.

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

#include "PrincipalComponentAnalysis.h"
#include "../../Matrix/Matrix.h"
#include "../../Matrix/MatrixTools.h"
#include "../../VectorTools.h"
#include "DualityDiagram.h"

#include <cmath>

using namespace bpp;
using namespace std;

PrincipalComponentAnalysis::PrincipalComponentAnalysis(
  const Matrix<double>& data,
  unsigned int nbAxes,
  const vector<double>& rowW,
  const vector<double>& colW,
  bool centered,
  bool scaled,
  double tol,
  bool verbose) throw (Exception) :
  DualityDiagram(),
  columnMeans_(),
  columnSd_()
{
  RowMatrix<double> tmpData = data;
  
  // Centering of data?
  if (centered)
  {
    center(tmpData, rowW);
  }
  
  // Scaling of data?
  if (scaled)
  {
    scale(tmpData, rowW);
  }

  setData(tmpData, rowW, colW, nbAxes, tol, verbose);
}

/******************************************************************************/

PrincipalComponentAnalysis::PrincipalComponentAnalysis(
  const Matrix<double>& data,
  unsigned int nbAxes,
  bool centered,
  bool scaled,
  double tol,
  bool verbose) throw (Exception) :
  DualityDiagram(),
  columnMeans_(),
  columnSd_()
{
  size_t nRow = data.getNumberOfRows();
  size_t nCol = data.getNumberOfColumns();

  vector<double> rowW(nRow);
  vector<double> colW(nCol);
  VectorTools::fill(rowW, 1. / static_cast<double>(nRow));
  VectorTools::fill(colW, 1.);

  RowMatrix<double> tmpData = data;

  // Centering of data?
  if (centered)
  {
    center(tmpData, rowW);
  }
  
  // Scaling of data?
  if (scaled)
  {
    scale(tmpData, rowW);
  }

  setData(tmpData, rowW, colW, nbAxes, tol, verbose);
}

/******************************************************************************/

void PrincipalComponentAnalysis::center(Matrix<double>& matrix, const vector<double>& rowW) throw (Exception)
{
  size_t nRow = matrix.getNumberOfRows();
  size_t nCol = matrix.getNumberOfColumns();
  if (nRow != rowW.size())
    throw Exception("PrincipalComponentAnalysis::center. The number of row weigths have to be equal to the number of rows!");

  double sumRowWeights = VectorTools::sum(rowW);

  vector<double> columnMeans(nCol);
  for (unsigned int i = 0; i < nCol; i++)
  {
    double tmp = 0.;
    for (unsigned int j = 0; j < nRow; j++)
    {
      tmp += matrix(j, i) * rowW[j];
    }
    columnMeans[i] = tmp / sumRowWeights;
  }

  for (unsigned int i = 0; i < nCol; i++)
  {
    for (unsigned int j = 0; j < nRow; j++)
    {
      matrix(j, i) -= columnMeans[i];
    }
  }
}

/******************************************************************************/

void PrincipalComponentAnalysis::scale(Matrix<double>& matrix, const vector<double>& rowW) throw (Exception)
{
  size_t nRow = matrix.getNumberOfRows();
  size_t nCol = matrix.getNumberOfColumns();
  if (nRow != rowW.size())
    throw Exception("PrincipalComponentAnalysis::scale. The number of row weigths have to be equal to the number of rows!");

  double sumRowWeights = VectorTools::sum(rowW);

  vector<double> columnSd(nCol);
  for (size_t i = 0; i < nCol; i++)
  {
    double tmp = 0.;
    for (unsigned int j = 0; j < nRow; j++)
    {
      tmp += pow(matrix(j, i), 2) * rowW[j];
    }
    columnSd[i] = sqrt(tmp / sumRowWeights);
  }

  for (size_t i = 0; i < nCol; i++)
  {
    for (unsigned int j = 0; j < nRow; j++)
    {
      if (columnSd[i] == 0.)
        matrix(j, i) = 0.;
      else
        matrix(j, i) /= columnSd[i];
    }
  }
}

