//
// File: FullHmmTransitionMatrix.cpp
// Created by: Laurent Guéguen
// Created on: samedi 21 septembre 2013, à 14h 43
//

/*
Copyright or © or Copr. Bio++Development Team, (November 16, 2004)

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

#include "FullHmmTransitionMatrix.h"

#include "../../Text/TextTools.h"

#include "../Matrix/MatrixTools.h"
#include "../VectorTools.h"

using namespace bpp;
using namespace std;

FullHmmTransitionMatrix::FullHmmTransitionMatrix(const HmmStateAlphabet* alph, const string& prefix) :
  AbstractHmmTransitionMatrix(alph),
  AbstractParametrizable(prefix),
  vSimplex_()
{
  size_t size=(size_t)getNumberOfStates();

  for (size_t i=0; i<size; i++)
    {
      vSimplex_.push_back(Simplex(size,1,false,prefix + TextTools::toString(i+1)+"."));
      addParameters_(vSimplex_[i].getParameters());
    }
}

FullHmmTransitionMatrix::FullHmmTransitionMatrix(const FullHmmTransitionMatrix& hptm) :
  AbstractHmmTransitionMatrix(hptm),
  AbstractParametrizable(hptm),
  vSimplex_(hptm.vSimplex_)
{
}

FullHmmTransitionMatrix& FullHmmTransitionMatrix::operator=(const FullHmmTransitionMatrix& hptm)
{
  AbstractHmmTransitionMatrix::operator=(hptm);
  AbstractParametrizable::operator=(hptm);
  
  return *this;
}

void FullHmmTransitionMatrix::setTransitionProbabilities(const Matrix<double>& mat)
{
  if (mat.getNumberOfRows()!=vSimplex_.size())
    throw BadSizeException("FullHmmTransitionMatrix::setTransitionProbabilities: Wrong number of rows in given Matrix", mat.getNumberOfRows(), vSimplex_.size());
  
  ParameterList pl;
  
  for (size_t i=0; i<mat.getNumberOfRows();i++)
  {
    vSimplex_[i].setFrequencies(mat.row(i));
    ParameterList pls=vSimplex_[i].getParameters();
    for (size_t j=0; j<pls.size(); j++)
    {
      Parameter* p=pls[j].clone();
      p->setName(TextTools::toString(i+1)+"."+p->getName());
      pl.addParameter(p);
    }
  }
  
  matchParametersValues(pl);
}


const Matrix<double>& FullHmmTransitionMatrix::getPij() const
 {
   if (!upToDate_){
     for (size_t i=0; i<vSimplex_.size(); i++)
       for (size_t j=0; j<vSimplex_[i].dimension(); j++)
         pij_(i,j)=vSimplex_[i].prob(j);
     upToDate_=true;
   }
   
   return pij_;
 }

const std::vector<double>& FullHmmTransitionMatrix::getEquilibriumFrequencies() const
{
  size_t salph=getNumberOfStates();
  
  if (!upToDate_){
    pij_=getPij();

    MatrixTools::pow(pij_, 256, tmpmat_);

    for (size_t i = 0; i < salph; i++)
      eqFreq_[i] = tmpmat_(0,i);
    
    upToDate_=true;
  }

  return eqFreq_;
}

void FullHmmTransitionMatrix::fireParameterChanged(const ParameterList& parameters)
{
  size_t salph=getNumberOfStates();

  for (size_t i=0; i< salph; i++)
    vSimplex_[i].matchParametersValues(parameters);
  
  upToDate_=false;
}


