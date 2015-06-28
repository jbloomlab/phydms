//
// File: AutoCorrelationTransitionMatrix.cpp
// Created by: Laurent Guéguen
// Created on: lundi 10 février 2014, à 09h 56
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

#include "AutoCorrelationTransitionMatrix.h"

#include "../../Text/TextTools.h"

#include "../Matrix/MatrixTools.h"
#include "../VectorTools.h"

using namespace bpp;
using namespace std;

AutoCorrelationTransitionMatrix::AutoCorrelationTransitionMatrix(const HmmStateAlphabet* alph, const string& prefix) :
  AbstractHmmTransitionMatrix(alph),
  AbstractParametrizable(prefix),
  vAutocorrel_()
{
  size_t size=(size_t)getNumberOfStates();

  for (size_t i=0; i<size; i++)
    {
      vAutocorrel_.push_back(1./(double)size);
      addParameter_(new Parameter(prefix + "lambda"+TextTools::toString(i+1), 1./(double)size, &Parameter::PROP_CONSTRAINT_EX));
    }

  for (size_t i = 0; i < size; i++)
    eqFreq_[i] = 1./(double)size;
}


AutoCorrelationTransitionMatrix::AutoCorrelationTransitionMatrix(const AutoCorrelationTransitionMatrix& aptm) :
  AbstractHmmTransitionMatrix(aptm),
  AbstractParametrizable(aptm),
  vAutocorrel_(aptm.vAutocorrel_)
{
}

AutoCorrelationTransitionMatrix& AutoCorrelationTransitionMatrix::operator=(const AutoCorrelationTransitionMatrix& aptm)
{
  AbstractHmmTransitionMatrix::operator=(aptm);
  AbstractParametrizable::operator=(aptm);
  
  vAutocorrel_=aptm.vAutocorrel_;
  
  return *this;
}

const Matrix<double>& AutoCorrelationTransitionMatrix::getPij() const
 {
   if (!upToDate_){
     for (size_t i = 0; i < vAutocorrel_.size(); ++i)
       for (size_t j = 0; j < vAutocorrel_.size(); ++j)
         pij_(i,j) = (i==j) ? vAutocorrel_[i] : (1 - vAutocorrel_[i]) / static_cast<double>(getNumberOfStates()-1);

     upToDate_ = true;
   }
   
   return pij_;
 }

const std::vector<double>& AutoCorrelationTransitionMatrix::getEquilibriumFrequencies() const
{
  return eqFreq_;
}

void AutoCorrelationTransitionMatrix::fireParameterChanged(const ParameterList& parameters)
{
  size_t salph=getNumberOfStates();

  for (size_t i=0; i< salph; i++)
    vAutocorrel_[i]=getParameterValue("lambda"+TextTools::toString(i+1));
  
  upToDate_=false;
}


