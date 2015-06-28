//
// File: HmmLikelihood.cpp
// Created by: Laurent Guéguen
// Created on: jeudi 26 septembre 2013, à 13h 55
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

#include "HmmLikelihood.h"

using namespace bpp;
using namespace std;

AbstractHmmLikelihood::AbstractHmmLikelihood() :
  dLogLik_(0),
  dVariable_(""),
  d2LogLik_(0),
  d2Variable_("") {}

AbstractHmmLikelihood::AbstractHmmLikelihood(const AbstractHmmLikelihood& adhlik) :
  dLogLik_(adhlik.dLogLik_),
  dVariable_(adhlik.dVariable_),
  d2LogLik_(adhlik.d2LogLik_),
  d2Variable_(adhlik.d2Variable_)
{}

AbstractHmmLikelihood& AbstractHmmLikelihood::operator=(const AbstractHmmLikelihood& adhlik)
{
  dLogLik_=adhlik.dLogLik_;
  dVariable_=adhlik.dVariable_;
  d2LogLik_=adhlik.d2LogLik_;
  d2Variable_=adhlik.d2Variable_;

  return *this;
}

double AbstractHmmLikelihood::getFirstOrderDerivative(const std::string& variable) const throw (Exception)
{
  if (variable!=dVariable_){
    dVariable_=variable;
    
    getHmmEmissionProbabilities().computeDEmissionProbabilities(dVariable_);
    computeDLikelihood_();
  }
  return -dLogLik_;
  
}
    
double AbstractHmmLikelihood::getSecondOrderDerivative(const std::string& variable) const throw (Exception)
{
  if (variable!=d2Variable_){
    d2Variable_=variable;

    getHmmEmissionProbabilities().computeD2EmissionProbabilities(d2Variable_);
    computeD2Likelihood_();
  }
  return -d2LogLik_;
}

