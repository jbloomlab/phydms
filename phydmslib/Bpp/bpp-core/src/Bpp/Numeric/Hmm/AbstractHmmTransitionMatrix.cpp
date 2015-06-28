//
// File: AbstractHmmTransitionMatrix.cpp
// Created by: Laurent Guéguen
// Created on: lundi 10 février 2014, à 10h 59
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

#include "AbstractHmmTransitionMatrix.h"

#include "../Matrix/MatrixTools.h"
#include "../VectorTools.h"
#include "../Random/RandomTools.h"

using namespace bpp;
using namespace std;

AbstractHmmTransitionMatrix::AbstractHmmTransitionMatrix(const HmmStateAlphabet* alph, const string& prefix) :
  alph_(alph),
  pij_((size_t)alph->getNumberOfStates(), (size_t)alph->getNumberOfStates()),
  tmpmat_((size_t)alph->getNumberOfStates(), (size_t)alph->getNumberOfStates()),
  eqFreq_((size_t)alph->getNumberOfStates()),
  upToDate_(false)
{
}

AbstractHmmTransitionMatrix::AbstractHmmTransitionMatrix(const AbstractHmmTransitionMatrix& hptm) :
  alph_(hptm.alph_),
  pij_(hptm.pij_),
  tmpmat_(hptm.tmpmat_),
  eqFreq_(hptm.eqFreq_),
  upToDate_(hptm.upToDate_)
{
}

AbstractHmmTransitionMatrix& AbstractHmmTransitionMatrix::operator=(const AbstractHmmTransitionMatrix& hptm)
{
  alph_=hptm.alph_;
  pij_=hptm.pij_;
  tmpmat_=hptm.tmpmat_;
  eqFreq_=hptm.eqFreq_;
  upToDate_=hptm.upToDate_;
  
  return *this;
}

void AbstractHmmTransitionMatrix::setHmmStateAlphabet(const HmmStateAlphabet* stateAlphabet) throw (HmmUnvalidAlphabetException)
{
  if (stateAlphabet==NULL)
    throw HmmUnvalidAlphabetException("Null alphabet in AbstractHmmTransitionMatrix::setHmmStateAlphabet");

  alph_=stateAlphabet;
}

vector<size_t> AbstractHmmTransitionMatrix::sample(size_t size) const
{
  vector<size_t> vres;
  if (size==0)
    return vres;

  size_t nbStates=getHmmStateAlphabet()->getNumberOfStates();

  // update pij_
  getPij();
    
  size_t sta=0, stb;
  double prob = RandomTools::giveRandomNumberBetweenZeroAndEntry(1.0);

  for (size_t i = 0; i < nbStates; ++i) {
    prob-=eqFreq_[i];
    if (prob < 0) {
      sta=i;
      break;
    }
  }
        
  vres.push_back(sta);

  for (size_t pos=1;pos<size;pos++)
  {
    prob = RandomTools::giveRandomNumberBetweenZeroAndEntry(1.0);

    const vector<double>& row=pij_.getRow(sta);

    for (size_t i = 0; i < nbStates; ++i) {
      prob-=row[i];
      if (prob < 0) {
        stb=i;
        break;
      }
    }
    vres.push_back(stb);
    sta=stb;
  }
  return vres;
}

  
