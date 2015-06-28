//
// File: NucleicAcidsReplication.cpp
// Created by: Julien Dutheil
// Created on: Fri May 20 14:40 2005
//

/*
Copyright or © or Copr. Bio++ Development Team, (November 17, 2004)

This software is a computer program whose purpose is to provide classes
for sequences analysis.

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

#include "NucleicAcidsReplication.h"

using namespace bpp;

using namespace std;

NucleicAcidsReplication::NucleicAcidsReplication(const NucleicAlphabet* nuc1, const NucleicAlphabet* nuc2) :
  nuc1_(nuc1), nuc2_(nuc2), trans_()
{
  trans_[-1] = -1;
  trans_[0] = 3;
  trans_[1] = 2;
  trans_[2] = 1;
  trans_[3] = 0;

  trans_[4] = 9;
  trans_[5] = 8;
  trans_[6] = 6;
  trans_[7] = 7;
  trans_[8] = 5;
  trans_[9] = 4;

  trans_[10] = 13;
  trans_[11] = 12;
  trans_[12] = 11;
  trans_[13] = 10;

  trans_[14] = 14;
}

int NucleicAcidsReplication::translate(int state) const throw (BadIntException)
{
  nuc1_->intToChar(state);
  return trans_[state];
}

std::string NucleicAcidsReplication::translate(const std::string& state) const throw (BadCharException)
{
  int i = nuc1_->charToInt(state);
  return nuc2_->intToChar(trans_[i]);
}

Sequence* NucleicAcidsReplication::translate(const Sequence& sequence) const throw (AlphabetMismatchException)
{
  if (sequence.getAlphabet()->getAlphabetType() != getSourceAlphabet()->getAlphabetType())
    throw AlphabetMismatchException("NucleicAcidsReplication::translate", getSourceAlphabet(), getTargetAlphabet());
  BasicSequence* tSeq = new BasicSequence(sequence.getName(), "", sequence.getComments(), getTargetAlphabet());
  for (unsigned int i = 0; i < sequence.size(); i++)
  {
    tSeq->addElement(translate(sequence.getValue(i)));
  }
  //tSeq->setSense(!tSeq->getSense());
  return tSeq;
}


int NucleicAcidsReplication::reverse(int state) const throw (BadIntException) 
{
  nuc2_->intToChar(state);
  return trans_[state];
}

std::string NucleicAcidsReplication::reverse(const std::string& state) const throw (BadCharException)
{
  int i = nuc2_->charToInt(state);
  return nuc1_->intToChar(trans_[i]);
}

Sequence* NucleicAcidsReplication::reverse(const Sequence& sequence) const throw (AlphabetMismatchException, Exception)
{
  if (sequence.getAlphabet()->getAlphabetType() != getTargetAlphabet()->getAlphabetType())
    throw AlphabetMismatchException("NucleicAcidsReplication::reverse", getSourceAlphabet(), getTargetAlphabet());
  BasicSequence* rSeq = new BasicSequence(sequence.getName(), "", sequence.getComments(), getSourceAlphabet());
  for (unsigned int i = 0; i < sequence.size(); i++)
  {
    rSeq->addElement(reverse(sequence.getValue(i)));
  }
  //rSeq->setSense(! rSeq->getSense());
  return rSeq;
}

