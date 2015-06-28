//
// File: SequenceWalker.cpp
// Created by: Julien Dutheil
// Created on: Thu Nov 24 2011
//

/*
Copyright or © or Copr. Bio++ Development Team, (November 17, 2011)

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

#include <iostream>
using namespace std;

#include "SequenceWalker.h"

#include <Bpp/Text/TextTools.h>

using namespace bpp;

size_t SequenceWalker::getAlignmentPosition(size_t seqPos) throw (Exception)
{
  if (seqPos == seqPos_) return alnPos_;
  if (seqPos > seqPos_) {
    //Move forward
    while (alnPos_ < seq_->size() && seqPos_ < seqPos) {
      if (alnPos_ == seq_->size() - 1)
        throw Exception("SequenceWalker::getAlignmentPosition(). Forward1. Position out of bound.");
      ++alnPos_;

      if ((*seq_)[alnPos_] != gap_) {
        ++seqPos_;
      }
    }
    if (seqPos_ != seqPos)
      throw Exception("SequenceWalker::getAlignmentPosition(). Forward2. Position out of bound (" + TextTools::toString(alnPos_) + ")");
  } else {
    //Move backward
    if (alnPos_ == 0)
      throw Exception("SequenceWalker::getAlignmentPosition(). Backward. Position out of bound.");
    while (alnPos_ > 0 && seqPos_ > seqPos) {
      --alnPos_;
      if ((*seq_)[alnPos_] != gap_) {
        --seqPos_;
      }
    }
    if (seqPos_ != seqPos)
      throw Exception("SequenceWalker::getAlignmentPosition(). Position out of bound.");
  }
  return alnPos_;
}

size_t SequenceWalker::getSequencePosition(size_t alnPos) throw (Exception)
{
  if (alnPos == alnPos_) return seqPos_;
  if (alnPos >= seq_->size())
    throw Exception("SequenceWalker::getSequencePosition(). Position out of bound.");
  if (alnPos > alnPos_) {
    //Move forward
    while (alnPos_ < alnPos) {
      ++alnPos_;
      if ((*seq_)[alnPos_] != gap_) {
        ++seqPos_;
      }
    }
  } else {
    //Move backward
    while (alnPos_ > alnPos) {
      if (seqPos_ == 0)
        return 0;
      --alnPos_;
      if ((*seq_)[alnPos_ + 1] != gap_) {
        --seqPos_;
      }
    }
  }
  return seqPos_;
}

