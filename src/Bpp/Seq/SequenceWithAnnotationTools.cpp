//
// File:       SequenceWithAnnotationTools.cpp
// Authors:    Julien Dutheil
// Created on: 06 Sep 2010
//

/*
Copyright or © or Copr. Bio++ Development Team, (Sep 06, 2010)

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

#include "SequenceWithAnnotationTools.h"
#include "Alphabet/CaseMaskedAlphabet.h"

using namespace bpp;
using namespace std;

const string SequenceMask::MASK = "Boolean mask";

/******************************************************************************/

void SequenceMask::afterSequenceChanged(const SymbolListEditionEvent& event)
{
  mask_.clear();
  mask_.insert(mask_.begin(), event.getSymbolList()->size(), false);
}

/******************************************************************************/

void SequenceMask::afterSequenceInserted(const SymbolListInsertionEvent& event)
{
  mask_.insert(mask_.begin() + static_cast<ptrdiff_t>(event.getPosition()),
      event.getLength(), false);
}

/******************************************************************************/

void SequenceMask::afterSequenceDeleted(const SymbolListDeletionEvent& event)
{
  mask_.erase(mask_.begin() + static_cast<ptrdiff_t>(event.getPosition()),
      mask_.begin() + static_cast<ptrdiff_t>(event.getPosition() + event.getLength()));
}

/******************************************************************************/

SequenceWithAnnotation* SequenceWithAnnotationTools::createMaskAnnotation(const Sequence& seq) throw (AlphabetException)
{
  const CaseMaskedAlphabet* cma = dynamic_cast<const CaseMaskedAlphabet*>(seq.getAlphabet());
  if (cma) {
    SequenceWithAnnotation* seqa = new SequenceWithAnnotation(seq.getName(), seq.toString(), seq.getComments(), seq.getAlphabet());
    vector<bool> mask(seq.size());
    for (unsigned int i = 0; i < seq.size(); ++i) {
      mask[i] = cma->isMasked(seq[i]);
    }
    seqa->addAnnotation(new SequenceMask(mask));
    return seqa;
  } else {
    throw AlphabetException("SequenceWithAnnotationTools::createMaskAnnotation. Alphabet should be a CaseMaskedAlphabet.", seq.getAlphabet());
  }
}

/******************************************************************************/

