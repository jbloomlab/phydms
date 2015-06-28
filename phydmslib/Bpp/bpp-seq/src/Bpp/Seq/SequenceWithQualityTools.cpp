//
// File:       SequenceWithQualityTools.h
// Authors:    Vincent Cahais
//             Sylvain Gaillard
// Created on: 16 Apr 2010
//

/*
Copyright or © or Copr. Bio++ Development Team, (Apr 16, 2010)

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

#include "SequenceWithQualityTools.h"

using namespace bpp;
using namespace std;

DNA SequenceWithQualityTools::DNA_;
RNA SequenceWithQualityTools::RNA_;
NucleicAcidsReplication SequenceWithQualityTools::DNARep_(& DNA_, & DNA_);
NucleicAcidsReplication SequenceWithQualityTools::RNARep_(& RNA_, & RNA_);
NucleicAcidsReplication SequenceWithQualityTools::transc_(& DNA_, & RNA_);

/******************************************************************************/

SequenceWithQuality* SequenceWithQualityTools::concatenate(const SequenceWithQuality& seqwq1, const SequenceWithQuality& seqwq2) throw (AlphabetMismatchException, Exception)
{
	// Sequence's alphabets matching verification
	if ((seqwq1.getAlphabet()->getAlphabetType()) != (seqwq2.getAlphabet()->getAlphabetType()))
		throw AlphabetMismatchException("SequenceTools::concatenate : Sequence's alphabets don't match ", seqwq1.getAlphabet(), seqwq2.getAlphabet());

	// Sequence's names matching verification
	if (seqwq1.getName() != seqwq2.getName())
    throw Exception ("SequenceTools::concatenate : Sequence's names don't match");

	// Concatenate sequences and send result
	vector<int> sequence  = seqwq1.getContent();
	vector<int> sequence2 = seqwq2.getContent();
	vector<int> qualities = seqwq1.getQualities();
	vector<int> qualities2 = seqwq2.getQualities();

	sequence.insert(sequence.end(), sequence2.begin(), sequence2.end());
	qualities.insert(qualities.end(), qualities2.begin(), qualities2.end());
	return new SequenceWithQuality(seqwq1.getName(), sequence, qualities, seqwq1.getComments(), seqwq1.getAlphabet());
}

/******************************************************************************/

SequenceWithQuality* SequenceWithQualityTools::complement(const SequenceWithQuality& sequence) throw (AlphabetException)
{
  // Alphabet type checking
  NucleicAcidsReplication* NAR;
  if (sequence.getAlphabet()->getAlphabetType() == "DNA alphabet")
  {
    NAR = &DNARep_;
  }
  else if(sequence.getAlphabet()->getAlphabetType() == "RNA alphabet")
  {
    NAR = &RNARep_;
  }
  else
  {
    throw AlphabetException ("SequenceTools::complement : Sequence must be nucleic.", sequence.getAlphabet());
  }
  Sequence* seq = NAR->translate(sequence);
  SequenceWithQuality* seqwq = new SequenceWithQuality(*seq, sequence.getQualities());
  delete seq;
  return seqwq;
}

/******************************************************************************/

SequenceWithQuality* SequenceWithQualityTools::transcript(const SequenceWithQuality& sequence) throw (AlphabetException)
{
  // Alphabet type checking
  if (sequence.getAlphabet()->getAlphabetType() != "DNA alphabet")
  {
    throw AlphabetException ("SequenceTools::transcript : Sequence must be DNA", sequence.getAlphabet());
  }
  Sequence* seq = transc_.translate(sequence);
  SequenceWithQuality* seqwq = new SequenceWithQuality(*seq, sequence.getQualities());
  delete seq;
  return seqwq;
}

/******************************************************************************/

SequenceWithQuality* SequenceWithQualityTools::reverseTranscript(const SequenceWithQuality& sequence) throw (AlphabetException)
{
  // Alphabet type checking
  if (sequence.getAlphabet()->getAlphabetType() != "RNA alphabet")
  {
    throw AlphabetException ("SequenceTools::reverseTranscript : Sequence must be RNA", sequence.getAlphabet());
  }

  Sequence* seq = transc_.reverse(sequence);
  //Here we must also reverse the scores:
  vector<int> scores(sequence.getQualities().rbegin(), sequence.getQualities().rend());
  SequenceWithQuality* seqwq = new SequenceWithQuality(*seq, scores);
  delete seq;
  return seqwq;
}

/******************************************************************************/

SequenceWithQuality* SequenceWithQualityTools::invert(const SequenceWithQuality& sequence)
{
  vector<int> iContent(sequence.getContent().rbegin(),sequence.getContent().rend());
  vector<int> iQualities(sequence.getQualities().rbegin(),sequence.getQualities().rend());
  SequenceWithQuality* iSeq = sequence.clone();
  iSeq->setContent(iContent);
  iSeq->setQualities(iQualities);

  return iSeq;
}

/******************************************************************************/

SequenceWithQuality* SequenceWithQualityTools::removeGaps(const SequenceWithQuality& seq)
{
  vector<int> content;
  vector<int> qualities;
  const Alphabet * alpha = seq.getAlphabet();
  for(unsigned int i = 0; i < seq.size(); i++)
  {
    if(! alpha->isGap(seq[i]))
    {
    	content.push_back(seq[i]);
    	qualities.push_back(seq.getQualities()[i]);
    }
  }
  SequenceWithQuality * newSeq = dynamic_cast<SequenceWithQuality *>(seq.clone());
  newSeq->setContent(content);
  newSeq->setQualities(qualities);
  return newSeq;
}

/******************************************************************************/

SequenceWithQuality& SequenceWithQualityTools::trimLeft(SequenceWithQuality& seq) {
  bool badqual = false;
  while (badqual) {
  }
  return seq;
}
