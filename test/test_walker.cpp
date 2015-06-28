//
// File: test_walker.cpp
// Created by: Julien Dutheil
// Created on: Thu Nov 24 14:42 2011
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

#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Seq/Alphabet/RNA.h>
#include <Bpp/Seq/Sequence.h>
#include <Bpp/Seq/SequenceWalker.h>
#include <Bpp/Numeric/Random/RandomTools.h>
#include <iostream>

using namespace bpp;
using namespace std;

bool testSeq(SequenceWalker& walker, unsigned int pos, unsigned int truth) {
  cout << 0 << endl;
  cout << walker.getSequencePosition(0) << endl;
  cout << 46 << endl;
  cout << walker.getSequencePosition(46) << endl;;
  for (unsigned int i = 0; i < 1000; ++i) {
    ApplicationTools::displayGauge(i, 999, '=');
    size_t r = RandomTools::giveIntRandomNumberBetweenZeroAndEntry<size_t>(47);
    size_t x = walker.getSequencePosition(r);
    if (walker.getSequencePosition(pos) != truth) {
      cout << endl;
      cerr << r << "\t" << x << endl;
      cerr << walker.getSequencePosition(pos) << "<>" << truth << endl;
      return false;
    }
  }
  cout << endl;
  return true;
}

bool testAln(SequenceWalker& walker, unsigned int pos, unsigned int truth) {
  cout << 0 << endl;
  cout << walker.getAlignmentPosition(0) << endl;
  cout << 26 << endl;
  cout << walker.getAlignmentPosition(26) << endl;
  for (unsigned int i = 0; i < 1000; ++i) {
    ApplicationTools::displayGauge(i, 999, '=');
    size_t r = RandomTools::giveIntRandomNumberBetweenZeroAndEntry<size_t>(27);
    walker.getAlignmentPosition(r);
    if (walker.getAlignmentPosition(pos) != truth) {
      cout << endl;
      cerr << walker.getSequencePosition(pos) << "<>" << truth << endl;
      return false;
    }
  }
  cout << endl;
  return true;
}

int main() {
  RNA* alpha = new RNA();
  BasicSequence seq1("seq1", "----AUGCCG---GCGU----UUU----G--G-CCGACGUGUUUU--", alpha);
  SequenceWalker walker(seq1);

  for (unsigned int i = 0; i < 27; ++i) {
    size_t j = walker.getAlignmentPosition(i);
    cout << i << "\t" << seq1.getChar(j) << "\t" << j << endl;
  }

  cout << endl;

  if (!testAln(walker, 5, 9)) return 1;
  if (!testAln(walker, 10, 21)) return 1;
  if (!testAln(walker, 22, 40)) return 1;

  cout << "_________________________________________________" << endl;

  for (unsigned int i = 0; i < seq1.size(); ++i) {
    cout << i << "\t" << seq1.getChar(i) << "\t" << walker.getSequencePosition(i) << endl;
  }

 
  cout << endl;
 
  if (!testSeq(walker, 9, 5)) return 1;
  if (!testSeq(walker, 21, 10)) return 1;
  if (!testSeq(walker, 40, 22)) return 1;

  return 0;
}
