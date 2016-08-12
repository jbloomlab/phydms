//
// File: test_bowker.cpp
// Created by: Julien Dutheil
// Created on: Sat Apr 16 13:19 2009
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

#include <Bpp/Seq/Alphabet/RNA.h>
#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Seq/SequenceTools.h>
#include <iostream>

using namespace bpp;
using namespace std;

BasicSequence* getRandomSequence(const Alphabet* alphabet, unsigned int size) {
  string seq = "";
  for (unsigned int i = 0; i < size; ++i)
    seq += alphabet->intToChar(RandomTools::giveIntRandomNumberBetweenZeroAndEntry(static_cast<int>(alphabet->getSize())));
  return new BasicSequence("random seq", seq, alphabet);
}

int test(const Alphabet* alphabet, unsigned int size, unsigned int rep) {
  double n01 = 0;
  double n05 = 0;
  double n10 = 0;
  double n50 = 0;

  //ofstream out("pvalues.txt", ios::out);
  for (unsigned int i = 0; i < rep; ++i) {
    ApplicationTools::displayGauge(i, rep-1);
    unique_ptr<BasicSequence> seq1(getRandomSequence(alphabet, size));
    unique_ptr<BasicSequence> seq2(getRandomSequence(alphabet, size));
    unique_ptr<BowkerTest> test(SequenceTools::bowkerTest(*seq1, *seq2));
    double p = test->getPValue();
    if (p <= 0.01) n01++;
    if (p <= 0.05) n05++;
    if (p <= 0.10) n10++;
    if (p <= 0.50) n50++;
    //out << p << endl;
  }
  //out.close();

  cout << n01 / rep <<
  "\t" << n05 / rep <<
  "\t" << n10 / rep <<
  "\t" << n50 / rep << endl;

  if (abs(n01*100 / rep - 1 ) > 1) return 1;
  if (abs(n05*100 / rep - 5 ) > 1) return 1;
  if (abs(n10*100 / rep - 10) > 1) return 1;
  if (abs(n50*100 / rep - 50) > 1) return 1;
  return 0;
}

int main() {
  //ProteicAlphabet* alpha = new ProteicAlphabet;
  RNA* alpha = new RNA();
  BasicSequence seq1("seq1", "----AUGCCG---GCGU----UUU----G--G-CCGACGUGUUUU--", alpha);
  BasicSequence seq2("seq2", "---GAAGGCG---G-GU----UUU----GC-GACCGACG--UUUU--", alpha);
  unique_ptr<BowkerTest> btest(SequenceTools::bowkerTest(seq1, seq2));
  cout << btest->getStatistic() << "\t" << btest->getPValue() << endl;
  delete alpha;

  test(&AlphabetTools::DNA_ALPHABET, 1000, 5000);
  test(&AlphabetTools::PROTEIN_ALPHABET, 20000, 5000);
 
  return 0;
}
