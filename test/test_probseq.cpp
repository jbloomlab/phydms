//
// File: test_probseq.cpp
// Created by: Murray Patterson
// Created on: Fri Oct 9 2015
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

// from the STL
#include <iostream>
#include <string>
#include <vector>

// alphabets
#include <Bpp/Seq/Alphabet/BinaryAlphabet.h>
#include <Bpp/Seq/Alphabet/DNA.h>
#include <Bpp/Seq/Alphabet/AlphabetTools.h>

// symbol lists
#include <Bpp/Seq/SymbolList.h>
#include <Bpp/Seq/ProbabilisticSymbolList.h>
#include <Bpp/Numeric/DataTable.h>

// sequences
#include <Bpp/Seq/Sequence.h>
#include <Bpp/Seq/ProbabilisticSequence.h>

// sites
#include <Bpp/Seq/Site.h>
#include <Bpp/Seq/ProbabilisticSite.h>

// containers
#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <Bpp/Seq/Container/VectorProbabilisticSiteContainer.h>

// file formats
#include <Bpp/Seq/Io/Fasta.h>
#include <Bpp/Seq/Io/Pasta.h>

using namespace std;
using namespace bpp;

// convert a vector of string to a string
const string str(const vector<string> & v) {

  string s = "[";
  for(vector<string>::const_iterator i = v.begin(); i != v.end(); ++i)
    s += " " + *i;
  s += " ]";

  return s;
}

// basic info about an alphabet
void alphabet_info(const Alphabet * a) {

  cout << "alphabet is of size : " << a->getSize() << endl;
  cout << "supported chars : " << str(a->getSupportedChars()) << endl;
  cout << "resolved chars : " << str(a->getResolvedChars()) << endl;
}

// initialize empty (Probabilistic-) SymbolLists just to see
void init_empty_lists(const Alphabet * a) {
    BasicSymbolList list(a);
    BasicProbabilisticSymbolList p_list(a);
}
void init_empty_lists(const vector<Alphabet *> as) { // for several alphabets
  for(vector<Alphabet *>::const_iterator a = as.begin(); a != as.end(); ++a) { init_empty_lists(*a); }
}

/*
 * MAIN
 */
int main() {

  // initialize alphabets
  cout << endl << "init binary alphabet...";
  const BinaryAlphabet * a = new BinaryAlphabet();
  cout << "OK." << endl;
  alphabet_info(a);

  cout << endl << "init DNA alphabet...";
  const DNA * dna = new DNA();
  cout << "OK." << endl;
  alphabet_info(dna);

  // initialize empty (Probabilistic-) SymbolLists as a start
  cout << endl << "init (Probabilistic-) SymbolList with just alphabets...";
  init_empty_lists(a);
  init_empty_lists(dna);
  cout << "OK." << endl;

  /*
   * *** lists with binary content ***
   */
  cout << endl;
  cout << "***" << endl;
  cout << "*** lists with binary content ***" << endl;
  cout << "***" << endl;

  // *** the normal version ***
  string c[] = {"1","0","0"};
  vector<string> content(c,c+sizeof(c)/sizeof(c[0]));

  cout << endl << "init binary symbol list with : " << str(content) << " ...";
  BasicSymbolList list(content,a);
  cout << "OK." << endl;

  cout << "list contains : " << list.toString() << endl;

  // sequence
  cout << endl << "init binary sequence with : " << str(content) << " ...";
  BasicSequence seq("basic binary sequence",content,a);
  cout << "OK." << endl;

  // site
  cout << endl << "init binary site with : " << str(content) << " and position 42...";
  Site site(content,a,42);
  cout << "OK." << endl;

  cout << "site has position : " << site.getPosition() << endl;
  
  // *** the probabilistic version ***
  istringstream iss("0 1\n0.85 0.15\n0.99 0.01");
  DataTable * data = DataTable::read(iss, " ", false);

  cout << endl << "init probabilistic symbol list with : " << endl;
  DataTable::write(*data, cout);
  cout << "...";
  BasicProbabilisticSymbolList p_list(*data,a);
  cout << "OK." << endl;

  cout << "prob-list contains : ";
  DataTable::write(p_list.getContent(), cout);

  // sequence
  cout << endl << "init binary probabilistic sequence with its content...";
  BasicProbabilisticSequence p_seq("basic probabilistic binary sequence",*data,a);
  cout << "OK." << endl;

  // site
  cout << endl << "init binary probabilistic site with its content and position 3...";
  BasicProbabilisticSite p_site(*data,a,3);
  cout << "OK." << endl;

  cout << "site has position : " << p_site.getPosition() << endl;

  /*
   * *** lists with DNA content ***
   */
  cout << endl;
  cout << "***" << endl;
  cout << "*** lists with DNA content ***" << endl;
  cout << "***" << endl;

  // *** the normal version ***
  string cc[] = {"G", "T", "C"};
  vector<string> dna_content(cc,cc+sizeof(cc)/sizeof(cc[0]));

  cout << endl << "init DNA symbol list with : " << str(dna_content) << " ...";
  BasicSymbolList dna_list(dna_content,dna);
  cout << "OK." << endl;

  cout << "list contains : " << dna_list.toString() << endl;

  // sequence
  cout << endl << "init DNA sequence with : " << str(dna_content) << " ...";
  BasicSequence dna_seq("basic DNA sequence",dna_content,dna);
  cout << "OK." << endl;

  // site
  cout << endl << "init DNA site with : " << str(dna_content) << " and position 23...";
  Site dna_site(dna_content,dna,23);
  cout << "OK." << endl;

  cout << "site has position : " << dna_site.getPosition() << endl;

  // *** the probabilistic version ***
  istringstream isss("0 0 1 0\n0.05 0 0.05 0.9\n0.01 0.97 0 0.02");
  DataTable * dna_data = DataTable::read(isss, " ", false);

  cout << endl << "init probabilistic DNA symbol list with : " << endl;
  DataTable::write(*dna_data, cout);
  cout << "...";
  BasicProbabilisticSymbolList dna_p_list(*dna_data,dna);
  cout << "OK." << endl;

  cout << "probabilistic list contains : ";
  DataTable::write(dna_p_list.getContent(), cout);

  // sequence
  cout << endl << "init DNA probabilistic sequence with its content...";
  BasicProbabilisticSequence dna_p_seq("basic probabilistic binary sequence",*dna_data,dna);
  cout << "OK." << endl;

  // site
  cout << endl << "init DNA probabilistic site with its content...";
  BasicProbabilisticSite dna_p_site(*dna_data,dna);
  cout << "OK." << endl;

  cout << "site has position : " << dna_p_site.getPosition() << endl;

  /*
   * *** vector (probabilistic) site containers ***
   */
  cout << endl;
  cout << "***" << endl;
  cout << "*** vector site containers ***" << endl;
  cout << "***" << endl;

  // *** the normal version ***
  cout << endl << "init vector site containers with just alphabets...";
  VectorSiteContainer container(a);
  VectorSiteContainer dna_container(dna);
  cout << "OK." << endl;

  cout << endl << "add binary sequence " << seq.toString() << " to binary container...";
  container.addSequence(seq);
  cout << "OK." << endl;

  cout << "binary container's first element is : ";
  cout << container.getSequence(0).toString() << endl;

  // *** the probabilistic version ***
  cout << endl << "init vector probabilistic site containers with just alphabets...";
  VectorProbabilisticSiteContainer p_container(a);
  VectorProbabilisticSiteContainer dna_p_container(dna);
  cout << "OK." << endl;

  cout << endl << "add binary probabilistic sequence to binary probabilistic container...";
  p_container.addSequence(p_seq);
  cout << "OK." << endl;

  cout << "binary probabilistic container's first element is : ";  
  DataTable::write(p_container.getProbabilisticSequence(0).getContent(), cout);

  /*
   * *** Fasta (and Pasta) files ***
   */
  cout << endl;
  cout << "***" << endl;
  cout << "*** Fasta (and Pasta) files ***" << endl;
  cout << "***" << endl;

  // *** the normal version ***
  Fasta * fasta = new Fasta();
  cout << endl << "created a handler of type : " << fasta->getFormatName() << endl;

  string fasta_in = ">another binary sequence\n101\n";
  istringstream fasta_iss(fasta_in);
  cout << "read the following into binary container..." << endl;
  cout << endl << fasta_in << endl;
  fasta->appendSequencesFromStream(fasta_iss, container);
  cout << "OK." << endl;
  
  cout << "binary container contains : " << endl << endl;
  for(size_t i = 0; i < container.getNumberOfSequences(); ++i)
    cout << container.getSequence(i).toString() << endl;
  cout << endl;

  // DNA ...
  Fasta * dna_fasta = new Fasta();
  cout << endl << "created a handler of type : " << fasta->getFormatName() << endl;
  string dna_fasta_in = ">another dna sequence\nACG\n";
  istringstream dna_fasta_iss(dna_fasta_in);
  cout << "read the following into binary container..." << endl;
  cout << endl << dna_fasta_in << endl;
  dna_fasta->appendSequencesFromStream(dna_fasta_iss, dna_container);
  cout << "OK." << endl;

  cout << "DNA container contains : " << endl << endl;
  for(size_t i = 0; i < dna_container.getNumberOfSequences(); ++i)
    cout << dna_container.getSequence(i).toString() << endl;
  cout << endl;

  // *** the probabilistic version ***
  Pasta * pasta = new Pasta();
  cout << "created a handler of type : " << pasta->getFormatName() << endl;

  string pasta_in = "0 1\n>a binary probabilistic sequence\n0.64 0.36\n0 1\n0.3 0.7\n";
  istringstream pasta_iss(pasta_in);
  cout << "read the following into binary probabilistic container..." << endl;
  cout << endl << pasta_in << endl;
  pasta->appendSequencesFromStream(pasta_iss, p_container);
  cout << "OK." << endl;

  string pasta_in2 = ">another binary probabilistic sequence\n0.8 0.4 0.333\n";
  istringstream pasta_iss2(pasta_in2);
  cout << "read the following into binary probabilistic container (in fast-track way for binary alphabets)..." << endl;
  cout << endl << pasta_in2 << endl;
  pasta->appendSequencesFromStream(pasta_iss2, p_container);
  cout << "OK." << endl;

  cout << "binary probabilistic container contains : " << endl << endl;
  for(size_t i = 0; i < p_container.getNumberOfProbabilisticSequences(); ++i) {
    DataTable::write(p_container.getProbabilisticSequence(i).getContent(), cout);
    cout << endl;
  }

  // DNA...
  Pasta * dna_pasta = new Pasta();
  cout << "created a handler of type : " << pasta->getFormatName() << endl;

  string dna_pasta_in = "A C G T\n>a dna prob. sequence\n0.1834088 0.6140376 0.132227880 0.07032571\n0.4960896 0.0523049 0.123549944 0.32805560\n";
  istringstream dna_pasta_iss(dna_pasta_in);
  cout << "read the following into dna prob. container" << endl;
  cout << endl << dna_pasta_in << endl;
  dna_pasta->appendSequencesFromStream(dna_pasta_iss, dna_p_container);
  cout << "OK." << endl;

  string dna_pasta_in2 = "C T A G\n>another dna prob. sequence\n0.1885256 0.2023275 0.570924031 0.03822292\n0.1122945 0.2366416 0.004093129 0.64697079\n";
  istringstream dna_pasta_iss2(dna_pasta_in2);
  cout << "read the following (permuted) sequence into dna prob. container" << endl;
  cout << endl << dna_pasta_in2 << endl;
  dna_pasta->appendSequencesFromStream(dna_pasta_iss2, dna_p_container);
  cout << "OK." << endl;

  cout << "dna prob. container contains : " << endl << endl;
  for(size_t i = 0; i < dna_p_container.getNumberOfProbabilisticSequences(); ++i) {
    DataTable::write(dna_p_container.getProbabilisticSequence(i).getContent(), cout);
    cout << endl;
  }

  // the end
  return 0;
}
