//
// File: Pasta.cpp
// Authors: Murray Patterson
// Created: Tue Oct 20 2015
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

#include "Pasta.h"

#include <fstream>

#include "../StringSequenceTools.h"
#include <Bpp/Text/TextTools.h>
#include <Bpp/Text/StringTokenizer.h>
#include <Bpp/Io/FileTools.h>
#include <Bpp/Numeric/DataTable.h>

using namespace bpp;
using namespace std;

/********************************************************************************/

bool Pasta::nextSequence(istream & input, ProbabilisticSequence & seq) const throw (Exception) {

  if(!input)
    throw IOException("Pasta::nextSequence : can't read from istream input");

  string seqname = "";
  vector<string> tokens;
  Comments seqcmts;
  short seqcpt = 0;
  string linebuffer = "";
  char c;

  while(!input.eof()) {

    c = static_cast<char>(input.peek());
    if(input.eof())
      c = '\n';

    // detect the beginning of a sequence
    if(c == '>') {

      // stop if we find a new sequence
      if(seqcpt++)
	break;
    }

    getline(input, linebuffer);

    if(c == '>') {

      // get the sequence name line
      seqname = string(linebuffer.begin() + 1, linebuffer.end());
    }

    /*
      While dealing with all the formatting, names, etc., including
      strictNames and extended formats will all be the same as with
      the Fasta format, this is where it will be radically different :
      the handling of the sequence content

      For the moment, I'm keeping it simple : assuming that all of the
      'sequence' has a 'width' of 1, i.e., one can read it line by
      line, like we do below, that is, that it pertains to a binary
      alphabet (and moreover, the probability that the character is
      1).  However, one can make this much more general -- for general
      alphabets, where the 'width' would be the size S of the
      alphabet, i.e., we would read in S lines at a time -- something
      I may come back and do someday

      Note also, that the goal here is to feed this 'sequence' into a
      ProbabilisticSequence object, which contains a DataTable, where
      the intension was to read column-wise 'sequences' with
      DataTable's read function, yet I'm now feeding this row-wise
      sequence into (a column of the DataTable of)
      ProbabilisticSequence.  This is because I like the Pasta format
      as is, but in the future, another generalization could be to
      really have this column-wise type of Pasta format ... we'll
      cross that river when we get to it I suppose

      idea : a format that generalizes this one is
      p(0);p(1);p(2);...;p(s-1) where p(i) is the probability that
      site has state i, and s is the number of states (p(s) is not
      needed, because p(s) = 1 - (p(0)+...+p(s-1)) which is checked
      for in probabilistic symbol list.  Since binary alphabet has 2
      states, the format we handle here is exactly this ... something
      to think about for the future
    */

    if(c != '>' && !TextTools::isWhiteSpaceCharacter(c)) {

      // sequence content : probabilities for each site are space-separated
      StringTokenizer st(linebuffer, " \t\n", false, false);
      while(st.hasMoreToken()) {

	string s = st.nextToken();
	tokens.push_back(s);
      }
    }
  }

  bool res = (!input.eof());

  // Sequence name and comments isolation (identical to that of Fasta)
  if(strictNames_ || extended_) {

    size_t pos = seqname.find_first_of(" \t\n");
    string seqcmt;

    if(pos != string::npos) {
      seqcmt = seqname.substr(pos + 1);
      seqname = seqname.substr(0, pos);
    }

    if(extended_) {
      StringTokenizer st(seqcmt, " \\", true, false);
      while(st.hasMoreToken()) {
	seqcmts.push_back(st.nextToken());
      }
    }
    else {
      seqcmts.push_back(seqcmt);
    }
    seq.setComments(seqcmts);
  }

  seq.setName(seqname);

  // finally, add pairs of probabilities that (binary) character is 0, resp. 1

  // setup the names first
  string names[] = {"0","1"};
  vector<string> names_v(names,names+2);
  DataTable content(names_v); // a data table with 2 columns for p(0), p(1)

  // now add the rows
  for(vector<string>::const_iterator i = tokens.begin(); i != tokens.end(); ++i) {

    // the following will throw an exception if v[i] is not a properly
    // formatted double : a check that we want to have
    int precision = 10; // I think this will be more than enough
    string pair[] = {TextTools::toString(double(1) - TextTools::toDouble(*i,'.','E'), precision), *i};
    vector<string> pair_v(pair,pair+2);
    content.addRow(pair_v); // p(0), p(1)
  }

  // finally, we set the content of the sequence to the above.  Since
  // any number format exception was taken care of above, as well as
  // that each sequence must be of the same length by construction of
  // a (named) DataTable object, the only thing left that could go
  // wrong is that p(0) + p(1) != 1 : a check that is done in the call
  // of the function below
  seq.setContent(content);

  return res;
}

/********************************************************************************/

void Pasta::writeSequence(ostream & output, const ProbabilisticSequence & seq) const throw (Exception) {

  if(!output)
    throw IOException("Pasta::writeSequence : can't write to ostream output");

  // sequence name
  output << ">" << seq.getName();

  // sequence comments
  if(extended_) {
    for(unsigned int i = 0; i < seq.getComments().size(); ++i) {
      output << " \\" << seq.getComments()[i];
    }
  }
  output << endl;

  /*
    Again, here is where we radically diverge from the Fasta way : in
    the sequence content handling.  For now we assume that the
    DataTable (within the ProbabilisticSequence) contains exactly one
    column, pertaining to the probability that the binary character is
    1.  We output this column as a single (pasta-style) row

    again, extending this to other alphabets is future work
  */

  // sequence content
  vector<string> column = seq.getContent().getColumn("1");

  vector<string>::iterator i = column.begin();
  if(i != column.end()) // output the first element
    output << *i;
  for( ; i != column.end(); ++i) // output the rest preceeded by a space
    output << " " + *i;
  output << endl;
}

/********************************************************************************/

// note: because we did the quick hack of making
// VectorProbabilisticSiteContainer a modified copy of
// VectorSiteContainer, it, nor any object from which it inherits can
// take in a ProbabilisticSequence, so we can only read into
// VectorProbabilisticSiteContainers ... all we'll need for working
// with bppAncestor and bppMl, etc.  But a more general solution could
// be future work
void Pasta::appendSequencesFromStream(istream & input, VectorProbabilisticSiteContainer & container) const throw (Exception) {

  if(!input)
    throw IOException("Pasta::appendFromStream: can't read from istream input");

  char c = '\n';
  char last_c;
  bool header = false;
  bool hasSeq = true;
  string line = "";
  Comments cmts;

  while(!input.eof() && hasSeq) {

    last_c = c;
    input.get(c);

    // detect the header
    if(extended_ && c == '#') {
      header = true;
      continue;
    }

    // detect the end of the header
    if(c == '\n') {
      if(extended_ && header) {
	if(line[0] == '\\') {

	  line.erase(line.begin());
	  cmts.push_back(line);
	}
	line = "";
	header = false;
      }
      continue;
    }

    // capture the header
    if(header) {
      line.append(1,c);
    }

    // detect the sequence
    if(c == '>' && last_c == '\n') {

      input.putback(c);
      c = last_c;
      BasicProbabilisticSequence tmpseq(container.getAlphabet()); // add probabilistic sequences instead
      hasSeq = nextSequence(input, tmpseq);
      container.addSequence(tmpseq, checkNames_);
    }
  }
  if(extended_ && cmts.size()) {
    container.setGeneralComments(cmts);
  }
}
