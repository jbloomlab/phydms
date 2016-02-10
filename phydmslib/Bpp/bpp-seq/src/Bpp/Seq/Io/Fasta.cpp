//
// File: Fasta.cpp
// Authors: Guillaume Deuchst
//          Julien Dutheil
//          Sylvain Gaillard
// Created: Tue Aug 21 2003
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

#include "Fasta.h"

#include <fstream>

#include "../StringSequenceTools.h"
#include <Bpp/Text/TextTools.h>
#include <Bpp/Text/StringTokenizer.h>
#include <Bpp/Io/FileTools.h>

using namespace bpp;
using namespace std;

/******************************************************************************/

bool Fasta::nextSequence(istream& input, Sequence& seq) const throw (Exception) {
  if (!input)
    throw IOException("Fasta::nextSequence: can't read from istream input");
  string seqname = "";
  string content = "";
  Comments seqcmts;
  short seqcpt = 0;
  string linebuffer = "";
  char c;
  while (!input.eof())
  {
    c = static_cast<char>(input.peek());
    if (input.eof())
      c = '\n';

    // Sequence begining detection
    if (c == '>')
    {
      // Stop if find a new sequence
      if (seqcpt++)
        break;
    }
    getline(input, linebuffer);
    if (c == '>')
    {
      // Get the sequence name line
      seqname = string(linebuffer.begin() + 1, linebuffer.end());
    }
    if (c != '>' && !TextTools::isWhiteSpaceCharacter(c)) {
      // Sequence content
      content += TextTools::toUpper(TextTools::removeWhiteSpaces(linebuffer));
    }
  }

  bool res = (!input.eof());
  // Sequence name and comments isolation
  if (strictNames_ || extended_) {
    size_t pos = seqname.find_first_of(" \t\n");
    string seqcmt;
    if (pos != string::npos) {
      seqcmt = seqname.substr(pos + 1);
      seqname = seqname.substr(0, pos);
    }
    if (extended_) {
      StringTokenizer st(seqcmt, " \\", true, false);
      while (st.hasMoreToken()) {
        seqcmts.push_back(st.nextToken());
      }
    } else {
      seqcmts.push_back(seqcmt);
    }
    seq.setComments(seqcmts);
  }
  seq.setName(seqname);
  seq.setContent(content);
  return res;
}

/******************************************************************************/

void Fasta::writeSequence(ostream& output, const Sequence& seq) const throw (Exception)
{
  if (!output)
    throw IOException("Fasta::writeSequence: can't write to ostream output");
  // Sequence name
  output << ">" << seq.getName();
  // Sequence comments
  if (extended_)
  {
    for (unsigned int i = 0 ; i < seq.getComments().size() ; i++)
    {
      output << " \\" << seq.getComments()[i];
    }
  }
  output << endl;
  // Sequence content
  string buffer; // use a buffer to format sequence with states > 1 char
  for (size_t i = 0 ; i < seq.size() ; ++i)
  {
    buffer += seq.getChar(i);
    if (buffer.size() >= charsByLine_)
    {
      output << string(buffer.begin(), buffer.begin() + charsByLine_) << endl;
      buffer.erase(0, charsByLine_);
    }
  }
  output << string(buffer.begin(), buffer.end()) << endl;
}

/******************************************************************************/

void Fasta::appendSequencesFromStream(istream& input, SequenceContainer& vsc) const throw (Exception)
{
  if (!input)
    throw IOException("Fasta::appendFromStream: can't read from istream input");
  char c = '\n';
  char last_c;
  bool header = false;
  bool hasSeq = true;
  string line = "";
  Comments cmts;
  while (!input.eof() && hasSeq)
  {
    last_c = c;
    input.get(c);
    // Header detection
    if (extended_ && c == '#')
    {
      header = true;
      continue;
    }
    // Header end detection
    if (c == '\n')
    {
      if (extended_ && header)
      {
        if (line[0] == '\\')
        {
          line.erase(line.begin());
          cmts.push_back(line);
        }
        line = "";
        header = false;
      }
      continue;
    }
    // Header capture
    if (header)
    {
      line.append(1, c);
    }
    // Sequence detection
    if (c == '>' && last_c == '\n')
    {
      input.putback(c);
      c = last_c;
      BasicSequence tmpseq("", "", vsc.getAlphabet());
      hasSeq = nextSequence(input, tmpseq);
      vsc.addSequence(tmpseq, checkNames_);
    }
  }
  if (extended_ && cmts.size()) {
    vsc.setGeneralComments(cmts);
  }
}

/******************************************************************************/

void Fasta::writeSequences(ostream& output, const SequenceContainer& sc) const throw (Exception)
{
	if (!output)
    throw IOException("Fasta::write: can't write to ostream output");

  if (extended_)
  {
    // Loop for all general comments
    for (unsigned int i = 0 ; i < sc.getGeneralComments().size() ; i++)
    {
      output << "#\\" << sc.getGeneralComments()[i] << endl;
    }
    output << endl;
  }

	// Main loop : for all sequences in vector container
	vector<string> names = sc.getSequencesNames();
	for (size_t i = 0; i < names.size(); ++i)
  {
    writeSequence(output, sc.getSequence(names[i]));
	}
}

/******************************************************************************/

// FileIndex class

void Fasta::FileIndex::build(const std::string& path, const bool strictSequenceNames) throw (Exception) {
  // open the file
  std::ifstream f_in(path.c_str());
  // get the size of the file
  f_in.seekg(0, std::ios::end);
  fileSize_ = f_in.tellg();
  // feed the map
  f_in.seekg(0, std::ios::beg);
  streampos pos = f_in.tellg();
  char ch;
  std::string seq_id = "";
  while (f_in.get(ch)) {
    if (ch == '>') {
      pos = static_cast<int>(f_in.tellg()) - 1;
      std::getline(f_in, seq_id);
      if (strictSequenceNames) {
        seq_id = seq_id.substr(0, seq_id.find_first_of(" \t\n"));
      }
      index_[seq_id] = pos;
    }
  }
  f_in.close();
}

streampos Fasta::FileIndex::getSequencePosition(const std::string& id) const throw (Exception) {
  std::map<std::string, streampos>::const_iterator it = index_.find(id);
  if (it != index_.end()) {
    return it->second;
  }
  throw Exception("Sequence not found: " + id);
}

void Fasta::FileIndex::read(const std::string& path) throw (Exception) {
  std::ifstream f_in(path.c_str());
  std::string line_buffer = "";
  while (!f_in.eof()) {
    std::getline(f_in, line_buffer);
    if (bpp::TextTools::isEmpty(bpp::TextTools::removeSurroundingWhiteSpaces(line_buffer))) {
      continue;
    }
    bpp::StringTokenizer tk(line_buffer, "\t");
    index_[tk.getToken(0)] = bpp::TextTools::toInt(tk.getToken(1));
  }
  f_in.close();
}

void Fasta::FileIndex::write(const std::string& path) throw (Exception) {
  std::ofstream f_out(path.c_str());
  for (std::map<std::string, streampos>::const_iterator it = index_.begin() ; it != index_.end() ; ++it) {
    f_out << it->first << "\t" << bpp::TextTools::toString(it->second) << std::endl;
  }
  f_out.close();
}

void Fasta::FileIndex::getSequence(const std::string& seqid, Sequence& seq, const std::string& path) const {
  getSequence(seqid, seq, path, false);
}

void Fasta::FileIndex::getSequence(const std::string& seqid, Sequence& seq, const std::string& path, const bool strictSequenceNames) const {
  Fasta fs(60);
  fs.strictNames(strictSequenceNames);
  streampos seq_pos = this->getSequencePosition(seqid);
  std::ifstream fasta(path.c_str());
  fasta.seekg(seq_pos);
  fs.nextSequence(fasta, seq);
  fasta.close();
}

/******************************************************************************/

