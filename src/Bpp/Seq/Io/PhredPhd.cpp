//
// File: PhredPhd.cpp
// Created by: Sylvain Gaillard
// Created on: Wed Nov 5 2008
//

/*
Copyright or © or Copr. CNRS, (November 5, 2008)

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

#include "PhredPhd.h"
#include <Bpp/Text/TextTools.h>
#include <Bpp/Text/StringTokenizer.h>

using namespace bpp;

/******************************************************************************/

//PhredPhd::PhredPhd() {}

/******************************************************************************/

bool PhredPhd::nextSequence(std::istream& input, Sequence& seq) const throw (Exception) {
  std::vector<int> pos;
  return nextSequence(input, seq, pos);
}

/******************************************************************************/

bool PhredPhd::nextSequence(std::istream& input, Sequence& seq, std::vector<int>& pos) const throw (Exception) {
  if (!input) {
    throw IOException ("PhredPhd::read: fail to open stream");
  }

  bool flag = false;
  std::string name, sequence = "";  // Initialization
  std::vector<int> q, p;

  flag = parseFile_(input, name, sequence, q, p);
  // Sequence creation
  if(name == "")
    throw Exception("PhredPhd::read: sequence without name!");
  seq.setName(name);
  seq.setContent(sequence);
  try {
    SequenceWithQuality& sq = dynamic_cast<SequenceWithQuality&>(seq);
    sq.setQualities(q);
  } catch (...) {
  }
  return flag;
}

/******************************************************************************/

bool PhredPhd::parseFile_(std::istream& input, std::string& name, std::string& sequence, std::vector<int>& qual, std::vector<int>& pos) const {
  bool readSeqFlag = false;
  std::string temp;
  // Read sequence info
  // Main loop : for all lines
  while (!input.eof()) {
    std::getline(input, temp, '\n');  // Copy current line in temporary string
    StringTokenizer st(temp, " ");
    if (st.hasMoreToken()) {
      if (st.getToken(0) == "BEGIN_SEQUENCE") {
        name = st.getToken(1);
      }
      std::string flag = st.getToken(0);
      while (flag != "END_SEQUENCE" && !input.eof()) {
        getline(input, temp, '\n');
        StringTokenizer st2(temp, " ");
        if (st2.hasMoreToken()) {
          flag = st2.getToken(0);
        }
        if (flag == "BEGIN_DNA") {
          readSeqFlag = parseDNA_(input, sequence, qual, pos);
          break; // End the whole loop after parsing DNA
        }
      }
    }
  }
  return readSeqFlag;
}

/******************************************************************************/

bool PhredPhd::parseDNA_(std::istream& input, std::string& sequence, std::vector<int>& qual, std::vector<int>& pos) const {
  bool readSeqFlag = false;
  std::string line_buffer;
  std::string flag;
  sequence.clear();
  qual.clear();
  pos.clear();
  while (flag != "END_DNA" && !input.eof()) {
    std::getline(input, line_buffer, '\n');
    StringTokenizer st(line_buffer, " ");
    if (st.hasMoreToken()) {
      flag = TextTools::toUpper(st.getToken(0));
      if (st.numberOfRemainingTokens() == 3) {
        sequence += flag;
        qual.push_back(TextTools::toInt(st.getToken(1)));
        pos.push_back(TextTools::toInt(st.getToken(2)));
        readSeqFlag = true;
      }
    }
  }
  return readSeqFlag;
}

/******************************************************************************/
