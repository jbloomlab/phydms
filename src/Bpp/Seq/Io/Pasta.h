//
// File: Pasta.h
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

#ifndef _BPP_SEQ_IO_PASTA_H_
#define _BPP_SEQ_IO_PASTA_H_

#include "Fasta.h"

#include "../ProbabilisticSequence.h"
#include "../Container/VectorProbabilisticSiteContainer.h"

namespace bpp
{

/**
 * @brief The Pasta sequence file format.
 *
 * Read and write from/to Pasta files -- a format that is more general
 * than the Fasta format : while the Fasta format contains sequence
 * information in the form of character states at each site of the
 * sequence, the Pasta format contains sequence information in the
 * form of probability of presence of each character state at each
 * site.  See implementation of methods below for more details
 */
class Pasta :
  public Fasta
{

 public :

  /**
   * @brief Build a new Pasta object.
   *
   * @param charsByLine Number of characters per line when writing files.
   * @param checkSequenceNames Tell if the names in the file should be checked for unicity (slower, in o(n*n) where n is the number of sequences).
   * @param extended Tell if we should read general comments and sequence comments in HUPO-PSI format.
   * @param strictSequenceNames Tells if the sequence names should be restricted to the characters between '>' and the first blank one.
   */
  Pasta(unsigned int charsByLine = 100, bool checkSequenceNames = true, bool extended = false, bool strictSequenceNames = false) : Fasta(charsByLine, checkSequenceNames, extended, strictSequenceNames) {}

  // class destructor
  virtual ~Pasta() {}

 public :

  /**
   * @brief Get the format name
   *
   * @return format name
   */
  const std::string getFormatName() const { return "PASTA file"; }

  /**
   * @name The ISequenceStream interface
   *
   * @{
   */
  bool nextSequence(std::istream & input, ProbabilisticSequence & seq) const throw (Exception);
  /**
   * @}
   */

  /**
   * @name The OSequenceStream interface
   *
   * @{
   */
  void writeSequence(std::ostream & output, const ProbabilisticSequence & seq) const throw (Exception);
  /**
   * @}
   */

  /**
   * @name The AbstractISequence interface
   *
   */
  void appendSequencesFromStream(std::istream & input, VectorProbabilisticSiteContainer & container) const throw (Exception);
  /**
   * @}
   */

};

} // end of namespace bpp

#endif // _BPP_SEQ_IO_PASTA_H_
