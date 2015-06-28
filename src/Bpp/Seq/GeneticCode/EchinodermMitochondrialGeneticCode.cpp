//
// File: EchinodermMitochondrialGeneticCode.cpp
// Created by: Eric Bazin
// Created on: 14 11:31:27 CET 2005

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

#include "EchinodermMitochondrialGeneticCode.h"

using namespace bpp;

#include <iostream>

using namespace std;

EchinodermMitochondrialGeneticCode::EchinodermMitochondrialGeneticCode(const NucleicAlphabet* alphabet) :
  GeneticCode(alphabet) 
{
  tlnTable_[0] = 2; //AAA -> N
  tlnTable_[1] = 2; //AAC -> N
  tlnTable_[2] = 11; //AAG -> K
  tlnTable_[3] = 2; //AAT -> N
  tlnTable_[4] = 16; //ACA -> T
  tlnTable_[5] = 16; //ACC -> T
  tlnTable_[6] = 16; //ACG -> T
  tlnTable_[7] = 16; //ACT -> T
  tlnTable_[8] = 15; //AGA -> S
  tlnTable_[9] = 15; //AGC -> S
  tlnTable_[10] = 15; //AGG -> S
  tlnTable_[11] = 15; //AGT -> S
  tlnTable_[12] = 9; //ATA -> I
  tlnTable_[13] = 9; //ATC -> I
  tlnTable_[14] = 12; //ATG -> M
  tlnTable_[15] = 9; //ATT -> I
  tlnTable_[16] = 5; //CAA -> Q
  tlnTable_[17] = 8; //CAC -> H
  tlnTable_[18] = 5; //CAG -> Q
  tlnTable_[19] = 8; //CAT -> H
  tlnTable_[20] = 14; //CCA -> P
  tlnTable_[21] = 14; //CCC -> P
  tlnTable_[22] = 14; //CCG -> P
  tlnTable_[23] = 14; //CCT -> P
  tlnTable_[24] = 1; //CGA -> R
  tlnTable_[25] = 1; //CGC -> R
  tlnTable_[26] = 1; //CGG -> R
  tlnTable_[27] = 1; //CGT -> R
  tlnTable_[28] = 10; //CTA -> L
  tlnTable_[29] = 10; //CTC -> L
  tlnTable_[30] = 10; //CTG -> L
  tlnTable_[31] = 10; //CTT -> L
  tlnTable_[32] = 6; //GAA -> E
  tlnTable_[33] = 3; //GAC -> D
  tlnTable_[34] = 6; //GAG -> E
  tlnTable_[35] = 3; //GAT -> D
  tlnTable_[36] = 0; //GCA -> A
  tlnTable_[37] = 0; //GCC -> A
  tlnTable_[38] = 0; //GCG -> A
  tlnTable_[39] = 0; //GCT -> A
  tlnTable_[40] = 7; //GGA -> G
  tlnTable_[41] = 7; //GGC -> G
  tlnTable_[42] = 7; //GGG -> G
  tlnTable_[43] = 7; //GGT -> G
  tlnTable_[44] = 19; //GTA -> V
  tlnTable_[45] = 19; //GTC -> V
  tlnTable_[46] = 19; //GTG -> V
  tlnTable_[47] = 19; //GTT -> V
  tlnTable_[48] = -99; //TAA -> STOP
  tlnTable_[49] = 18; //TAC -> Y
  tlnTable_[50] = -99; //TAG -> STOP
  tlnTable_[51] = 18; //TAT -> Y
  tlnTable_[52] = 15; //TCA -> S
  tlnTable_[53] = 15; //TCC -> S
  tlnTable_[54] = 15; //TCG -> S
  tlnTable_[55] = 15; //TCT -> S
  tlnTable_[56] = 17; //TGA -> W
  tlnTable_[57] = 4; //TGC -> C
  tlnTable_[58] = 17; //TGG -> W
  tlnTable_[59] = 4; //TGT -> C
  tlnTable_[60] = 10; //TTA -> L
  tlnTable_[61] = 13; //TTC -> F
  tlnTable_[62] = 10; //TTG -> L
  tlnTable_[63] = 13; //TTT -> F
  tlnTable_[codonAlphabet_.getUnknownCharacterCode()] = proteicAlphabet_.getUnknownCharacterCode();
}

