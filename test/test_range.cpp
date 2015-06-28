//
// File: test_range.cpp
// Created by: Julien Dutheil
// Created on: Mon Nov 5 16:12 2011
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

#include <Bpp/Numeric/Range.h>
#include <Bpp/Numeric/VectorTools.h>
#include <iostream>

using namespace bpp;
using namespace std;

template<class T> void test_range(RangeCollection<T>& collection) {
  cout << collection.toString() << endl;
 
  Range<unsigned int> r(10, 13);
  cout << "Adding " << r.toString() << endl;
  collection.addRange(r);
  cout << collection.toString() << endl;

  r = Range<unsigned int>(18, 21);
  cout << "Adding " << r.toString() << endl;
  collection.addRange(r);
  cout << collection.toString() << endl;

  r = Range<unsigned int>(25, 23);
  cout << "Adding " << r.toString() << endl;
  collection.addRange(r);
  cout << collection.toString() << endl;

  r = Range<unsigned int>(7, 11);
  cout << "Adding " << r.toString() << endl;
  collection.addRange(r);
  cout << collection.toString() << endl;

  r = Range<unsigned int>(9, 24);
  cout << "Restricting to " << r.toString() << endl;
  collection.restrictTo(r);
  cout << collection.toString() << endl;

  r = Range<unsigned int>(5, 23);
  cout << "Restricting to " << r.toString() << endl;
  collection.restrictTo(r);
  cout << collection.toString() << endl;

  r = Range<unsigned int>(19, 24);
  cout << "Adding " << r.toString() << endl;
  collection.addRange(r);
  cout << collection.toString() << endl;

  r = Range<unsigned int>(0, 50);
  cout << "Adding " << r.toString() << endl;
  collection.addRange(r);
  cout << collection.toString() << endl;
}

int main() {
  
  cout << endl << "..:: Definitions ::.." << endl;
  Range<unsigned int> r1(3, 6);
  Range<unsigned int> r2(10, 12);
  Range<unsigned int> r3(0, 2);
  Range<unsigned int> r4(0, 5);
  Range<unsigned int> r5(4, 25);
  Range<unsigned int> r6(4, 5);
  Range<unsigned int> r7(1, 50);
  Range<unsigned int> r8(1, 3);
  Range<unsigned int> r9(6, 8);
  cout << "r1: " << r1.toString() << endl;
  cout << "r2: " << r2.toString() << endl;
  cout << "r3: " << r3.toString() << endl;
  cout << "r4: " << r4.toString() << endl;
  cout << "r5: " << r5.toString() << endl;
  cout << "r6: " << r6.toString() << endl;
  cout << "r7: " << r7.toString() << endl;
  cout << "r8: " << r8.toString() << endl;
  cout << "r9: " << r9.toString() << endl;
  
  cout << endl << "..:: Overlaps ::.." << endl;
  cout << r1.overlap(r2) << endl;
  if (r1.overlap(r2)) return 1;
  cout << r1.overlap(r3) << endl;
  if (r1.overlap(r3)) return 1;
  cout << r1.overlap(r4) << endl;
  if (!r1.overlap(r4)) return 1;
  cout << r1.overlap(r5) << endl;
  if (!r1.overlap(r5)) return 1;
  cout << r1.overlap(r6) << endl;
  if (!r1.overlap(r6)) return 1;
  cout << r1.overlap(r7) << endl;
  if (!r1.overlap(r7)) return 1;
  cout << r1.overlap(r8) << endl;
  if (r1.overlap(r8)) return 1;
  cout << r1.overlap(r9) << endl;
  if (r1.overlap(r9)) return 1;

  cout << endl << "..:: Contiguous ::.." << endl;
  cout << r1.isContiguous(r2) << endl;
  if (r1.isContiguous(r2)) return 1;
  cout << r1.isContiguous(r3) << endl;
  if (r1.isContiguous(r3)) return 1;
  cout << r1.isContiguous(r4) << endl;
  if (r1.isContiguous(r4)) return 1;
  cout << r1.isContiguous(r5) << endl;
  if (r1.isContiguous(r5)) return 1;
  cout << r1.isContiguous(r6) << endl;
  if (r1.isContiguous(r6)) return 1;
  cout << r1.isContiguous(r7) << endl;
  if (r1.isContiguous(r7)) return 1;
  cout << r1.isContiguous(r8) << endl;
  if (!r1.isContiguous(r8)) return 1;
  cout << r1.isContiguous(r9) << endl;
  if (!r1.isContiguous(r9)) return 1;

  Range<unsigned int> r;

  cout << endl << "..:: Expand ::.." << endl;
  r = r1; r.expandWith(r2); cout << "r1 \\/ r2: " << r.toString() << endl;
  if (r != r1) return 1;
  r = r1; r.expandWith(r3); cout << "r1 \\/ r3: " << r.toString() << endl;
  if (r != r1) return 1;
  r = r1; r.expandWith(r4); cout << "r1 \\/ r4: " << r.toString() << endl;
  if (r != Range<unsigned int>(0, 6)) return 1;
  r = r1; r.expandWith(r5); cout << "r1 \\/ r5: "  << r.toString() << endl;
  if (r != Range<unsigned int>(3, 25)) return 1;
  r = r1; r.expandWith(r6); cout << "r1 \\/ r6: "  << r.toString() << endl;
  if (r != r1) return 1;
  r = r1; r.expandWith(r7); cout << "r1 \\/ r7: "  << r.toString() << endl;
  if (r != r7) return 1;

  cout << endl << "..:: Slice ::.." << endl;
  r = r1; r.sliceWith(r2); cout << "r1 /\\ r2: " << r.toString() << endl;
  if (r != Range<unsigned int>(0, 0)) return 1;
  r = r1; r.sliceWith(r3); cout << "r1 /\\ r3: " << r.toString() << endl;
  if (r != Range<unsigned int>(0, 0)) return 1;
  r = r1; r.sliceWith(r4); cout << "r1 /\\ r4: " << r.toString() << endl;
  if (r != Range<unsigned int>(3, 5)) return 1;
  r = r1; r.sliceWith(r5); cout << "r1 /\\ r5: "  << r.toString() << endl;
  if (r != Range<unsigned int>(4, 6)) return 1;
  r = r1; r.sliceWith(r6); cout << "r1 /\\ r6: "  << r.toString() << endl;
  if (r != r6) return 1;
  r = r1; r.sliceWith(r7); cout << "r1 /\\ r7: "  << r.toString() << endl;
  if (r != r1) return 1;

  cout << endl << "..:: Operators ::.." << endl;
  r = r1; r += 1; cout << "r1 += 1: " << r.toString() << endl;
  if (r != Range<unsigned int>(4, 7)) return 1;
  r = r1 + 2; cout << "r1 + 2: " << r.toString() << endl;
  if (r != Range<unsigned int>(5, 8)) return 1;
  r = r1; r -= 1; cout << "r1 -= 1: " << r.toString() << endl;
  if (r != Range<unsigned int>(2, 5)) return 1;
  r = r1 - 2; cout << "r1 - 2: " << r.toString() << endl;
  if (r != Range<unsigned int>(1, 4)) return 1;

  cout << endl << "..:: MultiRange ::.." << endl;
  MultiRange<unsigned int> mr;
  test_range(mr);

  cout << endl << "..:: RangeSet ::.." << endl;
  RangeSet<unsigned int> rs;
  test_range(rs);

  return 0;
}
