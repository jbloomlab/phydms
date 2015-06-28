//
// File: XFigPostscriptFontManager.cpp
// Created by: Julien Dutheil
// Created on: Wed Dec 30 2009
// From file: FontManager.h
//

/*
Copyright or © or Copr. CNRS, (November 17, 2004)

This software is a computer program whose purpose is to provide utilitary
classes. This file belongs to the Bio++ Project.

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

#include "XFigPostscriptFontManager.h"

using namespace bpp;

XFigPostscriptFontManager::XFigPostscriptFontManager()
{
  // Add "official" font codes, from 0 to 34:
  registerFont_(Font("Default", Font::STYLE_NORMAL, Font::WEIGHT_NORMAL, 12), -1);
  registerFont_(Font("Times", Font::STYLE_NORMAL, Font::WEIGHT_NORMAL, 12), 0); //Roman
  registerFont_(Font("Times", Font::STYLE_ITALIC, Font::WEIGHT_NORMAL, 12), 1);
  registerFont_(Font("Times", Font::STYLE_NORMAL, Font::WEIGHT_BOLD, 12), 2);
  registerFont_(Font("Times", Font::STYLE_ITALIC, Font::WEIGHT_BOLD, 12), 3);
  registerFont_(Font("AvantGarde", Font::STYLE_NORMAL, Font::WEIGHT_NORMAL, 12), 4); //Book
  registerFont_(Font("AvantGarde", Font::STYLE_ITALIC, Font::WEIGHT_NORMAL, 12), 5); //Book Oblique
  registerFont_(Font("AvantGarde", Font::STYLE_NORMAL, Font::WEIGHT_BOLD, 12), 6); //Demi
  registerFont_(Font("AvantGarde", Font::STYLE_ITALIC, Font::WEIGHT_BOLD, 12), 7); //Demi Oblique
  registerFont_(Font("Bookman", Font::STYLE_NORMAL, Font::WEIGHT_NORMAL, 12), 8); //Light
  registerFont_(Font("Bookman", Font::STYLE_ITALIC, Font::WEIGHT_NORMAL, 12), 9); //Light Italic
  registerFont_(Font("Bookman", Font::STYLE_NORMAL, Font::WEIGHT_BOLD, 12), 10); //Demi
  registerFont_(Font("Bookman", Font::STYLE_ITALIC, Font::WEIGHT_BOLD, 12), 11); //Demi Italic
  registerFont_(Font("Courier", Font::STYLE_NORMAL, Font::WEIGHT_NORMAL, 12), 12);
  registerFont_(Font("Courier", Font::STYLE_ITALIC, Font::WEIGHT_NORMAL, 12), 13); //Oblique
  registerFont_(Font("Courier", Font::STYLE_NORMAL, Font::WEIGHT_BOLD, 12), 14); //Bold
  registerFont_(Font("Courier", Font::STYLE_ITALIC, Font::WEIGHT_BOLD, 12), 15); //Bold Oblique
  registerFont_(Font("Helvetica", Font::STYLE_NORMAL, Font::WEIGHT_NORMAL, 12), 16);
  registerFont_(Font("Helvetica", Font::STYLE_ITALIC, Font::WEIGHT_NORMAL, 12), 17); //Oblique
  registerFont_(Font("Helvetica", Font::STYLE_NORMAL, Font::WEIGHT_BOLD, 12), 18); //Bold
  registerFont_(Font("Helvetica", Font::STYLE_ITALIC, Font::WEIGHT_BOLD, 12), 19); //Bold Oblique
  registerFont_(Font("Helvetica", Font::STYLE_NORMAL, Font::WEIGHT_NORMAL, 12), 20); //Narrow
  registerFont_(Font("Helvetica", Font::STYLE_ITALIC, Font::WEIGHT_NORMAL, 12), 21); //Narrow Oblique
  registerFont_(Font("Helvetica", Font::STYLE_NORMAL, Font::WEIGHT_BOLD, 12), 22); //Narrow Bold
  registerFont_(Font("Helvetica", Font::STYLE_ITALIC, Font::WEIGHT_BOLD, 12), 23); //Narrow Bold Oblique
  registerFont_(Font("New Century Schoolbook", Font::STYLE_NORMAL, Font::WEIGHT_NORMAL, 12), 24); //Roman
  registerFont_(Font("New Century Schoolbook", Font::STYLE_ITALIC, Font::WEIGHT_NORMAL, 12), 25); //Italic
  registerFont_(Font("New Century Schoolbook", Font::STYLE_NORMAL, Font::WEIGHT_BOLD, 12), 26); //Bold
  registerFont_(Font("New Century Schoolbook", Font::STYLE_ITALIC, Font::WEIGHT_BOLD, 12), 27); //Bold Italic
  registerFont_(Font("Palatino", Font::STYLE_NORMAL, Font::WEIGHT_NORMAL, 12), 28); //Roman
  registerFont_(Font("Palatino", Font::STYLE_ITALIC, Font::WEIGHT_NORMAL, 12), 29); //Italic
  registerFont_(Font("Palatino", Font::STYLE_NORMAL, Font::WEIGHT_BOLD, 12), 30); //Bold
  registerFont_(Font("Palatino", Font::STYLE_ITALIC, Font::WEIGHT_BOLD, 12), 31); //Bold Italic
  registerFont_(Font("Symbol", Font::STYLE_NORMAL, Font::WEIGHT_NORMAL, 12), 32);
  registerFont_(Font("Zapf Chancery Medium", Font::STYLE_ITALIC, Font::WEIGHT_NORMAL, 12), 33); //Italic
  registerFont_(Font("Zapf Dingbats", Font::STYLE_NORMAL, Font::WEIGHT_NORMAL, 12), 34);
}

