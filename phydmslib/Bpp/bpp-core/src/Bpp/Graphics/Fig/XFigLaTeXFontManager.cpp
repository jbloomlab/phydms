//
// File: XFigLaTeXFontManager.cpp
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

#include "XFigLaTeXFontManager.h"

using namespace bpp;

XFigLaTeXFontManager::XFigLaTeXFontManager()
{
  // Add "official" font codes, from 0 to 5:
  registerFont_(Font("Default", Font::STYLE_NORMAL, Font::WEIGHT_NORMAL, 12), 0);
  registerFont_(Font("Roman", Font::STYLE_NORMAL, Font::WEIGHT_NORMAL, 12), 1);
  registerFont_(Font("Bold", Font::STYLE_NORMAL, Font::WEIGHT_NORMAL, 12), 2);
  registerFont_(Font("Italic", Font::STYLE_NORMAL, Font::WEIGHT_NORMAL, 12), 3);
  registerFont_(Font("Sans Serif", Font::STYLE_NORMAL, Font::WEIGHT_NORMAL, 12), 4);
  registerFont_(Font("Typewriter", Font::STYLE_NORMAL, Font::WEIGHT_NORMAL, 12), 5);
}

