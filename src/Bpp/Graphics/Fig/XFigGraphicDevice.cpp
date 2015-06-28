//
// File: XFigGraphicDevice.cpp
// Created by: Julien Dutheil
// Created on: Mon Mar 03 2008
//

/*
Copyright or © or Copr. Bio++ Development Team, (November 17, 2004)

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

#include "XFigGraphicDevice.h"

using namespace bpp;
using namespace std;

const unsigned int XFigGraphicDevice::FONTFLAG_LATEX = 0;
const unsigned int XFigGraphicDevice::FONTFLAG_POSTSCRIPT = 4;

void XFigGraphicDevice::begin()
{
  content_.clear();
}

void XFigGraphicDevice::end()
{
  //Print to file:
  //Header:
  out_ << "#FIG 3.2 Produced by the Bio++ Graphic Device System" << endl;
  out_ << "Portrait" << endl;
  out_ << "Flush left" << endl;
  out_ << "Metric" << endl;
  out_ << "A4" << endl;
  out_ << "100" << endl;
  out_ << "Single" << endl;
  out_ << "0" << endl;
  //out << "254 2" << endl; // 1fig unit = 0.1mm
  out_ << "72 2" << endl; // 1fig unit = 1pt

  // Color definitions:
  out_ << "#Color definitions:" << endl; 
  vector<unsigned int> codes = colorManager_.getCodes();
  vector<RGBColor> colors    = colorManager_.getColors();
  for(unsigned int i = 32; i < colorManager_.getNumberOfColors(); i++)
  {
    string hexCode = colors[i].toHex();
    out_ << "0 " << codes[i] << " " << hexCode << endl;
  }

  // Print content:
  out_ << "#Drawing content:" << endl; 
  for(unsigned int i = 0; i < content_.size(); i++)
  {
    out_ << content_[i] << endl;
  }
}

void XFigGraphicDevice::setCurrentForegroundColor(const RGBColor& color)
{
  fgColorCode_ = colorManager_.getCode(color);
  AbstractGraphicDevice::setCurrentForegroundColor(color);
}

void XFigGraphicDevice::setCurrentBackgroundColor(const RGBColor& color)
{
  bgColorCode_ = colorManager_.getCode(color);
  AbstractGraphicDevice::setCurrentBackgroundColor(color);
}

void XFigGraphicDevice::setCurrentFont(const Font& font)
{
  if (fontFlag_ == FONTFLAG_LATEX)
    fontCode_ = latexFontManager_.getCode(font);
  else if (fontFlag_ == FONTFLAG_POSTSCRIPT)
    fontCode_ = postscriptFontManager_.getCode(font);
  else
    fontCode_ = 0;
  fontSize_ = font.getSize();
  AbstractGraphicDevice::setCurrentFont(font); 
}

void XFigGraphicDevice::drawLine(double x1, double y1, double x2, double y2)
{
  ostringstream oss;
  oss << "2 1 " << lineTypeCode_ << " " << getCurrentPointSize()
      << " " << fgColorCode_
      << " " << bgColorCode_
      << " " << getCurrentLayer()
      << " " << "-1 -1 -1 0 0 0 0 0 2" << endl;
  oss << round(x_(x1)) << " " << round(y_(y1)) << endl;
  oss << round(x_(x2)) << " " << round(y_(y2));
  content_.push_back(oss.str());
}
 
void XFigGraphicDevice::drawRect(double x, double y, double width, double height, short fill)
{
  ostringstream oss;
  oss << "2 2 0 " << getCurrentPointSize()
      << " " << fgColorCode_
      << " " << bgColorCode_
      << " " << getCurrentLayer()
      << " " << "-1"
      << " " << getFillCode(fill) << " -1 0 0 0 0 0 5" << endl;
  oss << round(x) << " " << round(y) << endl;
  oss << round(x_(x + width)) << " " << round(y_(y)) << endl;
  oss << round(x_(x + width)) << " " << round(y_(y + height)) << endl;
  oss << round(x_(x)) << " " << round(y_(y + height)) << endl;
  oss << round(x_(x)) << " " << round(y_(y));
  content_.push_back(oss.str());
}

void XFigGraphicDevice::drawCircle(double x, double y, double radius, short fill)
{
  ostringstream oss;
  oss << "1 3 0 " << getCurrentPointSize()
      << " " << fgColorCode_
      << " " << bgColorCode_
      << " " << getCurrentLayer()
      << " " << "-1"
      << " " << getFillCode(fill) << " -1 1 0 "
      << round(x_(x)) << " " << round(y_(y)) << " "
      << round(x_(radius)) << " " << round(y_(radius)) << " "
      << round(x_(x + radius)) << " " << round(y_(y)) << " "
      << round(x_(x + radius)) << " " << round(y_(y)) << endl;
  content_.push_back(oss.str());
}

void XFigGraphicDevice::drawText(double x, double y, const std::string& text, short hpos, short vpos, double angle) throw (UnvalidFlagException)
{
  int xrel = static_cast<int>(round(x_(x)));
  short sub = 0;
  if (hpos == TEXT_HORIZONTAL_LEFT)
    sub = 0;
  else if (hpos == TEXT_HORIZONTAL_CENTER) 
    sub = 1;
  else if (hpos == TEXT_HORIZONTAL_RIGHT)
    sub = 2;
  else throw UnvalidFlagException("XFigGraphicDevice::drawText(). Bad horizontal text alignment flag: " + TextTools::toString(hpos));
  
  int yrel = 0;
  if (vpos == TEXT_VERTICAL_BOTTOM)
    yrel = static_cast<int>(round(y_(y - 1.)));
  else if (vpos == TEXT_VERTICAL_CENTER)
    yrel = static_cast<int>(round(y + fontSize_ / 2 - 1)); 
  else if (vpos == TEXT_VERTICAL_TOP)
    yrel = static_cast<int>(round(y - fontSize_)); 
  else throw UnvalidFlagException("XFigGraphicDevice::drawText(). Bad vertical text alignment flag: " + TextTools::toString(vpos));

  ostringstream oss;
  oss << "4 " << sub << " " << fgColorCode_ << " " << 50 << " " << -1 << " " << fontCode_ << " " << fontSize_ << " " 
      << angle << " " << fontFlag_ << " " << -1 << " " << -1 << " " << xrel << " " << yrel << " " << text << "\\001";
  content_.push_back(oss.str());
}

int XFigGraphicDevice::getFillCode(short fill)
{
  if (fill == FILL_EMPTY) return -1;
  if (fill == FILL_FILLED) return 20;
  if (fill == FILL_PATTERN)
  {
    //TODO: define a field names currentPattern_, etc.
  }
  //Temp:
  return 20;
}

