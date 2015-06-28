//
// File: SvgGraphicDevice.cpp
// Created by: Julien Dutheil
// Created on: Mon Mar 10 2008
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

#include "SvgGraphicDevice.h"

using namespace bpp;
using namespace std;

void SvgGraphicDevice::begin()
{
  layers_.clear();
  minx_ = maxx_ = miny_ = maxy_ = 0;
}

void SvgGraphicDevice::end()
{
  //Header:
  out_ << "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>" << endl;
  out_ << "<!DOCTYPE svg PUBLIC \"-//W3C//DTD Svg 1.1//EN\"" << endl;
  out_ << "\"http://www.w3.org/Graphics/Svg/1.1/DTD/svg11.dtd\">" << endl;
  out_ << "<svg width=\"" << (maxx_ - minx_) << "\" height=\"" << (maxy_ - miny_) << "\" version=\"1.1\"" << endl;
  out_ << " xmlns=\"http://www.w3.org/2000/svg\"" << endl;
  if (inkscapeEnabled_)
    out_ << " xmlns:inkscape=\"http://www.inkscape.org/namespaces/inkscape\"";
  out_ << " >" << endl;
  
  out_ << "<g transform=\"translate(" << (-minx_) << "," << (-miny_) << ")\">" << endl;

  for(map<int, vector<string> >::iterator it = layers_.begin(); it != layers_.end(); it++)
  {
    out_ << "<g id=\"layer" << it->first << "\"";
    if(inkscapeEnabled_)
    {
      out_ << " inkscape:groupmode=\"layer\"";
    }
    out_ << " >" << endl;
    vector<string> * v = &it->second;
    for(unsigned int i = 0; i < v->size(); i++)
    {
      out_ << (*v)[i] << endl;
    }
    out_ << "</g>" << endl;
  }
  out_ << "</g>" << endl;
  
  out_ << "</svg>" << endl;
}

void SvgGraphicDevice::drawLine(double x1, double y1, double x2, double y2)
{
  x1 = x_(x1);
  x2 = x_(x2);
  y1 = y_(y1);
  y2 = y_(y2);
  string style = "stroke:" + colorToText(getCurrentForegroundColor()) + ";stroke-width:" + TextTools::toString(getCurrentPointSize());
  if(getCurrentLineType() == LINE_DASHED)
    style += ";stroke-dasharray:4,4";
  else if(getCurrentLineType() == LINE_DOTTED)
    style += ";stroke-dasharray:1,2";
  ostringstream oss;
  oss << "<line x1=\"" << x1 << "\" y1=\"" << y1 << "\" x2=\"" << x2 << "\" y2=\"" << y2 << "\" style=\"" << style << "\" />";
  layers_[getCurrentLayer()].push_back(oss.str());
  if (x1 < minx_) minx_ = x1;
  if (x2 < minx_) minx_ = x2;
  if (y1 < miny_) miny_ = y1;
  if (y2 < miny_) miny_ = y2;
  if (x1 > maxx_) maxx_ = x1;
  if (x2 > maxx_) maxx_ = x2;
  if (y1 > maxy_) maxy_ = y1;
  if (y2 > maxy_) maxy_ = y2;
}
 
void SvgGraphicDevice::drawRect(double x, double y, double width, double height, short fill)
{
  x = x_(x);
  y = y_(y);
  width = x_(width);
  height = y_(height);
  string style = "stroke:" + colorToText(getCurrentForegroundColor()) + ";stroke-width:" + TextTools::toString(getCurrentPointSize());
  if(fill == FILL_FILLED)
  {
    style += ";fill:" + colorToText(getCurrentBackgroundColor());
  }
  ostringstream oss;
  oss << "<rect x=\"" << x << "\" y=\"" << y << "\" width=\"" << width << "\" height=\"" << height << "\" style=\"" << style << "\" />";
  layers_[getCurrentLayer()].push_back(oss.str());
  if (x < minx_) minx_ = x;
  if (y < miny_) miny_ = y;
  if (x + width > maxx_) maxx_ = x + width;
  if (y + height > maxy_) maxx_ = y + height;
}

void SvgGraphicDevice::drawCircle(double x, double y, double radius, short fill)
{
  x = x_(x);
  y = y_(y);
  radius = x_(radius);
  string style = "stroke:" + colorToText(getCurrentForegroundColor()) + ";stroke-width:" + TextTools::toString(getCurrentPointSize());
  if(fill == FILL_FILLED)
  {
    style += ";fill:" + colorToText(getCurrentBackgroundColor());
  }
  ostringstream oss;
  oss << "<rect cx=\"" << x << "\" cy=\"" << y << "\" cr=\"" << radius << "\" style=\"" << style << "\" />";
  layers_[getCurrentLayer()].push_back(oss.str());
}

void SvgGraphicDevice::drawText(double x, double y, const std::string & text, short hpos, short vpos, double angle) throw (UnvalidFlagException)
{
  x = x_(x);
  y = y_(y);
  string style = "font-family:" + getCurrentFont().getFamily() + ";font-style:" + fontStyles_[getCurrentFont().getStyle()] + ";font-weight:" + fontWeights_[getCurrentFont().getWeight()] + ";font-size:" + TextTools::toString(getCurrentFont().getSize()) + "px";
  style += ";dominant-baseline:";
  if (vpos == TEXT_VERTICAL_BOTTOM)
    style += "before-edge";
  else if (vpos == TEXT_VERTICAL_TOP)
    style += "after-edge";
  else if (vpos == TEXT_VERTICAL_CENTER)
    style += "middle";
  else throw UnvalidFlagException("SvgGraphicDevice::drawText. Invalid vertical alignment option.");
  style += ";text-anchor:";
  if (hpos == TEXT_HORIZONTAL_LEFT)
    style += "start";
  else if (hpos == TEXT_HORIZONTAL_RIGHT)
    style += "end";
  else if (hpos == TEXT_HORIZONTAL_CENTER)
    style += "middle";
  else throw UnvalidFlagException("SvgGraphicDevice::drawText. Invalid horizontal alignment option.");
  style += ";fill:" + colorToText(getCurrentForegroundColor());

  ostringstream oss;
  oss << "<text x=\"" << x << "\" y=\"" << y << "\" rotate=\"" << angle << "\" style=\"" << style << "\" >" << text << "</text>";
  layers_[getCurrentLayer()].push_back(oss.str());
}

