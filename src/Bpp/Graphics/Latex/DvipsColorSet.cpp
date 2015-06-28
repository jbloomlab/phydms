//
// File: DvipsColorSet.cpp
// Created by: Julien Dutheil
// Created on: Mon Apr 14 2008
//

/*
Copyright or © or Copr. CNRS, (November 17, 2008)

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

#include "DvipsColorSet.h"

using namespace bpp;

DvipsColorSet::DvipsColorSet()
{
  colors_["GreenYellow"] = RGBColor(217, 255, 79);
  colors_["Yellow"] = RGBColor(255, 255, 0);
  colors_["Goldenrod"] = RGBColor(255, 230, 41);
  colors_["Dandelion"] = RGBColor(255, 181, 41);
  colors_["Apricot"] = RGBColor(255, 173, 122);
  colors_["Peach"] = RGBColor(255, 128, 77);
  colors_["Melon"] = RGBColor(255, 138, 128);
  colors_["YellowOrange"] = RGBColor(255, 148, 0);
  colors_["Orange"] = RGBColor(255, 99, 33);
  colors_["BurntOrange"] = RGBColor(255, 125, 0);
  colors_["Bittersweet"] = RGBColor(194, 48, 0);
  colors_["RedOrange"] = RGBColor(255, 59, 33);
  colors_["Mahogany"] = RGBColor(166, 25, 22);
  colors_["Maroon"] = RGBColor(173, 23, 55);
  colors_["BrickRed"] = RGBColor(184, 20, 11);
  colors_["Red"] = RGBColor(255, 0, 0);
  colors_["OrangeRed"] = RGBColor(255, 0, 128);
  colors_["RubineRed"] = RGBColor(255, 0, 222);
  colors_["WildStrawberry"] = RGBColor(255, 10, 156);
  colors_["Salmon"] = RGBColor(255, 120, 158);
  colors_["CarnationPink"] = RGBColor(255, 94, 255);
  colors_["Magenta"] = RGBColor(255, 0, 255);
  colors_["VioletRed"] = RGBColor(255, 48, 255);
  colors_["Rhodamine"] = RGBColor(255, 46, 255);
  colors_["Mulberry"] = RGBColor(165, 25, 250);
  colors_["RedViolet"] = RGBColor(157, 17, 168);
  colors_["Fuchsia"] = RGBColor(124, 21, 235);
  colors_["Lavender"] = RGBColor(255, 133, 255);
  colors_["Thistle"] = RGBColor(224, 105, 255);
  colors_["Orchid"] = RGBColor(173, 92, 255);
  colors_["DarkOrchid"] = RGBColor(153, 51, 204);
  colors_["Purple"] = RGBColor(140, 36, 255);
  colors_["Plum"] = RGBColor(128, 0, 255);
  colors_["Violet"] = RGBColor(54, 31, 255);
  colors_["RoyalPurple"] = RGBColor(64, 25, 255);
  colors_["BlueViolet"] = RGBColor(34, 22, 245);
  colors_["Periwinkle"] = RGBColor(110, 115, 255);
  colors_["CadetBlue"] = RGBColor(97, 110, 196);
  colors_["CornflowerBlue"] = RGBColor(89, 222, 255);
  colors_["MidnightBlue"] = RGBColor(3, 126, 145);
  colors_["NavyBlue"] = RGBColor(15, 117, 255);
  colors_["RoyalBlue"] = RGBColor(0, 128, 255);
  colors_["Blue"] = RGBColor(0, 0, 255);
  colors_["Cerulean"] = RGBColor(15, 227, 255);
  colors_["Cyan"] = RGBColor(0, 255, 255);
  colors_["ProcessBlue"] = RGBColor(10, 255, 255);
  colors_["SkyBlue"] = RGBColor(97, 255, 224);
  colors_["Turquoise"] = RGBColor(38, 255, 204);
  colors_["TealBlue"] = RGBColor(35, 250, 165);
  colors_["Aquamarine"] = RGBColor(46, 255, 178);
  colors_["BlueGreen"] = RGBColor(38, 255, 171);
  colors_["Emerald"] = RGBColor(0, 255, 128);
  colors_["JungleGreen"] = RGBColor(3, 255, 122);
  colors_["SeaGreen"] = RGBColor(79, 255, 128);
  colors_["Green"] = RGBColor(0, 255, 0);
  colors_["ForestGreen"] = RGBColor(20, 224, 27);
  colors_["PineGreen"] = RGBColor(15, 191, 78);
  colors_["LimeGreen"] = RGBColor(128, 255, 0);
  colors_["YellowGreen"] = RGBColor(143, 255, 66);
  colors_["SpringGreen"] = RGBColor(189, 255, 61);
  colors_["OliveGreen"] = RGBColor(55, 153, 8);
  colors_["RawSienna"] = RGBColor(140, 39, 0);
  colors_["Sepia"] = RGBColor(77, 13, 0);
  colors_["Brown"] = RGBColor(102, 19, 0);
  colors_["Tan"] = RGBColor(219, 148, 112);
  colors_["Gray"] = RGBColor(128, 128, 128);
  colors_["Black"] = RGBColor(0, 0, 0);
  colors_["White"] = RGBColor(255, 255, 255);
}

