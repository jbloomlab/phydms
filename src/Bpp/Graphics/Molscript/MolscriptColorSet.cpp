//
// File: MolscriptColorSet.cpp
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

#include "MolscriptColorSet.h"

using namespace bpp;

MolscriptColorSet::MolscriptColorSet()
{
  colors_["aliceblue"] = RGBColor(240, 248, 255);
  colors_["antiquewhite"] = RGBColor(250, 235, 215);
  colors_["aquamarine"] = RGBColor(127, 255, 212);
  colors_["azure"] = RGBColor(240, 255, 255);
  colors_["beige"] = RGBColor(245, 245, 220);
  colors_["bisque"] = RGBColor(255, 228, 196);
  colors_["black"] = RGBColor(0, 0, 0);
  colors_["blanchedalmond"] = RGBColor(255, 235, 205);
  colors_["blue"] = RGBColor(0, 0, 255);
  colors_["blueviolet"] = RGBColor(138, 43, 226);
  colors_["brown"] = RGBColor(165, 42, 42);
  colors_["burlywood"] = RGBColor(222, 184, 135);
  colors_["cadetblue"] = RGBColor(95, 158, 160);
  colors_["chartreuse"] = RGBColor(127, 255, 0);
  colors_["chocolate"] = RGBColor(210, 105, 30);
  colors_["coral"] = RGBColor(255, 127, 80);
  colors_["cornflowerblue"] = RGBColor(100, 149, 237);
  colors_["cornsilk"] = RGBColor(255, 248, 220);
  colors_["crimson"] = RGBColor(220, 20, 60);
  colors_["cyan"] = RGBColor(0, 255, 255);
  colors_["darkblue"] = RGBColor(0, 0, 139);
  colors_["darkcyan"] = RGBColor(0, 139, 139);
  colors_["darkgoldenrod"] = RGBColor(184, 134, 11);
  colors_["darkgray"] = RGBColor(169, 169, 169);
  colors_["darkgreen"] = RGBColor(0, 100, 0);
  colors_["darkgrey"] = RGBColor(169, 169, 169);
  colors_["darkkhaki"] = RGBColor(189, 183, 107);
  colors_["darkmagenta"] = RGBColor(139, 0, 139);
  colors_["darkolivegreen"] = RGBColor(85, 107, 47);
  colors_["darkorange"] = RGBColor(255, 140, 0);
  colors_["darkorchid"] = RGBColor(153, 50, 204);
  colors_["darkred"] = RGBColor(139, 0, 0);
  colors_["darksalmon"] = RGBColor(233, 150, 122);
  colors_["darkseagreen"] = RGBColor(143, 188, 143);
  colors_["darkslateblue"] = RGBColor(72, 61, 139);
  colors_["darkslategray"] = RGBColor(47, 79, 79);
  colors_["darkslategrey"] = RGBColor(47, 79, 79);
  colors_["darkturquoise"] = RGBColor(0, 206, 209);
  colors_["darkviolet"] = RGBColor(148, 0, 211);
  colors_["deeppink"] = RGBColor(255, 20, 147);
  colors_["deepskyblue"] = RGBColor(0, 191, 255);
  colors_["dimgray"] = RGBColor(105, 105, 105);
  colors_["dimgrey"] = RGBColor(105, 105, 105);
  colors_["dodgerblue"] = RGBColor(30, 144, 255);
  colors_["firebrick"] = RGBColor(178, 34, 34);
  colors_["floralwhite"] = RGBColor(255, 250, 240);
  colors_["forestgreen"] = RGBColor(34, 139, 34);
  colors_["gainsboro"] = RGBColor(220, 220, 220);
  colors_["ghostwhite"] = RGBColor(248, 248, 255);
  colors_["gold"] = RGBColor(255, 215, 0);
  colors_["goldenrod"] = RGBColor(218, 165, 32);
  colors_["gray"] = RGBColor(190, 190, 190);
  colors_["green"] = RGBColor(0, 255, 0);
  colors_["greenyellow"] = RGBColor(173, 255, 47);
  colors_["grey"] = RGBColor(190, 190, 190);
  colors_["honeydew"] = RGBColor(240, 255, 240);
  colors_["hotpink"] = RGBColor(255, 105, 180);
  colors_["indianred"] = RGBColor(205, 92, 92);
  colors_["indigo"] = RGBColor(75, 0, 130);
  colors_["ivory"] = RGBColor(255, 255, 240);
  colors_["khaki"] = RGBColor(240, 230, 140);
  colors_["lavender"] = RGBColor(230, 230, 250);
  colors_["lavenderblush"] = RGBColor(255, 240, 245);
  colors_["lawngreen"] = RGBColor(124, 252, 0);
  colors_["lemonchiffon"] = RGBColor(255, 250, 205);
  colors_["lightblue"] = RGBColor(173, 216, 230);
  colors_["lightcoral"] = RGBColor(240, 128, 128);
  colors_["lightcyan"] = RGBColor(224, 255, 255);
  colors_["lightgoldenrod"] = RGBColor(238, 221, 130);
  colors_["lightgoldenrodyellow"] = RGBColor(250, 250, 210);
  colors_["lightgray"] = RGBColor(211, 211, 211);
  colors_["lightgreen"] = RGBColor(144, 238, 144);
  colors_["lightgrey"] = RGBColor(211, 211, 211);
  colors_["lightpink"] = RGBColor(255, 182, 193);
  colors_["lightsalmon"] = RGBColor(255, 160, 122);
  colors_["lightseagreen"] = RGBColor(32, 178, 170);
  colors_["lightskyblue"] = RGBColor(135, 206, 250);
  colors_["lightslateblue"] = RGBColor(132, 112, 255);
  colors_["lightslategray"] = RGBColor(119, 136, 153);
  colors_["lightslategrey"] = RGBColor(119, 136, 153);
  colors_["lightsteelblue"] = RGBColor(176, 196, 222);
  colors_["lightyellow"] = RGBColor(255, 255, 224);
  colors_["limegreen"] = RGBColor(50, 205, 50);
  colors_["linen"] = RGBColor(250, 240, 230);
  colors_["magenta"] = RGBColor(255, 0, 255);
  colors_["maroon"] = RGBColor(176, 48, 96);
  colors_["mediumaquamarine"] = RGBColor(102, 205, 170);
  colors_["mediumblue"] = RGBColor(0, 0, 205);
  colors_["mediumorchid"] = RGBColor(186, 85, 211);
  colors_["mediumpurple"] = RGBColor(147, 112, 219);
  colors_["mediumseagreen"] = RGBColor(60, 179, 113);
  colors_["mediumslateblue"] = RGBColor(123, 104, 238);
  colors_["mediumspringgreen"] = RGBColor(0, 250, 154);
  colors_["mediumturquoise"] = RGBColor(72, 209, 204);
  colors_["mediumvioletred"] = RGBColor(199, 21, 133);
  colors_["midnightblue"] = RGBColor(25, 25, 112);
  colors_["mintcream"] = RGBColor(245, 255, 250);
  colors_["mistyrose"] = RGBColor(255, 228, 225);
  colors_["moccasin"] = RGBColor(255, 228, 181);
  colors_["navajowhite"] = RGBColor(255, 222, 173);
  colors_["navy"] = RGBColor(0, 0, 128);
  colors_["navyblue"] = RGBColor(0, 0, 128);
  colors_["oldlace"] = RGBColor(253, 245, 230);
  colors_["olivedrab"] = RGBColor(107, 142, 35);
  colors_["orange"] = RGBColor(255, 165, 0);
  colors_["orangered"] = RGBColor(255, 69, 0);
  colors_["orchid"] = RGBColor(218, 112, 214);
  colors_["palegoldenrod"] = RGBColor(238, 232, 170);
  colors_["palegreen"] = RGBColor(152, 251, 152);
  colors_["paleturquoise"] = RGBColor(175, 238, 238);
  colors_["palevioletred"] = RGBColor(219, 112, 147);
  colors_["papayawhip"] = RGBColor(255, 239, 213);
  colors_["peachpuff"] = RGBColor(255, 218, 185);
  colors_["peru"] = RGBColor(205, 133, 63);
  colors_["pink"] = RGBColor(255, 192, 203);
  colors_["plum"] = RGBColor(221, 160, 221);
  colors_["powderblue"] = RGBColor(176, 224, 230);
  colors_["purple"] = RGBColor(160, 32, 240);
  colors_["red"] = RGBColor(255, 0, 0);
  colors_["rosybrown"] = RGBColor(188, 143, 143);
  colors_["royalblue"] = RGBColor(65, 105, 225);
  colors_["saddlebrown"] = RGBColor(139, 69, 19);
  colors_["salmon"] = RGBColor(250, 128, 114);
  colors_["sandybrown"] = RGBColor(244, 164, 96);
  colors_["seagreen"] = RGBColor(46, 139, 87);
  colors_["seashell"] = RGBColor(255, 245, 238);
  colors_["sgibeet"] = RGBColor(142, 56, 142);
  colors_["sgibrightgray"] = RGBColor(197, 193, 170);
  colors_["sgibrightgrey"] = RGBColor(197, 193, 170);
  colors_["sgichartreuse"] = RGBColor(113, 198, 113);
  colors_["sgidarkgray"] = RGBColor(85, 85, 85);
  colors_["sgidarkgrey"] = RGBColor(85, 85, 85);
  colors_["sgilightblue"] = RGBColor(125, 158, 192);
  colors_["sgilightgray"] = RGBColor(170, 170, 170);
  colors_["sgilightgrey"] = RGBColor(170, 170, 170);
  colors_["sgimediumgray"] = RGBColor(132, 132, 132);
  colors_["sgimediumgrey"] = RGBColor(132, 132, 132);
  colors_["sgiolivedrab"] = RGBColor(142, 142, 56);
  colors_["sgisalmon"] = RGBColor(198, 113, 113);
  colors_["sgislateblue"] = RGBColor(113, 113, 198);
  colors_["sgiteal"] = RGBColor(56, 142, 142);
  colors_["sgiverydarkgray"] = RGBColor(40, 40, 40);
  colors_["sgiverydarkgrey"] = RGBColor(40, 40, 40);
  colors_["sgiverylightgray"] = RGBColor(214, 214, 214);
  colors_["sgiverylightgrey"] = RGBColor(214, 214, 214);
  colors_["sienna"] = RGBColor(160, 82, 45);
  colors_["skyblue"] = RGBColor(135, 206, 235);
  colors_["slateblue"] = RGBColor(106, 90, 205);
  colors_["slategray"] = RGBColor(112, 128, 144);
  colors_["slategrey"] = RGBColor(112, 128, 144);
  colors_["snow"] = RGBColor(255, 250, 250);
  colors_["springgreen"] = RGBColor(0, 255, 127);
  colors_["steelblue"] = RGBColor(70, 130, 180);
  colors_["tan"] = RGBColor(210, 180, 140);
  colors_["thistle"] = RGBColor(216, 191, 216);
  colors_["tomato"] = RGBColor(255, 99, 71);
  colors_["turquoise"] = RGBColor(64, 224, 208);
  colors_["violet"] = RGBColor(238, 130, 238);
  colors_["violetred"] = RGBColor(208, 32, 144);
  colors_["wheat"] = RGBColor(245, 222, 179);
  colors_["white"] = RGBColor(255, 255, 255);
  colors_["whitesmoke"] = RGBColor(245, 245, 245);
  colors_["yellow"] = RGBColor(255, 255, 0);
  colors_["yellowgreen"] = RGBColor(154, 205, 50);
}
