//
// File CodonSiteTools.cpp
// Author : Sylvain Glémin
// Last modification : October 2004
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

#include "CodonSiteTools.h"
#include "Alphabet/CodonAlphabet.h"
#include "Alphabet/DNA.h"
#include "Alphabet/AlphabetTools.h"
#include "SiteTools.h"
#include "GeneticCode/GeneticCode.h"
#include "GeneticCode/StandardGeneticCode.h"
#include <Bpp/Utils/MapTools.h>
#include <Bpp/Numeric/NumTools.h>
#include <Bpp/Numeric/VectorTools.h>

using namespace bpp;

// From the STL:
#include <cmath>

using namespace std;

/******************************************************************************/

bool CodonSiteTools::hasGapOrStop(const Site& site, const GeneticCode& gCode) throw (AlphabetException)
{
  // Alphabet checking
  if (!AlphabetTools::isCodonAlphabet(site.getAlphabet()))
    throw AlphabetException("CodonSiteTools::hasGapOrStop: alphabet is not CodonAlphabet", site.getAlphabet());
  for (size_t i = 0; i < site.size(); i++)
  {
    if (site[i] < 0)
      return true;
  }
  return false;
}

/******************************************************************************/

bool CodonSiteTools::hasStop(const Site& site, const GeneticCode& gCode) throw (AlphabetException)
{
  // Alphabet checking
  if (!AlphabetTools::isCodonAlphabet(site.getAlphabet()))
    throw AlphabetException("CodonSiteTools::hasStop: alphabet is not CodonAlphabet", site.getAlphabet());
  for (size_t i = 0; i < site.size(); i++)
  {
    if (gCode.isStop(site[i]))
      return true;
  }
  return false;
}

/******************************************************************************/

bool CodonSiteTools::isMonoSitePolymorphic(const Site& site) throw (AlphabetException, EmptySiteException)
{
  // Alphabet checking
  if (!AlphabetTools::isCodonAlphabet(site.getAlphabet()))
    throw AlphabetException("CodonSiteTools::isMonoSitePolymorphic: alphabet is not CodonAlphabet", site.getAlphabet());
  // Empty site checking
  if (site.size() == 0)
    throw EmptySiteException("CodonSiteTools::isMonoSitePolymorphic: Incorrect specified site", &site);

  // Global polymorphism checking
  if (SiteTools::isConstant(site))
    return false;
  // initialisation of the 3 sub-sites ot the codon
  vector<int> pos1, pos2, pos3;
  const CodonAlphabet* ca = dynamic_cast<const CodonAlphabet*>(site.getAlphabet());
  for (size_t i = 0; i < site.size(); i++)
  {
    pos1.push_back(ca->getFirstPosition(site[i]));
    pos2.push_back(ca->getSecondPosition(site[i]));
    pos3.push_back(ca->getThirdPosition(site[i]));
  }
  const NucleicAlphabet* na = ca->getNucleicAlphabet();
  Site s1(pos1, na), s2(pos2, na), s3(pos3, na);
  // polymorphism checking for each sub-sites
  size_t nbpol = 0;
  if (!SiteTools::isConstant(s1))
    nbpol++;
  if (!SiteTools::isConstant(s2))
    nbpol++;
  if (!SiteTools::isConstant(s3))
    nbpol++;
  if (nbpol > 1)
    return false;
  return true;
}

/******************************************************************************/

bool CodonSiteTools::isSynonymousPolymorphic(const Site& site, const GeneticCode& gCode)
throw (AlphabetException, AlphabetMismatchException, EmptySiteException)
{
  // Alphabet checking
  if (!AlphabetTools::isCodonAlphabet(site.getAlphabet()))
    throw AlphabetException("CodonSiteTools::isSynonymousPolymorphic: alphabet is not CodonAlphabet", site.getAlphabet());
  if (!site.getAlphabet()->equals(*gCode.getSourceAlphabet()))
    throw AlphabetMismatchException("CodonSiteTools::isSynonymousPolymorphic: site and genetic code have not the same codon alphabet.", site.getAlphabet(), gCode.getSourceAlphabet());
  // Empty site checking
  if (site.size() == 0)
    throw EmptySiteException("CodonSiteTools::isSynonymousPolymorphic: Incorrect specified site", &site);

  // Global polymorphism checking
  if (SiteTools::isConstant(site))
    return false;

  // Synonymous polymorphism checking
  vector<int> prot;
  int first_aa = gCode.translate(site[0]);
  for (size_t i = 1; i < site.size(); i++)
  {
    int aa = gCode.translate(site[i]);
    if (aa != first_aa)
      return false;
  }
  return true;
}

/******************************************************************************/

Site* CodonSiteTools::generateCodonSiteWithoutRareVariant(const Site& site, const GeneticCode& gCode, double freqmin)
throw (AlphabetException, EmptySiteException)
{
  // Alphabet checking
  if (!AlphabetTools::isCodonAlphabet(site.getAlphabet()))
    throw AlphabetException("CodonSiteTools::generateCodonSiteWithoutRareVariant: alphabet is not CodonAlphabet", site.getAlphabet());
  // Empty site checking
  if (site.size() == 0)
    throw EmptySiteException("CodonSiteTools::generateCodonSiteWithoutRareVariant: Incorrect specified site", &site);

  if (SiteTools::isConstant(site))
  {
    Site* noRareVariant = new Site(site);
    return noRareVariant;
  }
  else
  {
    // Computation
    map<int, double> freqcodon;
    SiteTools::getFrequencies(site, freqcodon);
    const CodonAlphabet* ca = dynamic_cast<const CodonAlphabet*>(site.getAlphabet());
    const NucleicAlphabet* na = ca->getNucleicAlphabet();
    int newcodon = -1;
    for (map<int, double>::iterator it = freqcodon.begin(); it != freqcodon.end(); it++)
    {
      if (it->second > freqmin && !gCode.isStop(it->first))
      {
        newcodon = it->first;
        break;
      }
    }
    vector<int> pos1, pos2, pos3;
    for (size_t i = 0; i < site.size(); i++)
    {
      pos1.push_back(ca->getFirstPosition(site[i]));
      pos2.push_back(ca->getSecondPosition(site[i]));
      pos3.push_back(ca->getThirdPosition(site[i]));
    }
    Site s1(pos1, na), s2(pos2, na), s3(pos3, na);
    map<int, double> freq1;
    SiteTools::getFrequencies(s1, freq1);
    map<int, double> freq2;
    SiteTools::getFrequencies(s2, freq2);
    map<int, double> freq3;
    SiteTools::getFrequencies(s3, freq3);
    vector<int> codon;
    for (size_t i = 0; i < site.size(); i++)
    {
      if (freq1[s1.getValue(i)] > freqmin && freq2[s2.getValue(i)] > freqmin && freq3[s3.getValue(i)] > freqmin)
      {
        codon.push_back(site.getValue(i));
      }
      else
        codon.push_back(newcodon);
    }
    Site* noRareVariant = new Site(codon, ca);
    return noRareVariant;
  }
}

/******************************************************************************/

size_t CodonSiteTools::numberOfDifferences(int i, int j, const CodonAlphabet& ca)
{
  size_t nbdif = 0;
  if (ca.getFirstPosition(i) != ca.getFirstPosition(j))
    nbdif++;
  if (ca.getSecondPosition(i) != ca.getSecondPosition(j))
    nbdif++;
  if (ca.getThirdPosition(i) != ca.getThirdPosition(j))
    nbdif++;
  return nbdif;
}

/******************************************************************************/

double CodonSiteTools::numberOfSynonymousDifferences(int i, int j, const GeneticCode& gCode, bool minchange)
{
  const CodonAlphabet* ca = dynamic_cast<const CodonAlphabet*>(gCode.getSourceAlphabet());
  vector<int> ci = ca->getPositions(i);
  vector<int> cj = ca->getPositions(j);

  switch (numberOfDifferences(i, j, *ca))
  {
  case 0: return 0;
  case 1:
  {
    if (gCode.areSynonymous(i, j))
      return 1;
    return 0;
  }
  case 2:
  {
    if (gCode.areSynonymous(i, j))
      return 2;
    vector<double> path(2, 0); // Vector of number of synonymous changes per path (2 here)
    vector<double> weight(2, 1); // Weight to exclude path through stop codon
    if (ci[0] == cj[0])
    {
      int trans1 = ca->getCodon(ci[0], cj[1], ci[2]); // transitory codon between NcNiNi et NcNjNj: NcNjNi, Nc = identical site
      int trans2 = ca->getCodon(ci[0], ci[1], cj[2]); // transitory codon between NcNiNi et NcNjNj: NcNiNj, Nc = identical site
      if (!gCode.isStop(trans1))
      {
        if (gCode.areSynonymous(i, trans1))
          path[0]++;
        if (gCode.areSynonymous(trans1, j))
          path[0]++;
      }
      else
        weight[0] = 0;
      if (!gCode.isStop(trans2))
      {
        if (gCode.areSynonymous(i, trans2))
          path[1]++;
        if (gCode.areSynonymous(trans2, j))
          path[1]++;
      }
      else
        weight[1] = 0;
    }
    if (ci[1] == cj[1])
    {
      int trans1 = ca->getCodon(cj[0], ci[1], ci[2]); // transitory codon between NiNcNi et NjNcNj: NjNcNi, Nc = identical site
      int trans2 = ca->getCodon(ci[0], ci[1], cj[2]); // transitory codon between NiNcNi et NjNcNj: NiNcNj, Nc = identical site
      if (!gCode.isStop(trans1))
      {
        if (gCode.areSynonymous(i, trans1))
          path[0]++;
        if (gCode.areSynonymous(trans1, j))
          path[0]++;
      }
      else
        weight[0] = 0;
      if (!gCode.isStop(trans2))
      {
        if (gCode.areSynonymous(i, trans2))
          path[1]++;
        if (gCode.areSynonymous(trans2, j))
          path[1]++;
      }
      else
        weight[1] = 0;
    }
    if (ci[2] == cj[2])
    {
      int trans1 = ca->getCodon(cj[0], ci[1], ci[2]); // transitory codon between NiNiNc et NjNjNc: NjNiNc, Nc = identical site
      int trans2 = ca->getCodon(ci[0], cj[1], ci[2]); // transitory codon between NiNiNc et NjNjNc: NiNjNc, Nc = identical site
      if (!gCode.isStop(trans1))
      {
        if (gCode.areSynonymous(i, trans1))
          path[0]++;
        if (gCode.areSynonymous(trans1, j))
          path[0]++;
      }
      else
        weight[0] = 0;
      if (!gCode.isStop(trans2))
      {
        if (gCode.areSynonymous(i, trans2))
          path[1]++;
        if (gCode.areSynonymous(trans2, j))
          path[1]++;
      }
      else
        weight[1] = 0;
    }
    if (minchange)
      return VectorTools::max(path);

    double nbdif = 0;
    for (size_t k = 0; k < 2; k++)
    {
      nbdif += path[k] * weight[k];
    }

    return nbdif / VectorTools::sum(weight);
  }
  case 3:
  {
    vector<double> path(6, 0);
    vector<double> weight(6, 1);
    // First transitory codons
    int trans100 = ca->getCodon(cj[0], ci[1], ci[2]);
    int trans010 = ca->getCodon(ci[0], cj[1], ci[2]);
    int trans001 = ca->getCodon(ci[0], ci[1], cj[2]);
    // Second transitory codons
    int trans110 = ca->getCodon(cj[0], cj[1], ci[2]);
    int trans101 = ca->getCodon(cj[0], ci[1], cj[2]);
    int trans011 = ca->getCodon(ci[0], cj[1], cj[2]);
    // Paths
    if (!gCode.isStop(trans100))
    {
      if (gCode.areSynonymous(i, trans100))
      {
        path[0]++; path[1]++;
      }
      if (!gCode.isStop(trans110))
      {
        if (gCode.areSynonymous(trans100, trans110))
          path[0]++;
        if (gCode.areSynonymous(trans110, j))
          path[0]++;
      }
      else
        weight[0] = 0;
      if (!gCode.isStop(trans101))
      {
        if (gCode.areSynonymous(trans100, trans101))
          path[1]++;
        if (gCode.areSynonymous(trans101, j))
          path[1]++;
      }
      else
        weight[1] = 0;
    }
    else
    {
      weight[0] = 0; weight[1] = 0;
    }
    if (!gCode.isStop(trans010))
    {
      if (gCode.areSynonymous(i, trans010))
      {
        path[2]++; path[3]++;
      }
      if (!gCode.isStop(trans110))
      {
        if (gCode.areSynonymous(trans010, trans110))
          path[2]++;
        if (gCode.areSynonymous(trans110, j))
          path[2]++;
      }
      else
        weight[2] = 0;
      if (!gCode.isStop(trans011))
      {
        if (gCode.areSynonymous(trans010, trans011))
          path[3]++;
        if (gCode.areSynonymous(trans011, j))
          path[3]++;
      }
      else
        weight[3] = 0;
    }
    else
    {
      weight[2] = 0; weight[3] = 0;
    }
    if (!gCode.isStop(trans001))
    {
      if (gCode.areSynonymous(i, trans001))
      {
        path[4]++; path[5]++;
      }
      if (!gCode.isStop(trans101))
      {
        if (gCode.areSynonymous(trans001, trans101))
          path[4]++;
        if (gCode.areSynonymous(trans101, j))
          path[4]++;
      }
      else
        weight[4] = 0;
      if (!gCode.isStop(trans011))
      {
        if (gCode.areSynonymous(trans001, trans011))
          path[5]++;
        if (gCode.areSynonymous(trans011, j))
          path[5]++;
      }
      else
        weight[5] = 0;
    }
    else
    {
      weight[4] = 0; weight[5] = 0;
    }
    if (minchange)
      return VectorTools::max(path);

    double nbdif = 0;
    for (size_t k = 0; k < 6; k++)
    {
      nbdif += path[k] * weight[k];
    }

    return nbdif / VectorTools::sum(weight);
  }
  }
  // This line is never reached but sends a warning if not there:
  return 0.;
}

/******************************************************************************/

double CodonSiteTools::piSynonymous(const Site& site, const GeneticCode& gCode, bool minchange)
throw (AlphabetException, AlphabetMismatchException, EmptySiteException)
{
  // Alphabet checking
  if (!AlphabetTools::isCodonAlphabet(site.getAlphabet()))
    throw AlphabetException("CodonSiteTools::piSynonymous: alphabet is not CodonAlphabet", site.getAlphabet());
  if (!site.getAlphabet()->equals(*gCode.getSourceAlphabet()))
    throw AlphabetMismatchException("CodonSiteTools::piSynonymous: site and genetic code have not the same codon alphabet.", site.getAlphabet(), gCode.getSourceAlphabet());
  // Empty site checking
  if (site.size() == 0)
    throw EmptySiteException("CodonSiteTools::piSynonymous: Incorrect specified site", &site);

  // General polymorphism checking
  if (SiteTools::isConstant(site))
    return 0;
  // Computation
  map<int, double> freq;
  SiteTools::getFrequencies(site, freq);
  double pi = 0;
  for (map<int, double>::iterator it1 = freq.begin(); it1 != freq.end(); it1++)
  {
    for (map<int, double>::iterator it2 = freq.begin(); it2 != freq.end(); it2++)
    {
      pi += (it1->second) * (it2->second) * (numberOfSynonymousDifferences(it1->first, it2->first, gCode, minchange));
    }
  }
  size_t n = site.size();
  return pi * static_cast<double>(n / (n - 1));
}

/******************************************************************************/

double CodonSiteTools::piNonSynonymous(const Site& site, const GeneticCode& gCode, bool minchange)
throw (AlphabetException, AlphabetMismatchException, EmptySiteException)
{
  // Alphabet checking
  if (!AlphabetTools::isCodonAlphabet(site.getAlphabet()))
    throw AlphabetException("CodonSiteTools::piNonSynonymous: alphabet is not CodonAlphabet", site.getAlphabet());
  if (!site.getAlphabet()->equals(*gCode.getSourceAlphabet()))
    throw AlphabetMismatchException("CodonSiteTools::piNonSynonymous: site and genetic code have not the same codon alphabet.", site.getAlphabet(), gCode.getSourceAlphabet());
  // Empty site checking
  if (site.size() == 0)
    throw EmptySiteException("CodonSiteTools::piSynonymous: Incorrect specified site", &site);

  // General polymorphism checking
  if (SiteTools::isConstant(site))
    return 0;
  if (isSynonymousPolymorphic(site, gCode))
    return 0;
  // Computation
  map<int, double> freq;
  SiteTools::getFrequencies(site, freq);
  const CodonAlphabet* ca = dynamic_cast<const CodonAlphabet*>(site.getAlphabet());
  double pi = 0;
  for (map<int, double>::iterator it1 = freq.begin(); it1 != freq.end(); it1++)
  {
    for (map<int, double>::iterator it2 = freq.begin(); it2 != freq.end(); it2++)
    {
      double nbtot = static_cast<double>(numberOfDifferences(it1->first, it2->first, *ca));
      double nbsyn = numberOfSynonymousDifferences(it1->first, it2->first, gCode, minchange);
      pi += (it1->second) * (it2->second) * (nbtot - nbsyn);
    }
  }
  size_t n = site.size();
  return pi * static_cast<double>(n / (n - 1));
}

/******************************************************************************/

double CodonSiteTools::numberOfSynonymousPositions(int i, const GeneticCode& gCode, double ratio) throw (Exception)
{
  const CodonAlphabet* ca = dynamic_cast<const CodonAlphabet*>(gCode.getSourceAlphabet());
  if (gCode.isStop(i))
    return 0;
  if (ca->isUnresolved(i))
    return 0;
  double nbsynpos = 0.0;
  vector<int> codon = ca->getPositions(i);
  int acid = gCode.translate(i);
  for (size_t pos = 0; pos < 3; ++pos)
  {
    for (int an = 0; an < 4; ++an)
    {
      if (an == codon[pos])
        continue;
      vector<int> mutcodon = codon;
      mutcodon[pos] = an;
      int intcodon = ca->getCodon(mutcodon[0], mutcodon[1], mutcodon[2]);
      if (gCode.isStop(intcodon))
        continue;
      int altacid = gCode.translate(intcodon);
      if (altacid == acid)   // if synonymous
      {
        if (((codon[pos] == 0 || codon[pos] == 2) && (mutcodon[pos] == 1 || mutcodon[pos] == 3)) ||
            ((codon[pos] == 1 || codon[pos] == 3) && (mutcodon[pos] == 0 || mutcodon[pos] == 2)))   // if it is a transversion
        {
          nbsynpos = nbsynpos + 1 / (ratio + 2);
        }
        else   // if transition
        {
          nbsynpos = nbsynpos + ratio / (ratio + 2);
        }
      }
    }
  }
  return nbsynpos;
}

/******************************************************************************/

double CodonSiteTools::meanNumberOfSynonymousPositions(const Site& site, const GeneticCode& gCode, double ratio)
throw (AlphabetException, AlphabetMismatchException, EmptySiteException)
{
  // Alphabet checking
  if (!AlphabetTools::isCodonAlphabet(site.getAlphabet()))
    throw AlphabetException("CodonSiteTools::meanNumberOfSynonymousPositions: alphabet is not CodonAlphabet", site.getAlphabet());
  if (!site.getAlphabet()->equals(*gCode.getSourceAlphabet()))
    throw AlphabetMismatchException("CodonSiteTools::meanNumberOfSynonymousPositions: site and genetic code have not the same codon alphabet.", site.getAlphabet(), gCode.getSourceAlphabet());
  // Empty site checking
  if (site.size() == 0)
    throw EmptySiteException("CodonSiteTools::meanNumberOfSynonymousPositions: Incorrect specified site", &site);

  // Computation
  double NbSyn = 0;
  map<int, double> freq;
  SiteTools::getFrequencies(site, freq);
  for (map<int, double>::iterator it = freq.begin(); it != freq.end(); it++)
  {
    NbSyn += (it->second) * numberOfSynonymousPositions(it->first, gCode, ratio);
  }
  return NbSyn;
}

/******************************************************************************/

size_t CodonSiteTools::numberOfSubsitutions(const Site& site, const GeneticCode& gCode, double freqmin)
throw (AlphabetException, EmptySiteException)
{
  // Alphabet checking
  if (!AlphabetTools::isCodonAlphabet(site.getAlphabet()))
    throw AlphabetException("CodonSiteTools::numberOfSubsitutions: alphabet is not CodonAlphabet", site.getAlphabet());
  // Empty site checking
  if (site.size() == 0)
    throw EmptySiteException("CodonSiteTools::numberOfSubsitutions: Incorrect specified site", &site);

  if (SiteTools::isConstant(site))
    return 0;
  Site* newsite;
  if (freqmin > 1. / static_cast<double>(site.size()))
    newsite = CodonSiteTools::generateCodonSiteWithoutRareVariant(site, gCode, freqmin);
  else
    newsite = new Site(site);
  // Computation
  if (SiteTools::hasGap(*newsite))
    return 0;
  vector<int> pos1, pos2, pos3;

  const CodonAlphabet* ca = dynamic_cast<const CodonAlphabet*>(site.getAlphabet());

  for (size_t i = 0; i < newsite->size(); i++)
  {
    pos1.push_back(ca->getFirstPosition(newsite->getValue(i)));
    pos2.push_back(ca->getSecondPosition(newsite->getValue(i)));
    pos3.push_back(ca->getThirdPosition(newsite->getValue(i)));
  }

  const NucleicAlphabet* na = ca->getNucleicAlphabet();

  Site s1(pos1, na), s2(pos2, na), s3(pos3, na);
  size_t Scodon = SiteTools::getNumberOfDistinctCharacters(*newsite) - 1;
  size_t Sbase = SiteTools::getNumberOfDistinctCharacters(s1) + SiteTools::getNumberOfDistinctCharacters(s2) + SiteTools::getNumberOfDistinctCharacters(s3) - 3;
  delete newsite;
  if (Scodon >= Sbase)
    return Scodon;
  else
    return Sbase;
}

/******************************************************************************/

size_t CodonSiteTools::numberOfNonSynonymousSubstitutions(const Site& site, const GeneticCode& gCode, double freqmin)
throw (AlphabetException, AlphabetMismatchException, EmptySiteException)
{
  // Alphabet checking
  if (!AlphabetTools::isCodonAlphabet(site.getAlphabet()))
    throw AlphabetException("CodonSiteTools::numberOfNonSynonymousSubstitutions: alphabet is not CodonAlphabet", site.getAlphabet());
  if (!site.getAlphabet()->equals(*gCode.getSourceAlphabet()))
    throw AlphabetMismatchException("CodonSiteTools::numberOfNonSynonymousSubstitutions: site and genetic code have not the same codon alphabet.", site.getAlphabet(), gCode.getSourceAlphabet());
  // Empty site checking
  if (site.size() == 0)
    throw EmptySiteException("CodonSiteTools::numberOfNonSynonymousSubstitutions: Incorrect specified site", &site);

  if (SiteTools::isConstant(site))
    return 0;
  Site* newsite;
  if (freqmin > 1. / static_cast<double>(site.size()))
    newsite = generateCodonSiteWithoutRareVariant(site, gCode, freqmin);
  else
    newsite = new Site(site);
  if (SiteTools::hasGap(*newsite))
    return 0;
  // computation
  map<int, size_t> count;
  SiteTools::getCounts(*newsite, count);
  size_t NaSup = 0;
  size_t Nminmin = 10;

  const CodonAlphabet* ca = dynamic_cast<const CodonAlphabet*>(site.getAlphabet());

  for (map<int, size_t>::iterator it1 = count.begin(); it1 != count.end(); it1++)
  {
    size_t Nmin = 10;
    for (map<int, size_t>::iterator it2 = count.begin(); it2 != count.end(); it2++)
    {
      size_t Ntot = numberOfDifferences(it1->first, it2->first, *ca);
      size_t Ns = (size_t)numberOfSynonymousDifferences(it1->first, it2->first, gCode, true);
      if (Nmin > Ntot - Ns && it1->first != it2->first)
        Nmin = Ntot - Ns;
    }
    NaSup += Nmin;
    if (Nmin < Nminmin)
      Nminmin = Nmin;
  }
  delete newsite;
  return NaSup - Nminmin;
}

/******************************************************************************/

vector<size_t> CodonSiteTools::fixedDifferences(const Site& siteIn, const Site& siteOut, int i, int j, const GeneticCode& gCode)
throw (AlphabetException, AlphabetMismatchException, EmptySiteException)
{
  // Alphabet checking
  if (!AlphabetTools::isCodonAlphabet(siteIn.getAlphabet()))
    throw AlphabetException("CodonSiteTools::fixedDifferences: alphabet is not CodonAlphabet (siteIn)", siteIn.getAlphabet());
  if (!AlphabetTools::isCodonAlphabet(siteOut.getAlphabet()))
    throw AlphabetException("CodonSiteTools::fixedDifferences: alphabet is not CodonAlphabet (siteOut)", siteOut.getAlphabet());
  if (!siteIn.getAlphabet()->equals(*gCode.getSourceAlphabet()))
    throw AlphabetMismatchException("CodonSiteTools::fixedDifferences: siteIn and genetic code have not the same codon alphabet.", siteIn.getAlphabet(), gCode.getSourceAlphabet());
  if (!siteOut.getAlphabet()->equals(*gCode.getSourceAlphabet()))
    throw AlphabetMismatchException("CodonSiteTools::fixedDifferences: siteOut and genetic code have not the same codon alphabet.", siteOut.getAlphabet(), gCode.getSourceAlphabet());
  // Empty site checking
  if (siteIn.size() == 0)
    throw EmptySiteException("CodonSiteTools::getFixedDifferences Incorrect specified site", &siteIn);
  if (siteOut.size() == 0)
    throw EmptySiteException("CodonSiteTools::getFixedDifferences Incorrect specified site", &siteOut);

  const CodonAlphabet* ca = dynamic_cast<const CodonAlphabet*>(gCode.getSourceAlphabet());

  size_t Ntot = numberOfDifferences(i, j, *ca);
  size_t Ns = (size_t) numberOfSynonymousDifferences(i, j, gCode, true);
  size_t Na = Ntot - Ns;
  size_t Nfix = Ntot;
  vector<int> pos1in, pos2in, pos3in, pos1out, pos2out, pos3out;

  for (size_t k = 0; k < siteIn.size(); k++)
  {
    pos1in.push_back(ca->getFirstPosition(siteIn[k]));
    pos2in.push_back(ca->getSecondPosition(siteIn[k]));
    pos3in.push_back(ca->getThirdPosition(siteIn[k]));
    pos1out.push_back(ca->getFirstPosition(siteOut[k]));
    pos2out.push_back(ca->getSecondPosition(siteOut[k]));
    pos3out.push_back(ca->getThirdPosition(siteOut[k]));
  }
  const NucleicAlphabet* na = ca->getNucleicAlphabet();

  Site s1in(pos1in, na), s2in(pos2in, na), s3in(pos3in, na);
  Site s1out(pos1out, na), s2out(pos2out, na), s3out(pos3out, na);
  bool test1 = false;
  bool test2 = false;
  bool test3 = false;
  if ( (!SiteTools::isConstant(s1in) || !SiteTools::isConstant(s1out)) && ca->getFirstPosition(i) != ca->getFirstPosition(j) )
  {
    test1 = true;
    Nfix--;
  }
  if ( (!SiteTools::isConstant(s2in) || !SiteTools::isConstant(s2out)) && ca->getSecondPosition(i) != ca->getSecondPosition(j) )
  {
    test2 = true;
    Nfix--;
  }
  if ( (!SiteTools::isConstant(s3in) || !SiteTools::isConstant(s3out)) && ca->getThirdPosition(i) != ca->getThirdPosition(j) )
  {
    test3 = true;
    Nfix--;
  }
  // Suppression of differences when not fixed
  vector<size_t> v(2);
  if (Nfix == 0)
  {
    v[0] = 0;
    v[1] = 0;
    return v;
  }
  if (Nfix < Ntot)
  {
    if (Na == 0)
      Ns = Nfix;
    if (Ns == 0)
      Na = Nfix;
    else
    {
      if (Ntot == 3)
      {
        if (Nfix == 1)
        {
          if (test1 && test2)
          {
            Na = 0; Ns = 1;
          }
          if (test1 && test3)
          {
            Na = 1; Ns = 0;
          }
          if (test2 && test3)
          {
            Na--; Ns--;
          }
        }
      }
      if (Nfix == 2)
      {
        if (test1)
        {
          Na = 1; Ns = 1;
        }
        if (test2)
          Na--;
        if (test3)
          Ns--;
      }
    }
    if (Ntot == 2)
    {
      if (test1)
      {
        if (ca->getSecondPosition(i) == ca->getSecondPosition(j))
          Na--;
        else
          Ns--;
      }
      if (test2)
        Na--;
      if (test3)
        Ns--;
    }
  }
  v[0] = Ns;
  v[1] = Na;
  return v;
}

/******************************************************************************/

bool CodonSiteTools::isFourFoldDegenerated(const Site& site, const GeneticCode& gCode)
{
  if (!SiteTools::isConstant(site, true))
  {
    /** If non-synonymous mutation **/
    if (!(CodonSiteTools::isSynonymousPolymorphic(site, gCode)))
      return false;

    for (size_t i = 0; i < site.size(); i++)
    {
      if (!(gCode.isFourFoldDegenerated(site.getValue(i))))
      {
        return false;
      }
    }
  }
  else
  {
    for (size_t i = 0; i < site.size(); i++)
    {
      if (!(gCode.isFourFoldDegenerated(site.getValue(i))))
      {
        return false;
      }
    }
  }
  return true;
}

/******************************************************************************/

