//
// File: YN98WithRateParameter.cpp
// Created by: Jesse Bloom
// Created on: July 2015
//

/*
   Copyright or © or Copr. Bio++ Development Team, (November 16, 2004)
   This software is a computer program whose purpose is to provide classes
   for phylogenetic data analysis.

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

#include "YN98WithRateParameter.h"
#include <Bpp/Phyl/Model/Nucleotide/K80.h>
#include <Bpp/Phyl/Model/FrequenciesSet/CodonFrequenciesSet.h>
#include <Bpp/Numeric/NumConstants.h>


bppextensions::YN98WithRateParameter::YN98WithRateParameter(const bpp::GeneticCode* gc, bpp::FrequenciesSet* codonFreqs) :
  AbstractBiblioSubstitutionModel("YN98WithRateParameter."),
  pmodel_(new bpp::CodonDistanceFrequenciesSubstitutionModel(gc, new bpp::K80(dynamic_cast<const bpp::CodonAlphabet*>(gc->getSourceAlphabet())->getNucleicAlphabet()), codonFreqs)),
  rateparameter_(1)
{
  addParameter_(new bpp::Parameter("YN98WithRateParameter.kappa", 1, &bpp::Parameter::R_PLUS_STAR));
  addParameter_(new bpp::Parameter("YN98WithRateParameter.omega", 1, new bpp::IntervalConstraint(bpp::NumConstants::MILLI(), 999, true, true), true));
  addParameter_(new bpp::Parameter("YN98WithRateParameter.rateparameter", 1, new bpp::IntervalConstraint(0.0001, 1000, true, true), true));

  pmodel_->setNamespace("YN98WithRateParameter.");
  addParameters_(codonFreqs->getParameters());

  lParPmodel_.addParameters(pmodel_->getParameters());

  std::vector<std::string> v = pmodel_->getFrequenciesSet()->getParameters().getParameterNames();

  for (size_t i = 0; i < v.size(); i++)
  {
    mapParNamesFromPmodel_[v[i]] = getParameterNameWithoutNamespace(v[i]);
  }
  mapParNamesFromPmodel_["YN98WithRateParameter.123_K80.kappa"] = "kappa";
  mapParNamesFromPmodel_["YN98WithRateParameter.beta"] = "omega";

  updateMatrices();
}


bppextensions::YN98WithRateParameter::YN98WithRateParameter(const bppextensions::YN98WithRateParameter& model) : AbstractBiblioSubstitutionModel(model),
  pmodel_(new bpp::CodonDistanceFrequenciesSubstitutionModel(*model.pmodel_)),
  rateparameter_(model.rateparameter_)
{}

bppextensions::YN98WithRateParameter& bppextensions::YN98WithRateParameter::operator=(const bppextensions::YN98WithRateParameter& model)
{
  AbstractBiblioSubstitutionModel::operator=(model);
  pmodel_.reset(new bpp::CodonDistanceFrequenciesSubstitutionModel(*model.pmodel_));
  rateparameter_ = model.rateparameter_;
  return *this;
}

void bppextensions::YN98WithRateParameter::fireParameterChanged(const bpp::ParameterList& parameterlist) {
  rateparameter_ = getParameterValue("rateparameter");
  pmodel_->fireParameterChanged(parameterlist);
  AbstractBiblioSubstitutionModel::fireParameterChanged(parameterlist);
}

double bppextensions::YN98WithRateParameter::getCodonsMulRate(size_t i, size_t j) const {
    return rateparameter_ * pmodel_->getCodonsMulRate(i, j);
}
