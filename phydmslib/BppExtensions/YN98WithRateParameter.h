//
// File: YN98WithRateParemeter.h
// Created by: Jesse Bloom
// Created on: July 2015
//

/*
  This file was created by modifying the YN98.h file distributed with Bio++
*/

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

#include <Bpp/Phyl/Model/AbstractBiblioSubstitutionModel.h>
#include <Bpp/Phyl/Model/Codon/CodonDistanceFrequenciesSubstitutionModel.h>

namespace bppextensions
{
/**
 * @brief Extension of the YN98 (Yang and Nielsen, 1998) substitution model for codons that adds a rate parameter.
 *
 * @author Jesse Bloom
 *
 * This model adds a single parameter corresponding to the rate parameter that scales
 * all substitution rates. You would use this over the standard YN98 if you have fixed
 * branch lengths and effectively want to scale these.
 */
class YN98WithRateParameter :
    public bpp::AbstractBiblioSubstitutionModel,
    public virtual bpp::CodonReversibleSubstitutionModel
{
private:
  std::auto_ptr<bpp::CodonDistanceFrequenciesSubstitutionModel> pmodel_;
  double rateparameter_;

public:
  YN98WithRateParameter(const bpp::GeneticCode* gc, bpp::FrequenciesSet* codonFreqs);

  YN98WithRateParameter(const YN98WithRateParameter& model);

  YN98WithRateParameter& operator=(const YN98WithRateParameter&);

  virtual ~YN98WithRateParameter() {}

  YN98WithRateParameter* clone() const { return new YN98WithRateParameter(*this); }

public:
  void fireParameterChanged(const bpp::ParameterList& parameterlist);

  std::string getName() const { return "YN98WithRateParameter"; }

  const SubstitutionModel& getModel() const { return *pmodel_.get(); }

  const bpp::GeneticCode* getGeneticCode() const { return pmodel_->getGeneticCode(); }
  
  double getCodonsMulRate(size_t i, size_t j) const;

private:
  SubstitutionModel& getModel() { return *pmodel_.get(); }

};

} // end of namespace bppextensions.
