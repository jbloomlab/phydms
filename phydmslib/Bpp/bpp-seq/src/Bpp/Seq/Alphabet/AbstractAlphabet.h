//
// File: AbstractAlphabet.h
// Authors: Guillaume Deuchst
//          Julien Dutheil
//          Sylvain Gaillard
// Created on: Tue Jul 22 2003
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

#ifndef _ABSTRACTALPHABET_H_
#define _ABSTRACTALPHABET_H_

#include "Alphabet.h"
#include "AlphabetState.h"
#include <Bpp/Exceptions.h>

// From the STL:
#include <string>
#include <vector>
#include <map>

namespace bpp
{

  /**
   * @brief A partial implementation of the Alphabet interface.
   *
   * It contains a vector of AlphabetState.
   * All methods are based uppon this vector
   * but do not provide any method to initialize it.
   * This is up to each constructor of the derived classes.
   *
   * @see Alphabet
   */
  class AbstractAlphabet:
    public Alphabet
  {
  private:
		
    /**
     * @brief Alphabet: vector of AlphabetState.
     */
    std::vector<AlphabetState*> alphabet_;
    /**
     * @name maps used to quick search for letter and num.
     * @{
     */
    std::map<std::string, size_t> letters_;
    std::map<int, size_t> nums_;
    /** @} */
    /**
     * @brief Update the private maps letters_ and nums_ when adding a state.
     *
     * @param pos The index of the state in the alphabet_ vector.
     * @param st The state that have been added or modified
     */
    void updateMaps_(size_t pos, const AlphabetState& st);

  protected:
    /**
     * @name Available codes
     *
     * These vectors will be computed the first time you call the getAvailableInts or getAvailableChars method.
     *
     * @{
     */
    mutable std::vector<std::string> charList_;
    mutable std::vector<int> intList_;
    /** @} */

  public:
		
    AbstractAlphabet(): alphabet_(), letters_(), nums_(), charList_(), intList_() {}

    AbstractAlphabet(const AbstractAlphabet& alph) : alphabet_(), letters_(alph.letters_), nums_(alph.nums_), charList_(alph.charList_), intList_(alph.intList_)
    {
      for (size_t i = 0; i < alph.alphabet_.size(); ++i)
        alphabet_.push_back(new AlphabetState(*alph.alphabet_[i]));
    }

    AbstractAlphabet& operator=(const AbstractAlphabet& alph)
    {
      for (size_t i = 0 ; i < alphabet_.size() ; ++i)
        delete alphabet_[i];

      for (size_t i = 0; i < alph.alphabet_.size(); ++i)
        alphabet_.push_back(new AlphabetState(*alph.alphabet_[i]));

      letters_  = alph.letters_;
      nums_     = alph.nums_;
      charList_ = alph.charList_;
      intList_  = alph.intList_;

      return *this;
    }

#ifndef NO_VIRTUAL_COV
    virtual AbstractAlphabet* clone() const = 0;
#endif

    virtual ~AbstractAlphabet()
    {
      for (size_t i = 0 ; i < alphabet_.size() ; ++i)
        delete alphabet_[i];
    }
	
  public:
    /**
     * @name Implement these methods from the Alphabet interface.
     *
     * @{
     */
    size_t getNumberOfStates() const { return alphabet_.size(); }
    unsigned int getNumberOfChars() const { return static_cast<unsigned int>(alphabet_.size()); }
    std::string getName(const std::string& state) const throw (BadCharException);
    std::string getName(int state) const throw (BadIntException);
    int charToInt(const std::string& state) const throw (BadCharException);
    std::string intToChar(int state) const throw (BadIntException);
    bool isIntInAlphabet(int state) const;
    bool isCharInAlphabet(const std::string& state) const;
    std::vector<int> getAlias(int state) const throw (BadIntException);
    std::vector<std::string> getAlias(const std::string& state) const throw (BadCharException);
    int getGeneric(const std::vector<int>& states) const throw (BadIntException);
    std::string getGeneric(const std::vector<std::string>& states) const throw (AlphabetException);
    const std::vector<int>& getSupportedInts() const;
    const std::vector<std::string>& getSupportedChars() const;
    const std::vector<std::string> & getResolvedChars() const;
    int getGapCharacterCode() const { return -1; }
    bool isGap(int state) const { return state == -1; }
    bool isGap(const std::string& state) const { return charToInt(state) == -1; }
    /** @} */

    /**
     * @name Specific methods to access AlphabetState
     * @{
     */
    /**
     * @brief Get a state at a position in the alphabet_ vector.
     *
     * This method must be overloaded in specialized classes to send back
     * a reference of the corect type.
     *
     * @param stateIndex The index of the state in the alphabet_ vector.
     * @throw IndexOutOfBoundsException If the index is invalid.
     */
    virtual AlphabetState& getStateAt(size_t stateIndex) throw (IndexOutOfBoundsException);
    
    /**
     * @brief Get a state at a position in the alphabet_ vector.
     *
     * This method must be overloaded in specialized classes to send back
     * a reference of the corect type.
     *
     * @param stateIndex The index of the state in the alphabet_ vector.
     * @throw IndexOutOfBoundsException If the index is invalid.
     */
    virtual const AlphabetState& getStateAt(size_t stateIndex) const throw (IndexOutOfBoundsException);

    /**
     * @brief Get a state by its letter.
     *
     * This method must be overloaded in specialized classes to send back
     * a reference of the corect type.
     *
     * @param letter The letter of the state to find.
     * @throw BadCharException If the letter is not in the Alphabet.
     */
    const AlphabetState& getState(const std::string& letter) const throw (BadCharException);

    AlphabetState& getState(const std::string& letter) throw (BadCharException);
    
    /**
     * @brief Get a state by its num.
     *
     * This method must be overloaded in specialized classes to send back
     * a reference of the corect type.
     *
     * @param num The num of the state to find.
     * @throw BadIntException If the num is not in the Alphabet.
     */
    const AlphabetState& getState(int num) const throw (BadIntException);

    AlphabetState& getState(int num) throw (BadIntException);

    int getIntCodeAt(size_t stateIndex) const throw (IndexOutOfBoundsException) {
      return getStateAt(stateIndex).getNum();
    }

    const std::string& getCharCodeAt(size_t stateIndex) const throw (IndexOutOfBoundsException) {
      return getStateAt(stateIndex).getLetter();
    }

    size_t getStateIndex(int state) const throw (BadIntException);
    
    size_t getStateIndex(const std::string& state) const throw (BadCharException);
    /** @} */

  protected:
    /**
     * @brief Add a state to the Alphabet.
     *
     * @param st The state to add.
     * @throw Exception If a wrong alphabet state is provided.
     */
    virtual void registerState(AlphabetState* st) throw (Exception);
    
    /**
     * @brief Set a state in the Alphabet.
     *
     * @param pos The index of the state in the alphabet_ vector.
     * @param st The new state to put in the Alphabet.
     * @throw Exception If a wrong alphabet state is provided.
     * @throw IndexOutOfBoundsException If an incorrect index is provided.
     */
    virtual void setState(size_t pos, AlphabetState* st) throw (Exception, IndexOutOfBoundsException);
    
    /**
     * @brief Resize the private alphabet_ vector.
     *
     * @param size The new size of the Alphabet.
     */
    void resize(size_t size) { alphabet_.resize(size); }
    
    /**
     * @brief Re-update the maps using the alphabet_ vector content.
     */
    void remap() {
      letters_.clear();
      nums_.clear();
      for (size_t i = 0; i < alphabet_.size(); ++i) {
        updateMaps_(i, *alphabet_[i]);
      }
    }

    unsigned int getStateCodingSize() const { return 1; }

    bool equals(const Alphabet& alphabet) const {
      return getAlphabetType() == alphabet.getAlphabetType();
    }
  };

} //end of namespace bpp.

#endif // _ABSTRACTALPHABET_H_

