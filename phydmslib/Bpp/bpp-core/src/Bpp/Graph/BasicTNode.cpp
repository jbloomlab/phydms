// 
// File:    BasicTNode.cpp
// Author:  Sylvain Gaillard
// Created: 14/01/2011 14:59:07
// 

/*
Copyright or © or Copr. Bio++ Development Team, (January 12, 2011)

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

#include "BasicTNode.h"

using namespace bpp;

#include <cstddef>

using namespace std;

BasicTNode::~BasicTNode() {
  if (father_) {
    father_->removeSon(this);
  }
  for (size_t i = 0 ; i < sons_.size() ; i++) {
    sons_[i]->removeFather();
  }
}

BasicTNode::BasicTNode(const BasicTNode& node):
  sons_(node.sons_),
  father_(node.father_)
{};

BasicTNode& BasicTNode::operator=(const BasicTNode& node) {
  sons_ = node.sons_;
  father_ = node.father_;
  return * this;
}

// Neighbors

const BasicTNode* BasicTNode::getNeighbor(int pos) const {
  if (pos < 0 || pos > static_cast<int>(sons_.size())) {
    throw IndexOutOfBoundsException("BasicTNode::getNeighbor() pos is out of bounds", static_cast<size_t>(pos), 0, sons_.size());
  }
  if (pos == 0)
    return father_;
  else
    return sons_[static_cast<size_t>(pos - 1)];
}

BasicTNode* BasicTNode::getNeighbor(int pos) {
  if (pos < 0 || pos > static_cast<int>(sons_.size())) {
    throw IndexOutOfBoundsException("BasicTNode::getNeighbor() pos is out of bounds", static_cast<size_t>(pos), 0, sons_.size());
  }
  if (pos == 0)
    return father_;
  else
    return sons_[static_cast<size_t>(pos - 1)];
}

const BasicTNode* BasicTNode::operator[](int i) const {
  if (i < 0) {
    return father_;
  } else {
    return sons_[static_cast<size_t>(i)];
  }
}

BasicTNode* BasicTNode::operator[](int i) {
  if (i < 0) {
    return father_;
  } else {
    return sons_[static_cast<size_t>(i)];
  }
}

// Fathers

const BasicTNode* BasicTNode::getFather(int pos) const {
  if (pos != 0) {
    throw IndexOutOfBoundsException("BasicTNode::getFather() pos must be 0 for TNode", static_cast<size_t>(pos), 0, 0);
  }
  return getFather();
}

BasicTNode* BasicTNode::getFather(int pos) {
  if (pos != 0) {
    throw IndexOutOfBoundsException("BasicTNode::getFather() pos must be 0 for TNode", static_cast<size_t>(pos), 0, 0);
  }
  return getFather();
}

const BasicTNode* BasicTNode::getFather() const {
  return father_;
}

BasicTNode* BasicTNode::getFather() {
  return father_;
}

bool BasicTNode::isFather(const BasicTNode* node) const {
  if (father_ == node)
    return true;
  return false;
}

void BasicTNode::addFather(BasicTNode* node) {
  if (!node)
    throw NullPointerException("BasicTNode::addFather() Empty node given as input");
  if (father_)
    throw Exception("BasicTNode::addFather() This node already has a father.");
  if (!isFather(node))
    father_ = node;
  if (!node->isSon(this))
    node->addSon(this);
}

BasicTNode* BasicTNode::removeFather() {
  if (hasFathers()) {
    BasicTNode* father = father_;
    father_ = 0;
    father->removeSon(this);
    return father;
  }
  return 0;
}

// Sons

const BasicTNode* BasicTNode::getSon(int pos) const {
  if (pos < 0 || pos > static_cast<int>(sons_.size()) - 1) {
    throw IndexOutOfBoundsException("BasicTNode::getSon() pos out of range", static_cast<size_t>(pos), 0, sons_.size() - 1);
  }
  return sons_[static_cast<size_t>(pos)];
}

BasicTNode* BasicTNode::getSon(int pos) {
  if (pos < 0 || pos > static_cast<int>(sons_.size()) - 1) {
    throw IndexOutOfBoundsException("BasicTNode::getSon() pos out of range", static_cast<size_t>(pos), 0, sons_.size() - 1);
  }
  return sons_[static_cast<size_t>(pos)];
}

bool BasicTNode::isSon(const BasicTNode* node) const {
  for (size_t i = 0 ; i < sons_.size() ; i++) {
    if (sons_[i] == node)
      return true;
  }
  return false;
}

void BasicTNode::addSon(BasicTNode* node) {
  if (!node)
    throw NullPointerException("BasicTNode::addSon() Empty node given as input.");
  if (!isSon(node))
    sons_.push_back(node);
  if (!node->isFather(this))
    node->addFather(this);
}

void BasicTNode::removeSon(BasicTNode* node) {
  if (!node)
    throw NullPointerException("BasicTNode::removeSon() Empty node given as input.");
  for (size_t i = 0 ; i < sons_.size() ; i++) {
    if (sons_[i] == node) {
      sons_.erase(sons_.begin() + static_cast<ptrdiff_t>(i));
      node->removeFather();
    }
  }
}

BasicTNode* BasicTNode::removeSon(int pos) {
  if (pos < 0 || pos > static_cast<int>(sons_.size() - 1))
    throw IndexOutOfBoundsException("BasicTNode::removeSon() pos out of bound", static_cast<size_t>(pos), 0, sons_.size() - 1);
  BasicTNode* node = sons_[static_cast<size_t>(pos)];
  sons_.erase(sons_.begin() + pos);
  node->removeFather();
  return node;
}
