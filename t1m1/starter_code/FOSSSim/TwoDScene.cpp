#include "TwoDScene.h"
#include <iostream>

scalar TwoDScene::computeKineticEnergy() const
{
  VectorXs vmp =  m_v.transpose().cwiseProduct(m_m);
  scalar ke = 0.5 * vmp.cwiseProduct(getV())(0);

  return ke;
}

