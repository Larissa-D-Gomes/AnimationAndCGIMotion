#include "SymplecticEuler.h"

bool SymplecticEuler::stepScene( TwoDScene& scene, scalar dt )
{
  // Your code goes here!
  
  VectorXs& x = scene.getX();
  VectorXs& v = scene.getV();
    
  int imax = scene.getNumParticles();
  int xsize = x.size();
  
  VectorXs F(xsize);
  F.setZero(xsize);
  scene.accumulateGradU(F);

  for(int i = 0; i < imax; i++)
  {
     if(scene.isFixed(i))
     {
         F(2*i) = 0;
         F(2*i + 1) = 0;
     }
  }
  
  // New velocity
  v = v + dt * scene.getM().cwiseInverse().cwiseProduct(-F);
  
  // New configuration
  x = x + dt * v;
  
  return true;
}






