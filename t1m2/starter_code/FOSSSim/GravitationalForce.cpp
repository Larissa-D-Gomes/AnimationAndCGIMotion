#include "GravitationalForce.h"

void GravitationalForce::addEnergyToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, scalar& E )
{
  assert( x.size() == v.size() );
  assert( x.size() == m.size() );
  assert( x.size()%2 == 0 );
  assert( m_particles.first >= 0 );  assert( m_particles.first < x.size()/2 );
  assert( m_particles.second >= 0 ); assert( m_particles.second < x.size()/2 );

  // Your code goes here!
  VectorXs posp1 = x.segment<2>(2 * m_particles.first);
  VectorXs posp2 = x.segment<2>(2 * m_particles.second);
  
  scalar l = (posp1 - posp2).norm();
  
  // U = -G*m1*m2/l
  E -= m_G * m(2 * m_particles.first) * m(2 * m_particles.second) / l; 
}

void GravitationalForce::addGradEToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, VectorXs& gradE )
{
  assert( x.size() == v.size() );
  assert( x.size() == m.size() );
  assert( x.size() == gradE.size() );
  assert( x.size()%2 == 0 );
  assert( m_particles.first >= 0 );  assert( m_particles.first < x.size()/2 );
  assert( m_particles.second >= 0 ); assert( m_particles.second < x.size()/2 );

  // Your code goes here!
  VectorXs posp1 = x.segment<2>(2 * m_particles.first);
  VectorXs posp2 = x.segment<2>(2 * m_particles.second);
  scalar l = (posp1 - posp2).norm();
  
  // unit vector
  VectorXs n = (posp1 - posp2) / l;
  
  // gradE = (G*m1*m2/l^2)*n^
  VectorXs grad = (m_G * m(2 * m_particles.first) * m(2 * m_particles.second) * n / (l * l)); 
  gradE.segment<2>(2 * m_particles.first) += grad; 
  gradE.segment<2>(2 * m_particles.second) -= grad;
  
}
