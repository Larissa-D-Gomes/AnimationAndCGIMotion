#include "SpringForce.h"

void SpringForce::addEnergyToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, scalar& E )
{
  assert( x.size() == v.size() );
  assert( x.size()%2 == 0 );
  assert( m_endpoints.first >= 0 );  assert( m_endpoints.first < x.size()/2 );
  assert( m_endpoints.second >= 0 ); assert( m_endpoints.second < x.size()/2 );

  // Your code goes here!
  
  // Calculate l
  VectorXs posp1 = x.segment<2>(2 * m_endpoints.first);
  VectorXs posp2 = x.segment<2>(2 * m_endpoints.second);
  scalar l = (posp1 - posp2).norm();
  
  // U = 1/2 * k (l - l0)^2
  E += 0.5 * m_k * (l - m_l0) * (l - m_l0); 
  
}


void SpringForce::addGradEToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, VectorXs& gradE )
{
  assert( x.size() == v.size() );
  assert( x.size() == gradE.size() );
  assert( x.size()%2 == 0 );
  assert( m_endpoints.first >= 0 );  assert( m_endpoints.first < x.size()/2 );
  assert( m_endpoints.second >= 0 ); assert( m_endpoints.second < x.size()/2 );

  // Your code goes here!
  
  // Calculate l
  VectorXs posp1 = x.segment<2>(2 * m_endpoints.first);
  VectorXs posp2 = x.segment<2>(2 * m_endpoints.second);
  scalar l = (posp1 - posp2).norm();
  
  // unit vector
  VectorXs n = (posp1 - posp2) / l;
  
  // gradE = k(l - l0) n^
  VectorXs grad = m_k * (l - m_l0) * n;
  
  // Damping gradE = m_b n^ . (vi - vj) n^
  VectorXs sub = v.segment<2>(2 * m_endpoints.first) - v.segment<2>(2 * m_endpoints.second);
  grad += (m_b * n).dot(sub) * n;
  
  gradE.segment<2>(2 * m_endpoints.first) += grad; 
  gradE.segment<2>(2 * m_endpoints.second) -= grad; 
}


