#include "SimpleGravityForce.h"

void SimpleGravityForce::addEnergyToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, scalar& E )
{
    assert( x.size() == v.size() );
    assert( x.size() == m.size() );
    assert( x.size()%2 == 0 );

    // Your code goes here!
    scalar imax = x.size()/2;

    for(int i = 0; i < imax; i++)
    {
        // E = -m*g*x
        E -= m(2*i) * m_gravity.dot(x.segment<2>(2 * i));
    }
}

void SimpleGravityForce::addGradEToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, VectorXs& gradE )
{
    assert( x.size() == v.size() );
    assert( x.size() == m.size() );
    assert( x.size() == gradE.size() );
    assert( x.size()%2 == 0 );
    
    // Your code goes here!
    scalar imax = x.size()/2;
    scalar g = m_gravity(1);
    
    for(int i = 0; i < imax; i++)
    {
        // gradE = -m*g
        gradE(2*i + 1) -= m(2*i) * g;
    }
}

