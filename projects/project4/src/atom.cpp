#include "atom.h"
#include "math/random.h"
#include <cmath>

Atom::Atom(double mass) :
    m_mass(mass)
{
    
}

void Atom::resetForce()
{
    force.zeros();
}

void Atom::resetVelocityMaxwellian(double temperature)
{
    // Resetting the velocity according to a Maxwell-Boltzmann distribution (see http://en.wikipedia.org/wiki/Maxwell%E2%80%93Boltzmann_distribution )
    double boltzmannConstant = 1.0; // In these units, the boltzmann constant equals 1
    double standardDeviation = sqrt(boltzmannConstant*temperature/m_mass);
    velocity.randomGaussian(0, standardDeviation);
}

double Atom::distance(Atom *otherAtom)
{
  double x1,y1,z1,x2,y2,z2,xx,yy,zz;

  x1 = this->position[0];
  y1 = this->position[1];
  z1 = this->position[2];

  x2 = otherAtom->position[0];
  y2 = otherAtom->position[1];
  z2 = otherAtom->position[2];

  xx = x1-x2;
  yy = y1-y2;
  zz = z1-z2;

  return sqrt(xx*xx + yy*yy + zz*zz);

}
