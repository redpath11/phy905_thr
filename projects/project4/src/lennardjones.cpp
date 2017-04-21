#include "lennardjones.h"
#include "system.h"
#include <math.h>
#include <iostream>

using namespace std;

double LennardJones::potentialEnergy() const
{
    return m_potentialEnergy;
}

double LennardJones::sigma() const
{
    return m_sigma;
}

void LennardJones::setSigma(double sigma)
{
    m_sigma = sigma;
}

double LennardJones::epsilon() const
{
    return m_epsilon;
}

void LennardJones::setEpsilon(double epsilon)
{
    m_epsilon = epsilon;
}

void LennardJones::calculateForces(System &system)
{
    m_potentialEnergy = 0; // Remember to compute this in the loop
    // take care of some coefficients
    double sigma6 = pow(m_sigma,6.);
    double FourEpsilon = 4.*m_epsilon;
    double coefficient = 24.*m_epsilon*sigma6;

    int natoms = system.getNumberOfAtoms();

    for(int i=0;i<natoms;i++)
    {
    
      Atom *current = system.m_atoms[i];
      cout << current->position[0] << endl;

      for(int j=i+1;j<natoms;j++)
      {
        Atom *other = system.m_atoms[j];
	double r        = current->distance(other);
//	cout << r << endl;
/**/
	for(int idx=0;idx<3;idx++)
	{
	  double relative = current->position[idx] - other->position[idx];
	  double r8       = pow(r,8.);
	  current->force[idx] -= coefficient*relative*(1.-(2.*sigma6/pow(r,6.)))/pow(r,8.);
	}
	double sigOverR = sigma6/pow(r,6.);
	m_potentialEnergy += FourEpsilon*sigOverR*(sigOverR - 1.);

      }

    }

}
