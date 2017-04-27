#include "lennardjones.h"
#include "system.h"
#include <math.h>
#include <iostream>
#include <iomanip>

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
    double epsilon24 = 24.*m_epsilon;
    double rcut = 2.5*m_sigma;
    double rcut2= rcut*rcut;
    double sigma6 = pow(m_sigma,6.);
    double sigma12= sigma6*sigma6;
    double potentialAtRcut = 4*m_epsilon*(sigma12*pow(rcut,-12) - sigma6*pow(rcut, -6));
    double sizex = system.systemSize()[0];
    double sizey = system.systemSize()[1];
    double sizez = system.systemSize()[2];


/**/
    int natoms = system.getNumberOfAtoms();
    int i=0;

      for(Atom *atom : system.atoms()) {
        double fx = 0.; double fy = 0.; double fz = 0.; double pe = 0.;
	double x = atom->position[0];
	double y = atom->position[1];
	double z = atom->position[2];

      for(int j=i+1;j<natoms;j++) {
        Atom *other = system.atoms()[j];

	double dx = x - other->position[0];
	double dy = y - other->position[1];
	double dz = z - other->position[2];

	if ((dx) <= -sizex/2.) dx += sizex;
	if ((dx) >   sizex/2.) dx -= sizex;

	if ((dy) <= -sizey/2.) dy += sizey;
	if ((dy) >   sizey/2.) dy -= sizey;

	if ((dz) <= -sizez/2.) dz += sizez;
	if ((dz) >   sizez/2.) dz -= sizez;

	double dr2 = dx*dx + dy*dy + dz*dz;

        if(dr2 < rcut2) {
	  double oneOverDr2 = 1./dr2;
	  double oneOverDr6 = oneOverDr2*oneOverDr2*oneOverDr2;
	  double oneOverDr12= oneOverDr6*oneOverDr6;
	  double fr = epsilon24*oneOverDr2*(2.*sigma12*oneOverDr12 - sigma6*oneOverDr6);

	  fx += fr*dx; fy += fr*dy; fz += fr*dz;

	  other->force[0] -= fr*dx;
	  other->force[1] -= fr*dy;
	  other->force[2] -= fr*dz;

	  pe += 4.*m_epsilon*(sigma12*oneOverDr12 - sigma6*oneOverDr6) - potentialAtRcut;
	}
      }

	atom->force[0] += fx; atom->force[1] += fy; atom->force[2] += fz;
	m_potentialEnergy += pe;
      i++;
    }
}
