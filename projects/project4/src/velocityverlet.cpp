#include "velocityverlet.h"
#include "system.h"
#include "atom.h"
#include <iostream>
#include <iomanip>

using namespace std;

void VelocityVerlet::integrate(System &system, double dt)
{
    // precalculate stuff
    double hOver2  = 0.5*dt;
    double hhOver2 = 0.5*dt*dt;
    // save old forces
//    double oldForces[3];
    vec3 oldForces;
/*
    cout << "Atom 2 force before VV:" << endl;
    cout << "  " << setw(14) << system.m_atoms[1]->force[0] <<
         setw(14) << system.m_atoms[1]->force[1] <<
         setw(14) << system.m_atoms[1]->force[2] << endl;
*/
    if(m_firstStep) {
        system.calculateForces();
        m_firstStep = false;
    }

    for(Atom *atom : system.atoms()) {
//	for(int idx=0;idx<3;idx++){
//	  oldForces[idx]=atom->force[idx];
//          atom->position[idx] += atom->velocity[idx]*dt + hhOver2*atom->force[idx]/atom->mass();
//        atom->velocity += atom->force*0.5*dt/atom->mass();
//        }

      atom->velocity += atom->force*0.5*dt/atom->mass();
      atom->position += atom->velocity*dt;
//	oldForces.zeros();
//	oldForces=atom->force;
//        atom->position += atom->velocity*dt + hhOver2*atom->force/atom->mass();

    }

/*
    cout << "Atom 2 updated forces VV:" << endl;
    cout << "  " << setw(14) << system.m_atoms[1]->force[0] <<
         setw(14) << system.m_atoms[1]->force[1] <<
         setw(14) << system.m_atoms[1]->force[2] << endl;
*/

    system.applyPeriodicBoundaryConditions();
    system.calculateForces(); // New positions, recompute forces

    for(Atom *atom : system.atoms()) {
//      for(int idx=0;idx<3;idx++){
//        atom->velocity[idx] += hOver2*(oldForces[idx] + atom->force[idx])/atom->mass();
//      }
//      atom->velocity += hOver2*(oldForces + atom->force)/atom->mass();
      atom->velocity += atom->force*0.5*dt/atom->mass();
    }
}
