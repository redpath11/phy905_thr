#include "system.h"
#include "velocityverlet.h"
#include "lennardjones.h"
#include "statisticssampler.h"
#include "unitconverter.h"
#include "math/random.h"

System::System()
{

}

System::~System()
{
    for(Atom *atom : m_atoms) {
        delete atom;
    }
    m_atoms.clear();
}

void System::applyPeriodicBoundaryConditions() {
  // Read here: http://en.wikipedia.org/wiki/Periodic_boundary_conditions#Practical_implementation:_continuity_and_the_minimum_image_convention
  // Loop over atoms and check positions relative to system center
  for(Atom *atom : m_atoms)
  {
    for(int j=0;j<3;j++)
    {
      if(atom->position[j]<-0.5*m_systemSize[j])
      {atom->position[j]=atom->position[j] + m_systemSize[j];}
      if(atom->position[j]>=0.5*m_systemSize[j])
      {atom->position[j]=atom->position[j] - m_systemSize[j];}
    }
  }
}

void System::removeTotalMomentum() {
    // Find the total momentum and remove momentum equally on each atom so the total momentum becomes zero.
  double Ptotal[3]={0.,0.,0.};
  for(Atom *atom : m_atoms)
  {
    for(int j=0;j<3;j++)
    { Ptotal[j] += atom->mass()*atom->velocity[j]; }
  }

  for(Atom *atom : m_atoms)
  {
    for(int j=0;j<3;j++)
    { atom->velocity[j]
      
    }
  }

}

void System::createFCCLattice(int numberOfUnitCellsEachDimension, double latticeConstant, double temperature) {
    // You should implement this function properly. Right now, 100 atoms are created uniformly placed in the system of size (10, 10, 10).

/*
    for(int i=0; i<100; i++) {
        Atom *atom = new Atom(UnitConverter::massFromSI(6.63352088e-26));
        double x = Random::nextDouble(0, 10); // random number in the interval [0,10]
        double y = Random::nextDouble(0, 10);
        double z = Random::nextDouble(0, 10);
        atom->position.set(x,y,z);
        atom->resetVelocityMaxwellian(temperature);
        m_atoms.push_back(atom);
    }
    setSystemSize(vec3(10, 10, 10)); // Remember to set the correct system size!
*/


}

void System::calculateForces() {
    for(Atom *atom : m_atoms) {
        atom->resetForce();
    }
    m_potential.calculateForces(*this); // this is a pointer, *this is a reference to this object
}

void System::calculateMass() {
  for(Atom *atom : m_atoms)
  { m_mass += atom->mass();}
}

void System::step(double dt) {
    m_integrator.integrate(*this, dt);
    m_steps++;
    m_time += dt;
}
