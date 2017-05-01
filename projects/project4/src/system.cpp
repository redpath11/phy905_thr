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
  /**/
  for(Atom *atom : m_atoms)
  {
    for(int j=0;j<3;j++)
    {
      if(atom->position[j]<0.) {
//        double multiplier = floor(-1.*atom->position[j] / m_systemSize[j]);
//        atom->position[j] += m_systemSize[j] * multiplier;
//        atom->initialPosition[j] += m_systemSize[j]*multiplier;
        atom->position[j] += m_systemSize[j];
        atom->initialPosition[j] += m_systemSize[j];
      }
      if(atom->position[j]>=m_systemSize[j]) {
//        double multiplier = floor(atom->position[j] / m_systemSize[j]);
//        atom->position[j]-= m_systemSize[j] * multiplier;
//        atom->initialPosition[j] -= m_systemSize[j] * multiplier;
        atom->position[j]-= m_systemSize[j];
        atom->initialPosition[j] -= m_systemSize[j];
      }
    }
  }

}

void System::removeTotalMomentum() {
    // Find the total momentum and remove momentum equally on each atom so the total momentum becomes zero.
  vec3 Ptotal = vec3(0,0,0);
  double Mtotal = 0.;
  for(Atom *atom : m_atoms)
  { 
    Mtotal += atom->mass(); 
    Ptotal += atom->mass()*atom->velocity;
  }

  for(Atom *atom : m_atoms)
  {
    atom->velocity -= Ptotal/Mtotal;
  }

}

void System::removeTotalVelocity() {
  vec3 Vtotal = vec3(0,0,0);
//  double natoms = this->getNumberOfAtoms();
  for(Atom *atom : m_atoms) {Vtotal += atom->velocity;}
//  Vtotal /= natoms;
  for(Atom *atom : m_atoms) {atom->velocity -= Vtotal/m_atoms.size();}
}

void System::createFCCLattice(int numberOfUnitCellsEachDimension, double latticeConstant, double temperature) {
    // You should implement this function properly. Right now, 100 atoms are created uniformly placed in the system of size (10, 10, 10).

/*
//    for(int i=0; i<100; i++) {
        Atom *atom = new Atom(UnitConverter::massFromSI(6.63352088e-26));
        double x = Random::nextDouble(0, 10); // random number in the interval [0,10]
        double y = Random::nextDouble(0, 10);
        double z = Random::nextDouble(0, 10);
        atom->position.set(x,y,z);
        atom->resetVelocityMaxwellian(temperature);
	atom->velocity.set(0.,0.,100.);
        m_atoms.push_back(atom);
//    }
*/
    double ssize = ((double) numberOfUnitCellsEachDimension) * latticeConstant;
    double offset = 0.;
//    double offset = ceil(ssize) - ssize;
//    ssize = ceil(ssize);
    setSystemSize(vec3(ssize, ssize, ssize)); // Remember to set the correct system size!
//    setSystemSize(vec3(10, 10, 10)); // Remember to set the correct system size!
    // First, set up a single cell of four atoms
    const double b = latticeConstant/2.;
    const double x[4] = {0.,b ,0.,b };
    const double y[4] = {0.,b ,b ,0.};
    const double z[4] = {0.,0.,b ,b };
    for(int i=0;i<numberOfUnitCellsEachDimension;i++)
    {// place unit cells along x
      double bx = ((double) i)*b*2. + offset;
      for(int j=0;j<numberOfUnitCellsEachDimension;j++)
      {// place unit cells along y
        double by = ((double) j)*b*2. + offset;
	  for(int k=0;k<numberOfUnitCellsEachDimension;k++)
	  {// place unit cells along z
	    double bz = ((double) k)*b*2. + offset;
//            Atom *atom[4];
            for(int a=0;a<4;a++)
            {// place 4 atoms
              Atom *atom;
              atom = new Atom(UnitConverter::massFromSI(6.63352088e-26));
              atom->position.set(x[a]+bx,y[a]+by,z[a]+bz);
	      atom->initialPosition.set(x[a]+bx,y[a]+by,z[a]+bz);
              atom->resetVelocityMaxwellian(temperature);
              m_atoms.push_back(atom);
            }// END place 4 atoms
          }// END place unit cells along z
      }// END place unit cells along y
    }// END place unit cells along x




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
