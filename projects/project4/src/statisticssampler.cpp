#include "system.h"
#include "statisticssampler.h"
#include "lennardjones.h"
#include "unitconverter.h"
#include <iostream>
#include <iomanip>

using namespace std;

using std::ofstream; using std::cout; using std::endl;

StatisticsSampler::StatisticsSampler()
{

}

void StatisticsSampler::saveToFile(System &system)
{
    // Save the statistical properties for each timestep for plotting etc.
    // First, open the file if it's not open already
    if(!m_file.good()) {
        m_file.open("statistics.txt", ofstream::out);
        // If it's still not open, something bad happened...
        if(!m_file.good()) {
            cout << "Error, could not open statistics.txt" << endl;
            exit(1);
        }
    }
  /* Column ordering:
   * Timestep, Time, Temp, KE, PE, Total Energy, diffusion const
  */

  // Print out values here
  m_file << setw(20) << system.steps() <<
  	    setw(20) << UnitConverter::timeToSI(system.time()) <<
	    setw(20) << UnitConverter::temperatureToSI(m_temperature) <<
	    setw(20) << m_kineticEnergy <<
	    setw(20) << m_potentialEnergy <<
	    setw(20) << m_kineticEnergy+m_potentialEnergy <<
	    setw(20) << m_diffusion << endl;
}

void StatisticsSampler::sample(System &system)
{
    // Here you should measure different kinds of statistical properties and save it to a file.
    sampleKineticEnergy(system);
    samplePotentialEnergy(system);
    sampleTemperature(system);
    sampleDensity(system);
    sampleDiffusion(system);
//    saveToFile(system);
}

void StatisticsSampler::sampleKineticEnergy(System &system)
{
  // All in MD units
  m_kineticEnergy = 0; // Remember to reset the value from the previous timestep
  for(Atom *atom : system.atoms()) {
    m_kineticEnergy += 0.5*atom->mass()*atom->velocity.lengthSquared();
  }
}

void StatisticsSampler::samplePotentialEnergy(System &system)
{
    m_potentialEnergy = system.potential().potentialEnergy();
}

void StatisticsSampler::sampleTemperature(System &system)
{
  // instantaneous temperature
  // In MD units; k=1
  // Hint: reuse the kinetic energy that we already calculated
  double frac = 2./3.;
  double natoms = ((double) system.getNumberOfAtoms());
  m_temperature = frac*m_kineticEnergy/natoms;
}

void StatisticsSampler::sampleDensity(System &system)
{
  m_density = UnitConverter::massToSI(system.mass())/system.volume();
  

}

void StatisticsSampler::sampleDiffusion(System &system)
{
  m_diffusion = 0.;
  double displacement = 0.;
  for(Atom *atom : system.atoms()) {
    vec3 dr = vec3(0,0,0);
    dr      = atom->position - atom->initialPosition;
    displacement += dr.lengthSquared();
  }
  m_diffusion = displacement/system.atoms().size();
  m_diffusion /= (6.*system.time());

}
