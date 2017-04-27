#include "math/random.h"
#include "lennardjones.h"
#include "velocityverlet.h"
#include "system.h"
#include "statisticssampler.h"
#include "atom.h"
#include "io.h"
#include "unitconverter.h"
#include <iostream>
#include <iomanip>

using namespace std;

int main(int numberOfArguments, char **argumentList)
{
    int numberOfUnitCells = 5;
    double initialTemperature = UnitConverter::temperatureFromSI(300.0); // measured in Kelvin
    double latticeConstant = UnitConverter::lengthFromAngstroms(5.26); // measured in angstroms

    // If a first argument is provided, it is the number of unit cells
    if(numberOfArguments > 1) numberOfUnitCells = atoi(argumentList[1]);
    // If a second argument is provided, it is the initial temperature (measured in kelvin)
    if(numberOfArguments > 2) initialTemperature = UnitConverter::temperatureFromSI(atof(argumentList[2]));
    // If a third argument is provided, it is the lattice constant determining the density (measured in angstroms)
    if(numberOfArguments > 3) latticeConstant = UnitConverter::lengthFromAngstroms(atof(argumentList[3]));

    double dt = UnitConverter::timeFromSI(1e-15); // Measured in seconds.

    cout << "One unit of length is " << UnitConverter::lengthToSI(1.0) << " meters" << endl;
    cout << "One unit of velocity is " << UnitConverter::velocityToSI(1.0) << " meters/second" << endl;
    cout << "One unit of time is " << UnitConverter::timeToSI(1.0) << " seconds" << endl;
    cout << "One unit of mass is " << UnitConverter::massToSI(1.0) << " kg" << endl;
    cout << "One unit of temperature is " << UnitConverter::temperatureToSI(1.0) << " K" << endl;

    System system;
    system.createFCCLattice(numberOfUnitCells, latticeConstant, initialTemperature);
    // Set for argon
//    system.potential().setEpsilon(UnitConverter::temperatureFromSI(119.8));
    system.potential().setEpsilon(UnitConverter::energyFromSI(UnitConverter::kb*119.8));
    system.potential().setSigma(UnitConverter::lengthFromAngstroms(3.405));
    cout << "epsilon = " << system.potential().epsilon() << endl;
    cout << "sigma = " << system.potential().sigma() << endl;

    cout << "System of " << system.getNumberOfAtoms() << " atoms initialized." << endl;

    /* Remove net velocity */
//    system.removeTotalMomentum();// use this if the system contains different types of atoms
    system.removeTotalVelocity();// use this if the system contains just one type of atoms

    StatisticsSampler statisticsSampler;
    IO movie("movie.xyz"); // To write the state to file

    cout << setw(20) << "Timestep" <<
            setw(20) << "Time" <<
            setw(20) << "Temperature" <<
            setw(20) << "KineticEnergy" <<
            setw(20) << "PotentialEnergy" <<
            setw(20) << "TotalEnergy" << endl;

/**/
    for(int timestep=0; timestep<500; timestep++) {
//        cout << " ----- Step: " << timestep << " ----- " << endl << endl;
        system.step(dt);
	system.applyPeriodicBoundaryConditions();
        statisticsSampler.sample(system);
	statisticsSampler.saveToFile(system);
        if( timestep % 100 == 0 ) {
            // Print the timestep every 100 timesteps
            cout << setw(20) << system.steps() <<
                    setw(20) << system.time() <<
//                    setw(20) << statisticsSampler.temperature() <<
                    setw(20) << UnitConverter::temperatureToSI(statisticsSampler.temperature()) <<
                    setw(20) << statisticsSampler.kineticEnergy() <<
                    setw(20) << statisticsSampler.potentialEnergy() <<
                    setw(20) << statisticsSampler.totalEnergy() << endl;
        }
        movie.saveState(system);
    }


    movie.close();

/*
    for(int i=0;i<system.getNumberOfAtoms();i++)
    { cout << system.m_atoms[i]->position[0] << endl; }
*/

    return 0;
}
