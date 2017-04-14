#ifndef SSYSTEM_H
#define SSYSTEM_H
#include "planet.h"
#include <vector>
#include <fstream>
using std::vector;

class ssystem
{
public:
    friend class planet;

    // properties
    double radius,total_mass,G,c;
    int total_planets;
    vector<planet> all_planets;
    double totalKinetic;
    double totalPotential;
    double totalAngularMomentum;

    // constants

    // initializers
    ssystem();
    ssystem(double radi);

    // functions
    void add(planet newplanet);
    void Euler(int nsteps, double years, int printNsteps, double epsilon);
    void VVerlet(int nsteps, double years, int printNsteps, double epsilon);
    void HgVerlet(int nsteps, double years, int printNsteps, double epsilon);
//    void GravAccel(planet &current,planet &other, double &ax,double &ay,double &az,double epsilon);
    void GravAccel(planet &current,planet &other, double *a,double epsilon);
    void HgGravAccel(planet &current,planet &other,double *a,double l);
    void GravForce(planet &current, planet &other, double *f, double epsilon);
    void KineticEnergySystem();
    void PotentialEnergySystem(double epsilon);
    void PotentialEnergy2body(double epsilon);
    void AngularMomentumSystem();
    bool Bound(planet p);

};

#endif // SSYSTEM_H
