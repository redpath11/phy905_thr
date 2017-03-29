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
    double radius,total_mass,G;
    int total_planets;
    vector<planet> all_planets;
    double totalKinetic;
    double totalPotential;

    // constants

    // initializers
    ssystem();
    ssystem(double radi);

    // functions
    void add(planet newplanet);
//    void addM(planet newplanet);
//    void GravitationalConstant();
//    void print_position(std::ofstream &output, int dimension, double time, int number);
//    void print_energy(std::ofstream &output, double time, double epsilon);
    void VVerlet(int nsteps, double years, int printNsteps, double epsilon);
//    void GravAccel(planet &current,planet &other, double &ax,double &ay,double &az,double epsilon);
    void GravAccel(planet &current,planet &other, double *a,double epsilon);
//    double **setup_matrix(int height, int width);
//    void delete_matrix(double **matrix);
//    void GravitationalForce(planet &current, planet &other, double &Fx, double &Fy, double &Fz, double epsilon);
//    void GravitationalForce_RK(double x_rel, double y_rel, double z_rel, double &Fx, double &Fy, double &Fz, double mass1, double mass2);
//    void KineticEnergySystem();
//    void PotentialEnergySystem(double epsilon);
//    double EnergyLoss();
//    bool Bound(planet OnePlanet);

};

#endif // SSYSTEM_H
