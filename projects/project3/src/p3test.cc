#include<iostream>
#include<fstream>
#include<cmath>
#include<string>

#include"MatrixLib.hh"
#include"vec3.h"
#include"planet.h"
#include"ssystem.h"

using namespace std;

/* Global variables */
//const double GM = 4.*3.14159265*3.14159265;
const double GM = 39.4784176;

/* Function Declarations */
void Euler(int n,double years,double **);
void Verlet(int n,double years,double **);
void SunEarth();
void SunJupiterEarth();
int main();

/* Function Definitions */

int main()
{
  SunEarth();

  return 1;
}

/** Solve the equations of motion for the Sun-Earth
 ** system using two different methods: Euler and Verlet
**/
void SunEarth()
{
  int n=10000; double yr = 10.;
  // initial positions and velocities
  double **PV = AllocateMatrix(4,n);
  double **VV = AllocateMatrix(4,n);
//  PV[0][0]=VV[0][0]=1.; PV[1][0]=VV[1][0] = 0.;
//  PV[2][0]=VV[2][0] = 365.25*(-3.158e-3); PV[3][0]=VV[3][0] = 365.25*(-1.701e-2);//[AU/yr]
  PV[0][0]=VV[0][0]=0.; PV[1][0]=VV[1][0] = 1.;
  PV[2][0]=VV[2][0] = 365.25*(-0.017301); PV[3][0]=VV[3][0] = 365.25*(0.);//[AU/yr]

  Euler(n,yr,PV);
  Verlet(n,yr,VV);

  FileMatrix("../Benchmark/SEsystemEuler.out",PV,4,n);
  FileMatrix("../Benchmark/SEsystemVerlet.out",VV,4,n);

  DeallocateMatrix(PV,4,n);
  DeallocateMatrix(VV,4,n);

}


void Euler(int n,double years,double ** xyvxvy)
{
/* The matrix xyvxvy will hold the positons and velocities
 * at each time step (i) with the following convention:
 * x -> xyvxvy[0][i]
 * y -> xyvxvy[1][i]
 * vx-> xyvxvy[2][i]
 * vy-> xyvxvy[3][i]
*/
  double h = years/((double) n);
  double NumeratorCoefficient = GM*h;
  for(int i=0;i<n-1;i++)
  {
    double ri3 = sqrt(xyvxvy[0][i]*xyvxvy[0][i] + xyvxvy[1][i]*xyvxvy[1][i]);
    ri3 = ri3*ri3*ri3;
    xyvxvy[0][i+1] = xyvxvy[0][i] + xyvxvy[2][i]*h;
    xyvxvy[1][i+1] = xyvxvy[1][i] + xyvxvy[3][i]*h;
    xyvxvy[2][i+1] = xyvxvy[2][i] - NumeratorCoefficient*xyvxvy[0][i]/ ri3;
    xyvxvy[3][i+1] = xyvxvy[3][i] - NumeratorCoefficient*xyvxvy[1][i]/ ri3;
  }

}

void Verlet(int n,double years,double ** xyvxvy)
{
/* The matrix xyvxvy will hold the positons and velocities
 * at each time step (i) with the following convention:
 * x -> xyvxvy[0][i]
 * y -> xyvxvy[1][i]
 * vx-> xyvxvy[2][i]
 * vy-> xyvxvy[3][i]
*/
  double h = years/((double) n);
  double hGMOver2 = h*GM/2.; double hhGMOver2 = h*h*GM/2.;
  for(int i=0;i<n-1;i++)
  {
    double ri  = sqrt(xyvxvy[0][i]*xyvxvy[0][i] + xyvxvy[1][i]*xyvxvy[1][i]);
    double xri = xyvxvy[0][i] / ri;
    double yri = xyvxvy[1][i] / ri;
//    printf("%8.2f %8.2f %8.2f\n",ri,xri,yri);
    // positions
    xyvxvy[0][i+1] = xyvxvy[0][i] + h*xyvxvy[2][i] - hhGMOver2*xri;// should be minus on last term??
    xyvxvy[1][i+1] = xyvxvy[1][i] + h*xyvxvy[3][i] - hhGMOver2*yri;// should be minus on last term??
    double rii = sqrt(xyvxvy[0][i+1]*xyvxvy[0][i+1] + xyvxvy[1][i+1]*xyvxvy[1][i+1]);
    double xrii = xyvxvy[0][i+1] / rii;
    double yrii = xyvxvy[1][i+1] / rii;
//    printf("%8.2f %8.2f %8.2f\n",rii,xrii,yrii);
    xyvxvy[2][i+1] = xyvxvy[2][i] - hGMOver2*(xrii + xri);// should be minus on last term??
    xyvxvy[3][i+1] = xyvxvy[3][i] - hGMOver2*(yrii + yri);// should be minus on last term??
  }

}


