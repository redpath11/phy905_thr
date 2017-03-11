#include<iostream>
#include<fstream>
#include<cmath>
#include<string>

#include"MatrixLib.hh"

using namespace std;

/* Global variables */
const double GM = 4.*3.14159265*3.14159265;

/* Function Declarations */
void Euler(int n,double years,double **);
void SunEarth();
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
  int n=100; double yr = 1.;
  // initial positions and velocities
  double **PV = AllocateMatrix(4,n);
  PV[0][0]=1.; PV[1][0] = 0.;
  PV[2][0] = 365.25*(-3.158e-3); PV[3][0] = 365.25*(-1.701e-2);//[AU/yr]

  Euler(n,yr,PV);

  FileMatrix("../Benchmark/SEsystemEuler.out",PV,4,n);

  DeallocateMatrix(PV,4,n);

}

void Euler(int n,double years,double ** xyvxvy)
{
/* The matrix xyvxvy will hold the positons and velocities
 * at each time step with the following convention:
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
