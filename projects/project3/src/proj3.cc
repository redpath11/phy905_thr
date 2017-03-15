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
void Verlet(int n,double years,double **);
void ThreeBodyEuler(int n,double years,double **);
void SunEarth();
void SunJupiterEarth();
int main();


/* Function Definitions */

int main()
{
//  SunEarth();
  SunJupiterEarth();

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

void SunJupiterEarth()
{
  int n=1000; double yr = 4.;
  double ** PV = AllocateMatrix(8,n);
  PV[0][0]=0.; PV[1][0]=1.;// Earth
  PV[2][0]=365.25*(-0.017301); PV[3][0]=365.25*(0.);//[AU/yr]

  PV[4][0]=0.; PV[5][0]=5.2;//Jupiter
  PV[6][0]=365.25*(-0.007179); PV[7][0]=365.25*(0.);//[AU/yr]

  ThreeBodyEuler(n,yr,PV);

  FileMatrix("../Benchmark/SJEsystemEuler.out",PV,8,n);

  DeallocateMatrix(PV,8,n);

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

void ThreeBodyEuler(int n,double years,double ** xyvxvy)
{
/* The matrix xyvxvy will hold the positons and velocities
 * at each time step (i) with the following convention:
 * ----- Earth -----
 * x -> xyvxvy[0][i]
 * y -> xyvxvy[1][i]
 * vx-> xyvxvy[2][i]
 * vy-> xyvxvy[3][i]
 * ---- Jupiter ----
 * x -> xyvxvy[4][i]
 * y -> xyvxvy[5][i]
 * vx-> xyvxvy[6][i]
 * vy-> xyvxvy[7][i]
*/
  double h = years/((double) n);
  double NumeratorCoefficient = GM*h;
  double MJMSunRatio = 9.5458e-4; double MEMSunRatio = 3.0e-6;
  for(int i=0;i<n-1;i++)
  {
    // Earth-Sun distance
    double ri3 = sqrt(xyvxvy[0][i]*xyvxvy[0][i] + xyvxvy[1][i]*xyvxvy[1][i]);
    ri3 = ri3*ri3*ri3;
    double rJi3 = sqrt(xyvxvy[4][i]*xyvxvy[4][i] + xyvxvy[5][i]*xyvxvy[5][i]);
    rJi3 = rJi3*rJi3*rJi3;
    // Earth-Jupiter distance
    double rEJx = xyvxvy[0][i] - xyvxvy[4][i];
    double rEJy = xyvxvy[1][i] - xyvxvy[5][i];
    double rEJ3 = sqrt(rEJx*rEJx + rEJy*rEJy);
    rEJ3 = rEJ3*rEJ3*rEJ3;
    xyvxvy[0][i+1] = xyvxvy[0][i] + xyvxvy[2][i]*h;
    xyvxvy[1][i+1] = xyvxvy[1][i] + xyvxvy[3][i]*h;
    xyvxvy[2][i+1] = xyvxvy[2][i] - NumeratorCoefficient*(xyvxvy[0][i]/ ri3 + MJMSunRatio*rEJx/rEJ3);
    xyvxvy[3][i+1] = xyvxvy[3][i] - NumeratorCoefficient*(xyvxvy[1][i]/ ri3 + MJMSunRatio*rEJy/rEJ3);

    xyvxvy[4][i+1] = xyvxvy[4][i] + xyvxvy[6][i]*h;
    xyvxvy[5][i+1] = xyvxvy[5][i] + xyvxvy[7][i]*h;
    xyvxvy[6][i+1] = xyvxvy[6][i] - NumeratorCoefficient*(xyvxvy[4][i]/rJi3 + MEMSunRatio*rEJx/rEJ3);
    xyvxvy[7][i+1] = xyvxvy[7][i] - NumeratorCoefficient*(xyvxvy[5][i]/rJi3 + MEMSunRatio*rEJx/rEJ3);
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
    // positions
    xyvxvy[0][i+1] = xyvxvy[0][i] + h*xyvxvy[2][i] - hhGMOver2*xri;// should be minus on last term??
    xyvxvy[1][i+1] = xyvxvy[1][i] + h*xyvxvy[3][i] - hhGMOver2*yri;// should be minus on last term??
    double rii = sqrt(xyvxvy[0][i+1]*xyvxvy[0][i+1] + xyvxvy[1][i+1]*xyvxvy[1][i+1]);
    double xrii = xyvxvy[0][i+1] / rii;
    double yrii = xyvxvy[1][i+1] / rii;
    xyvxvy[2][i+1] = xyvxvy[2][i] - hGMOver2*(xrii + xri);// should be minus on last term??
    xyvxvy[3][i+1] = xyvxvy[3][i] - hGMOver2*(yrii + yri);// should be minus on last term??
  }

}
