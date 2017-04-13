#include "ssystem.h"
#include "planet.h"
#include "MatrixLib.hh"
//#include <iostream>
#include <cmath>
#include "time.h"

ssystem::ssystem()
{
  total_planets  = 0;
  radius         = 100.;
  total_mass     = 0.;
  G              = 4.*M_PI*M_PI;
  totalKinetic   = 0.;
  totalPotential = 0.;
}


ssystem::ssystem(double rradius)
{
  total_planets  = 0;
  radius         = rradius;
  total_mass     = 0.;
  G              = 4.*M_PI*M_PI;
  totalKinetic   = 0.;
  totalPotential = 0.;
}


void ssystem::add(planet newplanet)
{
  total_planets += 1;
  total_mass    += newplanet.mass;
  all_planets.push_back(newplanet);
}

/*
void ssystem::GravitationalConstant()
{
}
*/

void ssystem::VVerlet(int nsteps,double years,int printNsteps,double epsilon)
{
  double h = years/((double) nsteps);
  double hOver2 = 0.5*h; double hhOver2 = 0.5*h*h;
  double **Positions = AllocateMatrix(nsteps,3*total_planets);

  // add initial positions to positions array
  if(total_planets>0)
  {
    for(int pi=0;pi<total_planets;pi++)
    {
      planet &p = all_planets[pi];
      for(int j=0;j<3;j++){ Positions[0][pi*3+j]=p.position[j]; }
    }
  }
  else cout << "No planets!!" << endl;

  double ** acceleration     = AllocateMatrix(total_planets,3);
  double ** acceleration_new = AllocateMatrix(total_planets,3);

  for(int n=1;n<nsteps;n++)
  {
    for(int pi=0;pi<total_planets;pi++)
    {
      planet &current = all_planets[pi];// load current planet
      for(int i=0;i<3;i++)
      {// reset accelerations
        acceleration[pi][i]=0.;
	acceleration_new[pi][i]=0.;
      }
//      ax = ay = az = aax = aay = aaz = 0.;// reset accelerations

      // Forces on current
      for(int pii=pi+1;pii<total_planets;pii++)
      {
        planet &other = all_planets[pii];
	GravAccel(current,other,acceleration[pi],epsilon);
//	GravAccel(current,other,ax,ay,az,epsilon);
      }

      // New position for current
      for(int i=0;i<3;i++)
      { current.position[i] += current.velocity[i]*h + hhOver2*acceleration[pi][i]; }

      // New forces on current
      for(int pii=pi+1;pii<total_planets;pii++)
      {
        planet &other = all_planets[pii];
	GravAccel(current,other,acceleration_new[pi],epsilon);
//	GravAccel(current,other,ax,ay,az,epsilon);
      }

      // New velocity for current
      for(int i=0;i<3;i++)
      { current.velocity[i] += hOver2*(acceleration[pi][i] + acceleration_new[pi][i]); }

      for(int j=0;j<3;j++){Positions[n][pi*3+j]=current.position[j];}
    }

//    FileMatrix("../Benchmark/ooSystem.out",Positions,printNsteps,3*total_planets);
//    FileMatrix("../../Benchmark/ooSystem.out",Positions,printNsteps,3*total_planets);

  }

}

void ssystem::GravAccel(planet &current,planet &other, double *a,double epsilon)
{
  double relative=0.;
  double r=0.;
  double smoothing = epsilon*epsilon*epsilon;

  for(int i=0;i<3;i++)
  {
    relative = current.position[i]-other.position[i];
    r        = current.distance(other);
    a[i] -= this->G*other.mass*relative/((r*r*r) + smoothing);
  }

}
/*
void GravAccel(planet &current,planet &other, double &ax,double &ay,double &az,double epsilon)
{
  double rdist[3];

  for(int i=0;i<3;i++){rdist[i] = current.position[i]-other.position[i];}
  double r = current.distance(other);

  ax -= this->G*other.mass*rdist[0]/((r*r*r) + smoothing);
  ay -= this->G*other.mass*rdist[1]/((r*r*r) + smoothing);
  az -= this->G*other.mass*rdist[2]/((r*r*r) + smoothing);

}
*/
