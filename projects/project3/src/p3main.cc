#include <iostream>
#include <cmath>
#include <cstring>
#include <fstream>
#include "planet.h"
#include "ssystem.h"

using namespace std;

void SunEarth(int Nsteps,double time);
void FindEscapeVel();
bool EscapeVelocity(double vescape);
void SunJupiterEarth(int Nsteps,double time);
void SunMercury(int Nsteps,double time);
void TheSolSystem(int Nsteps,double time);

int main()
{
//SunEarth(10000,10.);
//  FindEscapeVel();
//  SunJupiterEarth(20000,20.);
  SunMercury(1010000000,101.);
//  SunMercury(10,0.01);
//  SunMercury(101000,101.);
//  TheSolSystem(30000,30.);
  return 1;
}

void SunEarth(int Nsteps,double time)
{
//  planet earth(3.0e-6,-0.97615,0.17058,-1.5302e-4,365.25*(-3.1581e-3),365.25*(-1.7014e-2),365.25*(8.5947e-8));
//  planet sun(1.,3.1815e-3,4.3976e-3,-1.5181e-4,365.25*(-3.3248e-6),365.25*(6.6209e-6),365.25*(7.1779e-8));

  planet earth(3.0e-6,0.,1.0,0.,365.25*(-0.017301),0.,0.);
//  planet jupiter(9.5458e-4,0.,5.45,0.,365.25*(-7.179e-3),0.,0.);
  planet sun(1.,0.,0.,0.,0.,0.,0.);
  ssystem se(10.);

  se.add(earth);
  se.add(sun);
//  se.Euler(Nsteps,time,Nsteps,0);
  se.VVerlet(Nsteps,time,Nsteps,0);
}


void FindEscapeVel()
{
  bool e=true;
  double testv=-8.0;
  while(e)
  {
    e=EscapeVelocity(testv);
    testv-=0.1;
  }
}


bool EscapeVelocity(double vescape)
{
  planet p1(3.0e-2,0.,1.,0.,vescape,0.,0.);
  planet sun(1.,0.,0.,0.,0.,0.,0.);

  ssystem escape(100.);
  escape.add(p1);
  escape.add(sun);
//  escape.VVerlet(1000,10.,1000,0);
  escape.KineticEnergySystem();
  escape.PotentialEnergySystem(0.);
  bool b = escape.Bound(escape.all_planets[0]);
  if(b){printf("Planet with v=%g [AU/yr] is bound\n",abs(vescape));}
  else {printf("Planet with v=%g [AU/yr] escapes\n",abs(vescape));}
  return b;
}


void SunJupiterEarth(int Nsteps,double time)
{
/**/
  planet earth(3.0e-6,0.,1.,0.,365.25*(-0.017301),0.,0.);
  planet jupiter(9.5458e-4,0.,5.45,0.,365.25*(-7.179e-3),0.,0.);
  planet sun(1.,0.,0.005428,0.,365.25*(6.9e-6),0.,0.);


  ssystem sej(10.);

  sej.add(earth);
  sej.add(jupiter);
  sej.add(sun);
  sej.VVerlet(Nsteps,time,Nsteps,0);

}


void SunMercury(int Nsteps,double time)
{
  planet mercury(1.66e-7,0.3075,0.,0.,0.,12.44,0.);
  planet sun(1.,0.,0.,0.,0.,0.,0.);

  ssystem sm(10.);

  sm.add(mercury);
  sm.add(sun);
  sm.HgVerlet(Nsteps,time,Nsteps,0);

}


void TheSolSystem(int Nsteps,double time)
{
  planet mercury(1.66e-7,3.5028e-01,3.3641e-02,-2.9606e-02,365.25*(-7.7990e-03),365.25*(2.9301e-02),365.25*(3.1089e-03));
  planet venus(2.4478e-6,-6.7689e-01,2.3416e-01,4.2243e-02,365.25*(-6.5732e-03),365.25*(-1.9248e-02),365.25*(1.1514e-04));
  planet earth(3.0e-6,-0.97615,0.17058,-1.5302e-4,365.25*(-3.1581e-3),365.25*(-1.7014e-2),365.25*(8.5947e-8));
  planet mars(3.227e-7,0.83434,1.2360,5.2585e-3,365.25*(-1.1071e-2),365.25*(9.0283e-3),365.25*(4.6076e-4));
  planet jupiter(9.5458e-4,-5.2429,-1.4907,1.2345e-1,365.25*(1.9762e-3),365.25*(-6.9012e-3),365.25*(-1.5539e-5));
  planet saturn(2.8582e-4,-1.5068,-9.9316,2.3265e-1,365.25*(5.3203e-3),365.25*(-8.5501e-4),365.25*(-1.9230e-4));
  planet uranus(4.3658e-5,1.8232e+1,8.0664,-2.0625e-1,365.25*(-1.62e-3),365.25*(3.4135e-3),365.25*(3.3653e-5));
  planet neptune(5.1503e-5,2.8407e+1,-9.4830,-4.5939e-1,365.25*(9.7305e-4),365.25*(2.9966e-3),365.25*(-8.3749e-5));
  planet pluto(6.583e-9,9.8750e+0,-3.1780e+1,5.4423e-1,365.25*(3.0732e-3),365.25*(2.8367e-4),365.25*(-9.2485e-4));
  planet sun(1.,3.1815e-3,4.3976e-3,-1.5181e-4,365.25*(-3.3248e-6),365.25*(6.6209e-6),365.25*(7.1779e-8));


  ssystem TheSystem(100.);
  TheSystem.add(mercury);
  TheSystem.add(venus);
  TheSystem.add(earth);
  TheSystem.add(mars);
  TheSystem.add(jupiter);
  TheSystem.add(saturn);
  TheSystem.add(neptune);
  TheSystem.add(uranus);
  TheSystem.add(pluto);
  TheSystem.add(sun);
  TheSystem.VVerlet(Nsteps,time,Nsteps,0);

}
