#include <iostream>
#include <cmath>
#include <cstring>
#include <fstream>
#include "planet.h"
#include "ssystem.h"

using namespace std;

int main()
{
  planet venus(2.4478e-6,-6.7689e-01,2.3416e-01,4.2243e-02,365.25*(-6.5732e-03),365.25*(-1.9248e-02),365.25*(1.1514e-04));
  planet earth(3.0e-6,-0.97615,0.17058,-1.5302e-4,365.25*(-3.1581e-3),365.25*(-1.7014e-2),365.25*(8.5947e-8));
  planet mars(3.227e-7,0.83434,1.2360,5.2585e-3,365.25*(-1.1071e-2),365.25*(9.0283e-3),365.25*(4.6076e-4));
  planet jupiter(9.5458e-4,-5.2429,-1.4907,1.2345e-1,365.25*(1.9762e-3),365.25*(-6.9012e-3),365.25*(-1.5539e-5));
  planet saturn(2.8582e-4,-1.5068,-9.9316,2.3265e-1,365.25*(5.3203e-3),365.25*(-8.5501e-4),365.25*(-1.9230e-4));
  planet sun(1.,3.1815e-3,4.3976e-3,-1.5181e-4,365.25*(-3.3248e-6),365.25*(6.6209e-6),365.25*(7.1779e-8));

/*
  ssystem SunEarth(5.0);
  SunEarth.add(earth);
  SunEarth.add(sun);
  SunEarth.VVerlet(100,1.,100,0);
*/

  ssystem TheSystem(100.);
  TheSystem.add(venus);
  TheSystem.add(earth);
  TheSystem.add(mars);
  TheSystem.add(jupiter);
  TheSystem.add(saturn);
  TheSystem.add(sun);
  TheSystem.VVerlet(4000,40.,4000,0);

  return 1;
}
