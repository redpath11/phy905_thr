#include "planet.h"

planet::planet()
{
    mass = 1.;
    position[0] = 0.;
    position[1] = 0.;
    position[2] = 0.;
    velocity[0] = 0.;
    velocity[1] = 0.;
    velocity[2] = 0.;
    am[0]       = 0.;
    am[1]       = 0.;
    am[2]       = 0.;
    potential = 0.;
    kinetic = 0.;
}

planet::planet(double M, double x, double y, double z, double vx, double vy, double vz)
{
    mass = M;
    position[0] = x;
    position[1] = y;
    position[2] = z;
    velocity[0] = vx;
    velocity[1] = vy;
    velocity[2] = vz;
    am[0]       = 0.;
    am[1]       = 0.;
    am[2]       = 0.;
    potential = 0.;
    kinetic = 0.;
}


double planet::rmag()
{
  double xx,yy,zz;
  xx = this->position[0]*this->position[0];
  yy = this->position[1]*this->position[1];
  zz = this->position[2]*this->position[2];

  return sqrt(xx + yy + zz);

}

double planet::distance(planet otherPlanet)
{
    double x1,y1,z1,x2,y2,z2,xx,yy,zz;

    x1 = this->position[0];
    y1 = this->position[1];
    z1 = this->position[2];

    x2 = otherPlanet.position[0];
    y2 = otherPlanet.position[1];
    z2 = otherPlanet.position[2];

    xx = x1-x2;
    yy = y1-y2;
    zz = z1-z2;

    return sqrt(xx*xx + yy*yy + zz*zz);
 }

double planet::GravitationalForce(planet otherPlanet,double Gconst)
{
    double r = this->distance(otherPlanet);
    if(r!=0) return Gconst*this->mass*otherPlanet.mass/(r*r);
    else return 0;
}

double planet::Acceleration(planet otherPlanet, double Gconst)
{
    double r = this->distance(otherPlanet);
    if(r!=0) return this->GravitationalForce(otherPlanet,Gconst)/(this->mass*r);
    else return 0;
}

double planet::KineticEnergy()
{
    double velocity2 = (this->velocity[0]*this->velocity[0]) + (this->velocity[1]*this->velocity[1]) + (this->velocity[2]*this->velocity[2]);
    return 0.5*this->mass*velocity2;
}

double planet::PotentialEnergy(planet &otherPlanet, double Gconst, double epsilon)
{
    if(epsilon==0.0) return -Gconst*this->mass*otherPlanet.mass/this->distance(otherPlanet);
    else return (Gconst*this->mass*otherPlanet.mass/epsilon)*(atan(this->distance(otherPlanet)/epsilon) - (0.5*M_PI));
}

double planet::AngularMomentum()
{
  double yz,zy,zx,xz,xy,yx;

  yz = this->position[1]*this->velocity[2];
  zy = this->position[2]*this->velocity[1];
  zx = this->position[2]*this->velocity[0];
  xz = this->position[0]*this->velocity[2];
  xy = this->position[0]*this->velocity[1];
  yx = this->position[1]*this->velocity[0];

  am[0] = yz-zy; am[1] = zx - xz; am[2] = xy - yx;

  return sqrt(am[0]*am[0] + am[1]*am[1] + am[2]*am[2]);

}
