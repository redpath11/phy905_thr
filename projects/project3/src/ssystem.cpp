#include "ssystem.h"
#include "planet.h"
#include "MatrixLib.hh"
#include <iostream>
#include <cmath>
#include "time.h"

ssystem::ssystem()
{
  total_planets  = 0;
  radius         = 100.;
  total_mass     = 0.;
  G              = 4.*M_PI*M_PI;
  c		 = 63198.;// [AU/year]
  totalKinetic   = 0.;
  totalPotential = 0.;
}


ssystem::ssystem(double rradius)
{
  total_planets  = 0;
  radius         = rradius;
  total_mass     = 0.;
  G              = 4.*M_PI*M_PI;
  c		 = 63198.;// [AU/year]
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

void ssystem::Euler(int nsteps, double years, int printNsteps, double epsilon)
{
  double h = years/((double) nsteps);
  double **Positions = AllocateMatrix(nsteps,3*total_planets);
  double deltaT=0.,deltaU=0.,deltaL=0.;
//  double EarthAngularMomentum;

  // add initial positions to positions array
  if(total_planets>0)
  {
    for(int pi=0;pi<total_planets;pi++)
    {
      planet &p = all_planets[pi];
      // Hardcode L magnitude calculation for earth
//      if(pi==0){EarthAngularMomentum = p.AngularMomentum();}
      for(int j=0;j<3;j++){ Positions[0][pi*3+j]=p.position[j]; }
    }
  }
  else cout << "No planets!!" << endl;

  double ** force_components = AllocateMatrix(total_planets,3);

  /* Initial energy */
  this->KineticEnergySystem();
  this->PotentialEnergySystem(0.);
  this->AngularMomentumSystem();
  deltaT=totalKinetic; deltaU=abs(totalPotential); deltaL=totalAngularMomentum;
  /*
  cout << "Initial KE: " << totalKinetic << endl;
  cout << "Initial PE: " << totalPotential << endl;
  cout << "Initial L: "  << totalAngularMomentum << endl;
  */
//  cout << "Initial Earth L: " << EarthAngularMomentum << endl;

  clock_t time1_E,time2_E;
  time1_E = clock();

  for(int n=1;n<nsteps;n++)
  {// Loop over steps
    for(int pi=0;pi<total_planets;pi++)
    {// Loop over planets
      planet &current = all_planets[pi];// load current planet
      for(int i=0;i<3;i++){force_components[pi][i]=0.;}
      // Forces on current
      for(int pii=pi+1;pii<total_planets;pii++)
      {
        planet &other = all_planets[pii];
	GravForce(current,other,force_components[pi],epsilon);
      }
      
      // New position and velocity for current
      for(int j=0;j<3;j++)
      {
//	cout << "  " << current.velocity[j];
        current.position[j] += current.velocity[j]*h;
	current.velocity[j] += force_components[pi][j]*h;
	Positions[n][pi*3+j] = current.position[j];
      }
//      cout << endl;

    }
  }

  time2_E = clock();
  cout << "Euler time elapsed: " << "\t" << ((float)(time2_E - time1_E)/CLOCKS_PER_SEC) << " seconds" << endl;

  /* Final energy and angular momentum */
  this->KineticEnergySystem();
  this->PotentialEnergySystem(0.);
  this->AngularMomentumSystem();
//  EarthAngularMomentum = all_planets[0].AngularMomentum();
  deltaT -= totalKinetic; deltaU -= abs(totalPotential); deltaL -= totalAngularMomentum;
  cout << "initial - final KE: " << deltaT << endl;
  cout << "initial - final PE: " << deltaU << endl;
  cout << "initial - final L: "  << deltaL << endl;
  /*
  cout << "Final KE: " << totalKinetic << endl;
  cout << "Final PE: " << totalPotential << endl;
  cout << "Final L: "  << totalAngularMomentum << endl;
  */
//  cout << "Final Earth L: " << EarthAngularMomentum << endl;

  FileMatrix("../Benchmark/ooEuler.out",Positions,printNsteps,3*total_planets);

  DeallocateMatrix(force_components,total_planets,3);
  DeallocateMatrix(Positions,nsteps,3*total_planets);


}

void ssystem::VVerlet(int nsteps,double years,int printNsteps,double epsilon)
{
  double h = years/((double) nsteps);
  double hOver2 = 0.5*h; double hhOver2 = 0.5*h*h;
  double **Positions = AllocateMatrix(nsteps,3*total_planets);
  double deltaT=0.,deltaU=0.,deltaL=0.;
//  double **Energy    = AllocateMatrix(nsteps,2);
//  int LostPlanets=0;

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

//  if(total_planets>2){
  this->KineticEnergySystem();
  this->PotentialEnergySystem(0.);
  this->AngularMomentumSystem();
  deltaT=totalKinetic; deltaU=abs(totalPotential); deltaL=totalAngularMomentum;
  /*
  cout << "Initial KE: " << totalKinetic << endl;
  cout << "Initial PE: " << totalPotential << endl;
  cout << "Initial L: " << totalAngularMomentum << endl;
  */
//    Energy[0][0] = totalKinetic; Energy[0][1] = totalPotential;
//  }

  double ** acceleration     = AllocateMatrix(total_planets,3);
  double ** acceleration_new = AllocateMatrix(total_planets,3);

  clock_t time1_VV,time2_VV;
  time1_VV = clock();

  for(int n=1;n<nsteps;n++)
//  while(n<nsteps && LostPlanets==0)
  {
//    n+=1;
    for(int pi=0;pi<total_planets;pi++)
    {
      planet &current = all_planets[pi];// load current planet
      for(int i=0;i<3;i++)
      {// reset accelerations
        acceleration[pi][i]=0.;
	acceleration_new[pi][i]=0.;
      }

      // Forces on current
      for(int pii=pi+1;pii<total_planets;pii++)
      {
        planet &other = all_planets[pii];
	GravAccel(current,other,acceleration[pi],epsilon);
      }
      
      // New position for current
      for(int i=0;i<3;i++)
      { current.position[i] += current.velocity[i]*h + hhOver2*acceleration[pi][i]; }

      // New forces on current
      for(int pii=pi+1;pii<total_planets;pii++)
      {
        planet &other = all_planets[pii];
	GravAccel(current,other,acceleration_new[pi],epsilon);
      }

      // New velocity for current
      for(int i=0;i<3;i++)
      { current.velocity[i] += hOver2*(acceleration[pi][i] + acceleration_new[pi][i]); }

      for(int j=0;j<3;j++){Positions[n][pi*3+j]=current.position[j];}
    }


/*    for(int np=0;np<total_planets;np++)
    {
      planet &current = all_planets[0];
      if(!(this->Bound(current)))
      {
        LostPlanets += 1;
	printf("Planet %i Lost!\n KE: %g\n PE: %g\n time: %g\n",current.kinetic,
								current.potential,
								n*h);
	printf("Planet %i Lost!\n last position: %g %g %g \n last velocity: %g %g %g \n time: %g \n",np,
	  current.position[0],current.position[1],current.position[2],
	  current.velocity[0],current.velocity[1],current.velocity[2],
	  n*h);
      }
    }
*/


  }// END loop over nsteps

  // Stopwatch for VV algorithm
  time2_VV = clock();
  cout << "VV time elapsed: " << "\t" << ((float)(time2_VV - time1_VV)/CLOCKS_PER_SEC) << " seconds" << endl;

//    if(total_planets>2){
  this->KineticEnergySystem();
//  this->PotentialEnergySystem(0.);
  this->PotentialEnergy2body(0.);
  this->AngularMomentumSystem();
  deltaT -= totalKinetic; deltaU -= abs(totalPotential); deltaL -= totalAngularMomentum;
  cout << "initial - final KE: " << deltaT << endl;
  cout << "initial - final PE: " << deltaU << endl;
  cout << "initial - final L: "  << deltaL << endl;
  /*
  cout << "Final KE: " << totalKinetic << endl;
  cout << "Final PE: " << totalPotential << endl;
  cout << "Final L: "  << totalAngularMomentum << endl;
  */
//      Energy[n][0] = totalKinetic; Energy[n][1] = totalPotential;
//    }

  FileMatrix("../Benchmark/ooSystem.out",Positions,printNsteps,3*total_planets);
//  FileMatrix("../Benchmark/ooSystemEnergy.out",Energy,printNsteps,2);

  DeallocateMatrix(acceleration,total_planets,3);
  DeallocateMatrix(acceleration_new,total_planets,3);
  DeallocateMatrix(Positions,nsteps,3*total_planets);
//  DeallocateMatrix(Energy,nsteps,2);

}

void ssystem::GravAccel(planet &current,planet &other, double *a,double epsilon)
{// Calculate the acceleration due to gravity between the two planets
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


void ssystem::HgGravAccel(planet &current,planet &other,double *a,double l)
{// Calculate the acceleration due to gravity between Mercury and the Sun
  double relative=0.;
  double r=current.distance(other);
  double coeff = this->G*other.mass/(r*r*r);

  for(int i=0;i<3;i++)
  {
    relative = current.position[i]-other.position[i];

    a[i] -= coeff*(relative + ((3.*l*l)/(r*r*c*c)));
  }
}


void ssystem::GravForce(planet &current, planet &other, double *f, double epsilon)
{// Calculate the gravitational force between two planets
  double relative=0.;
  double r=current.distance(other);
  double smoothing = epsilon*epsilon*epsilon;

  for(int i=0;i<3;i++)
  {
    relative = current.position[i]-other.position[i];
    f[i]    -= this->G*relative/((r*r*r) + smoothing);
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


void ssystem::AngularMomentumSystem()
{
  totalAngularMomentum = 0.;
  double totalLi[3] = {0.,0.,0.};
  for(int np=0;np<total_planets;np++)
  {
    planet &p = all_planets[np];
    double lmag=p.AngularMomentum();
    for(int j=0;j<3;j++)
    {totalLi[j] += p.am[j];}
  }

  for(int j=0;j<3;j++)
  {totalAngularMomentum += totalLi[j]*totalLi[j];}
  totalAngularMomentum = sqrt(totalAngularMomentum);

}


void ssystem::KineticEnergySystem()
{
  totalKinetic = 0.;
  for(int np=0;np<total_planets;np++)
  {
    planet &current = all_planets[np];
    current.kinetic = current.KineticEnergy();
    totalKinetic += current.kinetic;
  }
}


void ssystem::PotentialEnergySystem(double epsilon)
{
  totalPotential = 0.;
  // zero potential for all planets
  for(int np=0;np<total_planets;np++)
  {
    planet &current = all_planets[np];
    current.potential = 0.;
  }

  for(int pi=0;pi<total_planets;pi++)
  {
    planet &Current = all_planets[pi];
    for(int pj=pi+1;pj<total_planets;pj++)
    {
      planet &Other = all_planets[pj];
      Current.potential += Current.PotentialEnergy(Other,G,epsilon);
      Other.potential += Other.PotentialEnergy(Current,G,epsilon);
      totalPotential += Current.potential + Other.potential;
    }

  }

}


void ssystem::PotentialEnergy2body(double epsilon)
{
  totalPotential = 0.;
  // zero potential for all planets
  for(int np=0;np<total_planets;np++)
  {
    planet &current = all_planets[np];
    current.potential = 0.;
  }

  planet &p1 = all_planets[0];
  planet &p2 = all_planets[1];
  p1.potential += p1.PotentialEnergy(p2,G,epsilon);
  p2.potential += p2.PotentialEnergy(p1,G,epsilon);
  totalPotential += p1.potential + p2.potential;

}

bool ssystem::Bound(planet p)
{
  return ((p.kinetic + p.potential) < 0.);
}


/* Specilized VVerlet solver for the Mercury sun system. All output to files
 * has been eliminated to limit runtime and disk space usage */
void ssystem::HgVerlet(int nsteps,double years,int printNsteps,double epsilon)
{
  double h = years/((double) nsteps);
  double hOver2 = 0.5*h; double hhOver2 = 0.5*h*h;
  double perihelion=0.;
  double dist=0.,prevDist=0.3075;
  double xp=0.,yp=0.;
  bool approach=false;
//  double **Positions = AllocateMatrix(nsteps,3*total_planets);
//  double **Energy    = AllocateMatrix(nsteps,2);

  // add initial positions to positions array
  if(total_planets>0)
  {
    planet &pHg = all_planets[0];
    planet &Sun = all_planets[1];
    perihelion = pHg.distance(Sun);
    cout << perihelion << endl;

/*
    for(int j=0;j<3;j++)
    {
      Positions[0][j]=pHg.position[j];
      Positions[0][3+j]=Sun.position[j];
    }
*/
  }
  else cout << "No planets!!" << endl;

  /* Energy Calculations 
  this->KineticEnergySystem();
  if(total_planets>2){this->PotentialEnergySystem(0.);}
  else{this->PotentialEnergy2body(0.);}
//  Energy[0][0] = totalKinetic; Energy[0][1] = totalPotential;
*/
  double ** acceleration     = AllocateMatrix(total_planets,3);
  double ** acceleration_new = AllocateMatrix(total_planets,3);
//  double * acceleration = new double[3];
//  double * acceleration_new = new double[3];

  clock_t time1_VV,time2_VV;
  time1_VV = clock();

  for(int n=1;n<nsteps;n++)
  {
    double xp = all_planets[0].position[0];
    double yp = all_planets[0].position[1];
/*
    planet &Hg = all_planets[0];
    planet &S  = all_planets[1];

    double HgAngularMomentum = Hg.AngularMomentum();

    for(int i=0;i<3;i++)
    {// reset accelerations
      acceleration[i]=0.;
      acceleration_new[i]=0.;
    }
    // Force on Hg
    GravAccel(Hg,S,acceleration,epsilon);
    // new position for Hg
    for(int i=0;i<3;i++)
    {Hg.position[i] += Hg.velocity[i]*h + hhOver2*acceleration[i];}
    // new force on Hg
    GravAccel(Hg,S,acceleration_new,epsilon);
    // new velocity for Hg
    for(int i=0;i<3;i++)
    {Hg.velocity[i] += hOver2*(acceleration[i] + acceleration_new[i]);}
*/
/**/
    for(int pi=0;pi<total_planets;pi++)
    {
      planet &current = all_planets[pi];// load current planet
      for(int i=0;i<3;i++)
      {// reset accelerations
        acceleration[pi][i]=0.;
	acceleration_new[pi][i]=0.;
      }

      // Forces on current
      for(int pii=pi+1;pii<total_planets;pii++)
      {
        planet &other = all_planets[pii];
	GravAccel(current,other,acceleration[pi],epsilon);
//	HgGravAccel(current,other,acceleration[pi],current.AngularMomentum());
      }
      
      // New position for current
      for(int i=0;i<3;i++)
      { current.position[i] += current.velocity[i]*h + hhOver2*acceleration[pi][i]; }

      // New forces on current
      for(int pii=pi+1;pii<total_planets;pii++)
      {
        planet &other = all_planets[pii];
	GravAccel(current,other,acceleration_new[pi],epsilon);
//	HgGravAccel(current,other,acceleration_new[pi],current.AngularMomentum());
      }

      // New velocity for current
      for(int i=0;i<3;i++)
      { current.velocity[i] += hOver2*(acceleration[pi][i] + acceleration_new[pi][i]); }

//      for(int j=0;j<3;j++){Positions[n][pi*3+j]=current.position[j];}
    }

    this->KineticEnergySystem();
    if(total_planets>2){this->PotentialEnergySystem(0.);}
    else{this->PotentialEnergy2body(0.);}
//    Energy[n][0] = totalKinetic; Energy[n][1] = totalPotential;

    /* Calculate distance */
    double time = ((double)n)*h;
    prevDist = dist;
    dist = all_planets[0].distance(all_planets[1]);
    if(dist<prevDist)// at aphelion
    {approach=true;}
    if(dist>prevDist && approach)// at perihelion
    {
      perihelion=prevDist;
      double Theta = atan2(yp,xp) * 180./M_PI;
      approach=false;
      cout << "Mercury at perihelion: " << perihelion << endl;
      cout << "Time: " << time << endl;
      cout << "xp = " << xp << " ; yp = " << yp << endl;
      cout << "Theta = " << Theta << endl;
    }

    /* check perihelion 

    
    if(((double)n)*h>100)
    {
    planet &hg = all_planets[0];
    planet &s  = all_planets[1];
    double dist = hg.distance(s);
    double tolerance = 13e-3;
    if(dist<=perihelion+tolerance)
    {
      cout << "Mercury at perihelion: " << hg.distance(s) << endl;
      double Theta = atan2(hg.position[1],hg.position[0]) * 180./3.14159265;
      cout << hg.position[0] << "   " << hg.position[1] << endl;
      cout << "Theta_p [deg]: " << Theta << endl;
    }
    }
*/


  }// END loop over nsteps

  // Stopwatch for VV algorithm
  time2_VV = clock();
  cout << "VV time elapsed: " << "\t" << ((float)(time2_VV - time1_VV)/CLOCKS_PER_SEC) << " seconds" << endl;

//  FileMatrix("../Benchmark/ooSystem.out",Positions,printNsteps,3*total_planets);
//  FileMatrix("../Benchmark/ooSystemEnergy.out",Energy,printNsteps,2);

  DeallocateMatrix(acceleration,total_planets,3);
  DeallocateMatrix(acceleration_new,total_planets,3);
//  DeallocateMatrix(Positions,nsteps,3*total_planets);
//  DeallocateMatrix(Energy,nsteps,2);

}

