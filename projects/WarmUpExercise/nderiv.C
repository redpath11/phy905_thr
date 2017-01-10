#include<stdio.h>
#include<math.h>


/* Function Declarations */
float  func(float x);
double func(double x);
double nDeriv1(double x,double h);
double nDeriv2(double x,double h);

/* Function Definitions */

int main()
{
  /* Test 1: function eval 
  float  xf=0.2;
  double xd=0.2;
  printf("The arc tangent (single precision) of %g is %.20e radians.\n",xf,func(xf));
  printf("The arc tangent (double precision) of %g is %.20e radians.\n",xd,func(xd));
  */

  /* Test 2: derivative eval 
  double step = 0.001;
  double x = sqrt(2.);
  double nd1x = nDeriv1(x,step);
  double nd2x = nDeriv2(x,step);

  printf("The first deriv of arctan(%g) is %g (method 1)\n",x,nd1x);
  printf("The first deriv of arctan(%g) is %g (method 2)\n",x,nd2x);
  */

  double x = sqrt(2.);
  double * h;
  h = new double[5];
  h[0]=0.1;
  for(int i=1;i<5;i++)
  { h[i] = h[i-1]*0.1; }

  double * n1;
  n1 = new double[5];
  double * n2;
  n2 = new double[5];

  for(int i=0;i<5;i++)
  {
    n1[i] = nDeriv1(x,h[i]);
    n2[i] = nDeriv2(x,h[i]);
    
    printf("For h = %g\n method 1 gives %g\n method 2 gives %g\n",h[i],n1[i],n2[i]);
  }

  return 0;
}


// Function to be differentiated
// result is in radians
float func(float x)
{ return atan(x); }

double func(double x)
{ return atan(x); }


// Numerical differentiation
// method 1
double nDeriv1(double x,double h)
{
  double fip1 = func(x+h);
  double fi   = func(x);

  return (fip1 - fi) / h;
}


// Numerical differentiation
// method 2
double nDeriv2(double x,double h)
{
  double fip1 = func(x+h);
  double fim1 = func(x-h);

  return (fip1 - fim1) / (2.0*h);
}
