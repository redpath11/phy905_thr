#include<iostream>
#include<fstream>
#include<iomanip>
#include<cmath>
#include<string>
/*
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
*/

using namespace std;

/* Function Declarations */
inline double f(double x);
inline double exact(double x);
int SolvePoisson1(int asize);// First step algorithm
int SolvePoisson2(int n);// Race Car algorithm
void FileArray(string filename,double* ary,int dim);
void FileSoltn(string filename,double* x,double* sol,double* ext,int dim);

/* Function Definitions */
int main()
{
//    SolvePoisson1(10);
    SolvePoisson2(10);
    return 0;

}

int SolvePoisson1(int asize)
{
    // Declare arrays
    double * x; x = new double[asize+1];// array of grid points
    double *ft; ft = new double[asize+1];// vector f
    double * d; d = new double[asize+1];// diagonal matrix elements
    double * e; e = new double[asize+1];// off diagonal matrix elements
    double * u; u = new double[asize+1];// solution vector
    double * ana; ana = new double[asize+1];// exact analytical result
 //   double * el;el= new double[asize];
    // Set step size
    double h = 1./((double) asize);
    double hh= h*h;

    // Initialize arrays
    u[0]=0.; u[asize]=0.;// boundary conditions
    d[0] = d[asize] = 2.;
    for(int i=0;i<asize+1;i++)
    {
        x[i] = i*h;
        e[i]=-1.;
        ft[i] = hh*f(i*h);
        ana[i] = exact(i*h);
    }
//    el[0]=-1.;
    // Write results/tests to file
    ofstream ofile; ofstream tfile;
    tfile.open("../Benchmark/test1.out");
    tfile << "x" << "   " << "d" << "   " << "e" << "   " << "f" << "   " << "ana" << endl;

    for(int i=0;i<asize;i++)
    {
        tfile << x[i] << "   " << d[i] << "   " << e[i] << "   " << ft[i] << "   " << ana[i] << endl;
    }
    tfile.close();


    for(int i=1;i<asize;i++)
    {
        d[i] = d[i] - (e[i-1]*e[i-1])/(d[i-1]);
        ft[i] = ft[i] - (e[i-1]*ft[i-1])/(d[i-1]);
    }

    u[asize-1] = ft[asize-1] / d[asize-1];
    for(int i=asize-1;i>=1;i--)
    {
        u[i] = (ft[i] - e[i]*u[i])/d[i];
    }

    // Test solver
    ofile.open("../Benchmark/result.out");
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    ofile << setw(15) << "x";
    ofile << setw(15) << "sol";
    ofile << setw(15) << "exact" << endl;
    for(int i=0;i<asize;i++)
    {
        ofile << setw(15) << setprecision(8) << x[i];
        ofile << setw(15) << setprecision(8) << u[i];
        ofile << setw(15) << setprecision(8) << ana[i] << endl;
    }
    delete [] x; delete [] ft; delete [] d; delete [] e; delete [] u;
}


int SolvePoisson2(int n)
{
    double * d = new double[n+1];// diagonal matrix elements
    double * ft= new double[n+1];// right hand side
    double * x = new double[n+1];// x-value
    double * u = new double[n+1];// solution
    double * ana=new double[n+1];// exact analytical result


    // Set step size
    double h = 1./((double) n);
    double hh= h*h;

    u[0] = u[n] = 0.0;//boundary conditions
    d[0] = d[n] = 2.0;//diagonal elements

    for(int i=0;i<=n;i++)
    {
        x[i] = i*h;
        ft[i] = hh*f(i*h);
        ana[i]= exact(i*h);
    }

    for(int i=1;i<n;i++)
    {
        d[i]=(i+1.0)/((double) i);
    }
    // Forward substitution
    for(int i=2;i<n;i++)
    {
        ft[i] += ft[i-1]/d[i-1];
    }
    // Backward substitution
    u[n-1] = ft[n-1]/d[n-1];
    for(int i=n-2;i>0;i--)
    {
        u[i] = (ft[i]+u[i+1])/d[i];
    }
    // Check shit
//    FileArray("solution.out",u,n+1);
    FileSoltn("solution.out",x,u,ana,n+1);
    delete [] x; delete [] ft; delete [] d; delete [] u;
}

void FileArray(string filename,double* ary,int dim)
{
    string outpath = "../Benchmark/";
    outpath.append(filename);
    ofstream ofile;
    ofile.open(outpath.c_str());
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    for(int i=0;i<dim;i++)
    {
        ofile << setw(15) << setprecision(8) << ary[i] << endl;
    }
    ofile.close();
}

// Write out x, solution, exact arrays to ../Benchmark/<filename>
// This function takes a filename, pointers to arrays and array
// dimension as input
void FileSoltn(string filename,double* x,double* sol,double* ext,int dim)
{
    string outpath = "../Benchmark/";
    outpath.append(filename);
    ofstream ofile;
    ofile.open(outpath.c_str());
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    ofile << setw(15) << "x";
    ofile << setw(15) << "solution";
    ofile << setw(15) << "exact";
    ofile << setw(15) << "RelError" << endl;
    for(int i=0;i<dim;i++)
    {
        double RelError=1.;
        if(ext[i]>0.)
        {
            RelError = fabs((ext[i]-sol[i])/ext[i]);
        }

        ofile << setw(15) << setprecision(8) << x[i];
        ofile << setw(15) << setprecision(8) << sol[i];
        ofile << setw(15) << setprecision(8) << ext[i];
        ofile << setw(15) << setprecision(8) << log10(RelError) << endl;
    }
}

inline double f(double x){return 100.0*exp(-10.0*x);}
inline double exact(double x){return 1.0-(1.0-exp(-10.))*x-exp(-10.*x);}
