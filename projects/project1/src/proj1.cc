#include<iostream>
#include<fstream>
#include<iomanip>
#include<cmath>
#include<string>
#include<stdio.h>
#include<math.h>
/*
#include<stdlib.h>
*/
#include "MatrixLib.hh"

using namespace std;

/* Function Declarations */
inline double f(double x);
inline double exact(double x);
void SolvePoisson1(int n);// First step algorithm
void SolvePoisson2(int n);// Race Car algorithm
void FileArray(string filename,double* ary,int dim);
void FileSoltn(string filename,double* x,double* sol,double* ext,int dim);

/* Function Definitions */
int main()
{
//    SolvePoisson1(10);
//    SolvePoisson2(10);
    int n=10;
    // test LUDecomposition with small dimensino
    double **M = AllocateMatrix(n,n);
    for(int j=0;j<n;j++)
    {
        M[j][j]=2.;
        M[j][j+1]=-1.;
        M[j+1][j]=-1.;
    }
    WriteMatrix(M,n);
    DeallocateMatrix(M,n,n);
    /*
    while(n<=1000)
    {
       SolvePoisson1(n);
       n*=10;
    }
    */
    return 0;

}

void SolvePoisson1(int n)
{
    // Declare arrays
    double * d = new double[n+1];// diagonal matrix elements
    double *ft = new double[n+1];// right hand side
    double *e1 = new double[n+1];// upper off-diagonal band
    double *e2 = new double[n+1];// lower off-diagonal band
    double * x = new double[n+1];// array of grid points
    double * u = new double[n+1];// solution vector
    double * ana = new double[n+1];// exact analytical result
    // Set step size
    double h = 1./((double) n);
    double hh= h*h;

    // Initialize arrays
    u[0]=0.; u[n]=0.;// boundary conditions
    for(int i=0;i<=n;i++)
    {
        x[i] = i*h;
	d[i] = 2.0;
        e1[i]=-1.; e2[i]=-1.;
        ft[i] = hh*f(i*h);
        ana[i] = exact(i*h);
    }

    // Forward substitution
    for(int i=2;i<n;i++)
    {
        d[i]  -= (e1[i]*e2[i])/d[i-1];
        ft[i] -= (e2[i-1]*ft[i-1])/(d[i-1]);
    }

    // Backward substitution
    u[n-1] = ft[n-1] / d[n-1];
    for(int i=n-2;i>0;i--)
    {
        u[i] = (ft[i] - e1[i]*u[i+1])/d[i];
    }

    /* Write results to file */
    char f[64];
    sprintf(f,"solution%i.out",n);
//    FileSoltn("solution1.out",x,u,ana,n+1);
    FileSoltn(f,x,u,ana,n+1);
    delete [] x; delete [] ft; delete [] d; delete [] e1; delete [] e2; delete [] u;
}


void SolvePoisson2(int n)
{
    // Declare arrays
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

    // Initialize x, ~f, analytic solution
    for(int i=0;i<=n;i++)
    {
        x[i] = i*h;
        ft[i] = hh*f(i*h);
        ana[i]= exact(i*h);
    }

    // Precalculate diagonal array
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
    /* Write results to file */
    char f[64];
    sprintf(f,"solution%i.out",n);
    FileSoltn(f,x,u,ana,n+1);
    delete [] x; delete [] ft; delete [] d; delete [] u;
}

// Write out the elements of the specified array to a file
// in ../Benchmark/<filename>
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
    string outpath = "../output/";
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
