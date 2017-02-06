#include<iostream>
#include<fstream>
#include<iomanip>
#include<cmath>
#include<string>
#include<stdio.h>
#include<math.h>
#include<ctime>
//#include<chrono>
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
void SolveLUDcmp(int n);// LU Decomposition method
void FileArray(string filename,double* ary,int dim);
void FileSoltn(string filename,double* x,double* sol,double* ext,int dim);
void SolveLUDcmp(int n);


/* Function Definitions */
int main()
{
//    SolvePoisson1(10);
//    SolvePoisson2(10);
    int n=10;
//    SolvePoisson2(n);
//    SolveLUDcmp(n);
    /**/
    while(n<=10000000)
    {
//       SolvePoisson1(n);
       SolvePoisson2(n);

       n*=10;
    }

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
    // For timing ...
    clock_t start,finish;

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

    start = clock();
//    auto t_start = std::chrono::high_resolution_clock::now();
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
    finish = clock();
//    auto t_end = std::chrono::high_resolution_clock::now();
    cout << fixed << setprecision(3) << "CPU time used: "
        <<( 1000000.0 * (finish - start)/CLOCKS_PER_SEC)
        << " us" << endl;
//        << "Wall clock time passed: "
//        << chrono::duration<double, milli>(t_end-t_start).count()
//        << " ms" << endl;

    //cout << ( (finish - start)) << endl;

    /* Write results to file */
    char f[64];
    sprintf(f,"solution%i.out",n);
//    FileSoltn("solution1.out",x,u,ana,n+1);
//    FileSoltn(f,x,u,ana,n+1);
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
//    FileArray("../Benchmark/P2RHS.out",ft,n);

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

void SolveLUDcmp(int nn)
{
    double h = 1./((double) nn);
    double hh= h*h;
    int n=nn;// to set up matrix with endpoints
    double **M = AllocateMatrix(n,n);
    int *idx = new int[n];
    double *b = new double[n];

    for(int i=0;i<n;i++)
    {
        M[i][i]=2.;
        if(i<n-1)
        {
            M[i][i+1]=-1.;
            M[i+1][i]=-1.;
        }
        /* the rhs doesn't include the endpoints
        if(i==0 || i==(n-1))
        {
            b[i]=0.;
        }
        else{b[i]=hh*f(i*h);}
        */
        b[i]=hh*f(double(i+1)*h);
    }
//    b[0]=b[n-1]=0.;
    FileArray("../Benchmark/RHS.out",b,n);
    FileMatrix("../Benchmark/M_init.out",M,n);
    LUDecomposition(M,n,idx);
    LUBackwardSubstitution(M,n,idx,b);
    FileArray("../Benchmark/LUsol.out",b,n);
    FileMatrix("../Benchmark/M_finl.out",M,n);
    DeallocateMatrix(M,n,n);
    delete [] idx; delete [] b; delete [] idx;

}


inline double f(double x){return 100.0*exp(-10.0*x);}
inline double exact(double x){return 1.0-(1.0-exp(-10.))*x-exp(-10.*x);}
