#include<iostream>
#include<fstream>
#include<iomanip>
#include<cmath>
#include<string>
#include<stdio.h>
#include<math.h>
#include<ctime>

#include "MatrixLib.hh"

using namespace std;

/* Global w^2 declaration */
const double w2 = 1.0;

/* Function Declarations */
inline double Potential(double );
void test0(int ,double );// test building matrix
void test1();// test FindMaxOffDiag finder
void test2();// test Jacobi algorithm matrix w/ known e.values/vectors
void test3(int ,double );// test with the noninteracting case
void jacobi(double **M,double **R,int n);
double FindMaxOffDiag(double **M, int *k, int *l, int n);
void Rotate(double **,double **,int ,int, int );
void TestOrthogonality(double **,int);
void PrintLowestEvalues(double ** ,int ,int );


/* Function Definitinos */
int main()
{
//    test0(100,15.);
//    test1();
//    test2();
    test3(500,5.);
    return 0;
}

// test setting up the matrix
void test0(int n,double rmax)
{
    // set step size
    double h  = rmax/n;
    double hh = h*h;
    /* set up matrix */
    // precalculate stuff
    double ei = (-1.)/hh;
    double d  = (2.)/hh;
    // allocate matrix and set matrix elements
    double **M = AllocateMatrix(n,n);
    for(int i=1;i<n;i++)
    {
        int idx   = i-1;
        M[idx][idx]   = d + Potential(i*h);
        M[idx+1][idx] = ei;
        M[idx][idx+1] = ei;
    }
    // handle the last diagonal element ... this might be wrong ...
    M[n-1][n-1] = d + Potential(rmax);

    // file matrix to test in Python
    FileMatrix("../Benchmark/test0.out",M,n);
    DeallocateMatrix(M,n,n);

}

// test FindMaxOffDiag finder
void test1()
{
    int n=4;
    int kk,ll;
    double **M = AllocateMatrix(n,n);
    M[2][1]=10.;M[1][3]=-16.4;M[1][3]=56.2;
    double max = FindMaxOffDiag(M,&kk,&ll,n);
    FileMatrix("../Benchmark/test1.out",M,n);
    DeallocateMatrix(M,n,n);

    cout << "Max element at " << kk <<"," << ll << endl;
    cout << max << endl;

}

void test2()
{
    int n=4;
    double ** M = AllocateMatrix(n,n);
    double ** E = AllocateMatrix(n,n);
    for(int i=0;i<n;i++)
    {
        E[i][i]=1.0;
    }
    ReadMatrix("../Benchmark/test2_4.in",M,n);// test 3x3
    //ReadMatrix("../Benchmark/test2_4.in",M,n);// test 4x4
    WriteMatrix(M,n);
    jacobi(M,E,n);
    WriteMatrix(E,n);
    TestOrthogonality(E,n);

}


// test with non interacting case
void test3(int n,double rmax)
{
    // set step size
    double h  = rmax/n;
    double hh = h*h;
    /* set up matrix */
    // precalculate stuff
    double ei = (-1.)/hh;
    double d  = (2.)/hh;
    // allocate matrix and set matrix elements
    // two matrices: one to diagonalize, one for eigenvectors
    double **M = AllocateMatrix(n,n);
    double **N = AllocateMatrix(n,n);
    // initialize eigenvectors
    for(int i=0;i<n;i++)
    {
        N[i][i] = 1.0;
    }
    // generate tridiagonal matrix
    for(int i=1;i<n;i++)
    {
        int idx   = i-1;
        M[idx][idx]   = d + Potential(i*h);
        M[idx+1][idx] = ei;
        M[idx][idx+1] = ei;
    }
    // handle the last diagonal element ... this might be wrong ...
    M[n-1][n-1] = d + Potential(rmax);

    jacobi(M,N,n);

    FileMatrix("efnts.out",N,n);
    DeallocateMatrix(M,n,n);
    DeallocateMatrix(N,n,n);

}


double FindMaxOffDiag(double **M, int *k, int *l, int n)
{
    double max=0.;
    double abselement=0.;
    for(int i=0;i<n;i++)
    {
        for(int j=i+1;j<n;j++)// symmetric matrix
        {
            abselement=abs(M[i][j]);
            if(abselement>max)
            {
                *k=i;
                *l=j;
                max=abselement;
            }
        }
    }

    return max;
}

// perform rotation
// **M - matrix pointer
// **V - eigenvectors
// k,l - indices for largest offdiagonal element
// n   - size of matrix
void Rotate(double **M,double **V,int p,int q, int n)
{
    double s,c;
    if( M[p][q] != 0.0)
    {
        double t,tau;
        tau = (M[q][q] - M[p][p])/(2*M[p][q]);
        if(tau>0)
        {
            t = 1./(tau+sqrt(1.+tau*tau));
        }
        else if(tau<0)
        {
            t = -1./(-tau+sqrt(1.+tau*tau));
        }
        else
        {
            t = 1.;
        }
        c = 1./sqrt(1.+t*t);
        s = c*t;
        /* Debugging
        cout << endl;
        cout << "tau = " << tau << endl;
        cout << "t = " << t << "  c = " << c << "  s = " << s << endl;
        */
    }
    else
    {
        c=1.0; s=0.0;
    }
    double m_pp, m_qq, m_ip, m_iq, v_ip, v_iq;
    m_pp = M[p][p];
    m_qq = M[q][q];
    M[p][p] = c*c*m_pp - 2.*c*s*M[p][q] + s*s*m_qq;
    M[q][q] = s*s*m_pp + 2.*c*s*M[p][q] + c*c*m_qq;
    M[p][q] = 0.; M[q][p] = 0.;
    for(int i=0;i<n;i++)
    {
        if(i!=p && i!=q)
        {
            m_ip = M[i][p]; m_iq = M[i][q];
            M[i][p] = c*m_ip - s*m_iq;
            M[p][i] = M[i][p];
            M[i][q] = c*m_iq + s*m_ip;
            M[q][i] = M[i][q];
        }
        // eigen vecotors
        v_ip = V[i][p]; v_iq = V[i][q];
        V[i][p] = c*v_ip - s*v_iq;
        V[i][q] = c*v_iq + s*v_ip;
    }

}

// Apply Jacobi method
// M** - matrix
// V** - eigenvalues
// n   - size of matrix
void jacobi(double **M,double **V,int n)
{
    int k,l;
    int iterations=0;
    double max_iterations = (double) n * (double) n * (double) n;
    double maxoffdiag=FindMaxOffDiag(M,&k,&l,n);
    double limit=1.0e-8;
    while(maxoffdiag>limit && (double) iterations < max_iterations)
    {
        //cout << "Iteration: " << iterations << endl;
        maxoffdiag = FindMaxOffDiag(M,&k,&l,n);
        if(iterations%1000==0){cout << "Max OffDiag: " << maxoffdiag << endl;}
        Rotate(M,V,k,l,n);
        iterations++;
    }
    cout << "Number of iterations: " << iterations << endl;
/*    cout << "Eigenvalues: " << endl;
    double min=100.;
    int minelement=0;
    for(int i=0;i<n;i++)
    {
        if(M[i][i]<min){minelement=i;min=M[i][i];}
        printf(" %g\n",M[i][i]);
    }
    cout << "Ground State Energy: " << min << endl;
    cout << "Index of ground state: " << minelement << endl;
*/
    PrintLowestEvalues(M,n,3);



}

// Find the lowest x eigenvalues
void PrintLowestEvalues(double ** M,int nn,int x)
{
    double * eigenvalues = new double[x];
    int * index = new int[x];
    for(int i=0;i<x;i++){eigenvalues[i]=10000.;}
    for(int i=0;i<nn;i++)
    {
        if(M[i][i]<eigenvalues[0])
        {
            eigenvalues[0]=M[i][i];
            index[0]=i;
        }
    }
    if(x>1)
    {
    for(int i=1;i<x;i++)
    {
        for(int j=0;j<nn;j++)
        {
            if(M[j][j]<eigenvalues[i] && M[j][j]>eigenvalues[i-1])
            {
                eigenvalues[i]=M[j][j];
                index[i]=j;
            }
        }
    }
    }
    printf("Lowest %i Eigenvalues\n",x);
    cout << "Eigenvalues   Index" << endl;
    for(int i=0;i<x;i++)
    {
        printf("%-14.6f %-14i\n",eigenvalues[i],index[i]);
    }
}

// Test orthogonality of eigen vectors
void TestOrthogonality(double **V,int n)
{
    for(int i=0;i<n;i++)
    {
        for(int j=i+1;j<n;j++)
        {
            double sum=0;
            printf("Dot product of eigenvectors %i and %i:\n",i+1,j+1);
            for(int k=0;k<n;k++)
            {
                sum+=V[k][i]*V[k][j];
            }
            cout << sum << endl;
        }
    }
}


//inline double Potential(double x){return x*x;}
inline double Potential(double x){return w2*x*x;}
