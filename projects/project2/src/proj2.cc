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

/* Function Declarations */
inline double V(double );
void test0(int ,double );// test building matrix
void test1();// test maxoffdiag finder
void jacobi(double **M,double **R,int n);
double maxoffdiag(double **M, int *k, int *l, int n);
void rotate(double **,double **,int ,int, int );

/* Function Definitinos */
int main()
{
//    test0(100,10.);
    test1();
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
        M[idx][idx]   = d + V(i*h);
        M[idx+1][idx] = ei;
        M[idx][idx+1] = ei;
    }
    // handle the last diagonal element ... this might be wrong ...
    M[n-1][n-1] = d + V(rmax);

    // file matrix to test in Python
    FileMatrix("../Benchmark/test0.out",M,n);
    DeallocateMatrix(M,n,n);

}

// test maxoffdiag finder
void test1()
{
    int n=4;
    int *kk,*ll;
    double **M = AllocateMatrix(n,n);
    M[2][1]=10.;
    double max = maxoffdiag(M,kk,ll,n);
    FileMatrix("../Benchmark/test1.out",M,n);
    DeallocateMatrix(M,n,n);

    cout << *kk << endl;
    cout << max << endl;

}

double maxoffdiag(double **M, int *k, int *l, int n)
{
    double max=0.;
    double abselement=0.;
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
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

/*
void jacobi(double **M,double **R,int n)
{

}

void rotate(double **M,double **R,int k,int l, int n)
{
    double s,c;
    if( M[k][l] != 0.0)
    {
        double t,tau;
        tau = (M[l][l] - M[k][k])/(2*M[k][l]);
        if(tau>0)
        {
            t = 1./(tau+sqrt(1.+tau*tau));
        }
        else
        {
            t = -1./(-tau+sqrt(1.+tau*tau));
        }
        c = 1./sqrt(1.+t*t);
        s = c*t;
        else
        {
            c=1.0; s=0.0;
        }
        double m_kk, m_ll, m_ik, m_il, m_ik, m_il;
        m_kk = M[k][k];
        m_ll = M[l][l];
        M[k][k] = c*c*m_kk - 2.*c*s*M[k][l] + s*s*m_ll;
        M[l][l] = s*s*m_kk + 2.*c*s*M[k][l] + c*c*m_ll;
        M[k][l] = 0.; M[l][k] = 0.;
    }
}
*/

inline double V(double x){return x*x;}
