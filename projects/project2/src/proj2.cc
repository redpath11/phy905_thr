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

/* Global omega^2 declaration */
const double w2[6] = {0.0001,0.25,1.0,25.0,0.0625,0.0025};
const int omega_index = 1;// set which omega^2 to use

/* Function Declarations */
inline double Potential(double );
void test0(int ,double );// test building matrix
void test1();// test FindMaxOffDiag finder
void test2();// test Jacobi algorithm matrix w/ known e.values/vectors
void test3(int ,double );// test with the noninteracting case
void SetMatrix(double **,int ,double);
void jacobi(double **M,double **R,int n);
double FindMaxOffDiag(double **M, int *k, int *l, int n);
void Rotate(double **,double **,int ,int, int );
void TestOrthogonality(double **,int);
int * PrintLowestEvalues(double ** ,int ,int );
void WriteEigenvector(double **,int ,int ,double rmax);


/* Function Definitinos */
int main()
{
//    test0(100,5.);
//    test1();
//    test2();
    test3(100,10.);
    return 0;
}

// test setting up the matrix
void test0(int n,double rmax)
{
    // allocate matrix and set matrix elements
    double **M = AllocateMatrix(n,n);
    SetMatrix(M,n,rmax);

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
//    M[2][1]=10.;M[2][3]=-16.4;M[1][3]=56.2;
    M[2][1]=10.;M[2][3]=-16.4;M[1][3]=-56.2;
    double max = FindMaxOffDiag(M,&kk,&ll,n);
    FileMatrix("../Benchmark/test1.out",M,n);
    DeallocateMatrix(M,n,n);

    cout << "Max element at " << kk <<"," << ll << endl;
    cout << max << endl;

}

void test2()
{
    int n=3;
    double ** M = AllocateMatrix(n,n);
    double ** E = AllocateMatrix(n,n);
    for(int i=0;i<n;i++)
    {
        E[i][i]=1.0;
    }
    ReadMatrix("../Benchmark/test2_3.in",M,n);// test 3x3
    //ReadMatrix("../Benchmark/test2_4.in",M,n);// test 4x4
    WriteMatrix(M,n);
    jacobi(M,E,n);
    WriteMatrix(E,n);
    TestOrthogonality(E,n);
    int * evidx = PrintLowestEvalues(M,n,3);
//    cout << endl;
//    for(int i=0;i<3;i++){cout << evidx[i] << endl;}
    int gsIndex = evidx[0];
//    cout << gsIndex << endl;
    WriteEigenvector(E,n,gsIndex,1.);

    DeallocateMatrix(M,n,n);
    DeallocateMatrix(E,n,n);

}


// test with non interacting case
void test3(int n,double rmax)
{
    // For timing ...
    clock_t start,finish;
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
    SetMatrix(M,n,rmax);

    start = clock();
    jacobi(M,N,n);
    finish = clock();
    cout << fixed << setprecision(3) << "CPU time used: "
    <<(finish - start)/((double)CLOCKS_PER_SEC)
    << " s" << endl;

//    TestOrthogonality(N,n);

//    FileMatrix("efnts.out",N,n);
    int * evidx = PrintLowestEvalues(M,n,3);
    int gsIndex = evidx[0];
    WriteEigenvector(N,n,gsIndex,rmax);
    DeallocateMatrix(M,n,n);
    DeallocateMatrix(N,n,n);

}

void SetMatrix(double **MM,int n,double rm)
{
    // set step size
    double h  = rm/n;
    double hh = h*h;
    // precalculate stuff
    double ei = (-1.)/hh;
    double d  = (2.)/hh;

    for(int i=1;i<n;i++)
    {
        int idx   = i-1;
        MM[idx][idx]   = d + Potential(i*h);
        MM[idx+1][idx] = ei;
        MM[idx][idx+1] = ei;
    }
    // handle the last diagonal element ...
    // in practice the boundary condition is
    // u(rmax+h)=0
    MM[n-1][n-1] = d + Potential(rm);
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

}

// Find the lowest x eigenvalues
int * PrintLowestEvalues(double ** M,int nn,int x)
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
    return index;
}

// Test orthogonality of eigen vectors
void TestOrthogonality(double **V,int n)
{
    for(int i=0;i<n;i++)
    {
        for(int j=i+1;j<n;j++)
        {
            double sum=0;
            for(int k=0;k<n;k++)
            {
                sum+=V[k][i]*V[k][j];
            }
            if(sum>10e-8)
            {
                printf("WARNING: Dot product of eigenvectors %i and %i:\n",i+1,j+1);
                cout << sum << endl;
            }
        }
    }
}

void WriteEigenvector(double **Ef,int n,int idx,double rmax)
{
    char fname[512];
    sprintf(fname,"../Benchmark/GSwavefunction%4.2f.out",w2[omega_index]);
//    sprintf(fname,"../Benchmark/wavefunction%i.out",0);
    FILE *ofile;
    ofile=fopen(fname,"w");

    double h = rmax/((double)n);


    for(int i=0;i<n;i++)
    {
        fprintf(ofile,"%16.4E %16.4E\n",((double)(i+1))*h,Ef[i][idx]*Ef[i][idx]);
    }
    fclose(ofile);

}


//inline double Potential(double x){return x*x;}
//inline double Potential(double x){return w2[omega_index]*x*x;}
inline double Potential(double x){return w2[omega_index]*x*x + 1./x;}
