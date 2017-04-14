#include "vec3.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <stdio.h>
using namespace std;

// Constructors
vec3::vec3()
{
    components[0] = 0;
    components[1] = 0;
    components[2] = 0;
}

vec3::vec3(const vec3 &copy)
{
    components[0] = copy.x();
    components[1] = copy.y();
    components[2] = copy.z();
}

vec3::vec3(const vec3 *copy)
{
    components[0] = copy->x();
    components[1] = copy->y();
    components[2] = copy->z();
}

vec3::vec3(double x, double y, double z)
{
    components[0] = x;
    components[1] = y;
    components[2] = z;
}

// Mathematical operations
double vec3::dot(vec3 otherVector)
{
    return otherVector[0]*components[0] + otherVector[1]*components[1] + otherVector[2]*components[2];
}

double vec3::lengthSquared()
{
    // Returns the square of the length (or norm) of the vector
    return components[0]*components[0]+components[1]*components[1]+components[2]*components[2];
}

double vec3::length()
{
    // Returns the length (or norm) of the vector
    return sqrt(lengthSquared());
}

// Operators
vec3 &vec3::operator+=(double rhs)
{
    components[0] += rhs;
    components[1] += rhs;
    components[2] += rhs;
    return *this;
}

vec3 &vec3::operator+=(vec3 rhs)
{
    components[0] += rhs[0];
    components[1] += rhs[1];
    components[2] += rhs[2];
    return *this;
}

vec3 &vec3::operator*=(double rhs)
{
    components[0] *= rhs;
    components[1] *= rhs;
    components[2] *= rhs;
    return *this;
}

vec3 &vec3::operator*=(vec3 rhs)
{
    components[0] *= rhs[0];
    components[1] *= rhs[1];
    components[2] *= rhs[2];
    return *this;
}

vec3 &vec3::operator-=(double rhs)
{
    components[0] -= rhs;
    components[1] -= rhs;
    components[2] -= rhs;
    return *this;
}

vec3 &vec3::operator-=(vec3 rhs)
{
    components[0] -= rhs[0];
    components[1] -= rhs[1];
    components[2] -= rhs[2];
    return *this;
}

vec3 &vec3::operator/=(double rhs)
{
    components[0] /= rhs;
    components[1] /= rhs;
    components[2] /= rhs;
    return *this;
}

vec3 &vec3::operator/=(vec3 rhs)
{
    components[0] /= rhs[0];
    components[1] /= rhs[1];
    components[2] /= rhs[2];
    return *this;
}

// Printing
void vec3::print(string file)
{
    ofstream ofile;
    ofile.open(file.c_str(),std::ofstream::out | std::ofstream::app);
    ofile << setw(12) << setprecision(4) << components[0];
    ofile << setw(12) << setprecision(4) << components[1];
    ofile << setw(12) << setprecision(4) << components[2] << endl;

    ofile.close();

}

void vec3::print()
{
    cout << "[" << components[0] << "," << components[1] << "," << components[2] << "]" << endl;
}
