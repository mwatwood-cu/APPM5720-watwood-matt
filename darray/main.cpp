#include <iostream>
#include <iomanip>
#include <cmath>
#include "D1array.h"
#include "D2array.h"
#include <chrono>

using namespace std::chrono;
using namespace std;

// Main program
int main()
{
    D2array u1(4, -15, 16, -31, 32);
    D2array u2(16, -31, 32, -31,32);
 
    /*
    D1array u(2,-1,3);
    double a;
    int count = 0;
    for (int i = -1; i<4; i++) {
        for (int c = 1; c<3; c++) {
            u(c,i) = count++;
        }
    }

    D1array v;
    v.copy(u);
       
    for (int i = -1; i<4; i++)
        for (int c = 1; c<3; c++)
            cout << " " << c << " " << i << " "<< v(c,i)<<" at index "<< v.m_base+c*v.m_offc+i*v.m_offi << "\n";


    D2array u2(2,-1,1,-1, 2);
    count = 0;
    for (int c = 1; c<3; c++) {
        for (int i = -1; i<2; i++) {
            for (int j = -1; j<3; j++) {
                u2(c,i,j) = count;
                count+=1;
            }
        }
    }

    D2array v2;
    v2.copy(u2);
    for (int c=1; c<3; c++) {
        for(int i=-1; i<2; i++) {
            for(int j=-1; j<3; j++) {
                cout << "(" << c << ", " << i << ", "<<j<<") = "<< v2(c,i,j)<<" at index "<< v2.m_base+c*v2.m_offc+i*v2.m_offi+j*v2.m_offj << "\n";
            }
        }
    }*/
       
    return 0;
}

