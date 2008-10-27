#include<iostream>
#include<limits>
#include<stdlib.h>
#include<iomanip>
#include<cmath>
#include<vector>
#include"mini_matrix.h"

using namespace std;




int main(){
    //    problem 5
    //    double d[] = {2,0,0, 1,-1,-2, -1,0,1};

    //    double d[] = {4,-6,9,5,8,-2,7,-3,1};


    double d[] = {3, -3, 4, -5};
    bob_mini_matrix house(2,2,d);
    house.print();
    cout<<endl;
    //    house.householder();
    //    house.QR_method();
    //    house.jacobi_method();
    house.power_method();
    cout<<endl;
    //    house.print();



    return 0;
}

    // householder example p 116 
    //    double d[]={1,-2,3,-2,4,1,3,1,2};

    // 3.8 on page 143
    //    double d[] = {1,2,2,1};


    //example power method 3.3.1 page 108
    //    double d[] = {-3.9,-0.6,-3.9,.4,  1,0,0,0,  0,1,0,0,   0,0,1,0};
