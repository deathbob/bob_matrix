#include<iostream>
#include<limits>
#include<stdlib.h>
#include<iomanip>
#include<cmath>
#include<vector>
#include"mini_matrix.h"

using namespace std;




int main(){
    double arr1[] = {0,1,2,3,4,5,6,7,8};
    double arr2[] = {1,0,0,0,1,0,0,0,1};

    bob_mini_matrix A(3,3,arr1);
    A.print();
    bob_mini_matrix B(3,3,arr2);
    B.print();
    bob_mini_matrix C;
    C = A * B;
    C.print();

    C = A - B;
    C.print();

    bob_mini_matrix X(3, 'r');
    X.print();
    bob_mini_matrix Y(4,1,'i');
    Y.print();


    return 0;
}
