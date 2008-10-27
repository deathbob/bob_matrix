#include "bobMatrix.h"
#include <iostream>
#include "mini_matrix.h"

using namespace std;

int main(int argc, char* argv[]){
    int N = atoi(argv[1]);
    bobMatrix bob(N);
    bob.genSol();

    double a[] = {1,0,0,.001,1,0};

    bob.prints();
    bob.findB();
    bob.printb();
    bob.makeC();
    bob.printC();
    
    if(bob.makeUT()){
	bob.printC();
	bob.backSub();
	bob.printx();
	bob.makeS_X();
	cout<<"Two Norm = "<<setprecision(20)<<bob.twoNorm()<<endl;
    }
    else{cout<<"couldn't make UT"<<endl;}


    

    /*
      bob.A[0][0] = 1;
      bob.A[0][1] = 0;
      bob.A[0][2] = 2;
      bob.A[1][0] = 2;
      bob.A[1][1] = -1;
      bob.A[1][2] = 3;
      bob.A[2][0] = 4;
      bob.A[2][1] = 1;
      bob.A[2][2] = 8;
      bob.s[0] = -9;
      bob.s[1] = -2;
      bob.s[2] = 5;
    */

    /*    
	  bob.A[0][0] = 0;
	  bob.A[0][1] = 1;
	  bob.A[0][2] = -1;
	  bob.A[1][0] = 2;
	  bob.A[1][1] = -2;
	  bob.A[1][2] = 1;
	  bob.A[2][0] = 1;
	  bob.A[2][1] = 2;
	  bob.A[2][2] = 0;
	  bob.printA();
	  bob.s[0] = -4;
	  bob.s[1] = 9;
	  bob.s[2] = 0;
    */

    /*  LU Decomp 
	bob.genSol();
	bob.A[0][0] = 2;
	bob.A[0][1] = 6;
	bob.A[0][2] = -2;
	bob.A[1][0] = 2;
	bob.A[1][1] = double(2.0/3.0);
	bob.A[1][2] = double(1.0/3.0);
	bob.A[2][0] = 1;
	bob.A[2][1] = 2;
	bob.A[2][2] = -1;
	bob.printA();
	bob.s[0] = 4;
	bob.s[1] = 2;
	bob.s[2] = 1;
    */
    /*
      bob.genSol();
      bob.A[0][0] = 5;
      bob.A[0][1] = 0;
      bob.A[0][2] = -2;
      bob.A[1][0] = 3;
      bob.A[1][1] = 5;
      bob.A[1][2] = 1;
      bob.A[2][0] = 0;
      bob.A[2][1] = -3;
      bob.A[2][2] = 4;
      bob.s[0] = 7;
      bob.s[1] = 2;
      bob.s[2] = -4;
      bob.printA();
      bob.prints();
      bob.jacobi_iter(10);
    */
    /*
      bob.genSol();
      bob.A[0][0] = 8.888;
      bob.A[0][1] = 3333.33;
      bob.A[0][2] = 15.002;
      bob.A[1][0] = 2007;
      bob.A[1][1] = 2008;
      bob.A[1][2] = .005;
      bob.A[2][0] = 1001.2;
      bob.A[2][1] = 2013;
      bob.A[2][2] = 20.425;
      bob.s[0] = 3200.12;
      bob.s[1] = 2753.001;
      bob.s[2] = 2043.403;
      bob.printA();
      bob.prints();
      bob.jacobi_iter(10);
    */

    //    bob.prints();
    //    bob.findB();
    //    bob.printb();
    //    bob.makeC();
    //    bob.printC();
    /*
      if(bob.makeUT()){
      bob.printC();
      bob.backSub();
      bob.printx();
      bob.makeS_X();
      cout<<"Two Norm = "<<setprecision(20)<<bob.twoNorm()<<endl;
      }
      else{cout<<"couldn't make UT"<<endl;}
    */
    /*
      if(bob.LU_factor()){
	
      bob.makeLittleD();
      bob.LU_decomp();
      //	bob.printD();
      }
      else{cout<<"couldn't Factor LU "<<endl;}
    */
    
    
    return 0;
}
