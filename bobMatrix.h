#include<iostream>
#include<limits>
#include<stdlib.h>
#include<iomanip>
#include<cmath>
#include<vector>

using namespace std;

double randomDouble(){
  return fmod(double(rand()), 1000000);
};

class bobMatrix{
 public:  // private:
  int size;
  vector< vector<double> > A;
  vector< vector<double> > C;
  vector< vector<double> > D;
  vector<double> b;
  vector<double> s;
  vector<double> x;
  vector<double> s_x;
  vector<double> d;
  vector<double> LUx;
  vector<double> LUy;
  double det;
  // public:
  bobMatrix(int N);
  void genSol();
  void findB();
  void printA();
  void makeC();
  bool makeUT();
  void backSub();
  void printx();
  void printb();
  void prints();
  void printC();
  void makeS_X();
  void printS_X();
  double twoNorm();
  void printD();

  bool LU_factor();
  void makeLittleD();
  void LU_decomp();
  void jacobi_iter(int iters);

};// end matrix

void bobMatrix::jacobi_iter(int iters){
  double temp = 0.0;
  vector<double> holder;
  holder.resize(size,0.0);
  x.clear();
  x.resize(size, 0.0);
  for(int i = 0; i < iters; i++){
    for(int row = 0; row < size; row++){
      for(int col = 0; col < size; col++){
	if(row == col){}
	else{
	  temp += A[row][col]*x[col];
	}
      }
      holder[row] = (s[row] - temp) * (1.0/A[row][row]);
      temp = 0.0;
    }
    x = holder;
    printx();
    temp = 0.0;
    double max = 0.0;
    for(int row = 0; row < size; row++){
      for(int col = 0; col < size; col++){
	temp += A[row][col] * x[col];
	//	cout<<"******"<<A[row][col]<<"  "<<x[col]<<"  "<<temp<<endl;
      }
      holder[row] = s[row] - temp;
      if(fabs(holder[row]) > max){
	max = fabs(holder[row]);
      }
      temp = 0.0;
    }
    cout<<"||r(x"<<i<<")|| = "<<max<<endl;
  }
}

void bobMatrix::LU_decomp(){
  LUx.resize(size, 0.0);
  LUy.resize(size, 0.0);
  for(int k = 0; k < size; k++){
    double temp = 0.0;
    for(int i = 0; i < k; i++){
      temp += D[k][i] * LUy[i];
    }
    LUy[k] = (d[k] - temp) / D[k][k];
    temp = 0.0;
    cout<<k<<" LUy  "<<LUy[k]<<endl;
  }
  for(int k = size - 1; k >= 0; k--){
    double temp = 0.0;
    for(int i = k+1; i < size; i++){
      temp += D[k][i] * LUx[i];
    }
    LUx[k] = (LUy[k] - temp);
    temp = 0.0;
    cout<<k<<" LUx  "<<LUx[k]<<endl;
  }
  cout<<" Calculating two-norm "<<endl;
  vector<double> err;
  err.resize(size, 0.0);
  for(int i = 0; i < size; i++){
    err[i] = s[i] - LUx[i];
  }
  double squareSol = 0.0;
  double squareErr = 0.0;
  for(int i = 0; i < size; i++){
    squareErr += err[i] * err[i];
    squareSol += s[i] * s[i];
  }
  cout<<"Two-Norm is "<<sqrt(squareErr) / sqrt(squareSol)<<endl;
  //NOT DONE
  //have to multiply A * LUx to get solution, then compare that computed solution 
  // with the known solution. TODO.
      
}

void bobMatrix::makeLittleD(){
  d.resize(size,0.0);
  prints();
  for(int row = 0; row < size; row++){
    double temp = 0.0;
    for(int col = size; col < size*2; col++){
      temp += (D[row][col] * s[col-size]);
      //      cout<<"D[row][col] "<<D[row][col]<<" s[col]  "<<s[col-size]<<" temp "<<temp<<endl;
    }
    d[row] = temp;
    cout<<d[row]<<endl;
  }
}


bool bobMatrix::LU_factor(){
  det = 1;
  D.resize(size);
  for(int i = 0; i < size; i++){
    D[i].resize(size*2,0.0);
  }

  for(int row = 0; row < size; row++){
    for(int col = 0; col < size; col++){
      D[row][col] = A[row][col];
    }
    D[row][row+size] = 1;
  }

  for(int j = 0; j < size; j++){
    // compute some crap 
    double foo = 0.0;
    for(int k = j; k < size; k++){
      for(int i = 0; i < j; i++){
	foo += D[k][i] * D[i][j];
	//	cout<<"Dki "<<D[k][i]<<" Dij "<<D[i][j]<<endl;
	//	cout<<foo<<endl;
      }
      //      cout<<"Dkj = "<<D[k][j]<<" - "<<foo<<endl;
      D[k][j] = D[k][j] - foo;
      foo = 0.0;
    }
    //find pivot
    int pivot = j;
    double max = fabs(D[j][j]);
    for(int k = j; k < size; k++){
      if(fabs(D[k][j]) > max){
	pivot = k;
	max = fabs(D[k][j]);
      }
    }
    // if pivot = 0, there is a zero on the diag.  
    // return false
    if(D[pivot][j] == 0){
      det = 0;
      cout<<" Exiting, 0 on diag "<<endl;
      return false;
    }
    //if necessary, swap rows
    if(pivot > j){
      det = -det;
      for(int swap = 0; swap < size * 2; swap++){
	double temp;
	temp = D[j][swap];
	D[j][swap] = D[pivot][swap];
	D[pivot][swap] = temp;
      }
    }
    for(int k = j+1; k < size; k++){
      double bar = 0.0;
      for(int i = 0; i < j ; i++){
	bar += D[j][i] * D[i][k];
      }
      D[j][k] = (D[j][k] - bar) / D[j][j];
      bar = 0.0;
    }
    det = det * D[j][j];
    cout<<"det "<<det<<"  End of run "<<" J "<<j<<endl;
    printD();
    cout<<endl;

  }// end j
  return true;
}

double bobMatrix::twoNorm(){
  double err = 0.0;
  double sol = 0.0;
  for(int i = 0; i < size; i++){
    err += s_x[i] * s_x[i];
    sol += s[i] * s[i];
  }
  return sqrt(err) / sqrt(sol);
}


void bobMatrix::makeS_X(){
  s_x.clear();
  for(int i = 0; i < size; i++){
    s_x.push_back(s[i] - x[i]);
  }
}

//precondition: makeUT run first
void bobMatrix::backSub(){
  x.resize(size,0.0);
  double temp = 0.0;
  for(int j = size-1; j >= 0; j--){
    for(int i = j+1; i < size; i++){
      temp += C[j][i] * x[i];
    }
    x[j] = (C[j][size] - temp)/C[j][j];
    temp = 0.0;
  }
}

// preconditions: C, A exist
bool bobMatrix::makeUT(){
  // find pivot
  for(int j = 0; j < size; j++){
    int pivot = j;
    double max = fabs(C[j][j]);
    for(int k = j; k < size; k++){
      if(fabs(C[k][j]) > max){
	pivot = k;
	max = fabs(C[k][j]);
      }
    }
    // if pivot = 0, there is a zero on the diag.  
    // return false
    if(C[pivot][j] == 0){return false;}
    //if necessary, swap rows
    if(pivot > j){
      for(int swap = 0; swap < size + 1; swap++){
	double temp;
	temp = C[j][swap];
	C[j][swap] = C[pivot][swap];
	C[pivot][swap] = temp;
      }
    }
    // multiply row j by scalar and subtract from other rows
    for(int row = j+1; row < size; row++){
      double scalar = C[row][j] / C[j][j];
      for(int col = j; col < size+1; col++){
	C[row][col] = C[row][col] - (scalar * C[j][col]);
      }
    }
  }
  return true;
}

//precondition: A and b exist
void bobMatrix::makeC(){
  C.resize(size);
  for(int i = 0; i < size; i++){
    C[i].resize(size + 1, 0.0);
  }
  for(int row = 0; row < size; row++){
    for(int col = 0; col < size; col++){
      C[row][col] = A[row][col];
    }
    C[row][size] = b[row];
  }
}
 
// precondition: A and s exist
void bobMatrix::findB(){
  double bTemp = 0.0;
  for(int row = 0; row < size; row++){
    for(int col = 0; col < size; col++){
      bTemp += (A[row][col] * s[col]);
    }
    b.push_back(bTemp);
    bTemp = 0.0;
  }
}

// precondition: none
void bobMatrix::genSol(){
  for(int i = 0; i < size; i++){
    s.push_back(randomDouble());
  }
}

//precondition: none
bobMatrix::bobMatrix(int N){
  size = N;
  A.resize(size);
  for(int i = 0; i < size; i++){
    A[i].resize(size,0.0);
  }
  for(int i = 0; i < size;i++){
    for(int j = 0; j < size; j++){
      A[i][j] = randomDouble();
    }
  }
  for(int row = 0; row < size; row++){
    double sum = 0.0;
    for(int col = 0; col < size; col++){
      if(row != col){
	sum += fabs(A[row][col]);
      }
    }
    A[row][row] = sum+100;
    sum = 0.0;
  }
}

void bobMatrix::printC(){
  cout<<"Matrix C  (A + b)"<<endl;
  for(int row = 0; row < size; row++){
    for(int col = 0; col < size + 1; col++){
      cout<<setw(15)<<setprecision(7)<<C[row][col];
    }
    cout<<endl;
  }
}
	    
void bobMatrix::prints(){
  cout<<"Vector S (generated)"<<endl;
  for(int i = 0; i < size; i++){
    cout<<setw(15)<<setprecision(20)<<s[i]<<endl;
  }
}

void bobMatrix::printb(){
  cout<<"Vector B (A*s)"<<endl;
  for(int i = 0; i < size; i++){
    cout<<setw(15)<<setprecision(7)<<b[i]<<endl;
  }
}

void bobMatrix::printx(){
  cout<<" x (Computed)"<<endl;
  for(int i = 0; i < size; i++){
    cout<<setw(15)<<setprecision(20)<<x[i]<<endl;
  }
}

// precondition: A exists
void bobMatrix::printA(){
  cout<<"   A  (generated)"<<endl;
  for(int row = 0; row< size; row++){
    for(int col = 0; col < size;col++){
      cout<<setw(15)<<setprecision(7)<<A[row][col];
    }
    cout<<endl;
  }
}

void bobMatrix::printS_X(){
  cout<<"S_X (error vector)"<<endl;
  for(int i = 0; i < size; i++){
    cout<<setw(15)<<setprecision(7)<<s_x[i]<<endl;
  }
}
union bDouble{
  int intD[2];
  double dubD;
};

void bobMatrix::printD(){
  cout<<" D (A & I) "<<endl;
  for(int row = 0; row < size; row++){
    for(int col = 0; col < size*2; col++){
      cout<<setw(15)<<D[row][col];
    }
    cout<<endl;
  }
}
