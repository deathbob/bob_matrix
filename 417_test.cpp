#include<iostream>
#include<limits>
#include<stdlib.h>
#include<iomanip>
#include<cmath>
#include<vector>

using namespace std;

union bDouble{
    int intD[2];
    double dubD;
};
    

double randomDouble(){
    bDouble bd;
    bd.intD[0] = rand();
    bd.intD[1] = rand();
    if(bd.dubD > 0){bd.dubD-=1;}
    else{bd.dubD += 1;}
    return bd.dubD;
};

class bobMatrix{
private:
    int size;
    vector< vector<double> > A;
    vector< vector<double> > C;
    vector<double> b;
    vector<double> s;
    vector<double> computedX;

public:
    bobMatrix(int N);
    void genSol();
    void findB();
    void printA();
    void makeC();
    bool makeUT();
    void backSub();
    void printcx();
    void printb();
    void prints();
    void printC();

};// end matrix

void printC(){
    cout<<"Matrix C"<<end;
    for(row = 0; row < size; row++){
	for(col = 0; col < size + 1; col++){
	    cout<<setw(15)<<C[row][col]<<endl;
	}
	cout<<endl;
    }
}
	    

void bobMatrix::prints(){
    cout<<"Vector S"<<endl;
    for(int i = 0; i < size; i++){
	cout<<setw(15)<<s[i]<<endl;
    }
}

void bobMatrix::printb(){
    cout<<"Vector B"<<endl;
    for(int i = 0; i < size; i++){
	cout<<setw(15)<<b[i]<<endl;
    }
}

void bobMatrix::printcx(){
    cout<<"Computed X"<<endl;
    for(int i = 0; i < size; i++){
	cout<<setw(15)<<computedX[i]<<endl;
    }
}

//precondition: makeUT run first
void bobMatrix::backSub(){
    computedX.resize(size,0.0);
    double temp = 0;
    for(int j = size; j > 0; j--){
	for(int i = j+1; i < size; i++){
	    temp += C[j][i] * computedX[i];
	}
	computedX[j] = (C[j][size+1] - temp)/C[j][j];
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

// precondition A exists
void bobMatrix::printA(){
    cout<<"   A"<<endl;
    for(int row = 0; row< size; row++){
	for(int col = 0; col < size;col++){
	    cout<<setw(15)<<A[row][col];
	}
	cout<<endl;
    }
}

// precondition: A and s exist
void bobMatrix::findB(){
    double bTemp = 0;
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
}










int main(){
    bobMatrix A(3);
    return 0;

}
