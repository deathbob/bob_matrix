#include <cassert>
#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>

using namespace std;

double random_double(){
  return fmod(double(rand()),1000000);
};

class bob_mini_matrix{
  int r;
  int c;
  vector< vector<double> > data;
  void fill(char fill_type);

 public:
  void print();
  bob_mini_matrix& operator=(const bob_mini_matrix& right);
  bool operator==(const bob_mini_matrix right);
  bob_mini_matrix operator*(const bob_mini_matrix right);
  bob_mini_matrix();
  bob_mini_matrix(const bob_mini_matrix& X);
  bob_mini_matrix(int N, char type);
  bob_mini_matrix(int rows, int columns, char type);
  bob_mini_matrix(int rows, int cols, double Tarr[]);
  bob_mini_matrix& operator-=(const bob_mini_matrix& right); 
};//end bob_mini_matrix



bob_mini_matrix operator-(bob_mini_matrix left, bob_mini_matrix right){
  bob_mini_matrix result = left;
  return left -= right;
}

bob_mini_matrix& bob_mini_matrix::operator-=(const bob_mini_matrix& right){
  assert(r = right.r);
  assert(c = right.c);
  for(int row = 0; row < r; row++){
    for(int col = 0; col < c; col++){
      data[row][col] = data[row][col] - right.data[row][col];
    }
  }
  return *this;
}
  


void bob_mini_matrix::print(){
  cout<<"r = "<<r<<"  c = "<<c<<endl;
  for(int row = 0;row < r;row++){
    for(int col = 0; col < c; col++){
      cout<<setw(15)<<setprecision(7)<<data[row][col];
    }
    cout<<endl;
  }
};

void bob_mini_matrix::fill(char fill_type){
  for(int row = 0; row < r; row++){
    for(int col = 0; col < c; col++){
      if(fill_type == 'i'){
	if(row == col){data[row][col] = 1;}
	else{data[row][col] = 0;}
      }
      if(fill_type == 'b'){data[row][col] = 0;}
      if(fill_type == 'r'){data[row][col] = random_double();}
    }
  }
};



bob_mini_matrix& bob_mini_matrix::operator=(const bob_mini_matrix& right){
  data = right.data;
  r = right.r;
  c = right.c;
  return *this;
};

bool bob_mini_matrix::operator==(const bob_mini_matrix right){
  if((c == right.c)&&(r == right.r)&&(data == right.data))
    return true;
  return false;
};
 
bob_mini_matrix bob_mini_matrix::operator*(const bob_mini_matrix right){
  assert(c == right.r);
  bob_mini_matrix result(r, right.c);
  double temp = 0.0;
  for(int row = 0; row < r; row++){
    //	for(int rowRight = 0; rowRight < right.r; rowRight++){
    for(int colRight = 0; colRight < right.c; colRight++){
      for(int col = 0; col < c; col++){
	temp += data[row][col] * right.data[col][colRight];
      }
      result.data[row][colRight] = temp;
      temp = 0.0;
      //	}
    }
  }
  return result;
};

bob_mini_matrix::bob_mini_matrix(){
  r = 0;
  c = 0;
};

bob_mini_matrix::bob_mini_matrix(const bob_mini_matrix& X){
  r = X.r;
  c = X.c;
  data = X.data;
};

bob_mini_matrix::bob_mini_matrix(int N, char fill_type){
  r = N;
  c = N; 
  data.resize(N);
  for(int i = 0; i < N; i++){
    data[i].resize(N,0.0);
  }
  fill(fill_type);
};

bob_mini_matrix::bob_mini_matrix(int rows, int cols, char fill_type){
  r = rows;
  c = cols;
  data.resize(rows);
  for(int i = 0; i < rows; i++){
    data[i].resize(cols, 0.0);
  }
  fill(fill_type);
};

bob_mini_matrix::bob_mini_matrix(int rows, int cols, double Tarr[]){
  r = rows;
  c = cols;
  data.resize(rows);
  for(int i = 0; i < rows; i++){
    data[i].resize(cols, 0.0);
  }
  int foo = 0;
  for(int i = 0; i < rows; i++){
    for(int j = 0; j < cols; j++){
      data[i][j] = Tarr[foo];
      foo++;
    }
  }
};
   
