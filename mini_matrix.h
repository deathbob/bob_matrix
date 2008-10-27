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
  void transpose();
  bob_mini_matrix& operator=(const bob_mini_matrix& right);
  bool operator==(const bob_mini_matrix right);
  bob_mini_matrix operator*(const bob_mini_matrix right);
  bob_mini_matrix();
  bob_mini_matrix(const bob_mini_matrix& X);
  bob_mini_matrix(int N, char type);
  bob_mini_matrix(int rows, int columns, char type);
  bob_mini_matrix(int rows, int cols, double Tarr[]);
  bob_mini_matrix& operator-=(const bob_mini_matrix& right); 
  bob_mini_matrix& scale(double d);

  void power_method();
  void jacobi_method();
  void householder();
  void QR_method();
  double row_sum_norm();

  
};//end bob_mini_matrix

double bob_mini_matrix::row_sum_norm(){
  double sum = 0.0;
  double max = 0.0;
  for(int row = 0; row < r; row++){
    for(int col = 0; col < c; col++){
      sum += fabs(data[row][col]);
    }
    if (sum > max){ max = sum; }
    sum = 0.0;
  }
  return max;
}
  

void bob_mini_matrix::power_method(){
  int iters = 100;
  double epsilon = 0.0000001;
  int k = 0;
  // y is eigenvector estimate
  bob_mini_matrix y(r,1,'i');
  y.data[1][0] = 1;
  cout<<" Y "<<endl;
  y.print();
  //mu is eigenvalue estimate
  double mu = 0.0;
  bob_mini_matrix x = *this * y;
  bob_mini_matrix resid_err(y);
  do{
    cout<<"******** Iteration "<<k<<" *********"<<endl;
    y = x.scale(1.0/x.row_sum_norm());
    x = *this * y;
    bob_mini_matrix yt(y);
    yt.transpose();
    bob_mini_matrix numer(yt);
    numer = numer * x;
    bob_mini_matrix denom(yt);
    denom = denom * y;
    double t1 = numer.data[0][0];
    double t2 = denom.data[0][0];
    mu = t1 / t2;
    resid_err = y;
    resid_err.scale(mu);
    resid_err-=x;
    k++;
    cout<<" mu (eigenvalue) "<<mu<<endl;
    cout<<" y  (eigenvector) "<<endl;
    y.print();
  }while((resid_err.row_sum_norm() > epsilon)&&(k < iters));
}


bob_mini_matrix& bob_mini_matrix::scale(double d){
  for(int row = 0; row < r; row++){
    for(int col = 0; col < c; col++){
      data[row][col] = data[row][col] * d;
    }
  }
  return *this;
}

int sgn(double a){
  if(a >= 0) return 1;
  else return -1;
}

void bob_mini_matrix::householder(){
  bob_mini_matrix Q(r, 'i');
  double temp = 0.0;
  double alpha = 0.0;
  double temparr[3] = {0.0};
  for(int k = 0; k < r-2; k++){
    for(int j = k+1; j < r;j++){
      temp += (data[j][k] * data[j][k]);
    }
    alpha = sgn(data[k+1][k]) * (sqrt(temp));
    for(int i = k+1; i < r; i++){
      if(i == k+1)temparr[i] = data[i][k] + alpha;
      else temparr[i] = data[i][k];
    }
    bob_mini_matrix u(r,1,temparr);
    bob_mini_matrix ut(u);
    ut.transpose();
    bob_mini_matrix P(r, 'i');
    bob_mini_matrix numer;
    bob_mini_matrix denom;
    numer = u * ut;
    numer.scale(2.0);
    denom = ut * u;
    double foo = 1.0 / denom.data[0][0];
    numer.scale(foo);
    P -= numer;
    Q = Q * P;
    *this = P * *this;
    *this = *this * P;
  }
}

void bob_mini_matrix::QR_method(){
  householder();
  //  print();
  double epsilon = 0.000001;
  double m = 10;
  int i = 0;
  do{
    bob_mini_matrix Qt(r,'i');
    for(int k = 0; k < r-1;k++){
      double c = data[k][k] / (sqrt((data[k][k] * data[k][k])+(data[k+1][k] * data[k+1][k])));
      double s = data[k+1][k] / (sqrt((data[k][k] * data[k][k])+(data[k+1][k] * data[k+1][k])));
      bob_mini_matrix P(r, 'i');
      P.data[k][k] = c;
      P.data[k+1][k+1] = c;
      P.data[k+1][k] = -s;
      P.data[k][k+1] = s;
      
      *this = P * *this;
      Qt = P * Qt;
    }
    bob_mini_matrix Q(Qt);
    Q.transpose();
    *this = *this * Q;
    i++;
  }while(i < m);
  cout<<"Eigenvalues are "<<endl;
  for(int row = 0; row < r; row++){
    cout<<data[row][row]<<", ";
  }
  cout<<endl;
}

bob_mini_matrix operator-(bob_mini_matrix left, bob_mini_matrix right){
  bob_mini_matrix result = left;
  return result -= right;
}

void bob_mini_matrix::transpose(){
  bob_mini_matrix temp(c, r, 'b');
  for(int row = 0; row < r; row++){
    for(int col = 0; col < c; col++){
      temp.data[col][row] = data[row][col];
    }
  }
  swap(r, c);
  data = temp.data;
}

void bob_mini_matrix::jacobi_method(){
  // only for symmetric matrices :(  
  assert(r == c);
  double epsilon = .000000001;
  bob_mini_matrix P(r, 'i');
  bob_mini_matrix A(*this);
  double phi = 0.0;
  double max = 0.0;
  int p = 0;
  int q = 0;
  for(int row = 0; row < r; row++){
    for(int col = row+1; col < c; col++){
      if(fabs(data[row][col]) > max){
	max = data[row][col];
	p = row;
	q = col;
      }
    }
  }
  int counter = 0;
  while(fabs(max) > epsilon){
    cout<<"Iteration "<<counter<<endl;
    counter++;
    double t1 = 2.0 * data[p][q];
    double t2 = data[p][p] - data[q][q];
    double tphi = .5 * atan(t1/t2);
    //    phi = .5 * atan2(2.0 * data[p][q], (data[p][p] - data[q][q]));
    //    phi = .5 * atan2(t1,t2);
    phi = tphi;
    bob_mini_matrix R(r, 'i');
    R.data[q][q] = cos(phi);
    R.data[p][p] = R.data[q][q];
    R.data[p][q] = -sin(phi);
    R.data[q][p] = -R.data[p][q];
    P = P * R;
    A = A * R;
    R.transpose();
    A = R * A;
    cout<<endl;
    *this = A;
    max = 0.0;
    for(int row = 0; row < r; row++){
      for(int col = row+1; col < c; col++){
	if(fabs(data[row][col]) > max){
	  max = data[row][col];
	  p = row;
	  q = col;
	}
      }
    }
  }
  //  print();
  cout<<"Eigenvalues are ";
  for(int i = 0; i < r; i++){
    cout<<data[i][i]<<", ";
  }
  cout<<endl;
  cout<<endl;
  //  P.print();
  cout<<"Eigenvectors are "<<endl;
    for(int row = 0; row < r; row++){
      for(int col = 0; col < c; col++){
	cout<<P.data[col][row]<<", ";
      }
      cout<<" transposed "<<endl;
    }
  
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
  //  cout<<"r = "<<r<<"  c = "<<c<<endl;
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
  bob_mini_matrix result(r, right.c, 'b');
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
   
