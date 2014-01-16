#include <iostream>
#include <cmath>

#include "RungeKutta.hh"

using namespace std;

const double RK_STEP_SIZE = 0.1;
double f_x(double t,double x){
  return -x;
}
double f_deriv_x(double t,double deriv_x){
  return deriv_x;
}
double RK_increment(double (*func)(double,double), double &t, double delta){
  double k1 = 0;
  double k2 = 0;
  double k3 = 0;
  double k4 = 0;
  k1 = func(t,delta);
  k2 = func(t+RK_STEP_SIZE/2,delta+RK_STEP_SIZE*k1/2);
  k3 = func(t+RK_STEP_SIZE/2,delta+RK_STEP_SIZE*k2/2);
  k4 = func(t+RK_STEP_SIZE,delta+RK_STEP_SIZE*k3);
  delta = RK_STEP_SIZE/6*(k1+2*k2+2*k3+k4);
  t += RK_STEP_SIZE;
  return delta;
}

vector<double> f_vector_x(double t, vector<double> x){
  vector<double> output;
  output.push_back(x[1]);
  output.push_back(-x[0]);
  return output;
}

int main(){
  double t_final = 100;
  int iterations = t_final/RK_STEP_SIZE;

  double t = 0;
  double x = 0;
  double deriv_x = 1;
  for(int i=0;i<iterations;i++){
    deriv_x += RK_increment(f_x,t,x);
    t -= RK_STEP_SIZE;
    x += RK_increment(f_deriv_x,t,deriv_x);
  }
  cout << t << " " << x << " " << sin(t)-x << endl;

  vector<double> state;
  state.push_back(0);
  state.push_back(1);
  RungeKutta rk = RungeKutta(f_vector_x,RK_STEP_SIZE,0,state);
  for(int i=0;i<iterations;i++){
    rk.PerformStep();
  }
  cout << rk.GetTime() << " " << rk.GetState(0) << " " << sin(rk.GetTime()) - rk.GetState(0) << endl;
}
