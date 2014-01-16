#ifndef RUNGEKUTTA_H
#define RUNGEKUTTA_H

#include <vector>
#include <iostream>
using namespace std;

class RungeKutta{
public:
  RungeKutta(vector<double> (*func)(double,vector<double>),double step_size,int dimension=1);
  RungeKutta(vector<double> (*func)(double,vector<double>),double step_size,double t, vector<double> state);
  void SetInitial(double t, vector<double> state);
  double GetTime(){return m_t;}
  vector<double>* GetState(){return &m_state;}
  double GetState(int dim){return m_state[dim];}
  void PerformStep();

private:
  vector<double> (*m_func)(double,vector<double>);
  double m_step_size;
  double m_t;
  vector<double> m_state;
};

template <class T>
vector<T> VectorAdd(vector<T> a, vector<T> b){
  if(a.size()!=b.size()){
    throw 1; //Should subclass exception, but I'm lazy.
  }
  vector<T> output;
  for(int i=0; i<a.size(); i++){
    output.push_back(a[i]+b[i]);
  }
  return output;
}

template <class T, class U>
vector<T> VectorScale(vector<T> a, U b){
  vector<T> output;
  for(int i=0; i<a.size(); i++){
    output.push_back(a[i]*b);
  }
  return output;
}

#endif
