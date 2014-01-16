#include "RungeKutta.hh"

RungeKutta::RungeKutta(vector<double> (*func)(double,vector<double>),double step_size,int dimension){
  m_func = func;
  m_step_size = step_size;
  m_t = 0;
  for(int i=0; i<dimension; i++){
    m_state.push_back(0);
  }
}


RungeKutta::RungeKutta(vector<double> (*func)(double,vector<double>),double step_size,
		       double t, vector<double> state){
  m_func = func;
  m_step_size = step_size;
  m_t = t;
  m_state = state;
}

void RungeKutta::PerformStep(){
  vector<double> k1 = m_func(m_t,m_state);
  vector<double> k2 = m_func(m_t+m_step_size/2,VectorAdd(m_state,VectorScale(k1,m_step_size/2)));
  vector<double> k3 = m_func(m_t+m_step_size/2,VectorAdd(m_state,VectorScale(k2,m_step_size/2)));
  vector<double> k4 = m_func(m_t+m_step_size,VectorAdd(m_state,VectorScale(k3,m_step_size)));
  m_t += m_step_size;
  vector<double> new_state;
  for(int i=0; i<m_state.size(); i++){
    new_state.push_back( m_state[i] + (k1[i]+2*k2[i]+2*k3[i]+k4[i])*m_step_size/6);
  }
  m_state = new_state;
}
